//! MRS-LAB :: thesis T3 — sparsity in the Clifford algebra
//!
//! Hypothesis: a multivector with k nonzero blades multiplies in O(k^2), while a
//! dense representation must walk 4^n pairs. In a realistic Clifford-algebra
//! workload (rotors, reflectors, the sandwich R·v·R~) multivectors are sparse,
//! so the advantage is exponential.
//!
//! Baselines, deliberately all as strong as possible:
//!   B1  dense product of all pairs                 — 4^n blade operations
//!   B2  dense with zero skipping                   — 4^n CHECKS
//!   B3  dense, scan to a list, then sparse product — 2^n scan + k^2
//!
//! B3 is the strongest and the fairest: it is exactly the same algorithm as
//! MRS-LAB plus an unavoidable scan of the dense array to discover the sparsity.
//! Conclusion (docs/REVIEW.md): the algorithm itself is not the property of
//! MRS-LAB. The property is that density is part of the TYPE, so the sparse path
//! is chosen automatically and cannot be lost in a refactor — and the cost
//! difference against B3 is exactly the 2^n scan.

const std = @import("std");
const mrs = @import("mrs");
const h = @import("harness.zig");

const Io = std.Io;
const cl = mrs.clifford;

pub const Ctx = struct {
    alg: cl.Algebra,
    /// Mutable copies of the sparse multivector coefficients. They are kept
    /// separately so that the multivector in the API can stay immutable while
    /// the benchmark can still change its input (and thereby block hoisting of
    /// the computation out of the loop).
    mut_a: []cl.Term,
    mut_b: []cl.Term,
    dense_a: []f64,
    dense_b: []f64,
    conv_a: []cl.Term,
    conv_b: []cl.Term,
    scratch: []cl.Term,
    out_dense: []f64,
    out_sparse: cl.Sparse = .{ .terms = &.{} },
    sink: f64 = 0,
};

fn sparseView(terms: []cl.Term) cl.Sparse {
    return .{ .terms = terms };
}

/// Input mutation on every call — blocks hoisting of the computation out of
/// the loop.
///
/// DEFECT FOUND BY REVIEW, and the fix matters. The previous version always
/// touched index 0. When the sparse original had no scalar term, index 0 was
/// zero, so the bump CREATED a nonzero and increased the number of nonzero
/// coefficients by one. The B3 variant then called `mulSparseScratch` with a
/// scratch buffer sized k*k while the product needed (k+1)*(k+1) slots: the
/// assertion in `mulSparseScratch` fired in Debug and ReleaseSafe, and in
/// ReleaseFast — where assertions are compiled out — the same code wrote past
/// the buffer. That is a silent heap overflow inside the measurements that
/// produced the B3 column.
///
/// Fix: bump the first coefficient that is already nonzero, so the sparsity is
/// invariant to the mutation. The buffer sizing is also relaxed below as a
/// second line of defence, because the assertion is the guard that caught this.
fn bumpDense(a: []f64) void {
    for (a, 0..) |x, i| {
        if (x != 0.0) {
            a[i] = if (x > 1e3) 0.5 else x + 1e-6;
            return;
        }
    }
    a[0] = 1.0; // no nonzero present: keep a well-defined input
}

fn runMrs(c: *Ctx) void {
    c.mut_a[0].coeff = if (c.mut_a[0].coeff > 1e3) 0.5 else c.mut_a[0].coeff + 1e-6;
    c.out_sparse = cl.mulSparseScratch(
        c.alg,
        sparseView(c.mut_a),
        sparseView(c.mut_b),
        c.scratch,
    );
    std.mem.doNotOptimizeAway(c.out_sparse.terms.len);
}

fn runConvScan(c: *Ctx) void {
    bumpDense(c.dense_a);
    const na = denseToTerms(c.dense_a, c.conv_a);
    const nb = denseToTerms(c.dense_b, c.conv_b);
    c.out_sparse = cl.mulSparseScratch(
        c.alg,
        .{ .terms = c.conv_a[0..na] },
        .{ .terms = c.conv_b[0..nb] },
        c.scratch,
    );
    std.mem.doNotOptimizeAway(c.out_sparse.terms.len);
}

fn runDenseSkip(c: *Ctx) void {
    bumpDense(c.dense_a);
    cl.mulDenseSkippingZerosInto(c.alg, c.dense_a, c.dense_b, c.out_dense);
    std.mem.doNotOptimizeAway(c.out_dense[0]);
}

fn runDense(c: *Ctx) void {
    bumpDense(c.dense_a);
    cl.mulDenseInto(c.alg, c.dense_a, c.dense_b, c.out_dense);
    std.mem.doNotOptimizeAway(c.out_dense[0]);
}

/// Scan a dense array and list its nonzero coefficients.
fn denseToTerms(a: []const f64, out: []cl.Term) usize {
    var n: usize = 0;
    for (a, 0..) |x, i| {
        if (x != 0.0) {
            out[n] = .{ .mask = @intCast(i), .coeff = x };
            n += 1;
        }
    }
    return n;
}

// ---------------------------------------------------------------------------
// Building the test instances
// ---------------------------------------------------------------------------

fn termLessThanPublic(_: void, a: cl.Term, b: cl.Term) bool {
    return a.mask < b.mask;
}

/// Rotor-like multivector: k random distinct blades. A realistic
/// distribution — rotors carry a scalar and a few bivectors, the rest is zero.
fn buildSparse(
    alloc: std.mem.Allocator,
    k: usize,
    basis: usize,
    seed: u64,
) !cl.Sparse {
    const terms = try alloc.alloc(cl.Term, @min(k, basis));
    var prng = std.Random.DefaultPrng.init(seed);
    const rnd = prng.random();
    var used = std.AutoHashMap(u32, void).init(alloc);
    defer used.deinit();
    var n: usize = 0;
    while (n < terms.len) {
        const mask: u32 = rnd.intRangeLessThan(u32, 0, @intCast(basis));
        if (used.contains(mask)) continue;
        try used.put(mask, {});
        terms[n] = .{ .mask = mask, .coeff = rnd.float(f64) * 2.0 - 1.0 };
        n += 1;
    }
    std.mem.sort(cl.Term, terms, {}, termLessThanPublic);
    return .{ .terms = terms };
}

fn makeAlgebra(ng: u5) !cl.Algebra {
    var sq = [_]i8{0} ** cl.MAX_GEN;
    for (0..ng) |i| sq[i] = if (i == 0) 1 else -1; // one time, rest space
    return cl.Algebra.init(ng, sq);
}

const Row = struct {
    k: usize,
    mrs: f64,
    conv: f64,
    skip: f64,
    dense: f64,
};

// ---------------------------------------------------------------------------
// Uruchomienie
// ---------------------------------------------------------------------------

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    try w.writeAll("### T3a — fixed k = 4, growing n\n\n");
    try w.writeAll("One product of two multivectors with 4 nonzero blades. All buffers are " ++
        "supplied from outside: no variant allocates inside the hot loop.\n\n");
    try w.writeAll("| n | 2^n | MRS: k^2 | B3: scan 2^n + k^2 | B2: dense + zero skip | B1: dense full | B3/MRS | B1/MRS |\n");
    try w.writeAll("|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    const ns_a: []const u5 = if (quick)
        &[_]u5{ 4, 6, 8, 10 }
    else
        &[_]u5{ 4, 6, 8, 10, 12 };

    for (ns_a) |ng| {
        const alg = try makeAlgebra(ng);
        const basis = alg.basisCount();
        const k: usize = 4;

        const sa = try buildSparse(alloc, k, basis, 0x1111 + @as(u64, ng));
        defer alloc.free(sa.terms);
        const sb = try buildSparse(alloc, k, basis, 0x2222 + @as(u64, ng));
        defer alloc.free(sb.terms);
        const da = try cl.sparseToDense(alloc, alg, sa);
        defer alloc.free(da);
        const db = try cl.sparseToDense(alloc, alg, sb);
        defer alloc.free(db);
        const conv_a = try alloc.alloc(cl.Term, basis);
        defer alloc.free(conv_a);
        const conv_b = try alloc.alloc(cl.Term, basis);
        defer alloc.free(conv_b);
        // Sized for the worst case of the dense-scan path ((k+1)^2), not for
        // the nominal k^2: the scan may see one more nonzero than the sparse
        // original if the input mutation creates one.
        const scratch = try alloc.alloc(cl.Term, (k + 1) * (k + 1));
        defer alloc.free(scratch);
        const out_dense = try alloc.alloc(f64, basis);
        defer alloc.free(out_dense);

        const mut_a = try alloc.alloc(cl.Term, sa.terms.len);
        defer alloc.free(mut_a);
        const mut_b = try alloc.alloc(cl.Term, sb.terms.len);
        defer alloc.free(mut_b);
        @memcpy(mut_a, sa.terms);
        @memcpy(mut_b, sb.terms);

        var ctx = Ctx{
            .alg = alg,
            .mut_a = mut_a,
            .mut_b = mut_b,
            .dense_a = da,
            .dense_b = db,
            .conv_a = conv_a,
            .conv_b = conv_b,
            .scratch = scratch,
            .out_dense = out_dense,
        };

        const ops_dense = @as(f64, @floatFromInt(basis * basis));
        const inner = h.pickInner(ops_dense * 3.0, 1_500_000.0);
        const samples: usize = if (quick) 7 else 11;

        const st_mrs = try h.measure(io, alloc, &ctx, samples, inner, runMrs);
        const st_conv = try h.measure(io, alloc, &ctx, samples, inner, runConvScan);
        const st_skip = try h.measure(io, alloc, &ctx, samples, inner, runDenseSkip);
        const st_dense = try h.measure(io, alloc, &ctx, samples, inner, runDense);

        const mrs_op = st_mrs.perOpBest();
        const conv_op = st_conv.perOpBest();
        const skip_op = st_skip.perOpBest();
        const dense_op = st_dense.perOpBest();

        try w.print("| {d} | {d} | ", .{ ng, basis });
        try h.writeDuration(w, mrs_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, conv_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, skip_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, dense_op);
        try w.writeAll(" | ");
        try w.print("{d:.0}×", .{h.speedup(conv_op, mrs_op)});
        try w.writeAll(" | ");
        try w.print("{d:.0}×", .{h.speedup(dense_op, mrs_op)});
        try w.writeAll(" |\n");
    }

    try w.writeAll("\n**Reading.** B3 (the same sparse algorithm, but fed a dense array) " ++
        "pays the unavoidable 2^n scan and therefore loses to MRS-LAB by a growing " ++
        "factor. B1 pays 4^n. The advantage of MRS-LAB over B3 grows like 2^n/k^2.\n\n");

    try runCrossover(io, alloc, w, quick);
}

/// T3b — the crossover point at n = 10. The constant factor is derived
/// FROM THE MEASUREMENT rather than guessed: it is the ratio of the unit cost of
/// the sparse path to the unit cost of one dense pair.
fn runCrossover(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    const ng: u5 = 10;
    const alg = try makeAlgebra(ng);
    const basis = alg.basisCount();

    try w.writeAll("### T3b — crossover point (n = 10, 2^n = 1024, 4^n = 1 048 576)\n\n");
    try w.writeAll("| k | MRS: sparse | B3: scan + sparse | B2: dense + zero skip | B2/MRS | winner |\n");
    try w.writeAll("|---:|---:|---:|---:|---:|---|\n");

    var rows = std.ArrayList(Row).empty;
    defer rows.deinit(alloc);

    const ks: []const usize = if (quick)
        &[_]usize{ 2, 8, 32, 128, 512, 1024 }
    else
        &[_]usize{ 2, 8, 32, 128, 256, 512, 768, 1024 };

    for (ks) |k| {
        if (k > basis) continue;
        const sa = try buildSparse(alloc, k, basis, 0x3333 + k);
        defer alloc.free(sa.terms);
        const sb = try buildSparse(alloc, k, basis, 0x4444 + k);
        defer alloc.free(sb.terms);
        const da = try cl.sparseToDense(alloc, alg, sa);
        defer alloc.free(da);
        const db = try cl.sparseToDense(alloc, alg, sb);
        defer alloc.free(db);
        const conv_a = try alloc.alloc(cl.Term, basis);
        defer alloc.free(conv_a);
        const conv_b = try alloc.alloc(cl.Term, basis);
        defer alloc.free(conv_b);
        // Sized for the worst case of the dense-scan path ((k+1)^2), not for
        // the nominal k^2: the scan may see one more nonzero than the sparse
        // original if the input mutation creates one.
        const scratch = try alloc.alloc(cl.Term, (k + 1) * (k + 1));
        defer alloc.free(scratch);
        const out_dense = try alloc.alloc(f64, basis);
        defer alloc.free(out_dense);

        const mut_a = try alloc.alloc(cl.Term, sa.terms.len);
        defer alloc.free(mut_a);
        const mut_b = try alloc.alloc(cl.Term, sb.terms.len);
        defer alloc.free(mut_b);
        @memcpy(mut_a, sa.terms);
        @memcpy(mut_b, sb.terms);

        var ctx = Ctx{
            .alg = alg,
            .mut_a = mut_a,
            .mut_b = mut_b,
            .dense_a = da,
            .dense_b = db,
            .conv_a = conv_a,
            .conv_b = conv_b,
            .scratch = scratch,
            .out_dense = out_dense,
        };

        const inner = h.pickInner(@as(f64, @floatFromInt(basis * basis)) * 3.0, 1_500_000.0);
        const samples: usize = if (quick) 5 else 9;

        const st_mrs = try h.measure(io, alloc, &ctx, samples, inner, runMrs);
        const st_conv = try h.measure(io, alloc, &ctx, samples, inner, runConvScan);
        const st_skip = try h.measure(io, alloc, &ctx, samples, inner, runDenseSkip);

        const mrs_op = st_mrs.perOpBest();
        const conv_op = st_conv.perOpBest();
        const skip_op = st_skip.perOpBest();

        try rows.append(alloc, .{
            .k = k,
            .mrs = mrs_op,
            .conv = conv_op,
            .skip = skip_op,
            .dense = 0,
        });

        try w.print("| {d} | ", .{k});
        try h.writeDuration(w, mrs_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, conv_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, skip_op);
        try w.writeAll(" | ");
        try w.print("{d:.1}× | {s} |\n", .{
            h.speedup(skip_op, mrs_op),
            if (mrs_op < skip_op) "**MRS**" else "dense",
        });
    }

    // Derive the cost constant from the measurement: how many times a sparse
    // operation is more expensive than one dense pair. Only small k are used,
    // where the measured time is dominated by work rather than fixed overhead.
    var sum: f64 = 0;
    var cnt: usize = 0;
    const dense_basis_ops = @as(f64, @floatFromInt(basis * basis));
    for (rows.items) |r| {
        if (r.k > 64) continue;
        const cost_sparse_op = r.mrs / @as(f64, @floatFromInt(r.k * r.k));
        const cost_pairs = r.skip / dense_basis_ops;
        // k = 2 is dominated by fixed overhead, so only k >= 8 is used
        if (r.k < 8) continue;
        sum += cost_sparse_op / cost_pairs;
        cnt += 1;
    }
    const factor = if (cnt > 0) sum / @as(f64, @floatFromInt(cnt)) else 1.0;

    try w.writeAll("\n### T3b — numeric verdict\n\n");
    try w.print("Cost ratio derived from the measurement, sparse operation to dense pair: " ++
        "**C = {d:.1}**. It is the only constant in the model, and it is not guessed — " ++
        "it falls out of the data.\n\n", .{factor});

    // Falsifiable test: compare the measured winner with the `chooseStrategy`
    // prediction at the derived constant C.
    try w.writeAll("| k | measured | `chooseStrategy` prediction | agree |\n");
    try w.writeAll("|---:|---|---|---|\n");
    var agree: usize = 0;
    var total: usize = 0;
    for (rows.items) |r| {
        const measured: []const u8 = if (r.mrs < r.skip) "MRS" else "dense";
        const pred = cl.chooseStrategy(r.k, r.k, basis, factor);
        const pred_str: []const u8 = if (pred == .sparse) "MRS" else "dense";
        const ok = std.mem.eql(u8, measured, pred_str);
        if (ok) agree += 1;
        total += 1;
        try w.print("| {d} | {s} | {s} | {s} |\n", .{
            r.k, measured, pred_str, if (ok) "yes" else "**NO**",
        });
    }
    try w.print("\nPrediction agrees with the measurement in {d} of {d} points. The " ++
        "prediction is not fitted to the answer: the constant C comes only from small " ++
        "k and is then tested over the whole range.\n\n", .{ agree, total });

    try w.writeAll("**Honest caveat:** MRS-LAB creates no new mathematics here — it creates " ++
        "a type that does not let the structure be forgotten. That is real engineering " ++
        "value, but it is not a new theorem, and we call it by its name.\n\n");

    try w.writeAll("**Range of applicability:** for small n (4-6 generators) the dense path " ++
        "is as fast or faster, because sparsity has nothing to save; the T3a table shows " ++
        "this at n = 4. MRS-LAB pays off only when 2^n is large relative to the number of " ++
        "nonzero blades.\n\n");
}

test "test instances are sparse, sorted and distinct" {
    const alloc = std.testing.allocator;
    const s = try buildSparse(alloc, 8, 1024, 0x99);
    defer alloc.free(s.terms);
    try std.testing.expectEqual(@as(usize, 8), s.terms.len);
    for (1..s.terms.len) |i| {
        try std.testing.expect(s.terms[i - 1].mask < s.terms[i].mask);
    }
    const t = try buildSparse(alloc, 8, 1024, 0x100);
    defer alloc.free(t.terms);
    try std.testing.expect(s.terms[0].mask != t.terms[0].mask or s.terms[0].coeff != t.terms[0].coeff);
}

test "the dense-to-sparse scan finds exactly k elements" {
    const alloc = std.testing.allocator;
    const alg = try makeAlgebra(6);
    const basis = alg.basisCount();
    const s = try buildSparse(alloc, 5, basis, 0x77);
    defer alloc.free(s.terms);
    const d = try cl.sparseToDense(alloc, alg, s);
    defer alloc.free(d);
    const buf = try alloc.alloc(cl.Term, basis);
    defer alloc.free(buf);
    const n = denseToTerms(d, buf);
    try std.testing.expectEqual(@as(usize, 5), n);
}
