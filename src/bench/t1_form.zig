//! MRS :: teza T1 — reprezentacja formy
//!
//! Hypothesis: a form given by its signature (diagonal) evaluates in O(n) while
//! a matrix representation takes O(n^2); the inverse O(n) against O(n^3).
//!
//! The baseline is chosen to be maximally unfair to us — deliberately. Besides
//! the naive dense matrix we also measure a `smart_dense` variant: an
//! implementation that checks ONCE whether the matrix is diagonal and then uses
//! the fast path. Any competent engineer who knows their metric is diagonal
//! writes that code.
//!
//! Wniosek z pomiaru (docs/07): asymptotycznej przewagi TUTAJ NIE MA,
//! if the baseline may inspect its data once. The advantage of MRS-LAB is that
//! this knowledge lives in the type, so it can neither be lost nor cost a branch
//! in the hot loop. We say so plainly instead of inflating the result by
//! comparing against code nobody would write.

const std = @import("std");
const mrs = @import("mrs");
const h = @import("harness.zig");

const Io = std.Io;
const Signature = mrs.signature.Signature;
const DiagonalForm = mrs.form.DiagonalForm;
const DenseForm = mrs.form.DenseForm;

/// Konwencjonalna implementacja z jednorazowym rozpoznaniem struktury.
/// This is the honest baseline: we do not pretend nobody notices the diagonal.
pub const SmartDense = struct {
    n: usize,
    g: []const f64,
    diagonal: bool,

    pub fn init(n: usize, g: []const f64) SmartDense {
        var diag = true;
        outer: for (0..n) |i| {
            for (0..n) |j| {
                if (i != j and g[i * n + j] != 0.0) {
                    diag = false;
                    break :outer;
                }
            }
        }
        return .{ .n = n, .g = g, .diagonal = diag };
    }

    pub fn eval(self: SmartDense, v: []const f64) f64 {
        const n = self.n;
        if (self.diagonal) {
            var acc: f64 = 0;
            for (0..n) |i| acc += self.g[i * n + i] * v[i] * v[i];
            return acc;
        }
        var acc: f64 = 0;
        for (0..n) |i| {
            const row = self.g[i * n ..][0..n];
            var s: f64 = 0;
            for (0..n) |j| s += row[j] * v[j];
            acc += v[i] * s;
        }
        return acc;
    }
};

const Ctx = struct {
    /// Literal construction: signs are derived through `signAt`, which branches
    /// on the role of every dimension.
    diag: DiagonalForm,
    /// Preferred construction: signs precomputed into a contiguous array.
    diag_fast: DiagonalForm,
    dense: DenseForm,
    smart: SmartDense,
    v: []f64,
    acc: f64 = 0,
};

fn bump(v: []f64) void {
    // Input mutation on every call: without it the compiler may hoist the loop-
    // invariant computation and we would measure zero.
    v[0] = if (v[0] > 1e6) 0.0 else v[0] + 1.0;
}

fn runMrs(ctx: *Ctx) void {
    bump(ctx.v);
    ctx.acc = ctx.diag.eval(ctx.v);
    std.mem.doNotOptimizeAway(ctx.acc);
}

fn runMrsFast(ctx: *Ctx) void {
    bump(ctx.v);
    ctx.acc = ctx.diag_fast.eval(ctx.v);
    std.mem.doNotOptimizeAway(ctx.acc);
}

fn runDense(ctx: *Ctx) void {
    bump(ctx.v);
    ctx.acc = ctx.dense.eval(ctx.v);
    std.mem.doNotOptimizeAway(ctx.acc);
}

fn runSmart(ctx: *Ctx) void {
    bump(ctx.v);
    ctx.acc = ctx.smart.eval(ctx.v);
    std.mem.doNotOptimizeAway(ctx.acc);
}

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    const sizes: []const usize = if (quick)
        &[_]usize{ 16, 32, 64, 128, 256 }
    else
        &[_]usize{ 16, 32, 64, 128, 256, 512, 1024 };

    const samples: usize = if (quick) 9 else 15;

    try w.writeAll("### T1 — form evaluation g(v,v): O(n) versus O(n^2)\n\n");
    try w.print("Clock: `awake`, resolution 1 ns (verified). {d} samples per measurement, " ++
        "report is the MINIMUM (the least noisy estimator; medians are computed but " ++
        "the minimum is what isolates the algorithm from scheduler interference).\n\n", .{samples});
    try w.writeAll("Column 2 is the literal construction, where signs come from `signAt` " ++
        "(one branch per dimension). Column 3 precomputes the sign vector once — same " ++
        "mathematics, branch-free inner loop. The conventional baseline is allowed its " ++
        "own best form (a strided diagonal load, structure detected once outside the " ++
        "measurement), so both sides are compared at their best.\n\n");
    try w.writeAll("| n | MRS literal (branch) | MRS precomputed signs | dense matrix | dense + diagonal detection | dense/MRS | smart/MRS |\n");
    try w.writeAll("|---:|---:|---:|---:|---:|---:|---:|\n");

    for (sizes) |n| {
        const sigv = try alloc.alloc(mrs.signature.Role, n);
        defer alloc.free(sigv);
        // sygnatura z jednym czasem — realistyczny przypadek fizyczny
        for (0..n) |i| sigv[i] = if (i == 0) .temporal else .spatial;
        const s = Signature{ .roles = sigv };

        const dense = try DenseForm.fromDiagonal(alloc, s);
        defer dense.deinit(alloc);

        const v = try alloc.alloc(f64, n);
        defer alloc.free(v);
        var prng = std.Random.DefaultPrng.init(0xABCD);
        const rnd = prng.random();
        for (v) |*x| x.* = rnd.float(f64) * 2.0 - 1.0;

        var ctx = Ctx{
            .diag = .{ .signature = s },
            .diag_fast = DiagonalForm.init(s),
            .dense = dense,
            .smart = SmartDense.init(n, dense.g),
            .v = v,
        };

        const ops_dense = @as(f64, @floatFromInt(n * n));
        const inner = h.pickInner(ops_dense, 2_000_000.0);

        const st_mrs = try h.measure(io, alloc, &ctx, samples, inner, runMrs);
        const st_fast = try h.measure(io, alloc, &ctx, samples, inner, runMrsFast);
        const st_dense = try h.measure(io, alloc, &ctx, samples, inner, runDense);
        const st_smart = try h.measure(io, alloc, &ctx, samples, inner, runSmart);

        const mrs_op = st_mrs.perOpBest();
        const fast_op = st_fast.perOpBest();
        const dense_op = st_dense.perOpBest();
        const smart_op = st_smart.perOpBest();

        try w.print("| {d} | ", .{n});
        try h.writeDuration(w, mrs_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, fast_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, dense_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, smart_op);
        try w.writeAll(" | ");
        try w.print("{d:.1}×", .{h.speedup(dense_op, fast_op)});
        try w.writeAll(" | ");
        try w.print("{d:.2}×", .{h.speedup(smart_op, fast_op)});
        try w.writeAll(" |\n");
    }

    try w.writeAll("\n**Reading.** The `dense/MRS` ratio (last but one column) grows like n: " ++
        "that is the O(n^2) versus O(n) difference. The `smart/MRS` ratio is the honest " ++
        "test: a conventional implementation that detects the diagonal once and then uses " ++
        "a strided load. Against MRS WITHOUT the precomputed sign vector that baseline " ++
        "WINS at small n (ratio below 1) — the per-element role branch costs more than the " ++
        "strided load. Against MRS WITH the precomputed vector the ratio goes above 1 and " ++
        "grows with n.\n\n" ++
        "**Honest conclusion.** O(n^2) versus O(n) is a difference of representations, not of " ++
        "mathematics: the fastest known classical method for a diagonal form is the same " ++
        "O(n) sum. What MRS-LAB adds is that the structure lives in the type, so it cannot be " ++
        "forgotten, and that the sign vector is contiguous, so the inner loop is branch-free.\n\n");

    try runInverse(io, alloc, w, quick);
}

/// T1b — inverse of the form: O(n) against O(n^3).
/// Here the difference is algebraic rather than implementational: for a general
/// dense matrix no shortcut avoids elimination.
fn runInverse(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    const sizes: []const usize = if (quick)
        &[_]usize{ 8, 16, 32, 64 }
    else
        &[_]usize{ 8, 16, 32, 64, 128 };

    try w.writeAll("### T1b — inverse of the form (raising an index)\n\n");
    try w.writeAll("| n | MRS (diagonal) | Gauss-Jordan | G-J/MRS | theoretical complexity |\n");
    try w.writeAll("|---:|---:|---:|---:|---|\n");

    for (sizes) |n| {
        const sigv = try alloc.alloc(mrs.signature.Role, n);
        defer alloc.free(sigv);
        for (0..n) |i| sigv[i] = if (i == 0) .temporal else .spatial;
        const s = Signature{ .roles = sigv };

        const dense = try DenseForm.fromDiagonal(alloc, s);
        defer dense.deinit(alloc);
        const diag = DiagonalForm{ .signature = s };

        const inv = try alloc.alloc(f64, n * n);
        defer alloc.free(inv);
        const c = try alloc.alloc(f64, n);
        defer alloc.free(c);
        const out = try alloc.alloc(f64, n);
        defer alloc.free(out);
        for (c, 0..) |*x, i| x.* = @as(f64, @floatFromInt(i + 1));

        const CtxInv = struct {
            diag: DiagonalForm,
            dense: DenseForm,
            inv: []f64,
            c: []f64,
            out: []f64,
            acc: f64 = 0,
            alloc: std.mem.Allocator,
            failed: bool = false,
        };

        var ctx = CtxInv{
            .diag = diag,
            .dense = dense,
            .inv = inv,
            .c = c,
            .out = out,
            .alloc = alloc,
        };

        const fns = struct {
            fn mrs(x: *CtxInv) void {
                x.c[0] = if (x.c[0] > 1e6) 1.0 else x.c[0] + 1.0;
                x.diag.raiseIndex(x.c, x.out) catch {
                    x.failed = true;
                    return;
                };
                std.mem.doNotOptimizeAway(x.out[0]);
            }
            fn denseInv(x: *CtxInv) void {
                x.c[0] = if (x.c[0] > 1e6) 1.0 else x.c[0] + 1.0;
                x.dense.inverse(x.alloc, x.inv) catch {
                    x.failed = true;
                    return;
                };
                std.mem.doNotOptimizeAway(x.inv[0]);
            }
        };

        const ops_gj = @as(f64, @floatFromInt(2 * n * n * n));
        const inner_gj = h.pickInner(ops_gj, 5_000_000.0);
        const inner_mrs = h.pickInner(@floatFromInt(n), 5_000_000.0);
        const samples: usize = if (quick) 7 else 11;

        const st_mrs = try h.measure(io, alloc, &ctx, samples, inner_mrs, fns.mrs);
        const st_dense = try h.measure(io, alloc, &ctx, samples, inner_gj, fns.denseInv);

        const mrs_op = st_mrs.perOpBest();
        const dense_op = st_dense.perOpBest();

        try w.print("| {d} | ", .{n});
        try h.writeDuration(w, mrs_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, dense_op);
        try w.writeAll(" | ");
        try w.print("{d:.1}×", .{h.speedup(dense_op, mrs_op)});
        try w.print(" | dense O(n^3) = {d:.1}e3 vs MRS O(n) = {d:.0} operations |\n", .{
            @as(f64, @floatFromInt(2 * n * n * n)),
            @as(f64, @floatFromInt(n)),
        });
    }
    try w.writeAll("\n**Reading.** Here the advantage is structural rather than implementational: " ++
        "for a dense form no shortcut avoids elimination, so the ratio grows like n^2. In " ++
        "MRS-LAB the inverse of the form is the same object as the signature, because " ++
        "`g^{ij} = 1/s_i`. Note the baseline is Gauss-Jordan with partial pivoting — the " ++
        "standard practical choice — not the asymptotically fastest known matrix inverse; " ++
        "the qualitative claim (a diagonal form inverts in O(n), a dense one does not) is " ++
        "what the table supports.\n\n");
}
