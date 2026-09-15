//! MRS-LAB :: thesis T6 — the exact mode, pushed to its limit
//!
//! THE MODE THIS MEASURES. No floating point, no tolerance, and zero is exactly
//! zero: every coefficient is an i64, every sign is ±1 or 0, and a decision
//! procedure answers a question about ALL real vectors by evaluating a
//! polynomial on a finite determining set. Nothing here can be "approximately
//! right", which is the whole point — but the naive way to evaluate it scans
//! all m = 2^n coefficients of both factors although the determining set has at
//! most two nonzero ones.
//!
//! WHAT WAS REMOVED, and nothing else:
//!
//!   dense (baseline)  scans m·m = 1024 coefficient slots per pair, of which at
//!                     most four do anything; recomputes the norm of grid[j]
//!                     once per i, i.e. k times too often; and pays a bit loop
//!                     inside `bladeMul` for every product.
//!   sparse            walks only the nonzero terms (Terms), looks the product
//!                     up in a precomputed table, and computes each norm once.
//!
//! The arithmetic is the same term for term — a test asserts the two produce
//! identical verdicts AND identical witnesses on every signature the engine can
//! build — so this is a pure representation change on the decision path, not a
//! change of method. That distinction is the reason the result is allowed to be
//! called a speedup at all: nothing was traded away for it.
//!
//! The honest baseline is the DENSE path, because that is what MRS-LAB itself
//! did until 0.1.5. There is no external baseline to hide behind here.

const std = @import("std");
const Io = std.Io;
const mrs = @import("mrs");
const h = @import("harness.zig");
const sigs = @import("../explore/signatures.zig");
const exact = @import("../explore/exact.zig");
const nrm = @import("../explore/norm_mult.zig");

const Ctx = struct {
    alg: mrs.clifford.Algebra,
    grid: *[nrm.MAX_GRID]exact.IntVec,
    sink: usize = 0,
    /// Set when the kernel REFUSED the algebra instead of deciding it. The
    /// harness calls a `fn (*Ctx) void`, so a refusal cannot be propagated as
    /// an error — it is recorded here and checked after the measurement, which
    /// is the next best thing to a `try` and does not hide it.
    refused: bool = false,
};

fn runSparse(ctx: *Ctx) void {
    const v = nrm.checkScalarMultiplicative(ctx.alg, ctx.grid) catch {
        ctx.refused = true;
        return;
    };
    ctx.sink += if (v.holds) 1 else 0;
    std.mem.doNotOptimizeAway(ctx.sink);
}

fn runDense(ctx: *Ctx) void {
    const v = nrm.checkScalarMultiplicativeBrute(ctx.alg, ctx.grid);
    ctx.sink += if (v.holds) 1 else 0;
    std.mem.doNotOptimizeAway(ctx.sink);
}

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    try w.writeAll("### T6 — the exact mode: no floats, no tolerance, 0 is exactly 0\n\n");
    try w.writeAll("P2a decided over the determining set. The baseline is the dense path " ++
        "this project shipped until 0.1.5, not an outside implementation; a test asserts " ++
        "the two agree on verdict AND witness for all 55 signatures, so the column " ++
        "measures representation and not method.\n\n");
    try w.writeAll("| (p,q,r) | n | grid | pairs | dense | sparse | gain |\n");
    try w.writeAll("|---|---:|---:|---:|---:|---:|---:|\n");

    const cases = if (quick)
        &[_]sigs.Triple{
            .{ .p = 1, .q = 3 },
            .{ .p = 0, .q = 4 },
            .{ .p = 0, .q = 0, .r = 4 },
        }
    else
        &[_]sigs.Triple{
            .{ .p = 1, .q = 3 },
            .{ .p = 0, .q = 4 },
            .{ .p = 2, .q = 2 },
            .{ .p = 0, .q = 1, .r = 3 },
            .{ .p = 0, .q = 0, .r = 4 },
            .{ .p = 1, .q = 4 },
            .{ .p = 0, .q = 5 },
            .{ .p = 0, .q = 0, .r = 5 },
        };

    const samples: usize = if (quick) 3 else 5;

    for (cases) |t| {
        const sb = sigs.SigBuf.build(t, .mostly_minus);
        const alg = try sb.algebra();
        const grid = try alloc.create([nrm.MAX_GRID]exact.IntVec);
        defer alloc.destroy(grid);

        var ctx = Ctx{ .alg = alg, .grid = grid };
        // One evaluation is already hundreds of thousands of pairs; the harness
        // repeats the whole decision rather than an inner step.
        const inner: usize = 1;

        const st_dense = try h.measure(io, alloc, &ctx, samples, inner, runDense);
        const st_sparse = try h.measure(io, alloc, &ctx, samples, inner, runSparse);
        if (ctx.refused) {
            std.debug.print("T6: the kernel refused ({d},{d},{d}); no timing is meaningful\n", .{ t.p, t.q, t.r });
            return error.TooManyBlades;
        }

        const dense = st_dense.perOpBest();
        const sparse = st_sparse.perOpBest();
        const gain = if (sparse > 0) dense / sparse else 0;
        const k = nrm.gridSize(alg.basisCount());

        try w.writeAll("| ");
        try t.writeTo(w);
        try w.print(" | {d} | {d} | {d} | ", .{ t.n(), k, k * k });
        try h.writeDuration(w, dense);
        try w.writeAll(" | ");
        try h.writeDuration(w, sparse);
        try w.print(" | **{d:.1}×** |\n", .{gain});
    }

    try w.writeAll("\n**Verdict.** The gain comes from not doing three things: scanning " ++
        "1024 coefficient slots per pair when at most four are nonzero, recomputing a " ++
        "norm once per pair when it does not depend on the other index, and paying a bit " ++
        "loop for a product that is one table lookup. What was NOT traded: the arithmetic " ++
        "is identical term for term, and the cross-check test requires identical verdicts " ++
        "and identical witnesses, so the exact mode still decides rather than estimates.\n\n");
    try w.writeAll("What this does not claim: that the deciding itself is now cheap. At " ++
        "n = 5 the count of pairs is fixed by the MATHEMATICS — |D|² = 279 841 — and no " ++
        "representation change moves that. The remaining cost is the bound worth " ++
        "attacking next, and widening MAX_BASIS to 64 for n = 6 multiplies it by sixteen.\n\n");
}
