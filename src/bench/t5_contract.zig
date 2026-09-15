//! MRS-LAB :: thesis T5 — the price of carrying an error radius
//!
//! The contract layer of 0.1.3 promises that the zero test stops being an
//! argument the caller supplies and becomes a consequence of the computation.
//! That promise is worth nothing if it costs throughput in the inner loop, so
//! this is a PRICE LIST, not a victory lap.
//!
//! What is compared, on the same form and the same vector:
//!
//!   plain    `f.eval(v)`                  — the 0.1.2 code path, no radius
//!   bounded  `contract.evalForm(f, v)`    — value plus propagated radius
//!
//! The third row is not a measurement and is not timed. For the contracts
//! `bit_exact` and `correctly_rounded` the layer returns a plain `f64` from
//! `evalFormWith`, whose body IS `f.eval(v)`. Those paths cannot cost anything
//! because there is no radius in the type — that is the death criterion of the
//! design discharged by construction rather than by measurement, and the test
//! `a contract that carries no radius produces the value bit for bit` guards it.
//!
//! The honest verdict is expected in the ratios column, and if the bounded path
//! is more than a small factor slower the layer should not be used in a hot
//! loop — the report says so with numbers instead of adjectives.

const std = @import("std");
const Io = std.Io;
const mrs = @import("mrs");
const h = @import("harness.zig");

const sig = mrs.signature;
const form = mrs.form;
const contract = mrs.contract;

const Ctx = struct {
    f: form.DiagonalForm,
    v: [1024]f64 = undefined,
    len: usize = 0,
    sink: f64 = 0,
};

fn runPlain(ctx: *Ctx) void {
    const g = ctx.f.eval(ctx.v[0..ctx.len]);
    ctx.sink += g;
    std.mem.doNotOptimizeAway(ctx.sink);
}

fn runBounded(ctx: *Ctx) void {
    const g = contract.evalForm(ctx.f, ctx.v[0..ctx.len]);
    ctx.sink += g.value + g.radius;
    std.mem.doNotOptimizeAway(ctx.sink);
}

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    const sizes: []const usize = if (quick)
        &[_]usize{ 16, 64, 256 }
    else
        &[_]usize{ 16, 64, 256, 1024 };
    const samples: usize = if (quick) 9 else 15;

    try w.writeAll("### T5 — the price of an error radius\n\n");
    try w.writeAll("Same form, same vector, same arithmetic; the only difference is " ++
        "that the bounded path also computes and propagates `|value − exact|`. " ++
        "Both columns include the loop-invariant form lookup, so the ratio isolates " ++
        "the radius and not the call.\n\n");
    try w.writeAll("| n | plain `eval` | bounded `evalForm` | ratio |\n");
    try w.writeAll("|---:|---:|---:|---:|\n");

    var worst: f64 = 0;
    for (sizes) |n| {
        // The signature must have EXACTLY n dimensions: `DiagonalForm` borrows
        // this slice, and `eval` reads one sign per component. Handing it a
        // four-dimensional signature and a kilobyte-long vector reads past the
        // sign cache — in ReleaseFast silently, in ReleaseSafe as a panic.
        // `bench-check` is what caught it.
        const roles = try alloc.alloc(sig.Role, n);
        defer alloc.free(roles);
        for (0..n) |i| roles[i] = if (i == 0) .temporal else .spatial;
        const s = sig.Signature{ .roles = roles };

        var ctx = Ctx{ .f = form.DiagonalForm.init(s), .len = n };
        for (0..n) |i| {
            // Deterministic, non-degenerate values: no zeros, no symmetry.
            ctx.v[i] = 1.0 + @as(f64, @floatFromInt(i)) * 0.125;
            if (i % 3 == 0) ctx.v[i] = -ctx.v[i];
        }
        const inner: usize = if (n >= 256) 2000 else 20000;

        const st_plain = try h.measure(io, alloc, &ctx, samples, inner, runPlain);
        const st_bounded = try h.measure(io, alloc, &ctx, samples, inner, runBounded);
        const plain_op = st_plain.perOpBest();
        const bounded_op = st_bounded.perOpBest();
        const ratio = if (plain_op > 0) bounded_op / plain_op else 0;
        if (ratio > worst) worst = ratio;

        try w.print("| {d} | ", .{n});
        try h.writeDuration(w, plain_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, bounded_op);
        try w.print(" | {d:.2}× |\n", .{ratio});
    }

    try w.print("\n**Verdict.** The bounded path costs {d:.2}× at worst over these sizes. " ++
        "That is the price of knowing how far the value may be from the truth, and it " ++
        "is paid only where a radius is asked for: `evalFormWith(.bit_exact, …)` returns " ++
        "a plain `f64` from the same body as `f.eval`, so the contracts that carry no " ++
        "radius cost nothing at all — by construction, not by an optimiser.\n\n", .{worst});

    try w.writeAll("What this does NOT say: that the bounded path is cheap in every " ++
        "setting. These are `n` up to 1024 on a diagonal form, one evaluation per " ++
        "iteration, no allocation. A longer per-value contract chain would multiply " ++
        "the constant, and the radius rules take absolute values, so a loop that " ++
        "branches on the radius pays a dependency the plain path does not have.\n\n");
}
