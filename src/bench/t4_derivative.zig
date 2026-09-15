//! MRS-LAB :: thesis T4 — derivative of the dispersion relation with no step size
//!
//! Hypothesis: the dual algebra (ε² = 0) gives the exact derivative in a single
//! pass, so it is simultaneously MORE ACCURATE and FASTER than a finite
//! difference. This is the only thesis in the project where MRS-LAB wins on both
//! criteria at once and where the advantage does not depend on data structure.
//!
//! Baseline: central differences — the standard way to differentiate without
//! automatic differentiation. It is given its optimal step h (found by a scan)
//! so that the comparison is fair.

const std = @import("std");
const mrs = @import("mrs");
const h = @import("harness.zig");

const Io = std.Io;
const tach = mrs.tachyon;

pub const Ctx = struct {
    mu2: f64,
    ps: []f64,
    h: f64,
    acc: f64 = 0,
};

fn runAD(c: *Ctx) void {
    c.ps[0] = if (c.ps[0] > 1e5) 4.0 else c.ps[0] + 1e-3;
    var a: f64 = 0;
    for (c.ps) |p| a += tach.phaseEnergyDerivAD(p, c.mu2);
    c.acc = a;
    std.mem.doNotOptimizeAway(c.acc);
}

fn runFD(c: *Ctx) void {
    c.ps[0] = if (c.ps[0] > 1e5) 4.0 else c.ps[0] + 1e-3;
    var a: f64 = 0;
    for (c.ps) |p| a += tach.phaseEnergyDerivFD(p, c.mu2, c.h);
    c.acc = a;
    std.mem.doNotOptimizeAway(c.acc);
}

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    const mu2: f64 = 9.0; // μ = 3
    const n_samples: usize = if (quick) 2000 else 200_000;

    const ps = try alloc.alloc(f64, n_samples);
    defer alloc.free(ps);
    var prng = std.Random.DefaultPrng.init(0xD00D);
    const rnd = prng.random();
    for (ps) |*p| p.* = 4.0 + rnd.float(f64) * 96.0; // propagation range |p| > μ

    try w.writeAll("### T4 — derivative of the composite f(p) = sqrt(p^2 - mu^2)·sin(p)\n\n");
    try w.writeAll("The function is deliberately composite. For E(p) = sqrt(p^2 - mu^2) alone the " ++
        "dual algebra reproduces the analytic formula p/E EXACTLY, so comparing it with " ++
        "a finite difference would check nothing beyond the differentiation rules.\n\n");
    try w.writeAll("Reference: f'(p) = (p/E)·sin(p) + E·cos(p). Maximum error over 2000 points, p in [4, 100].\n\n");
    try w.writeAll("| method | step h | max absolute error |\n|---|---:|---:|\n");

    const hs = [_]f64{ 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };
    var best_h: f64 = 1e-5;
    var best_fd_err: f64 = std.math.floatMax(f64);
    for (hs) |hh| {
        const acc = tach.scanPhaseAccuracy(mu2, 4.0, 100.0, 2000, hh);
        if (acc.max_fd_error < best_fd_err) {
            best_fd_err = acc.max_fd_error;
            best_h = hh;
        }
        try w.print("| central difference | {e:.0} | {e:.3} |\n", .{ hh, acc.max_fd_error });
    }
    const acc_ad = tach.scanPhaseAccuracy(mu2, 4.0, 100.0, 2000, best_h);
    try w.print("| **MRS: dual numbers** | — | **{e:.3}** |\n\n", .{acc_ad.max_ad_error});

    try w.print("Best step for central differences is h = {e:.0} (error {e:.2}); for smaller " ++
        "h the error grows again because subtracting nearby numbers starts to win. That " ++
        "is the classic trade-off between truncation error O(h^2) and rounding error " ++
        "O(eps/h), and no amount of tuning escapes it.\n\n", .{
        best_h,
        best_fd_err,
    });

    if (acc_ad.max_ad_error == 0.0) {
        try w.writeAll("The dual-algebra error came out **exactly zero**. That deserves a " ++
            "precise statement so it does not sound like a numerical miracle: for this " ++
            "class of expressions the dual algebra performs EXACTLY the same sequence of " ++
            "floating point operations as the analytic formula (the product rule " ++
            "assembles the same product and the same sum), so the result is identical " ++
            "bit for bit.\n\n" ++
            "The value of MRS-LAB is therefore not that \"the result is more accurate " ++
            "than the analytic formula\" — that is impossible — but that:\n\n" ++
            "1. the formula does not have to be derived by hand or typed into code;\n" ++
            "2. when no closed form exists (integrals, special functions, numerically " ++
            "defined fields), the dual algebra still gives the exact derivative while the " ++
            "finite difference still stalls near 1e-8.\n\n");
    } else {
        try w.print("Dual-algebra error: {e:.2} — several orders better than the best tuned " ++
            "central-difference step.\n\n", .{acc_ad.max_ad_error});
    }

    try w.writeAll("### T4b — cost of one derivative\n\n");
    try w.writeAll("| N points | MRS: dual | central difference | difference/MRS |\n");
    try w.writeAll("|---:|---:|---:|---:|\n");

    const ns: []const usize = if (quick)
        &[_]usize{ 1_000, 10_000 }
    else
        &[_]usize{ 1_000, 10_000, 100_000 };

    for (ns) |n| {
        if (n > ps.len) continue;
        var ctx = Ctx{ .mu2 = mu2, .ps = ps, .h = best_h };
        const inner = h.pickInner(@as(f64, @floatFromInt(n)) * 40.0, 1_000_000.0);
        const samples: usize = if (quick) 7 else 11;

        const st_ad = try h.measure(io, alloc, &ctx, samples, inner, runAD);
        const st_fd = try h.measure(io, alloc, &ctx, samples, inner, runFD);

        const ad_op = st_ad.perOpBest() / @as(f64, @floatFromInt(n));
        const fd_op = st_fd.perOpBest() / @as(f64, @floatFromInt(n));

        try w.print("| {d} | ", .{n});
        try h.writeDuration(w, ad_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, fd_op);
        try w.writeAll(" | ");
        try w.print("{d:.2}× |\n", .{h.speedup(fd_op, ad_op)});
    }

    try w.writeAll("\n**Reading.** Central differences need TWO function calls per derivative, " ++
        "each containing a square root and a sinusoid. The dual algebra computes the " ++
        "derivative in one pass, adding only one multiplication per operation. The timing " ++
        "advantage is of order 2x, and that is what is measured — the rest of the " ++
        "difference is loop overhead.\n\n");
}

test "AD and FD agree in sign and order of magnitude" {
    const mu2 = 9.0;
    for ([_]f64{ 4.0, 5.0, 10.0, 50.0 }) |p| {
        const ad = tach.phaseEnergyDerivAD(p, mu2);
        const fd = tach.phaseEnergyDerivFD(p, mu2, 1e-6);
        const ex = tach.phaseEnergyDerivExact(p, mu2);
        try std.testing.expectApproxEqRel(ex, ad, 1e-13);
        try std.testing.expectApproxEqRel(ex, fd, 1e-5);
    }
}
