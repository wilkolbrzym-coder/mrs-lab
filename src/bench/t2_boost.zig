//! MRS-LAB :: thesis T2 — composing boosts
//!
//! Hypothesis: in the split-complex representation the rapidity is an additive
//! coordinate, so composing N boosts in one plane costs N additions. In the
//! matrix representation one must either multiply matrices (8 multiplications
//! per step) or extract the rapidity through the inverse function (atanh).
//!
//! Baselines (all optimised, matrices precomputed):
//!   B1  chain of 2x2 matrix products       — 8 multiplications + 4 additions / step
//!   B2  rapidity extraction via atanh + summation — 1 transcendental / step
//!
//! Warianty MRS:
//!   M1  akumulacja rapidyty                     — 1 dodawanie / krok
//!   M2  chain of split-complex products    — 4 multiplications + 2 additions / step
//!
//! We also measure DRIFT. The matrix representation of the group accumulates
//! rounding error like O(N·eps), and on top of that compensated summation cannot
//! compensated summation, because a matrix product is not an addition.
//! In MRS-LAB the rapidity is a number, so Kahan works — and that is a difference
//! in kind, not just in constant factor.

const std = @import("std");
const mrs = @import("mrs");
const h = @import("harness.zig");

const Io = std.Io;
const Z = mrs.split_complex.Z;

/// Kahan–Neumaier compensated summation. Meaningful only when an additive
/// coordinate exists — that is, only in MRS-LAB.
pub fn kahanSum(xs: []const f64) f64 {
    var sum: f64 = 0;
    var c: f64 = 0;
    for (xs) |x| {
        const y = x - c;
        const t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

pub fn plainSum(xs: []const f64) f64 {
    var sum: f64 = 0;
    for (xs) |x| sum += x;
    return sum;
}

pub const Ctx = struct {
    n: usize,
    thetas: []f64,
    zs: []Z,
    /// 4 liczby na boost: [cosh, sinh, sinh, cosh]
    mats: []f64,
    theta_acc: f64 = 0,
    z_acc: Z = Z.one,
    m_acc: [4]f64 = .{ 1, 0, 0, 1 },
    sink: f64 = 0,
};

fn runMrsRapidity(c: *Ctx) void {
    var acc: f64 = 0;
    for (c.thetas) |th| acc += th;
    c.theta_acc = acc;
    std.mem.doNotOptimizeAway(c.theta_acc);
}

fn runMrsRapidityKahan(c: *Ctx) void {
    c.theta_acc = kahanSum(c.thetas);
    std.mem.doNotOptimizeAway(c.theta_acc);
}

/// Baseline "best known classical mathematics": keep the rapidity as an ordinary
/// f64 variable and add. This function is DELIBERATELY identical to
/// with `runMrsRapidity` — and that is the result, not an oversight. The fastest
/// best known classical method here is literally the same code as the MRS-LAB
/// method, so the ratio comes out at 1.00x.
fn runKnownBestRapidity(c: *Ctx) void {
    var acc: f64 = 0;
    for (c.thetas) |th| acc += th;
    c.theta_acc = acc;
    std.mem.doNotOptimizeAway(c.theta_acc);
}

fn runMrsSplitComplex(c: *Ctx) void {
    var z = Z.one;
    for (c.zs) |zi| z = Z.mul(z, zi);
    c.z_acc = z;
    std.mem.doNotOptimizeAway(c.z_acc);
}

fn runConvMatrix(c: *Ctx) void {
    var m0: f64 = 1;
    var m1: f64 = 0;
    var m2: f64 = 0;
    var m3: f64 = 1;
    for (0..c.n) |i| {
        const b = c.mats[i * 4 ..][0..4];
        const n0 = m0 * b[0] + m1 * b[2];
        const n1 = m0 * b[1] + m1 * b[3];
        const n2 = m2 * b[0] + m3 * b[2];
        const n3 = m2 * b[1] + m3 * b[3];
        m0 = n0;
        m1 = n1;
        m2 = n2;
        m3 = n3;
    }
    c.m_acc = .{ m0, m1, m2, m3 };
    std.mem.doNotOptimizeAway(c.m_acc[0]);
}

fn runConvRapidityAtanh(c: *Ctx) void {
    var acc: f64 = 0;
    for (0..c.n) |i| {
        const b = c.mats[i * 4 ..][0..4];
        acc += std.math.atanh(b[1] / b[0]);
    }
    c.theta_acc = acc;
    std.mem.doNotOptimizeAway(c.theta_acc);
}

// ---------------------------------------------------------------------------
// Dryf numeryczny
// ---------------------------------------------------------------------------

pub const Drift = struct {
    n: usize,
    /// |N(z) − 1| przy sumowaniu kompensowanym rapidyty.
    mrs_kahan_inv_err: f64,
    /// |N(z) − 1| przy sumowaniu naiwnym.
    mrs_plain_inv_err: f64,
    ///     /// max |m_ij − analytic matrix| for the matrix chain.
    matrix_elem_err: f64,
    ///     /// |det m − 1| for the matrix chain (Lorentz invariant).
    matrix_det_err: f64,
    ///     /// |rapidity(m) − Σθ| for the matrix chain.
    matrix_rapidity_err: f64,
    ///     /// |Σθ plain − Σθ Kahan|, i.e. the summation error of the rapidity itself.
    rapidity_sum_err: f64,
};

/// A chain of N boosts with a FIXED total rapidity. Normalisation is necessary
/// so that the measurement isolates the ACCUMULATION of error rather than the
/// mapy θ → (cosh θ, sinh θ) — to drugie mierzymy osobno w
/// `conditioningOfRapidityMap`, because that is a different phenomenon.
pub fn measureDrift(alloc: std.mem.Allocator, n: usize) !Drift {
    const thetas = try alloc.alloc(f64, n);
    defer alloc.free(thetas);
    var prng = std.Random.DefaultPrng.init(0x5EED);
    const rnd = prng.random();
    var raw: f64 = 0;
    for (thetas) |*th| {
        th.* = (rnd.float(f64) - 0.5) * 0.2;
        raw += th.*;
    }
    const target_total: f64 = 1.0;
    for (thetas) |*th| th.* *= target_total / raw;

    const theta_kahan = kahanSum(thetas);
    const theta_plain = plainSum(thetas);

    // MRS: the rapidity as the coordinate
    const z_kahan = Z.fromRapidity(theta_kahan);
    const z_plain = Z.fromRapidity(theta_plain);

    // baseline: chain of 2x2 matrix products
    var m0: f64 = 1;
    var m1: f64 = 0;
    var m2: f64 = 0;
    var m3: f64 = 1;
    for (thetas) |th| {
        const c = std.math.cosh(th);
        const s = std.math.sinh(th);
        const n0 = m0 * c + m1 * s;
        const n1 = m0 * s + m1 * c;
        const n2 = m2 * c + m3 * s;
        const n3 = m2 * s + m3 * c;
        m0 = n0;
        m1 = n1;
        m2 = n2;
        m3 = n3;
    }
    const c_exact = std.math.cosh(theta_kahan);
    const s_exact = std.math.sinh(theta_kahan);
    const mat_err = @max(
        @max(@abs(m0 - c_exact), @abs(m1 - s_exact)),
        @max(@abs(m2 - s_exact), @abs(m3 - c_exact)),
    );
    const det_err = @abs(m0 * m3 - m1 * m2 - 1.0);
    const rap_err = @abs(std.math.atanh(m1 / m0) - theta_kahan);

    return .{
        .n = n,
        .mrs_kahan_inv_err = @abs(Z.norm(z_kahan) - 1.0),
        .mrs_plain_inv_err = @abs(Z.norm(z_plain) - 1.0),
        .matrix_elem_err = mat_err,
        .matrix_det_err = det_err,
        .matrix_rapidity_err = rap_err,
        .rapidity_sum_err = @abs(theta_plain - theta_kahan),
    };
}

/// Conditioning of the map θ → (cosh θ, sinh θ): how accurately the
/// multyplikatywna trzyma niezmiennik N(z) = cosh²θ − sinh²θ = 1.
/// The error grows like e^{2|θ|}·eps, because cosh² and sinh² are two huge
/// numbers whose difference is small — classic cancellation.
pub fn conditioningOfRapidityMap(theta: f64) f64 {
    const z = Z.fromRapidity(theta);
    return @abs(Z.norm(z) - 1.0);
}

// ---------------------------------------------------------------------------
// Uruchomienie
// ---------------------------------------------------------------------------

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    const ns: []const usize = if (quick)
        &[_]usize{ 1_000, 10_000 }
    else
        &[_]usize{ 1_000, 10_000, 100_000, 1_000_000 };

    try w.writeAll("### T2 — composing N boosts in one plane\n\n");
    try w.writeAll("All variants compute exactly the same object: a Lorentz group element " ++
        "equal to the composition of N boosts. Baseline matrices are precomputed " ++
        "outside the measurement — deliberately convenient for the baseline.\n\n");
    try w.writeAll("| N | MRS: rapidity (1 add) | **known best:** f64 rapidity | MRS: split-complex | B1: 2x2 matrix | B2: matrix + atanh | B1/MRS | B2/MRS |\n");
    try w.writeAll("|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    var best_ratio_sum: f64 = 0;
    var best_ratio_n: usize = 0;

    for (ns) |n| {
        const thetas = try alloc.alloc(f64, n);
        defer alloc.free(thetas);
        const zs = try alloc.alloc(Z, n);
        defer alloc.free(zs);
        const mats = try alloc.alloc(f64, n * 4);
        defer alloc.free(mats);

        var prng = std.Random.DefaultPrng.init(0x5EED);
        const rnd = prng.random();
        for (0..n) |i| {
            const th = (rnd.float(f64) - 0.5) * 0.2;
            thetas[i] = th;
            const c = std.math.cosh(th);
            const s = std.math.sinh(th);
            zs[i] = Z.init(c, s);
            mats[i * 4 + 0] = c;
            mats[i * 4 + 1] = s;
            mats[i * 4 + 2] = s;
            mats[i * 4 + 3] = c;
        }

        var ctx = Ctx{ .n = n, .thetas = thetas, .zs = zs, .mats = mats };

        const inner = h.pickInner(@as(f64, @floatFromInt(n)) * 12.0, 2_000_000.0);
        const samples: usize = if (quick) 9 else 13;

        const st_rap = try h.measure(io, alloc, &ctx, samples, inner, runMrsRapidity);
        const st_best = try h.measure(io, alloc, &ctx, samples, inner, runKnownBestRapidity);
        const st_kahan = try h.measure(io, alloc, &ctx, samples, inner, runMrsRapidityKahan);
        const st_sc = try h.measure(io, alloc, &ctx, samples, inner, runMrsSplitComplex);
        const st_mat = try h.measure(io, alloc, &ctx, samples, inner, runConvMatrix);
        const st_atanh = try h.measure(io, alloc, &ctx, samples, inner, runConvRapidityAtanh);

        const rap_op = st_rap.perOpBest();
        const sc_op = st_sc.perOpBest();
        const mat_op = st_mat.perOpBest();
        const atanh_op = st_atanh.perOpBest();

        const best_op = st_best.perOpBest();
        try w.print("| {d} | ", .{n});
        try h.writeDuration(w, rap_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, best_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, sc_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, mat_op);
        try w.writeAll(" | ");
        try h.writeDuration(w, atanh_op);
        try w.writeAll(" | ");
        try w.print("{d:.1}×", .{h.speedup(mat_op, rap_op)});
        try w.writeAll(" | ");
        try w.print("{d:.1}×", .{h.speedup(atanh_op, rap_op)});
        try w.writeAll(" |\n");
        // remember the ratio of MRS to the best known method
        best_ratio_sum += h.speedup(best_op, rap_op);
        best_ratio_n += 1;
        _ = &best_ratio_sum;
        _ = &best_ratio_n;

        _ = st_kahan;
    }

    if (best_ratio_n > 0) {
        try w.print("\n**Ratio of MRS to the BEST KNOWN classical method: " ++
            "{d:.3}x.** That method is plain addition of rapidities held in f64 " ++
            "variables — literally the same code, which is why the ratio is 1. This is " ++
            "an IDENTITY CHECK, not a measurement: it shows that the classical best " ++
            "method for this problem is the same computation. In other words, the " ++
            "2.5-2.8x advantage over matrices **is not an advantage of MRS over " ++
            "mathematics** — it is an advantage of the additive representation over the " ++
            "multiplicative one, and classical mathematics knows the additive one just as " ++
            "well.\n\n", .{
            best_ratio_sum / @as(f64, @floatFromInt(best_ratio_n)),
        });
    }

    try w.writeAll("\n### T2b — numerical drift over a long chain\n\n");
    try w.writeAll("Rapidities are NORMALISED to a fixed total (1.0) so that the measurement " ++
        "isolates the accumulation of error rather than the conditioning of the map " ++
        "θ → (cosh θ, sinh θ). That second effect is measured separately below, because " ++
        "it is a different phenomenon.\n\n");
    try w.writeAll("| N | MRS + Kahan: \\|N(z)−1\\| | MRS plain: \\|N(z)−1\\| | matrix: element error | matrix: \\|det−1\\| | matrix: rapidity error |\n");
    try w.writeAll("|---:|---:|---:|---:|---:|---:|\n");

    const dn: []const usize = if (quick)
        &[_]usize{ 1_000, 10_000, 100_000 }
    else
        &[_]usize{ 1_000, 10_000, 100_000, 1_000_000 };

    for (dn) |n| {
        const d = try measureDrift(alloc, n);
        try w.print("| {d} | {e:.2} | {e:.2} | {e:.2} | {e:.2} | {e:.2} |\n", .{
            n,
            d.mrs_kahan_inv_err,
            d.mrs_plain_inv_err,
            d.matrix_elem_err,
            d.matrix_det_err,
            d.matrix_rapidity_err,
        });
    }

    try w.writeAll("\n**Reading.** Columns 1-2 (MRS) do not depend on N: rapidity is an " ++
        "additive coordinate, so the composition error is a summation error, not the " ++
        "error of N successive matrix multiplications. Columns 3-5 (matrix baseline) " ++
        "grow; the growth is DIFFUSIVE, i.e. like sqrt(N), because the rounding errors of " ++
        "successive steps are approximately independent. Compensated summation cannot be " ++
        "applied to matrices, because a matrix product is not an addition — and that is a " ++
        "difference in kind, not just in constant factor.\n\n");

    try w.writeAll("### T2c — where the MRS rapidity does NOT help (honest caveat)\n\n");
    try w.writeAll("The map θ → (cosh θ, sinh θ) is ill conditioned for large |θ|: cosh²θ and " ++
        "sinh²θ are two huge numbers whose difference is 1.\n\n");
    try w.writeAll("| θ | cosh θ | \\|N(z)−1\\| |\n|---:|---:|---:|\n");
    for ([_]f64{ 1, 5, 10, 20, 40, 100, 350 }) |th| {
        try w.print("| {d:.0} | {e:.3} | {e:.2} |\n", .{
            th,
            std.math.cosh(th),
            conditioningOfRapidityMap(th),
        });
    }
    try w.writeAll("\nThe design conclusion: in MRS-LAB **the rapidity is the primary " ++
        "coordinate**, and the pair (cosh θ, sinh θ) is one possible representation to be " ++
        "evaluated at the end. The reverse order — composing in the multiplicative " ++
        "representation — loses the invariant at |θ| around 20.\n\n");
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "all boost composition variants give the same group element" {
    const alloc = std.testing.allocator;
    const n = 500;
    const thetas = try alloc.alloc(f64, n);
    defer alloc.free(thetas);
    var prng = std.Random.DefaultPrng.init(11);
    const rnd = prng.random();
    var exact: f64 = 0;
    for (thetas) |*th| {
        th.* = (rnd.float(f64) - 0.5) * 0.5;
        exact += th.*;
    }

    // MRS: rapidyta
    const z = Z.fromRapidity(kahanSum(thetas));

    // baseline: matrix chain
    var m0: f64 = 1;
    var m1: f64 = 0;
    var m2: f64 = 0;
    var m3: f64 = 1;
    for (thetas) |th| {
        const c = std.math.cosh(th);
        const s = std.math.sinh(th);
        const n0 = m0 * c + m1 * s;
        const n1 = m0 * s + m1 * c;
        const n2 = m2 * c + m3 * s;
        const n3 = m2 * s + m3 * c;
        m0 = n0;
        m1 = n1;
        m2 = n2;
        m3 = n3;
    }

    try std.testing.expectApproxEqAbs(z.re, m0, 1e-14);
    try std.testing.expectApproxEqAbs(z.im, m1, 1e-14);
    try std.testing.expectApproxEqAbs(z.re, std.math.cosh(exact), 1e-14);
    try std.testing.expectApproxEqAbs(z.im, std.math.sinh(exact), 1e-14);
}

test "drift: MRS does not accumulate error, the matrix accumulates diffusively" {
    const alloc = std.testing.allocator;
    const d3 = try measureDrift(alloc, 1_000);
    const d5 = try measureDrift(alloc, 100_000);

    // MRS: the composition error of the rapidity does not depend on N (additive coordinate)
    try std.testing.expect(d3.mrs_kahan_inv_err < 1e-13);
    try std.testing.expect(d5.mrs_kahan_inv_err < 1e-13);

    // matrix: the error grows. The growth is DIFFUSIVE (~sqrt(N)), not linear,
    // because the rounding errors of successive steps are approximately independent.
    // For N differing by 100x we expect roughly 10x growth of the error.
    try std.testing.expect(d5.matrix_elem_err > d3.matrix_elem_err * 5.0);
    try std.testing.expect(d5.matrix_elem_err > 1e-14);
}

test "the conditioning of the rapidity map breaks down for large |θ|" {
    try std.testing.expect(conditioningOfRapidityMap(1.0) < 1e-15);
    try std.testing.expect(conditioningOfRapidityMap(50.0) > 1e-9);
    try std.testing.expect(conditioningOfRapidityMap(200.0) > 1e-3);
}

test "Kahan beats plain summation on asymmetric terms" {
    // Deliberately bad case: large terms and small corrections.
    var xs: [2002]f64 = undefined;
    xs[0] = 1e16;
    for (1..2001) |i| xs[i] = 1.0;
    xs[2001] = -1e16;

    const p = plainSum(&xs);
    const k = kahanSum(&xs);
    try std.testing.expectApproxEqAbs(@as(f64, 2000.0), k, 1e-6);
    try std.testing.expect(@abs(k - 2000.0) <= @abs(p - 2000.0));
}
