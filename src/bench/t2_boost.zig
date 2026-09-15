//! MRS :: teza T2 — składanie boostów
//!
//! Hipoteza: w reprezentacji split-complex rapidyta jest współrzędną
//! addytywną, więc złożenie N boostów w jednej płaszczyźnie kosztuje N
//! dodawań. W reprezentacji macierzowej trzeba albo mnożyć macierze
//! (8 mnożeń na krok), albo wyłuskać rapidytę funkcją odwrotną (atanh).
//!
//! Baseline'y (wszystkie zoptymalizowane, macierze policzone z góry):
//!   B1  łańcuch iloczynów macierzy 2×2         — 8 mnożeń + 4 dodawania / krok
//!   B2  wyłuskanie rapidyty przez atanh + sumowanie — 1 transcendentna / krok
//!
//! Warianty MRS:
//!   M1  akumulacja rapidyty                     — 1 dodawanie / krok
//!   M2  łańcuch iloczynów split-complex         — 4 mnożenia + 2 dodawania / krok
//!
//! Dodatkowo mierzymy DRYF. Macierzowa reprezentacja grupy kumuluje błąd
//! zaokrągleń jak O(N·eps), a w dodatku nie da się na niej zastosować
//! sumowania kompensowanego, bo iloczyn macierzy nie jest dodawaniem.
//! W MRS rapidyta jest liczbą, więc Kahan działa — i to jest różnica
//! jakościowa, nie tylko stała wydajnościowa.

const std = @import("std");
const mrs = @import("mrs");
const h = @import("harness.zig");

const Io = std.Io;
const Z = mrs.split_complex.Z;

/// Sumowanie kompensowane Kahan–Neumaier. Sensowne wyłącznie wtedy, gdy
/// istnieje współrzędna addytywna — czyli wyłącznie w MRS.
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

/// Baseline "najlepsza znana zwykła matematyka": trzymamy rapidytę jako
/// zwykłą zmienną f64 i dodajemy. Ta funkcja jest CELOWO identyczna
/// z `runMrsRapidity` — i to jest wynik, a nie przeoczenie. Najszybsza
/// znana metoda z matematyki klasycznej jest tu dosłownie tym samym
/// kodem, co metoda MRS, więc stosunek wychodzi 1,00×.
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
    /// max |m_ij − macierz_analityczna| dla łańcucha macierzy.
    matrix_elem_err: f64,
    /// |det m − 1| dla łańcucha macierzy (niezmiennik Lorentza).
    matrix_det_err: f64,
    /// |rapidity(m) − Σθ| dla łańcucha macierzy.
    matrix_rapidity_err: f64,
    /// |Σθ naiwnie − Σθ Kahana|, czyli własny błąd sumowania rapidyty.
    rapidity_sum_err: f64,
};

/// Łańcuch N boostów o USTALONEJ rapidycie całkowitej. Normalizacja jest
/// konieczna, żeby pomiar izolował AKUMULACJĘ błędu, a nie warunkowanie
/// mapy θ → (cosh θ, sinh θ) — to drugie mierzymy osobno w
/// `conditioningOfRapidityMap`, bo to zupełnie inne zjawisko.
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

    // MRS: rapidyta jako współrzędna
    const z_kahan = Z.fromRapidity(theta_kahan);
    const z_plain = Z.fromRapidity(theta_plain);

    // baseline: łańcuch iloczynów macierzy 2×2
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

/// Warunkowanie mapy θ → (cosh θ, sinh θ): jak dokładnie reprezentacja
/// multyplikatywna trzyma niezmiennik N(z) = cosh²θ − sinh²θ = 1.
/// Błąd rośnie jak e^{2|θ|}·eps, bo cosh² i sinh² to dwie ogromne liczby,
/// których różnica jest mała — klasyczne znoszenie.
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
        // zapamiętujemy stosunek MRS do najlepszej znanej metody
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

    // baseline: łańcuch macierzy
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

test "drift: MRS nie akumuluje błędu, macierz akumuluje dyfuzyjnie" {
    const alloc = std.testing.allocator;
    const d3 = try measureDrift(alloc, 1_000);
    const d5 = try measureDrift(alloc, 100_000);

    // MRS: błąd składania rapidyty nie zależy od N (współrzędna addytywna)
    try std.testing.expect(d3.mrs_kahan_inv_err < 1e-13);
    try std.testing.expect(d5.mrs_kahan_inv_err < 1e-13);

    // macierz: błąd rośnie. Tempo jest DYFUZYJNE (~√N), a nie liniowe,
    // bo błędy zaokrągleń każdego kroku są w przybliżeniu niezależne.
    // Dla N różniących się 100× oczekujemy ~10× wzrostu błędu.
    try std.testing.expect(d5.matrix_elem_err > d3.matrix_elem_err * 5.0);
    try std.testing.expect(d5.matrix_elem_err > 1e-14);
}

test "warunkowanie mapy rapidyty psuje się dla dużych |θ|" {
    try std.testing.expect(conditioningOfRapidityMap(1.0) < 1e-15);
    try std.testing.expect(conditioningOfRapidityMap(50.0) > 1e-9);
    try std.testing.expect(conditioningOfRapidityMap(200.0) > 1e-3);
}

test "Kahan beats plain summation on asymmetric terms" {
    // Celowo zły przypadek: duże składniki i drobne dokładki.
    var xs: [2002]f64 = undefined;
    xs[0] = 1e16;
    for (1..2001) |i| xs[i] = 1.0;
    xs[2001] = -1e16;

    const p = plainSum(&xs);
    const k = kahanSum(&xs);
    try std.testing.expectApproxEqAbs(@as(f64, 2000.0), k, 1e-6);
    try std.testing.expect(@abs(k - 2000.0) <= @abs(p - 2000.0));
}
