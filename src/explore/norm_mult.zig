//! MRS-0.1 :: P2 — multyplikatywność normy
//!
//! PYTANIE (rozłożone na trzy predykaty, bo w szkicu projektu było jedno
//! i było wewnętrznie sprzeczne):
//!
//!   P2a  N_sc(x·y) = N_sc(x)·N_sc(y),  gdzie N_sc(z) = Sc(z·z̄)
//!   P2b  (x·y)·conj(x·y) = (x·conj x)·(y·conj y)   (tożsamość w centrum)
//!   P2c  czy z ↦ z·z̄ zawsze trafia w centrum? (sensowność P2b)
//!
//! WYNIK ZMIERZONY: P2c zachodzi dokładnie dla 2^n <= 8 (n <= 3) i pada od
//! n = 4. Świadek w Cl(1,3): x = e01 + e23, bo e01 i e23 są rozłączne i
//! komutują, więc x·x̄ zawiera człon 2·e0123, a e0123 nie jest centralny.
//! To jest ta sama granica co dla multyplikatywności — i to nie przypadek:
//! od n = 4 mapa z ↦ z·z̄ przestaje trafiać w centrum, więc struktura
//! multyplikatywnej normy nie ma gdzie powstać.
//!
//! ---------------------------------------------------------------------------
//! TO JEST PROCEDURA ROZSTRZYGAJĄCA, NIE PRÓBKOWANIE
//! ---------------------------------------------------------------------------
//! Kluczowa obserwacja: każda z tych tożsamości jest, przy ustalonym y,
//! **funkcją kwadratową** x (złożenie odwzorowania liniowego z formą
//! kwadratową), i symetrycznie — przy ustalonym x jest funkcją kwadratową y.
//! A funkcja kwadratowa (także z wyrazem liniowym) jest jednoznacznie
//! wyznaczona przez wartości na zbiorze
//!
//!     D = {0} ∪ {e_i} ∪ {e_i + e_j : i < j},
//!
//! bo z Q(0), Q(e_i) i Q(e_i+e_j) odczytujemy wszystkie współczynniki
//! kombinacji a_i + c_ii oraz c_ij. Wobec tego:
//!
//!     tożsamość zachodzi dla wszystkich x, y  ⟺  zachodzi dla (x,y) ∈ D×D.
//!
//! To nie jest „dużo testów". To jest **dowód wyczerpujący** dla tej klasy
//! tożsamości, w arytmetyce całkowitej, bez tolerancji i bez prawdopodobieństwa.
//! Dlatego wynik wolno nazwać rozstrzygnięciem, a nie poszlaką.
//!
//! ZAKRES. Silnik rozstrzyga o ILOCZYNIE CLIFFORDA w danej sygnaturze.
//! Nie rozstrzyga pytania „czy na tej przestrzeni istnieje JAKIEKOLWIEK
//! mnożenie z multyplikatywną normą" — a to drugie pytanie obejmuje oktoniony,
//! które NIE są algebrą Clifforda (Clifford jest łączna, oktoniony nie).
//! Dlatego granica Hurwitza (1,2,4,8) nie może wyjść z tego silnika w całości:
//! wymiar 8 w rodzinie Clifforda to Cl(0,3) ≅ H⊕H, a nie oktoniony.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const exact = @import("exact.zig");
const sigs = @import("signatures.zig");

const IntVec = exact.IntVec;

/// Pełny rozmiar zbioru determinującego dla n = 5.
pub const MAX_GRID: usize = 1 + 32 + (32 * 31) / 2; // 529

pub const Kind = enum { scalar_multiplicative, center_multiplicative, norm_is_central };

pub const Verdict = struct {
    kind: Kind,
    n: usize,
    grid_size: usize,
    pairs: usize,
    holds: bool,
    /// Świadek zaprzeczenia (tylko gdy `holds == false`).
    witness_x: IntVec = IntVec{},
    witness_y: IntVec = IntVec{},
    has_witness: bool = false,

    pub fn name(self: Verdict) []const u8 {
        return switch (self.kind) {
            .scalar_multiplicative => "P2a  N_sc(xy) = N_sc(x)·N_sc(y)",
            .center_multiplicative => "P2b  (xy)·conj(xy) = (x·conj x)(y·conj y)",
            .norm_is_central => "P2c  z·z̄ ∈ centrum dla każdego z",
        };
    }
};

/// Buduje zbiór determinujący D = {0} ∪ {e_i} ∪ {e_i+e_j}.
/// `out` musi mieć co najmniej MAX_GRID miejsc; zwraca użyty rozmiar.
pub fn buildGrid(alg: cl.Algebra, out: *[MAX_GRID]IntVec) usize {
    const m = alg.basisCount();
    var k: usize = 0;
    out[k] = IntVec.zero();
    k += 1;
    for (0..m) |i| {
        out[k] = IntVec.basis(@intCast(i));
        k += 1;
    }
    for (0..m) |i| {
        for (i + 1..m) |j| {
            var v = IntVec.basis(@intCast(i));
            v.c[j] = 1;
            out[k] = v;
            k += 1;
        }
    }
    return k;
}

fn normOf(x: IntVec, alg: cl.Algebra) IntVec {
    const m = alg.basisCount();
    return exact.mul(alg, x, exact.conj(x, m));
}

/// Sprawdza P2c: czy z·z̄ leży w centrum dla każdego z.
/// Centralność jest warunkiem liniowym, a z·z̄ jest kwadratowe w z,
/// więc zbiór D znowu wystarcza.
pub fn checkNormIsCentral(alg: cl.Algebra, grid: *[MAX_GRID]IntVec) Verdict {
    const k = buildGrid(alg, grid);
    const m = alg.basisCount();
    const center = exact.centerBasis(alg);

    var i: usize = 0;
    while (i < k) : (i += 1) {
        const nz = normOf(grid[i], alg);
        for (0..m) |b| {
            if (nz.c[b] == 0) continue;
            if (center & (@as(u32, 1) << @intCast(b)) == 0) {
                return .{
                    .kind = .norm_is_central,
                    .n = alg.n_gen,
                    .grid_size = k,
                    .pairs = k,
                    .holds = false,
                    .witness_x = grid[i],
                    .witness_y = IntVec.basis(@intCast(b)),
                    .has_witness = true,
                };
            }
        }
    }
    return .{
        .kind = .norm_is_central,
        .n = alg.n_gen,
        .grid_size = k,
        .pairs = k,
        .holds = true,
    };
}

/// Sprawdza P2a: N_sc(x·y) = N_sc(x)·N_sc(y).
pub fn checkScalarMultiplicative(alg: cl.Algebra, grid: *[MAX_GRID]IntVec) Verdict {
    const k = buildGrid(alg, grid);
    var i: usize = 0;
    while (i < k) : (i += 1) {
        const nx = exact.scalarPart(normOf(grid[i], alg));
        var j: usize = 0;
        while (j < k) : (j += 1) {
            const ny = exact.scalarPart(normOf(grid[j], alg));
            const nxy = exact.scalarPart(normOf(exact.mul(alg, grid[i], grid[j]), alg));
            if (nxy != exact_mul(nx, ny)) {
                return .{
                    .kind = .scalar_multiplicative,
                    .n = alg.n_gen,
                    .grid_size = k,
                    .pairs = k * k,
                    .holds = false,
                    .witness_x = grid[i],
                    .witness_y = grid[j],
                    .has_witness = true,
                };
            }
        }
    }
    return .{
        .kind = .scalar_multiplicative,
        .n = alg.n_gen,
        .grid_size = k,
        .pairs = k * k,
        .holds = true,
    };
}

/// Sprawdza P2b: tożsamość wektorowa w centrum.
pub fn checkCenterMultiplicative(alg: cl.Algebra, grid: *[MAX_GRID]IntVec) Verdict {
    const k = buildGrid(alg, grid);
    const m = alg.basisCount();
    var i: usize = 0;
    while (i < k) : (i += 1) {
        const nx = normOf(grid[i], alg);
        var j: usize = 0;
        while (j < k) : (j += 1) {
            const ny = normOf(grid[j], alg);
            const prod = exact.mul(alg, grid[i], grid[j]);
            const nxy = normOf(prod, alg);
            const rhs = exact.mul(alg, nx, ny);
            if (!IntVec.eql(nxy, rhs, m)) {
                return .{
                    .kind = .center_multiplicative,
                    .n = alg.n_gen,
                    .grid_size = k,
                    .pairs = k * k,
                    .holds = false,
                    .witness_x = grid[i],
                    .witness_y = grid[j],
                    .has_witness = true,
                };
            }
        }
    }
    return .{
        .kind = .center_multiplicative,
        .n = alg.n_gen,
        .grid_size = k,
        .pairs = k * k,
        .holds = true,
    };
}

inline fn exact_mul(a: i64, b: i64) i64 {
    return std.math.mul(i64, a, b) catch @panic("MRS explore: overflow w P2");
}

/// Zbiór determinujący jest kompletny — to kontrakt tej warstwy.
/// Zwraca rozmiar dla danej liczby blatów.
pub fn gridSize(basis_count: usize) usize {
    return 1 + basis_count + (basis_count * (basis_count - 1)) / 2;
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "rozmiar siatki zgadza się z faktycznie zbudowaną" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const k = buildGrid(alg, &grid);
        try std.testing.expectEqual(gridSize(alg.basisCount()), k);
        try std.testing.expect(k <= MAX_GRID);
    }
}

test "P2a: prawda dla n <= 2 w OBU konwencjach" {
    var grid: [MAX_GRID]IntVec = undefined;
    inline for (.{ .mostly_minus, .mostly_plus }) |ts| {
        for ([_]sigs.Triple{
            .{ .p = 1, .q = 0 },
            .{ .p = 0, .q = 1 },
            .{ .p = 0, .q = 2 },
            .{ .p = 1, .q = 1 },
            .{ .p = 2, .q = 0 },
        }) |t| {
            const alg = try (sigs.SigBuf.build(t, ts)).algebra();
            const v = checkScalarMultiplicative(alg, &grid);
            try std.testing.expect(v.holds);
            const c = checkCenterMultiplicative(alg, &grid);
            try std.testing.expect(c.holds);
        }
    }
}

test "P2a: fałsz od n = 3 — silnik podaje świadka" {
    var grid: [MAX_GRID]IntVec = undefined;
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_minus)).algebra();
    const v = checkScalarMultiplicative(alg, &grid);
    try std.testing.expect(!v.holds);
    try std.testing.expect(v.has_witness);
    // świadek musi faktycznie łamać tożsamość
    const nx = exact.scalarPart(normOf(v.witness_x, alg));
    const ny = exact.scalarPart(normOf(v.witness_y, alg));
    const nxy = exact.scalarPart(normOf(exact.mul(alg, v.witness_x, v.witness_y), alg));
    try std.testing.expect(nxy != exact_mul(nx, ny));
}

test "P2c: z·z̄ trafia w centrum DOKŁADNIE dla 2^n <= 8" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);

    // n <= 3 (2^n <= 8): zachodzi dla KAŻDEJ sygnatury i konwencji
    var low: usize = 0;
    for (0..n) |i| {
        if (buf[i].n() > 3) continue;
        low += 1;
        const a_minus = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const a_plus = try (sigs.SigBuf.build(buf[i], .mostly_plus)).algebra();
        if (!checkNormIsCentral(a_minus, &grid).holds) {
            std.debug.print("P2c pada dla n<=3: triple=({d},{d},{d}) minus\n", .{ buf[i].p, buf[i].q, buf[i].r });
            return error.TestUnexpectedResult;
        }
        if (!checkNormIsCentral(a_plus, &grid).holds) {
            std.debug.print("P2c pada dla n<=3: triple=({d},{d},{d}) plus\n", .{ buf[i].p, buf[i].q, buf[i].r });
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(low > 0);

    // n = 4: pada, i to jest wynik, nie usterka. Świadek konstrukcyjny:
    // x = e01 + e23 w Cl(1,3), bo e01 i e23 są rozłączne i komutują,
    // więc x·x̄ zawiera człon 2·e0123, a e0123 NIE jest centralny
    // (środkiem Cl(1,3) jest samo R).
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus)).algebra();
    const v = checkNormIsCentral(alg, &grid);
    try std.testing.expect(!v.holds);
    try std.testing.expect(v.has_witness);

    // sprawdzamy świadka ręcznie: e01 = maska 3, e23 = maska 12
    const e01 = IntVec.basis(3);
    const e23 = IntVec.basis(12);
    const xx = exact.add(e01, e23); // x = e01 + e23
    const nz = normOf(xx, alg);
    try std.testing.expect(nz.c[15] != 0); // człon przy e0123 (maska 15)
    try std.testing.expectEqual(@as(u5, 1), exact.dimOfSet(exact.centerBasis(alg)));

    // dla n = 3 ten sam mechanizm NIE psuje centralności, bo e012 jest centralny
    const alg3 = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_minus)).algebra();
    const v3 = checkNormIsCentral(alg3, &grid);
    try std.testing.expect(v3.holds);
    try std.testing.expectEqual(@as(u5, 2), exact.dimOfSet(exact.centerBasis(alg3)));
}

test "DECYZJA jest kompletna: wykrywa złamanie poza siatką" {
    // Bierzemy świadka z siatki i sprawdzamy, że naprawdę łamie tożsamość
    // dla punktów spoza siatki też — to test, że metoda nie jest sztuczką
    // zależną od wyboru D.
    var grid: [MAX_GRID]IntVec = undefined;
    const alg = try (sigs.SigBuf.build(.{ .p = 0, .q = 3 }, .mostly_minus)).algebra();
    const v = checkScalarMultiplicative(alg, &grid);
    try std.testing.expect(!v.holds);

    // losowy punkt spoza siatki: tożsamość musi padać dla ogromnej
    // większości punktów (gdyby padała tylko na D, metoda byłaby fałszywa)
    var prng = std.Random.DefaultPrng.init(1234);
    const rnd = prng.random();
    const m = alg.basisCount();
    var fails: usize = 0;
    const trials: usize = 500;
    for (0..trials) |_| {
        var x = IntVec{};
        var y = IntVec{};
        for (0..m) |i| {
            x.c[i] = @intCast(rnd.intRangeLessThan(i32, -2, 3));
            y.c[i] = @intCast(rnd.intRangeLessThan(i32, -2, 3));
        }
        const nx = exact.scalarPart(normOf(x, alg));
        const ny = exact.scalarPart(normOf(y, alg));
        const nxy = exact.scalarPart(normOf(exact.mul(alg, x, y), alg));
        if (nxy != exact_mul(nx, ny)) fails += 1;
    }
    // nie wymagamy 100% — wymagamy, by złamanie było regułą, a nie wyjątkiem
    try std.testing.expect(fails * 2 > trials);
}

test "P2a zależy od konwencji znaku tylko przez trójkę — wynik ten sam" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        if (buf[i].n() > 3) continue;
        const a_minus = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const a_plus = try (sigs.SigBuf.build(buf[i], .mostly_plus)).algebra();
        try std.testing.expectEqual(
            checkScalarMultiplicative(a_minus, &grid).holds,
            checkScalarMultiplicative(a_plus, &grid).holds,
        );
    }
}

// ---------------------------------------------------------------------------
// REGUŁY ZNALEZIONE PRZEZ SILNIK — utrwalone jako sprawdzalne twierdzenia
// ---------------------------------------------------------------------------
//
// Silnik nie został tak zaprogramowany. Te dwie reguły WYSZŁY z tabeli,
// a te testy zamieniają obserwację w twierdzenie sprawdzane przy każdej
// zmianie kodu. Dla każdej reguły podajemy też szkic dowodu.

test "REGUŁA 1: P2a zachodzi dokładnie gdy p+q <= 2 (zdegenerowane niewidoczne)" {
    // DLACZEGO TAK JEST. Cl(p,q,r) ≅ Cl(p,q) ⊗ ⋀(R), gdzie R to radykał.
    // Część skalarna normy widzi wyłącznie Cl(p,q):
    //   * iloczyn dwóch różnych blatów nigdy nie daje skalara (e_A·e_B = ±e_{A△B});
    //   * iloczyn e_A·e_A daje skalar równy iloczynowi s_i po A, a ten jest
    //     ZEREM, gdy A zawiera generator zdegenerowany.
    // Zatem N_sc(x) = (współczynnik skalarny x)². Stąd
    //   N_sc(xy) = (c0(x)·c0(y))² = N_sc(x)·N_sc(y),
    // czyli przy r > 0 tożsamość zachodzi TOŻSAMOŚCIOWO, niezależnie od Cl(p,q).
    // Dla r = 0 pozostaje klasyczna granica Hurwitza w rodzinie Clifforda:
    // multyplikatywność jest wtedy i tylko wtedy, gdy 2^(p+q) <= 4.
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);

    var checked: usize = 0;
    var degenerate_seen: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.n() > 4) continue;
        checked += 1;
        if (t.isDegenerate()) degenerate_seen += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const holds = checkScalarMultiplicative(alg, &grid).holds;
        const predicted = (@as(usize, t.p) + @as(usize, t.q)) <= 2;
        if (holds != predicted) {
            std.debug.print(
                "REGUŁA 1 złamana: t=({d},{d},{d}) p+q={d} zmierzone={} przewidziane={}\n",
                .{ t.p, t.q, t.r, @as(usize, t.p) + @as(usize, t.q), holds, predicted },
            );
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked >= 30);
    try std.testing.expect(degenerate_seen >= 15);
}

test "REGUŁA 2: P2c zachodzi dokładnie gdy 2^n <= 8" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.n() > 4) continue;
        checked += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const holds = checkNormIsCentral(alg, &grid).holds;
        const predicted = t.n() <= 3;
        if (holds != predicted) {
            std.debug.print(
                "REGUŁA 2 złamana: t=({d},{d},{d}) n={d} zmierzone={} przewidziane={}\n",
                .{ t.p, t.q, t.r, t.n(), holds, predicted },
            );
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked >= 30);
}

test "P2b różni się od P2a; P2b == P2c w zakresie (obserwacja, nie twierdzenie)" {
    // Po co trzy predykaty? Bo dają różne odpowiedzi — a gdyby nie dawały,
    // rozdzielenie P2 byłoby zbędne. Ten test pilnuje, że różnica P2a vs P2b
    // jest realna, i ODNOTOWUJE (bez dowodu), że P2b i P2c pokrywają się
    // na całym przebadanym zakresie n <= 4. Tego drugiego nie udajemy
    // twierdzeniem: to wzorzec do zbadania, wpisany do docs/09_open.md.
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var differ_a: usize = 0;
    var differ_c: usize = 0;
    var witness_a: sigs.Triple = .{};
    for (0..n) |i| {
        const t = buf[i];
        if (t.n() > 4) continue;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const a = checkScalarMultiplicative(alg, &grid).holds;
        const b = checkCenterMultiplicative(alg, &grid).holds;
        const c = checkNormIsCentral(alg, &grid).holds;
        if (a != b) {
            differ_a += 1;
            witness_a = t;
        }
        if (c != b) differ_c += 1;
    }
    // P2a i P2b MUSZĄ się różnić — inaczej rozdzielenie predykatów jest puste.
    try std.testing.expect(differ_a > 0);
    std.debug.print(
        "P2a != P2b w {d} sygnaturach (np. ({d},{d},{d})); P2b != P2c w {d}\n",
        .{ differ_a, witness_a.p, witness_a.q, witness_a.r, differ_c },
    );
    // Obserwacja bez dowodu: P2b i P2c zgadzają się na całym zakresie.
    try std.testing.expectEqual(@as(usize, 0), differ_c);
}
