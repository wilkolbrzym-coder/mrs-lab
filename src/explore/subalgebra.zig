//! MRS-0.1 :: P4 — podalgebry, ideały, centrum
//!
//! PYTANIE. Które podprzestrzenie rozpięte na blatach są zamknięte na iloczyn,
//! które są unitarne, przemienne, które są dwustronnymi ideałami całej algebry?
//!
//! DLACZEGO TO JEST DOKŁADNE. Blaty są liniowo niezależne, a iloczyn dwóch
//! blatów jest — z dokładnością do znaku — pojedynczym blatem
//! (e_A·e_B = ±e_{A△B}) albo zerem. Zatem podprzestrzeń rozpięta na zbiorze
//! blatów B jest zamknięta na iloczyn **wtedy i tylko wtedy**, gdy dla każdej
//! pary blatów z B ich iloczyn (lub zero) leży w rozpięciu B. Sprawdzenie
//! sprowadza się więc do testu przynależności maski — bez arytmetyki
//! zmiennoprzecinkowej i bez tolerancji.
//!
//! WYNIKI ZMIERZONE (weryfikowane niżej, wyczerpująco dla n ≤ 4):
//!   * blaty zawierające generator zdegenerowany rozpinają właściwy,
//!     niezerowy, dwustronny ideał I (ideał radykału);
//!   * I jest nilpotentny o indeksie DOKŁADNIE r + 1. W szczególności
//!     I² = 0 tylko dla r = 1 — dla r ≥ 2 różne generatory zdegenerowane
//!     mnożą się niezerowo, np. ζ₁·ζ₂ ≠ 0.
//!   * dla r = 0 nie ma właściwych niezerowych ideałów rozpiętych na blatach.
//!
//! ZAKRES — czytaj uważnie, bo to granica, nie drobiazg. Enumerujemy wyłącznie
//! podprzestrzenie ROZPIĘTE NA BLATACH. To jest podkratа wszystkich podalgebr,
//! a nie wszystkie podalgebry. Skutek: w rozszczepionych algebrach (np.
//! Cl(1,0) ≅ R⊕R) istnieją właściwe ideały rozpięte na idempotentach
//! (1 ± ω)/2, których ten silnik **nie zobaczy**, bo idempotenty nie są
//! blatami. Silnik tego nie ukrywa — raportuje wynik razem z zakresem.
//!
//! Wyczerpująco dla n ≤ 4 (2^16 = 65 536 podzbiorów). Dla n = 5 podzbiorów
//! jest 2^32, więc enumeracja odmawia — jawnie, a nie po cichu.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const exact = @import("exact.zig");
const sigs = @import("signatures.zig");

/// Największa liczba blatów, przy której enumeracja podzbiorów jest wykonalna.
pub const MAX_EXHAUSTIVE_BASIS: usize = 16;

pub const EnumerateError = error{TooManyBlades};

pub const Info = struct {
    /// Podzbiór blatów, bit i = blat o masce i.
    set: u32 = 0,
    dim: u5 = 0,
    unital: bool = false,
    commutative: bool = false,
    /// Dwustronny ideał CAŁEJ algebry.
    ideal: bool = false,
    proper: bool = false,
    /// I·I = 0 (iloczyn znika w całości) — własność ideału nilpotentnego.
    square_zero: bool = false,

    pub fn writeTo(self: Info, alg: cl.Algebra, w: anytype) !void {
        try w.print("dim={d} {{", .{self.dim});
        var first = true;
        const m = alg.basisCount();
        for (0..m) |i| {
            if (self.set & (@as(u32, 1) << @intCast(i)) == 0) continue;
            if (!first) try w.writeAll(", ");
            try alg.writeBlade(w, @intCast(i));
            first = false;
        }
        if (first) try w.writeAll("0");
        try w.writeAll("}");
    }
};

/// Czy rozpięcie zbioru blatów jest zamknięte na iloczyn.
pub fn isClosed(alg: cl.Algebra, set: u32) bool {
    const m = alg.basisCount();
    for (0..m) |i| {
        if (set & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (set & (@as(u32, 1) << @intCast(j)) == 0) continue;
            const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue; // zero należy do każdej podprzestrzeni
            if (set & (@as(u32, 1) << @intCast(bp.mask)) == 0) return false;
        }
    }
    return true;
}

/// Czy podalgebra jest przemienna (blaty parami się komutują).
pub fn isCommutative(alg: cl.Algebra, set: u32) bool {
    const m = alg.basisCount();
    for (0..m) |i| {
        if (set & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (set & (@as(u32, 1) << @intCast(j)) == 0) continue;
            const ab = cl.bladeMul(alg, @intCast(i), @intCast(j));
            const ba = cl.bladeMul(alg, @intCast(j), @intCast(i));
            if (ab.mask != ba.mask or ab.sign != ba.sign) return false;
        }
    }
    return true;
}

/// Czy rozpięcie zbioru jest dwustronnym ideałem całej algebry.
pub fn isTwoSidedIdeal(alg: cl.Algebra, set: u32) bool {
    if (set == 0) return true; // ideał zerowy
    const m = alg.basisCount();
    for (0..m) |b| {
        if (set & (@as(u32, 1) << @intCast(b)) == 0) continue;
        for (0..m) |a| { // po wszystkich blatach całej algebry
            const ab = cl.bladeMul(alg, @intCast(a), @intCast(b));
            if (ab.sign != 0 and set & (@as(u32, 1) << @intCast(ab.mask)) == 0) return false;
            const ba = cl.bladeMul(alg, @intCast(b), @intCast(a));
            if (ba.sign != 0 and set & (@as(u32, 1) << @intCast(ba.mask)) == 0) return false;
        }
    }
    return true;
}

/// Czy I·I = 0 (na blatach — wystarcza, bo iloczyn jest dwuliniowy).
pub fn isSquareZero(alg: cl.Algebra, set: u32) bool {
    const m = alg.basisCount();
    for (0..m) |i| {
        if (set & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (set & (@as(u32, 1) << @intCast(j)) == 0) continue;
            if (cl.bladeMul(alg, @intCast(i), @intCast(j)).sign != 0) return false;
        }
    }
    return true;
}

/// Zbiera informacje o zbiorze blatów.
pub fn info(alg: cl.Algebra, set: u32) Info {
    const m = alg.basisCount();
    return .{
        .set = set,
        .dim = @intCast(@popCount(set)),
        .unital = (set & 1) != 0,
        .commutative = if (isClosed(alg, set)) isCommutative(alg, set) else false,
        .ideal = isTwoSidedIdeal(alg, set),
        .proper = @popCount(set) < m,
        .square_zero = isSquareZero(alg, set),
    };
}

/// Wyczerpująca enumeracja zamkniętych podzbiorów blatów.
/// Zapisuje do `out` i zwraca liczbę znalezionych. `out` musi być
/// wystarczająco duże; reszta jest cicho pomijana, więc do diagnozy
/// używaj `countClosed`.
pub fn enumerateClosed(alloc: std.mem.Allocator, alg: cl.Algebra) ![]Info {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;

    var list = std.ArrayList(Info).empty;
    errdefer list.deinit(alloc);

    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 0;
    while (s < total) : (s += 1) {
        const set: u32 = @intCast(s);
        if (!isClosed(alg, set)) continue;
        try list.append(alloc, info(alg, set));
    }
    return list.toOwnedSlice(alloc);
}

/// Liczba zamkniętych podzbiorów (bez alokacji wyniku).
pub fn countClosed(alg: cl.Algebra) EnumerateError!usize {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;
    var count: usize = 0;
    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 0;
    while (s < total) : (s += 1) {
        if (isClosed(alg, @intCast(s))) count += 1;
    }
    return count;
}

/// Liczba właściwych, niezerowych, dwustronnych ideałów rozpiętych na blatach.
pub fn countProperIdeals(alg: cl.Algebra) EnumerateError!usize {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;
    var count: usize = 0;
    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 1; // pomijamy zbiór pusty (ideał zerowy)
    while (s < total) : (s += 1) {
        const set: u32 = @intCast(s);
        const dim: u5 = @intCast(@popCount(set));
        if (dim == 0 or dim == m) continue; // trywialne
        if (isTwoSidedIdeal(alg, set)) count += 1;
    }
    return count;
}

// ---------------------------------------------------------------------------
// Potęgi ideału i nilpotentność
// ---------------------------------------------------------------------------

/// Iloczyn dwóch podprzestrzeni rozpiętych na blatach.
///
/// FAKT STRUKTURALNY, na którym to stoi: iloczyn dwóch blatów jest
/// z dokładnością do znaku pojedynczym blatem. Dlatego iloczyn podprzestrzeni
/// rozpiętych na blatach JEST znowu rozpięty na blatach i nie potrzeba do tego
/// żadnej algebry liniowej nad ciałem — wystarczy zbiór masek.
pub fn mulSets(alg: cl.Algebra, a: u32, b: u32) u32 {
    const m = alg.basisCount();
    var out: u32 = 0;
    for (0..m) |i| {
        if (a & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (b & (@as(u32, 1) << @intCast(j)) == 0) continue;
            const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue;
            out |= (@as(u32, 1) << @intCast(bp.mask));
        }
    }
    return out;
}

/// Indeks nilpotentności: najmniejsze k ≥ 1 takie, że I^k = {0}.
/// Zwraca null, gdy I^k ≠ 0 przez `max_steps` kroków.
pub fn nilpotencyIndex(alg: cl.Algebra, set: u32, max_steps: u32) ?u32 {
    if (set == 0) return 1;
    var cur = set;
    var k: u32 = 1;
    while (k <= max_steps) : (k += 1) {
        cur = mulSets(alg, cur, set);
        if (cur == 0) return k + 1;
    }
    return null;
}

// ---------------------------------------------------------------------------
// Testy — sprawdzają PRZEWIDYWANIA, nie tylko działanie kodu
// ---------------------------------------------------------------------------

test "blaty parzyste tworzą unitarną podalgebrę wymiaru 2^(n-1)" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        if (buf[i].n() > 4) continue;
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const ev = exact.evenBlades(alg);
        try std.testing.expect(isClosed(alg, ev));
        try std.testing.expect((ev & 1) != 0); // unitarna: zawiera skalar 1
        try std.testing.expectEqual(@as(u5, @intCast(alg.basisCount() / 2)), @as(u5, @intCast(@popCount(ev))));
    }
}

test "centrum zawsze jest podalgebrą" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        if (buf[i].n() > 4) continue;
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const ctr = exact.centerBasis(alg);
        try std.testing.expect(isClosed(alg, ctr));
        try std.testing.expect((ctr & 1) != 0);
    }
}

test "PRZEWIDYWANIE: ideał radykału jest nilpotentny o indeksie DOKŁADNIE r+1" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (!t.isDegenerate() or t.n() > 4) continue;
        checked += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const I = exact.bladesWithDegenerate(alg);

        const closed = isClosed(alg, I);
        const ideal = isTwoSidedIdeal(alg, I);
        const proper = @as(usize, @popCount(I)) < alg.basisCount();
        const idx = nilpotencyIndex(alg, I, 8);

        if (!closed or !ideal or !proper or idx == null or idx.? != @as(u32, t.r) + 1) {
            std.debug.print(
                "ideał radykału: t=({d},{d},{d}) dim={d}/{d} closed={} ideal={} proper={} idx={?}\n",
                .{ t.p, t.q, t.r, @popCount(I), alg.basisCount(), closed, ideal, proper, idx },
            );
            return error.TestUnexpectedResult;
        }

        // I² = 0 zachodzi DOKŁADNIE dla r = 1 — dla r >= 2 różne generatory
        // zdegenerowane mnożą się niezerowo (np. ζ1·ζ2 ≠ 0).
        try std.testing.expectEqual(t.r == 1, isSquareZero(alg, I));

        const cnt = try countProperIdeals(alg);
        if (cnt < 1) {
            std.debug.print("brak właściwego ideału: t=({d},{d},{d})\n", .{ t.p, t.q, t.r });
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked > 0);
}

test "PRZEWIDYWANIE: r = 0 nie daje właściwych ideałów z blatów" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.isDegenerate() or t.n() > 4) continue;
        checked += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        // Uwaga: to jest twierdzenie O ZAKRESIE blatowym. Rozszczepione
        // przypadki (np. Cl(1,0) ≅ R⊕R) mają ideały rozpięte na idempotentach,
        // których ten silnik nie widzi — i to jest świadome ograniczenie.
        try std.testing.expectEqual(@as(usize, 0), try countProperIdeals(alg));
    }
    try std.testing.expect(checked > 0);
}

test "zbiory trywialne są zawsze zamknięte: {0} i cała algebra" {
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus)).algebra();
    try std.testing.expect(isClosed(alg, 0));
    const full: u32 = @intCast((@as(u64, 1) << @intCast(alg.basisCount())) - 1);
    try std.testing.expect(isClosed(alg, full));
    try std.testing.expect(isTwoSidedIdeal(alg, full));
    // pełna algebra nie jest właściwa
    try std.testing.expect(!info(alg, full).proper);
}

test "enumeracja odmawia dla n = 5 zamiast po cichu liczyć pół godziny" {
    const alg = try (sigs.SigBuf.build(.{ .p = 0, .q = 5 }, .mostly_minus)).algebra();
    try std.testing.expectEqual(@as(usize, 32), alg.basisCount());
    try std.testing.expectError(error.TooManyBlades, countClosed(alg));
    try std.testing.expectError(error.TooManyBlades, countProperIdeals(alg));
}

test "enumeracja jest deterministyczna i pokrywa się z licznikiem" {
    const alloc = std.testing.allocator;
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_minus)).algebra();
    const a = try enumerateClosed(alloc, alg);
    defer alloc.free(a);
    const b = try enumerateClosed(alloc, alg);
    defer alloc.free(b);
    try std.testing.expectEqual(a.len, b.len);
    try std.testing.expectEqual(try countClosed(alg), a.len);
    for (a, b) |x, y| {
        try std.testing.expectEqual(x.set, y.set);
        try std.testing.expectEqual(x.dim, y.dim);
        try std.testing.expectEqual(x.ideal, y.ideal);
    }
    // wyniki są niezależne od konwencji znaku (własność algebraiczna)
    const alg_plus = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_plus)).algebra();
    try std.testing.expectEqual(try countClosed(alg), try countClosed(alg_plus));
}
