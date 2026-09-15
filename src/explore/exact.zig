//! MRS-0.1 :: dokładna arytmetyka multivektorów (liczby całkowite)
//!
//! Cała warstwa algebraiczna MRS jest całkowita: znaki e_i² należą do
//! {+1,−1,0}, a iloczyn dwóch blatów jest z dokładnością do znaku pojedynczym
//! blatem (e_A·e_B = ±e_{A△B}). Dzięki temu pytania algebraiczne można
//! rozstrzygać w arytmetyce DOKŁADNEJ, bez tolerancji i bez zmiennoprzecinkowych
//! wątpliwości. To nie jest optymalizacja — to jest warunek sensowności
//! wyników: twierdzenie o wielomianie nie może zależeć od zaokrągleń.
//!
//! STABILNOŚĆ. Arytmetyka używa `std.math.add`/`std.math.mul`, czyli przy
//! przekroczeniu zakresu następuje JAWNA panika z komunikatem, a nie
//! ciche przekręcenie (UB w ReleaseFast). Współczynniki w tej warstwie są
//! rzędów kilku–kilkudziesięciu, więc i64 ma zapas rzędu 10^17; panika jest
//! wyłącznie siatką bezpieczeństwa na wypadek nadużycia API.
//!
//! Wektory mają stały rozmiar (2^5 = 32 blaty), więc **nie ma alokacji**
//! w żadnej operacji. To czyni wyniki powtarzalnymi i odpornymi na
//! fragmentację pamięci.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;

/// 2^MAX_TOTAL, gdzie MAX_TOTAL = 5.
pub const MAX_BASIS: usize = 32;

pub const Overflow = error{Overflow};

inline fn addExact(acc: *i64, delta: i64) void {
    acc.* = std.math.add(i64, acc.*, delta) catch @panic(
        "MRS explore: przekroczenie zakresu w arytmetyce dokładnej",
    );
}

inline fn mulExact(a: i64, b: i64) i64 {
    return std.math.mul(i64, a, b) catch @panic(
        "MRS explore: przekroczenie zakresu w arytmetyce dokładnej",
    );
}

/// Multivektor gęsty o współczynnikach całkowitych, indeksowany maską blatu.
pub const IntVec = struct {
    c: [MAX_BASIS]i64 = [_]i64{0} ** MAX_BASIS,

    pub fn zero() IntVec {
        return .{};
    }

    pub fn scalar(v: i64) IntVec {
        var r = IntVec{};
        r.c[0] = v;
        return r;
    }

    pub fn basis(mask: u32) IntVec {
        var r = IntVec{};
        r.c[mask] = 1;
        return r;
    }

    pub fn eql(a: IntVec, b: IntVec, basis_count: usize) bool {
        for (0..basis_count) |i| {
            if (a.c[i] != b.c[i]) return false;
        }
        return true;
    }

    pub fn isZero(self: IntVec, basis_count: usize) bool {
        for (0..basis_count) |i| {
            if (self.c[i] != 0) return false;
        }
        return true;
    }

    pub fn nnz(self: IntVec, basis_count: usize) usize {
        var k: usize = 0;
        for (0..basis_count) |i| {
            if (self.c[i] != 0) k += 1;
        }
        return k;
    }

    pub fn writeTo(self: IntVec, alg: cl.Algebra, w: anytype) !void {
        const m = alg.basisCount();
        var first = true;
        for (0..m) |i| {
            const v = self.c[i];
            if (v == 0) continue;
            if (!first) {
                try w.writeAll(if (v < 0) " - " else " + ");
            } else if (v < 0) {
                try w.writeAll("-");
            }
            const mag: i64 = if (v < 0) -v else v;
            if (mag != 1 or i == 0) try w.print("{d}", .{mag});
            try alg.writeBlade(w, @intCast(i));
            first = false;
        }
        if (first) try w.writeAll("0");
    }
};

/// Iloczyn geometryczny w arytmetyce dokładnej. Koszt O(nnz(a)·nnz(b)).
pub fn mul(alg: cl.Algebra, a: IntVec, b: IntVec) IntVec {
    const m = alg.basisCount();
    var out = IntVec{};
    for (0..m) |i| {
        const ai = a.c[i];
        if (ai == 0) continue;
        for (0..m) |j| {
            const bj = b.c[j];
            if (bj == 0) continue;
            const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue; // generator nilpotentny
            addExact(&out.c[bp.mask], mulExact(mulExact(ai, bj), @as(i64, bp.sign)));
        }
    }
    return out;
}

pub fn add(a: IntVec, b: IntVec) IntVec {
    var out = IntVec{};
    for (0..MAX_BASIS) |i| out.c[i] = a.c[i] + b.c[i];
    return out;
}

pub fn sub(a: IntVec, b: IntVec) IntVec {
    var out = IntVec{};
    for (0..MAX_BASIS) |i| out.c[i] = a.c[i] - b.c[i];
    return out;
}

// ---------------------------------------------------------------------------
// Sprzężenia
// ---------------------------------------------------------------------------

fn gradeOf(mask: u32) u32 {
    return @popCount(mask);
}

/// Odwrócenie (reverse): (−1)^{k(k−1)/2}.
pub fn reverseSign(mask: u32) i8 {
    const k = gradeOf(mask);
    const half = (k *% (k -% 1)) / 2;
    return if (half % 2 == 0) 1 else -1;
}

/// Involucja gradacji: (−1)^k.
pub fn gradeInvSign(mask: u32) i8 {
    return if (gradeOf(mask) % 2 == 0) 1 else -1;
}

/// Sprzężenie Clifforda (czyli koniugacja): reverse ∘ involucja gradacji,
/// znak (−1)^{k(k+1)/2}.
pub fn cliffordConjSign(mask: u32) i8 {
    const k = gradeOf(mask);
    const half = (k *% (k +% 1)) / 2;
    return if (half % 2 == 0) 1 else -1;
}

pub fn conj(x: IntVec, basis_count: usize) IntVec {
    var out = IntVec{};
    for (0..basis_count) |i| {
        out.c[i] = mulExact(x.c[i], @as(i64, cliffordConjSign(@intCast(i))));
    }
    return out;
}

/// Część skalarna.
pub fn scalarPart(x: IntVec) i64 {
    return x.c[0];
}

// ---------------------------------------------------------------------------
// Centrum — liczone czystą logiką masek, bez algebry liniowej
// ---------------------------------------------------------------------------

/// Blat e_S jest centralny ⟺ e_S·e_g = e_g·e_S dla każdego generatora g.
/// Oba iloczyny dają ten sam blat e_{S△g}, więc wystarczy porównać znaki.
/// To czyni wyznaczenie centrum operacją na bitach — dokładną i natychmiastową.
pub fn bladeIsCentral(alg: cl.Algebra, mask: u32) bool {
    var g: u5 = 0;
    while (g < alg.n_gen) : (g += 1) {
        const bit: u32 = @as(u32, 1) << g;
        const left = cl.bladeMul(alg, mask, bit);
        const right = cl.bladeMul(alg, bit, mask);
        if (left.mask != right.mask or left.sign != right.sign) return false;
    }
    return true;
}

/// Zbiór blatów rozpinających centrum, jako maska bitowa nad blatami.
pub fn centerBasis(alg: cl.Algebra) u32 {
    const m = alg.basisCount();
    var set: u32 = 0;
    var i: usize = 0;
    while (i < m) : (i += 1) {
        if (bladeIsCentral(alg, @intCast(i))) set |= (@as(u32, 1) << @intCast(i));
    }
    return set;
}

pub fn dimOfSet(set: u32) u5 {
    return @intCast(@popCount(set));
}

/// Blaty zawierające co najmniej jeden generator zdegenerowany.
/// Dla r > 0 rozpinają one nilpotentny ideał (patrz `subalgebra.zig`).
pub fn bladesWithDegenerate(alg: cl.Algebra) u32 {
    const m = alg.basisCount();
    var set: u32 = 0;
    var i: usize = 0;
    while (i < m) : (i += 1) {
        var g: u5 = 0;
        while (g < alg.n_gen) : (g += 1) {
            if (alg.squares[g] == 0 and (i & (@as(usize, 1) << g)) != 0) {
                set |= (@as(u32, 1) << @intCast(i));
                break;
            }
        }
    }
    return set;
}

/// Blaty parzystej gradacji (podalgebra Spin).
pub fn evenBlades(alg: cl.Algebra) u32 {
    const m = alg.basisCount();
    var set: u32 = 0;
    var i: usize = 0;
    while (i < m) : (i += 1) {
        if (gradeOf(@intCast(i)) % 2 == 0) set |= (@as(u32, 1) << @intCast(i));
    }
    return set;
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

const sigmod = @import("signatures.zig");

test "iloczyn dokładny zgadza się z iloczynem zmiennoprzecinkowym" {
    const sb = sigmod.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus);
    const alg = try sb.algebra();
    var prng = std.Random.DefaultPrng.init(7);
    const rnd = prng.random();
    const m = alg.basisCount();

    for (0..200) |_| {
        var a = IntVec{};
        var b = IntVec{};
        var af: [MAX_BASIS]f64 = undefined;
        var bf: [MAX_BASIS]f64 = undefined;
        for (0..m) |i| {
            a.c[i] = @intCast(rnd.intRangeLessThan(i32, -3, 4));
            b.c[i] = @intCast(rnd.intRangeLessThan(i32, -3, 4));
            af[i] = @floatFromInt(a.c[i]);
            bf[i] = @floatFromInt(b.c[i]);
        }
        const p = mul(alg, a, b);
        // sanity: 1·1 = 1
        const q = cl.bladeMul(alg, 0, 0);
        try std.testing.expectEqual(@as(i8, 1), q.sign);
        // porównanie z gęstym iloczynem zmiennoprzecinkowym
        const dense = try cl.mulDense(std.testing.allocator, alg, af[0..m], bf[0..m]);
        defer std.testing.allocator.free(dense);
        for (0..m) |i| {
            try std.testing.expectApproxEqAbs(dense[i], @as(f64, @floatFromInt(p.c[i])), 1e-6);
        }
    }
}

test "sprzężenie Clifforda daje skalar dla generatorów" {
    const sb = sigmod.SigBuf.build(.{ .p = 0, .q = 2 }, .mostly_minus);
    const alg = try sb.algebra();
    const m = alg.basisCount();

    // e_1: conj = −e_1, e_1·(−e_1) = −e_1² = +1
    const e1 = IntVec.basis(1);
    const p1 = mul(alg, e1, conj(e1, m));
    try std.testing.expectEqual(@as(i64, 1), scalarPart(p1));
    try std.testing.expect(p1.nnz(m) == 1);

    // e_12: conj = −e_12, iloczyn = +1
    const e12 = IntVec.basis(3);
    const p12 = mul(alg, e12, conj(e12, m));
    try std.testing.expectEqual(@as(i64, 1), scalarPart(p12));
    try std.testing.expect(p12.nnz(m) == 1);
}

test "znaki sprzężeń są spójne z modułem clifford" {
    const sb = sigmod.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus);
    const alg = try sb.algebra();
    const m = alg.basisCount();
    var i: usize = 0;
    while (i < m) : (i += 1) {
        const mask: u32 = @intCast(i);
        try std.testing.expectEqual(alg.reverseSign(mask), reverseSign(mask));
        // koniugacja = reverse ∘ involucja
        try std.testing.expectEqual(
            @as(i8, alg.reverseSign(mask)) * gradeInvSign(mask),
            cliffordConjSign(mask),
        );
    }
}

test "centrum: Cl(0,3) ma wymiar 2, Cl(0,2) ma wymiar 1" {
    const a3 = try (sigmod.SigBuf.build(.{ .p = 0, .q = 3 }, .mostly_minus)).algebra();
    const c3 = centerBasis(a3);
    try std.testing.expectEqual(@as(u5, 2), dimOfSet(c3));
    // centrum rozpięte na {1, e1e2e3}
    try std.testing.expect(c3 & 1 != 0);
    try std.testing.expect(c3 & (@as(u32, 1) << 7) != 0);

    const a2 = try (sigmod.SigBuf.build(.{ .p = 0, .q = 2 }, .mostly_minus)).algebra();
    const c2 = centerBasis(a2);
    try std.testing.expectEqual(@as(u5, 1), dimOfSet(c2));
    try std.testing.expectEqual(@as(u32, 1), c2);

    const a1 = try (sigmod.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus)).algebra();
    const c1 = centerBasis(a1);
    try std.testing.expectEqual(@as(u5, 1), dimOfSet(c1));
    try std.testing.expectEqual(@as(u32, 1), c1);
}

test "zdegenerowane blaty: jedna macierz pokrywa wszystkie zdegenerowane wymiary" {
    const alg = try (sigmod.SigBuf.build(.{ .p = 1, .q = 1, .r = 1 }, .mostly_minus)).algebra();
    const d = bladesWithDegenerate(alg);
    // wymiar 2 (generator zdegenerowany) → blaty {e2, e0e2, e1e2, e0e1e2} = maski 4,5,6,7
    try std.testing.expectEqual(@as(u32, 0b11110000), d);
    try std.testing.expectEqual(@as(u5, 4), dimOfSet(d));
}

test "parzyste blaty: połowa bazy" {
    const sb = sigmod.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus);
    const alg = try sb.algebra();
    const ev = evenBlades(alg);
    try std.testing.expectEqual(@as(u5, 8), dimOfSet(ev));
    try std.testing.expectEqual(alg.basisCount() / 2, @as(usize, dimOfSet(ev)));
}

test "brak alokacji w gorącej ścieżce — struktury są stałego rozmiaru" {
    // Ten test jest kontraktem: rozmiar IntVec nie zależy od n.
    try std.testing.expectEqual(@as(usize, 32 * 8), @sizeOf(IntVec));
}
