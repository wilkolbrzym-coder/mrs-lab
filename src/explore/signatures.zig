//! MRS-0.1 :: enumerator sygnatur (warstwa bazowa pytań P1–P4)
//!
//! Pytania P2 i P4 zależą wyłącznie od trójki (p,q,r), a nie od kolejności ról
//! w wektorze. Dlatego enumerujemy TRÓJKI, a sygnatury budujemy z nich
//! w ustalonej kolejności kanonicznej: najpierw wymiary czasowe, potem
//! przestrzenne, potem zdegenerowane.
//!
//! Kolejność jest częścią kontraktu: raport musi być identyczny między
//! uruchomieniami, inaczej nie da się go porównywać ani testować.
//!
//! Bufor `SigBuf` trzyma role w tablicy o stałym rozmiarze — dzięki temu
//! sygnatura nie wymaga alokacji i nie może wyciekać ani się przesunąć.

const std = @import("std");
const mrs = @import("mrs");
const sig = mrs.signature;

/// Największy wymiar obsługiwany przez wyrocznię wyczerpującą.
/// 2^5 = 32 blaty; dla n = 5 podzbiory blatów są już nieprzeliczalne
/// (2^32), więc enumeracja podalgebr jest ograniczona do n ≤ 4 (patrz
/// `subalgebra.zig`). Pytania o normę i centrum działają do n = 5.
pub const MAX_TOTAL: u5 = 5;
pub const MAX_DIM: usize = 8;

pub const Triple = struct {
    p: u5 = 0,
    q: u5 = 0,
    r: u5 = 0,

    pub fn n(self: Triple) usize {
        return @as(usize, self.p) + @as(usize, self.q) + @as(usize, self.r);
    }

    /// Sygnatura lorentzowska: dokładnie jeden wymiar czasowy, bez jądra.
    /// To jest dokładnie warunek na porządek częściowy (Twierdzenie 5.2).
    pub fn isLorentzian(self: Triple) bool {
        return self.p == 1 and self.r == 0;
    }

    pub fn isDegenerate(self: Triple) bool {
        return self.r > 0;
    }

    pub fn eql(a: Triple, b: Triple) bool {
        return a.p == b.p and a.q == b.q and a.r == b.r;
    }

    pub fn writeTo(self: Triple, w: anytype) !void {
        try w.print("({d},{d},{d})", .{ self.p, self.q, self.r });
    }
};

/// Sygnatura bez alokacji: role w tablicy stałego rozmiaru.
pub const SigBuf = struct {
    roles: [MAX_DIM]sig.Role = undefined,
    len: usize = 0,
    time_sign: sig.TimeSign = .mostly_minus,

    pub fn build(t: Triple, ts: sig.TimeSign) SigBuf {
        var b = SigBuf{ .time_sign = ts };
        var i: usize = 0;
        var k: u5 = 0;
        while (k < t.p) : (k += 1) {
            b.roles[i] = .temporal;
            i += 1;
        }
        k = 0;
        while (k < t.q) : (k += 1) {
            b.roles[i] = .spatial;
            i += 1;
        }
        k = 0;
        while (k < t.r) : (k += 1) {
            b.roles[i] = .degenerate;
            i += 1;
        }
        b.len = i;
        return b;
    }

    pub fn signature(self: *const SigBuf) sig.Signature {
        return .{ .roles = self.roles[0..self.len], .time_sign = self.time_sign };
    }

    pub fn algebra(self: *const SigBuf) mrs.clifford.AlgebraError!mrs.clifford.Algebra {
        return mrs.clifford.Algebra.fromSignature(self.signature());
    }
};

/// Liczba trójek o sumie od 1 do `max_total`.
pub fn tripleCount(max_total: u5) usize {
    var total: usize = 0;
    var n: u5 = 1;
    while (n <= max_total) : (n += 1) {
        var p: u5 = 0;
        while (p <= n) : (p += 1) {
            total += @as(usize, n - p) + 1;
        }
    }
    return total;
}

/// Wypełnia `out` wszystkimi trójkami o sumie 1..max_total.
/// Porządek jest deterministyczny: rosnące n, potem p, potem q, potem r.
/// Zwraca liczbę zapisanych elementów; jeśli `out` jest za mały, przerywa.
pub fn enumerateTriples(out: []Triple, max_total: u5) usize {
    var k: usize = 0;
    var n: u5 = 1;
    while (n <= max_total) : (n += 1) {
        var p: u5 = 0;
        while (p <= n) : (p += 1) {
            var q: u5 = 0;
            while (q <= n - p) : (q += 1) {
                const r: u5 = n - p - q;
                if (k >= out.len) return k;
                out[k] = .{ .p = p, .q = q, .r = r };
                k += 1;
            }
        }
    }
    return k;
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "liczba trójek zgadza się z symbolem Newtona" {
    // trójki o sumie <= max, bez (0,0,0): C(max+3,3) - 1
    try std.testing.expectEqual(@as(usize, 3), tripleCount(1)); // (1,0,0),(0,1,0),(0,0,1)
    try std.testing.expectEqual(@as(usize, 9), tripleCount(2)); // C(5,3)-1 = 9
    try std.testing.expectEqual(@as(usize, 19), tripleCount(3)); // C(6,3)-1 = 19
    try std.testing.expectEqual(@as(usize, 34), tripleCount(4)); // C(7,3)-1 = 34
    try std.testing.expectEqual(@as(usize, 55), tripleCount(5)); // C(8,3)-1 = 55
}

test "enumeracja jest deterministyczna i bez powtórzeń" {
    var buf_a: [64]Triple = undefined;
    var buf_b: [64]Triple = undefined;
    const na = enumerateTriples(&buf_a, MAX_TOTAL);
    const nb = enumerateTriples(&buf_b, MAX_TOTAL);
    try std.testing.expectEqual(na, nb);
    try std.testing.expectEqual(tripleCount(MAX_TOTAL), na);

    for (0..na) |i| {
        try std.testing.expect(buf_a[i].eql(buf_b[i]));
        // brak powtórzeń
        for (i + 1..na) |j| {
            try std.testing.expect(!buf_a[i].eql(buf_a[j]));
        }
        // zakres
        try std.testing.expect(buf_a[i].n() >= 1);
        try std.testing.expect(buf_a[i].n() <= MAX_TOTAL);
    }
}

test "SigBuf buduje poprawne sygnatury i algebry" {
    const alloc = std.testing.allocator;
    _ = alloc;
    var buf: [64]Triple = undefined;
    const n = enumerateTriples(&buf, MAX_TOTAL);
    for (0..n) |i| {
        const t = buf[i];
        const sb = SigBuf.build(t, .mostly_minus);
        const s = sb.signature();
        try std.testing.expectEqual(t.n(), s.n());
        try std.testing.expectEqual(@as(usize, t.p), s.p());
        try std.testing.expectEqual(@as(usize, t.q), s.q());
        try std.testing.expectEqual(@as(usize, t.r), s.r());
        const alg = try sb.algebra();
        try std.testing.expectEqual(@as(usize, @as(usize, 1) << @intCast(t.n())), alg.basisCount());
    }
}

test "konwencja nie zmienia trójki ani rozkładu ról" {
    var buf: [64]Triple = undefined;
    const n = enumerateTriples(&buf, MAX_TOTAL);
    for (0..n) |i| {
        const sb_minus = SigBuf.build(buf[i], .mostly_minus);
        const sb_plus = SigBuf.build(buf[i], .mostly_plus);
        try std.testing.expectEqual(sb_minus.len, sb_plus.len);
        for (0..sb_minus.len) |j| {
            try std.testing.expectEqual(sb_minus.roles[j], sb_plus.roles[j]);
        }
        // znaki formy są przeciwne, rola ta sama
        try std.testing.expectApproxEqAbs(
            sb_minus.signature().signAt(0),
            -sb_plus.signature().signAt(0),
            0,
        );
    }
}

test "algebra zdegenerowana odrzuca podniesienie wskaźnika" {
    const sb = SigBuf.build(.{ .p = 2, .q = 1, .r = 1 }, .mostly_minus);
    const f = mrs.form.DiagonalForm{ .signature = sb.signature() };
    const c = [_]f64{ 1, 1, 1, 1 };
    var out: [4]f64 = undefined;
    try std.testing.expectError(error.DegenerateForm, f.raiseIndex(&c, &out));
}
