//! MRS-LAB :: signature enumerator (base layer for questions P2 and P4)
//!
//! Questions P2 and P4 depend only on the triple (p,q,r), not on the order of
//! the roles in the vector. So we enumerate TRIPLES and build signatures from
//! them in a fixed canonical order: temporal dimensions first, then spatial,
//! then degenerate.
//!
//! The order is part of the contract: the report must be identical between
//! runs, otherwise it cannot be compared or tested.
//!
//! `SigBuf` keeps the roles in a fixed-size array, so a signature needs no
//! allocation and cannot leak or be moved.

const std = @import("std");
const mrs = @import("mrs");
const sig = mrs.signature;

/// Largest dimension the exhaustive oracle supports.
/// 2^5 = 32 blades; for n = 5 the blade subsets are already uncountable
/// (2^32), so subalgebra enumeration is limited to n <= 4 (see
/// `subalgebra.zig`). The norm and centre questions work up to n = 5.
pub const MAX_TOTAL: u5 = 5;
pub const MAX_DIM: usize = 8;

pub const Triple = struct {
    p: u5 = 0,
    q: u5 = 0,
    r: u5 = 0,

    pub fn n(self: Triple) usize {
        return @as(usize, self.p) + @as(usize, self.q) + @as(usize, self.r);
    }

    /// Lorentzian signature: exactly one temporal dimension, no radical.
    /// This is precisely the condition for a partial order (Theorem 5.2).
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

/// Allocation-free signature: roles in a fixed-size array.
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

/// Number of triples with sum from 1 to `max_total`.
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

/// Fills `out` with all triples of sum 1..max_total.
/// The order is deterministic: increasing n, then p, then q, then r.
/// Returns the number written; if `out` is too small it stops early.
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

test "the number of triples matches the binomial coefficient" {
    // triples of sum <= max, excluding (0,0,0): C(max+3,3) - 1
    try std.testing.expectEqual(@as(usize, 3), tripleCount(1)); // (1,0,0),(0,1,0),(0,0,1)
    try std.testing.expectEqual(@as(usize, 9), tripleCount(2)); // C(5,3)-1 = 9
    try std.testing.expectEqual(@as(usize, 19), tripleCount(3)); // C(6,3)-1 = 19
    try std.testing.expectEqual(@as(usize, 34), tripleCount(4)); // C(7,3)-1 = 34
    try std.testing.expectEqual(@as(usize, 55), tripleCount(5)); // C(8,3)-1 = 55
}

test "enumeration is deterministic and free of duplicates" {
    var buf_a: [64]Triple = undefined;
    var buf_b: [64]Triple = undefined;
    const na = enumerateTriples(&buf_a, MAX_TOTAL);
    const nb = enumerateTriples(&buf_b, MAX_TOTAL);
    try std.testing.expectEqual(na, nb);
    try std.testing.expectEqual(tripleCount(MAX_TOTAL), na);

    for (0..na) |i| {
        try std.testing.expect(buf_a[i].eql(buf_b[i]));
        // no duplicates
        for (i + 1..na) |j| {
            try std.testing.expect(!buf_a[i].eql(buf_a[j]));
        }
        // in range
        try std.testing.expect(buf_a[i].n() >= 1);
        try std.testing.expect(buf_a[i].n() <= MAX_TOTAL);
    }
}

test "SigBuf builds valid signatures and algebras" {
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

test "the convention changes neither the triple nor the roles" {
    var buf: [64]Triple = undefined;
    const n = enumerateTriples(&buf, MAX_TOTAL);
    for (0..n) |i| {
        const sb_minus = SigBuf.build(buf[i], .mostly_minus);
        const sb_plus = SigBuf.build(buf[i], .mostly_plus);
        try std.testing.expectEqual(sb_minus.len, sb_plus.len);
        for (0..sb_minus.len) |j| {
            try std.testing.expectEqual(sb_minus.roles[j], sb_plus.roles[j]);
        }
        // the form signs are opposite, the role is the same
        try std.testing.expectApproxEqAbs(
            sb_minus.signature().signAt(0),
            -sb_plus.signature().signAt(0),
            0,
        );
    }
}

test "a degenerate algebra rejects raising an index" {
    const sb = SigBuf.build(.{ .p = 2, .q = 1, .r = 1 }, .mostly_minus);
    const f = mrs.form.DiagonalForm{ .signature = sb.signature() };
    const c = [_]f64{ 1, 1, 1, 1 };
    var out: [4]f64 = undefined;
    try std.testing.expectError(error.DegenerateForm, f.raiseIndex(&c, &out));
}
