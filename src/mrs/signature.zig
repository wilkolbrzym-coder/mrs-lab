//! MRS-LAB :: signature module
//!
//! The signature is a first-class type in MRS-LAB, not incidental information
//! about a metric matrix. Two things follow:
//!
//!   1. The shape of the quadratic form is known statically (diagonal), so
//!      evaluation `g(v,v)` is O(n), not O(n^2).
//!   2. The sign convention (`time_sign`) is an explicit parameter, so the same
//!      mathematical object can be described in (+,-,-,-) and (-,+,+,+) without
//!      renaming theorems.
//!
//! Terminology note: the word "tachyon" is NOT invariant under the choice of
//! convention, so we define it relative to the declared time cone rather than
//! relative to the raw sign of g(v,v).

const std = @import("std");

/// Largest dimension handled by the fixed-size buffers used across the library
/// and by the property engine. Chosen so that 2^8 = 256 blades stay cheap.
pub const MAX_DIM: usize = 8;

/// Sign of the temporal dimensions. The enum VALUE is the sign of the temporal
/// dimension, while the name says which metric convention is meant — that
/// distinction is the source of the most common error in projects of this kind,
/// so it is resolved in one place.
pub const TimeSign = enum(i8) {
    /// Convention (+,-,-,-): temporal dimensions have sign +1,
    /// spatial −1. A temporal vector has g > 0.
    mostly_minus = 1,
    /// Convention (-,+,+,+): temporal dimensions have sign -1,
    /// spatial +1. A temporal vector has g < 0.
    mostly_plus = -1,

    pub fn f(self: TimeSign) f64 {
        return @floatFromInt(@intFromEnum(self));
    }
};

/// Role of a dimension. This is the content of the triple (p,q,r): p = number of
/// `temporal`, q = number of `spatial`, r = number of `degenerate` (the radical).
pub const Role = enum {
    temporal,
    spatial,
    degenerate,
};

pub const Signature = struct {
    /// Length = dim X. The order of dimensions is part of the representation.
    roles: []const Role,
    time_sign: TimeSign = .mostly_minus,

    pub fn n(self: Signature) usize {
        return self.roles.len;
    }

    pub fn count(self: Signature, role: Role) usize {
        var c: usize = 0;
        for (self.roles) |ri| {
            if (ri == role) c += 1;
        }
        return c;
    }

    /// p — number of temporal dimensions.
    pub fn p(self: Signature) usize {
        return self.count(.temporal);
    }

    /// q — number of spatial dimensions.
    pub fn q(self: Signature) usize {
        return self.count(.spatial);
    }

    /// r — number of degenerate dimensions (the radical of the form).
    pub fn r(self: Signature) usize {
        return self.count(.degenerate);
    }

    /// Sign of the i-th basis vector: +1, −1 or 0.
    pub fn signAt(self: Signature, i: usize) f64 {
        return switch (self.roles[i]) {
            .temporal => self.time_sign.f(),
            .spatial => -self.time_sign.f(),
            .degenerate => 0.0,
        };
    }

    pub fn isNondegenerate(self: Signature) bool {
        return self.r() == 0;
    }

    /// Lorentzian = exactly one temporal dimension and no radical.
    /// This is the condition under which the causal relation is transitive
    /// (Theorem 5.1 in `mrs/causal.zig`).
    pub fn isLorentzian(self: Signature) bool {
        return self.p() == 1 and self.isNondegenerate();
    }

    /// The signature has more than one temporal dimension — the zone where the
    /// naive causality derived from the form stops being transitive.
    pub fn hasMultipleTimes(self: Signature) bool {
        return self.p() > 1;
    }

    /// Prints the signature in (p,q,r) notation together with the convention.
    pub fn writeTo(self: Signature, w: anytype) !void {
        try w.print("({d},{d},{d}) [{s}]", .{
            self.p(), self.q(), self.r(),
            if (self.time_sign == .mostly_minus) "+,−,−,…" else "−,+,+,…",
        });
    }
};

// ---------------------------------------------------------------------------
// Canonical signatures used in tests and benchmarks.
// ---------------------------------------------------------------------------

const t = Role.temporal;
const s = Role.spatial;
const d = Role.degenerate;

/// Minkowski 3+1, convention (+,−,−,−). The validation target of MRS-0.
pub const minkowski_3_1: Signature = .{
    .roles = &.{ t, s, s, s },
    .time_sign = .mostly_minus,
};

/// Minkowski 3+1 in the (−,+,+,+) convention. The same object, another description.
pub const minkowski_3_1_flipped: Signature = .{
    .roles = &.{ t, s, s, s },
    .time_sign = .mostly_plus,
};

/// Minkowski 1+1 — the arena of the split-complex algebra.
pub const minkowski_1_1: Signature = .{ .roles = &.{ t, s } };

/// Signature (2,1): two times, one space. The transitivity witness of the
/// causal relation — three points `u ⪯ w ⪯ v` with `u ⋠ v`, built in
/// `mrs/causal.zig :: canonicalWitness`.
pub const two_times_2_1: Signature = .{ .roles = &.{ t, t, s } };

/// (2,1,1) — with a degenerate dimension. The cone is not proper and the form
/// has no inverse on the whole space.
pub const degenerate_2_1_1: Signature = .{ .roles = &.{ t, t, s, d } };

/// (0,4): Euclidean — no temporal dimension, no cone.
pub const euclidean_4: Signature = .{ .roles = &.{ s, s, s, s } };

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "counting (p,q,r)" {
    try std.testing.expectEqual(@as(usize, 1), minkowski_3_1.p());
    try std.testing.expectEqual(@as(usize, 3), minkowski_3_1.q());
    try std.testing.expectEqual(@as(usize, 0), minkowski_3_1.r());
    try std.testing.expectEqual(@as(usize, 4), minkowski_3_1.n());
    try std.testing.expect(minkowski_3_1.isLorentzian());

    try std.testing.expectEqual(@as(usize, 2), two_times_2_1.p());
    try std.testing.expect(two_times_2_1.hasMultipleTimes());
    try std.testing.expect(!two_times_2_1.isLorentzian());

    try std.testing.expectEqual(@as(usize, 1), degenerate_2_1_1.r());
    try std.testing.expect(!degenerate_2_1_1.isNondegenerate());

    try std.testing.expectEqual(@as(usize, 0), euclidean_4.p());
}

test "the sign convention does not change the role of a dimension" {
    // The same vector under two conventions: the form sign flips, the classification does not
    // (temporal/spatial) stays the same.
    const mm = minkowski_3_1.signAt(0);
    const mp = minkowski_3_1_flipped.signAt(0);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), mm, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -1.0), mp, 1e-15);
}
