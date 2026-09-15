//! MRS-LAB :: form module
//!
//! Dwie reprezentacje tej samej formy kwadratowej g: X × X → K:
//!
//!   * `DiagonalForm` — the MRS-native representation: the form is given
//!     by the signature, so `g(v,v) = Σ s_i v_i²`. Cost O(n); the inverse
//!     (raising an index) is also O(n), because the matrix is diagonal.
//!
//!   * `DenseForm` — the "textbook" representation: a full n×n matrix.
//!     Evaluation cost O(n^2), inverse O(n^3) (Gauss-Jordan).
//!
//! Both describe THE SAME mathematical object. The cost difference is purely
//! representational — and that is falsifiable thesis T1.

const std = @import("std");
const sig = @import("signature.zig");

const Signature = sig.Signature;
const TimeSign = sig.TimeSign;

/// Classification of an element with respect to the declared time cone.
/// Note: this is NOT an "absolute" classification — it depends on `time_sign`.
pub const Class = enum {
    /// Consistent with the time cone (in 3+1: temporal).
    temporal,
    /// Self-orthogonal, g(v,v) = 0 (in 3+1: lightlike).
    null_like,
    /// Opposite to the cone (in 3+1: spatial). Elements with nonzero norm here
    /// normie w tej klasie to kandydaci na tachiony.
    spatial,

    pub fn isCausalLike(self: Class) bool {
        return self == .temporal or self == .null_like;
    }

    pub fn label(self: Class) []const u8 {
        return switch (self) {
            .temporal => "czasowy",
            .null_like => "null (lightlike)",
            .spatial => "przestrzenny (tachionowy)",
        };
    }
};

// ---------------------------------------------------------------------------
// MRS-native representation: the diagonal form given by the signature
// ---------------------------------------------------------------------------

pub const DiagonalForm = struct {
    signature: Signature,
    /// Precomputed signs, used when `use_signs` is set. Reason: `signAt`
    /// branches on the role of every dimension, and MEASUREMENT (see the
    /// benchmark review) showed that branch is 74% of the evaluation cost at
    /// n = 128: 1709 ns with the switch against 439 ns reading a contiguous
    /// sign array. The field has a default, so existing struct literals keep
    /// their exact previous behaviour.
    signs: [sig.MAX_DIM]f64 = [_]f64{0} ** sig.MAX_DIM,
    use_signs: bool = false,

    /// Preferred constructor: precomputes the sign vector once, so the hot loop
    /// is branch-free and reads contiguous memory.
    ///
    /// BOUNDS. The cache is `sig.MAX_DIM` entries. Beyond that dimension the
    /// constructor deliberately falls back to the general path instead of
    /// overflowing: an earlier version wrote past the array, which in ReleaseFast
    /// is a segfault rather than a panic. The guard is the fix and is covered by
    /// a test with n > MAX_DIM.
    pub fn init(s: Signature) DiagonalForm {
        if (s.n() > sig.MAX_DIM) return .{ .signature = s, .use_signs = false };
        var f = DiagonalForm{ .signature = s, .use_signs = true };
        for (0..s.n()) |i| f.signs[i] = s.signAt(i);
        return f;
    }

    pub fn n(self: DiagonalForm) usize {
        return self.signature.n();
    }

    /// g(v,v) = Σ_i s_i v_i². Cost: n multiply-adds.
    pub fn eval(self: DiagonalForm, v: []const f64) f64 {
        std.debug.assert(v.len == self.n());
        var acc: f64 = 0;
        if (self.use_signs) {
            for (v, 0..) |vi, i| acc += self.signs[i] * vi * vi;
        } else {
            for (v, 0..) |vi, i| acc += self.signature.signAt(i) * vi * vi;
        }
        return acc;
    }

    /// Bilinear form: g(u,v) = Σ_i s_i u_i v_i. Cost O(n).
    pub fn bilinear(self: DiagonalForm, u: []const f64, v: []const f64) f64 {
        std.debug.assert(u.len == v.len and u.len == self.n());
        var acc: f64 = 0;
        if (self.use_signs) {
            for (u, v, 0..) |ui, vi, i| acc += self.signs[i] * ui * vi;
        } else {
            for (u, v, 0..) |ui, vi, i| acc += self.signature.signAt(i) * ui * vi;
        }
        return acc;
    }

    /// Classification with respect to the declared convention.
    pub fn classify(self: DiagonalForm, v: []const f64) Class {
        return self.classifyTol(v, 0.0);
    }

    pub fn classifyTol(self: DiagonalForm, v: []const f64, tol: f64) Class {
        const g = self.eval(v);
        if (@abs(g) <= tol) return .null_like;
        const ts = self.signature.time_sign.f();
        return if (g * ts > 0) .temporal else .spatial;
    }

    /// Is v a tachyon: opposite cone and nonzero norm.
    /// This definition is INVARIANT under the choice of sign convention.
    pub fn isTachyon(self: DiagonalForm, v: []const f64, tol: f64) bool {
        return self.classifyTol(v, tol) == .spatial;
    }

    /// Raising an index: out = g^{-1} c, i.e. out_i = c_i / s_i.
    /// Cost O(n). For degenerate dimensions (s_i = 0) it returns an error,
    /// because the form has no inverse there.
    pub fn raiseIndex(self: DiagonalForm, c: []const f64, out: []f64) error{DegenerateForm}!void {
        std.debug.assert(c.len == out.len and c.len == self.n());
        for (c, 0..) |ci, i| {
            const s_i = if (self.use_signs) self.signs[i] else self.signature.signAt(i);
            if (s_i == 0.0) return error.DegenerateForm;
            out[i] = ci / s_i;
        }
    }

    /// Operacje zmiennoprzecinkowe potrzebne do jednej ewaluacji formy.
    pub fn mulCountEval(self: DiagonalForm) usize {
        return self.n();
    }

    /// Operations needed for the inverse of the form.
    pub fn mulCountInverse(self: DiagonalForm) usize {
        return self.n();
    }
};

// ---------------------------------------------------------------------------
// Reference representation: the full dense matrix
// ---------------------------------------------------------------------------

pub const DenseForm = struct {
    dim: usize,
    /// Wiersz po wierszu (row-major), dim × dim.
    g: []const f64,

    /// Builds the diagonal matrix corresponding to the signature. This is the same
    /// form as `DiagonalForm`, only written densely.
    pub fn fromDiagonal(alloc: std.mem.Allocator, s: Signature) !DenseForm {
        const nn = s.n();
        const g = try alloc.alloc(f64, nn * nn);
        @memset(g, 0.0);
        for (0..nn) |i| g[i * nn + i] = s.signAt(i);
        return .{ .dim = nn, .g = g };
    }

    /// g(v,v) = vᵀ G v. Cost: n^2 multiplications (no structure exploited).
    pub fn eval(self: DenseForm, v: []const f64) f64 {
        const nn = self.dim;
        std.debug.assert(v.len == nn);
        var acc: f64 = 0;
        for (0..nn) |i| {
            const gi = self.g[i * nn ..][0..nn];
            var row: f64 = 0;
            for (0..nn) |j| row += gi[j] * v[j];
            acc += v[i] * row;
        }
        return acc;
    }

    /// Inverse by Gauss-Jordan elimination with partial pivoting.
    /// Cost O(n^3). `out` must have n*n elements.
    pub fn inverse(
        self: DenseForm,
        alloc: std.mem.Allocator,
        out: []f64,
    ) !void {
        const nn = self.dim;
        std.debug.assert(out.len == nn * nn);

        const aug = try alloc.alloc(f64, nn * 2 * nn);
        defer alloc.free(aug);

        for (0..nn) |i| {
            for (0..nn) |j| {
                aug[i * 2 * nn + j] = self.g[i * nn + j];
                aug[i * 2 * nn + nn + j] = if (i == j) 1.0 else 0.0;
            }
        }

        for (0..nn) |col| {
            // pivot
            var piv = col;
            var best = @abs(aug[col * 2 * nn + col]);
            for (col + 1..nn) |row| {
                const cand = @abs(aug[row * 2 * nn + col]);
                if (cand > best) {
                    best = cand;
                    piv = row;
                }
            }
            if (best < 1e-300) return error.SingularForm;
            if (piv != col) {
                const w = 2 * nn;
                for (0..w) |k| {
                    const tmp = aug[col * w + k];
                    aug[col * w + k] = aug[piv * w + k];
                    aug[piv * w + k] = tmp;
                }
            }
            const d = aug[col * 2 * nn + col];
            for (0..2 * nn) |k| aug[col * 2 * nn + k] /= d;
            for (0..nn) |row| {
                if (row == col) continue;
                const f = aug[row * 2 * nn + col];
                if (f == 0.0) continue;
                for (0..2 * nn) |k| {
                    aug[row * 2 * nn + k] -= f * aug[col * 2 * nn + k];
                }
            }
        }

        for (0..nn) |i| {
            for (0..nn) |j| out[i * nn + j] = aug[i * 2 * nn + nn + j];
        }
    }

    pub fn mulCountEval(self: DenseForm) usize {
        return self.dim * self.dim;
    }

    /// Rough operation count for Gauss-Jordan: O(n^3).
    pub fn mulCountInverse(self: DenseForm) usize {
        return 2 * self.dim * self.dim * self.dim;
    }

    pub fn deinit(self: DenseForm, alloc: std.mem.Allocator) void {
        alloc.free(@constCast(self.g));
    }
};

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "both representations give the identical form" {
    const alloc = std.testing.allocator;
    const s = sig.minkowski_3_1;
    const diag = DiagonalForm{ .signature = s };
    const dense = try DenseForm.fromDiagonal(alloc, s);
    defer dense.deinit(alloc);

    const probes = [_][4]f64{
        .{ 1, 0, 0, 0 }, // czasowy
        .{ 0, 1, 0, 0 }, // przestrzenny
        .{ 1, 1, 0, 0 }, // zerowy
        .{ 0.3, -1.7, 2.2, 0.5 },
        .{ 1e3, -2e3, 3e3, -4e3 },
    };
    for (probes) |p| {
        try std.testing.expectApproxEqRel(diag.eval(&p), dense.eval(&p), 1e-12);
    }
}

test "dense inverse equals diagonal inverse" {
    const alloc = std.testing.allocator;
    const s = sig.minkowski_3_1;
    const diag = DiagonalForm{ .signature = s };
    const dense = try DenseForm.fromDiagonal(alloc, s);
    defer dense.deinit(alloc);

    const nn = s.n();
    const inv = try alloc.alloc(f64, nn * nn);
    defer alloc.free(inv);
    try dense.inverse(alloc, inv);

    // g · g^{-1} = I
    for (0..nn) |i| {
        for (0..nn) |j| {
            var acc: f64 = 0;
            for (0..nn) |k| acc += dense.g[i * nn + k] * inv[k * nn + j];
            try std.testing.expectApproxEqAbs(if (i == j) @as(f64, 1) else @as(f64, 0), acc, 1e-12);
        }
    }

    // raiseIndex agrees between the diagonal and the dense representation
    const c = [_]f64{ 2.0, -4.0, 6.0, -8.0 };
    var raised: [4]f64 = undefined;
    try diag.raiseIndex(&c, &raised);
    for (0..nn) |i| {
        var acc: f64 = 0;
        for (0..nn) |j| acc += inv[i * nn + j] * c[j];
        try std.testing.expectApproxEqAbs(raised[i], acc, 1e-12);
    }
}

test "classification depends on the convention, not on the object" {
    const mm = DiagonalForm{ .signature = sig.minkowski_3_1 };
    const mp = DiagonalForm{ .signature = sig.minkowski_3_1_flipped };

    const v_time = [_]f64{ 2, 0, 0, 0 };
    const v_space = [_]f64{ 0, 2, 0, 0 };
    const v_null = [_]f64{ 1, 1, 0, 0 };

    for ([_]DiagonalForm{ mm, mp }) |f| {
        try std.testing.expectEqual(Class.temporal, f.classify(&v_time));
        try std.testing.expectEqual(Class.spatial, f.classify(&v_space));
        try std.testing.expectEqual(Class.null_like, f.classify(&v_null));
        // tachyon = spatial with nonzero norm, independent of the convention
        try std.testing.expect(f.isTachyon(&v_space, 1e-12));
        try std.testing.expect(!f.isTachyon(&v_time, 1e-12));
        try std.testing.expect(!f.isTachyon(&v_null, 1e-12));
    }
}

test "a degenerate dimension blocks the inverse of the form" {
    const f = DiagonalForm{ .signature = sig.degenerate_2_1_1 };
    const c = [_]f64{ 1, 1, 1, 1 };
    var out: [4]f64 = undefined;
    try std.testing.expectError(error.DegenerateForm, f.raiseIndex(&c, &out));
}

test "init above MAX_DIM falls back instead of overflowing" {
    // Witness of a real defect: before the guard this wrote past `signs`
    // and segfaulted in ReleaseFast for n = 1024.
    var roles: [64]sig.Role = undefined;
    for (&roles) |*r| r.* = .spatial;
    const s = Signature{ .roles = &roles };
    const f = DiagonalForm.init(s);
    try std.testing.expect(!f.use_signs);
    var v: [64]f64 = undefined;
    for (&v) |*x| x.* = 1.0;
    try std.testing.expectApproxEqAbs(@as(f64, -64.0), f.eval(&v), 1e-12);

    // and at MAX_DIM the fast path is used and agrees with the general path
    var roles8: [sig.MAX_DIM]sig.Role = undefined;
    for (&roles8) |*r| r.* = .spatial;
    roles8[0] = .temporal;
    const s8 = Signature{ .roles = &roles8 };
    const f8 = DiagonalForm.init(s8);
    const g8 = DiagonalForm{ .signature = s8 };
    try std.testing.expect(f8.use_signs);
    var v8: [sig.MAX_DIM]f64 = undefined;
    for (&v8, 0..) |*x, i| x.* = @as(f64, @floatFromInt(i)) - 3.0;
    try std.testing.expectApproxEqAbs(g8.eval(&v8), f8.eval(&v8), 0.0);
}
