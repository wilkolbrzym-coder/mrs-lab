//! MRS-0 :: Liczby split-complex (hiperboliczne)
//!
//! z = re + j·im,  j² = +1.
//!
//! This is not decoration: the split-complex algebra IS the algebra of Lorentz
//! boosts in the 1+1 plane. Three facts on which the whole module rests:
//!
//!   L1. The norm N(z) = re² - im² coincides with the Minkowski form in 1+1:
//!       for v = (t, x) we have g(v,v) = N(t + j·x).
//!
//!   L2. Zero divisors (N(z) = 0, i.e. re = ±im) correspond to lightlike
//!       vectors. The light cone is the set of zero divisors of the ring.
//!       That explains why a lightlike vector "has no inverse":
//!       it is not a numerical singularity, it is ring structure.
//!
//!   L3. Elements with N = 1 (i.e. e^{jθ} = cosh θ + j·sinh θ) form an
//!       abelian group — they are Lorentz boosts. Composing boosts is
//!       split-complex multiplication, and since the group is abelian,
//!       composition is COMMUTATIVE. That gives thesis T2.

const std = @import("std");

pub const Z = struct {
    re: f64 = 0,
    im: f64 = 0,

    pub const zero: Z = .{ .re = 0, .im = 0 };
    pub const one: Z = .{ .re = 1, .im = 0 };
    /// Jednostka hiperboliczna j, j² = +1.
    pub const j: Z = .{ .re = 0, .im = 1 };

    pub fn init(re: f64, im: f64) Z {
        return .{ .re = re, .im = im };
    }

    pub fn add(a: Z, b: Z) Z {
        return .{ .re = a.re + b.re, .im = a.im + b.im };
    }

    pub fn sub(a: Z, b: Z) Z {
        return .{ .re = a.re - b.re, .im = a.im - b.im };
    }

    /// Naive product: 4 multiplications, 2 additions.
    pub fn mul(a: Z, b: Z) Z {
        return .{
            .re = a.re * b.re + a.im * b.im,
            .im = a.re * b.im + a.im * b.re,
        };
    }

    /// Karatsuba product: 3 multiplications, 5 additions.
    /// Trading one multiplication for three additions pays off when multiplication
    /// is substantially more expensive than addition (fp, extended precision, GPU).
    pub fn mulFast(a: Z, b: Z) Z {
        const k1 = a.re * b.re;
        const k2 = a.im * b.im;
        const k3 = (a.re + a.im) * (b.re + b.im);
        return .{ .re = k1 + k2, .im = k3 - k1 - k2 };
    }

    /// Conjugate with respect to j: re − j·im.
    pub fn conj(a: Z) Z {
        return .{ .re = a.re, .im = -a.im };
    }

    /// Norma formy hiperbolicznej: N(z) = re² − im².
    /// NOT nonnegative — and that is the entire usefulness of this algebra.
    pub fn norm(a: Z) f64 {
        return a.re * a.re - a.im * a.im;
    }

    /// Inverse z^{-1} = conj(z) / N(z). Exists ⟺ N(z) ≠ 0.
    ///
    /// THE COMPUTATION IS SCALED, and that is not decoration. The naive version
    /// `conj(z)/N(z)` with N(z) = re²−im² breaks at both ends of the range:
    ///
    ///   * re = 1e-200 → N underflows to 0 → a false `NotInvertible`,
    ///     although the inverse 1e200 is fully representable;
    ///   * re = 1e200 → N overflows to +inf → `conj/inf = 0`, i.e. the function
    ///     returned the ZERO element as the inverse, silently, and z·z⁻¹ ≠ 1.
    ///
    /// Dlatego najpierw skalujemy: z = s·(r + j·i) przy s = |re|+|im|,
    /// then N(z) = s²·Ñ, where Ñ = r²−i² is computed on order-one numbers, and the
    /// order one, so it neither underflows nor overflows. The result:
    /// conj(r+ji) / (s·Ñ).
    pub fn inv(a: Z) error{NotInvertible}!Z {
        const s = @abs(a.re) + @abs(a.im);
        if (s == 0.0) return error.NotInvertible; // z = 0
        const r = a.re / s;
        const i = a.im / s;
        const nn_scaled = r * r - i * i;
        if (nn_scaled == 0.0) return error.NotInvertible; // prawdziwy dzielnik zera
        const denom = s * nn_scaled;
        return .{ .re = r / denom, .im = -i / denom };
    }

    /// Zero divisor: N(z) = 0. Equivalently z ∈ j·(1+j) ∪ j·(1−j),
    /// i.e. z lies on the light cone (Lemma L2).
    pub fn isZeroDivisor(a: Z, tol: f64) bool {
        return @abs(a.norm()) <= tol;
    }

    /// Boost of rapidity θ: e^{jθ} = cosh θ + j·sinh θ. N = 1 exactly
    /// in exact arithmetic; in fp up to rounding error.
    pub fn fromRapidity(theta: f64) Z {
        return .{ .re = std.math.cosh(theta), .im = std.math.sinh(theta) };
    }

    /// Rapidyta elementu o N = 1: θ = atanh(im/re).
    /// Cost: one transcendental. In MRS-LAB that is a signal that the
    /// matrix representation is the wrong one — see thesis T2.
    pub fn rapidity(a: Z) f64 {
        return std.math.atanh(a.im / a.re);
    }

    /// Macierz boostu w bazie (t, x): [[cosh, sinh], [sinh, cosh]].
    /// Serves only as the reference representation in benchmark T2.
    pub fn boostMatrix(a: Z) [4]f64 {
        return .{ a.re, a.im, a.im, a.re };
    }

    /// Zastosowanie boostu do wektora (t, x).
    pub fn apply(a: Z, t: f64, x: f64) [2]f64 {
        return .{ a.re * t + a.im * x, a.im * t + a.re * x };
    }

    /// Number of multiplications in the operations.
    pub const mulCountNaive: usize = 4;
    pub const mulCountKaratsuba: usize = 3;
};

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "j² = +1" {
    const jj = Z.mul(Z.j, Z.j);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), jj.re, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), jj.im, 1e-15);
}

test "L1: the norm coincides with the Minkowski form in 1+1" {
    const sig = @import("signature.zig");
    const form = @import("form.zig");
    const f = form.DiagonalForm{ .signature = sig.minkowski_1_1 };

    const probes = [_][2]f64{ .{ 1, 0 }, .{ 0, 1 }, .{ 1, 1 }, .{ 1, -1 }, .{ 2.5, -0.75 } };
    for (probes) |p| {
        const z = Z.init(p[0], p[1]);
        try std.testing.expectApproxEqAbs(f.eval(&p), Z.norm(z), 1e-13);
    }
}

test "L2: zero divisor equals lightlike vector" {
    // (1+j)(1−j) = 1 − j² = 0
    const a = Z.init(1, 1);
    const b = Z.init(1, -1);
    const p = Z.mul(a, b);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), p.re, 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), p.im, 1e-15);
    try std.testing.expect(Z.isZeroDivisor(a, 1e-15));
    try std.testing.expect(Z.isZeroDivisor(b, 1e-15));
    try std.testing.expectError(error.NotInvertible, Z.inv(a));

    // a plain representative of the cone: (t,x) = (3,3)
    const light = Z.init(3, 3);
    try std.testing.expect(Z.isZeroDivisor(light, 1e-15));
    try std.testing.expectError(error.NotInvertible, Z.inv(light));
}

test "L3: boosts form an abelian group" {
    var prng = std.Random.DefaultPrng.init(0xC0FFEE);
    const rnd = prng.random();

    for (0..200) |_| {
        const th1 = (rnd.float(f64) - 0.5) * 6.0;
        const th2 = (rnd.float(f64) - 0.5) * 6.0;
        const b1 = Z.fromRapidity(th1);
        const b2 = Z.fromRapidity(th2);

        // N = 1 (closure on the group)
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), Z.norm(b1), 1e-12);

        // commutativity
        const p1 = Z.mul(b1, b2);
        const p2 = Z.mul(b2, b1);
        try std.testing.expectApproxEqAbs(p1.re, p2.re, 1e-14);
        try std.testing.expectApproxEqAbs(p1.im, p2.im, 1e-14);

        // additivity of rapidity: θ1+θ2. The tolerance is looser than elsewhere and
        // elsewhere and that is DELIBERATE: recovering the rapidity through atanh
        // that is DELIBERATE: recovering the rapidity through atanh loses precision
        // for large |θ| (atanh is steep near argument 1). That is exactly the cost
        // MRS-LAB avoids by keeping the rapidity as the coordinate.
        const ths = Z.rapidity(p1);
        try std.testing.expectApproxEqAbs(th1 + th2, ths, 1e-8);

        // agreement with the boost matrix applied to a vector
        const mv = Z.boostMatrix(p1);
        const viaMatrix = [2]f64{
            mv[0] * 0.7 + mv[1] * -1.3,
            mv[2] * 0.7 + mv[3] * -1.3,
        };
        const viaZ = Z.apply(p1, 0.7, -1.3);
        try std.testing.expectApproxEqAbs(viaMatrix[0], viaZ[0], 1e-13);
        try std.testing.expectApproxEqAbs(viaMatrix[1], viaZ[1], 1e-13);
    }
}

test "Karatsuba agrees with the naive product" {
    var prng = std.Random.DefaultPrng.init(7);
    const rnd = prng.random();
    for (0..1000) |_| {
        const a = Z.init(rnd.float(f64) * 4 - 2, rnd.float(f64) * 4 - 2);
        const b = Z.init(rnd.float(f64) * 4 - 2, rnd.float(f64) * 4 - 2);
        const p = Z.mul(a, b);
        const q = Z.mulFast(a, b);
        try std.testing.expectApproxEqRel(p.re, q.re, 1e-12);
        try std.testing.expectApproxEqRel(p.im, q.im, 1e-12);
    }
}

test "the norm is multiplicative: N(z1·z2) = N(z1)·N(z2)" {
    // This is why split-complex is the right language for the Lorentz GROUP and
    // not for tachyons themselves: multiplicativity of the norm preserves the
    // class only because the sign of the product is the product of the signs.
    // The set N = 1 is a subgroup; the set N < 0 is not.
    var prng = std.Random.DefaultPrng.init(1234);
    const rnd = prng.random();
    for (0..2000) |_| {
        const a = Z.init(rnd.float(f64) * 6 - 3, rnd.float(f64) * 6 - 3);
        const b = Z.init(rnd.float(f64) * 6 - 3, rnd.float(f64) * 6 - 3);
        // ABSOLUTE tolerance, not relative: the product of norms can be close to
        // zero (when one factor is near the null cone), and then the relative error
        // grows without bound. The same trap is noted in check A3 for the interval
        // invariant.
        try std.testing.expectApproxEqAbs(Z.norm(Z.mul(a, b)), Z.norm(a) * Z.norm(b), 1e-12);
    }
}

test "tachyons do NOT form a subgroup: two spatial elements multiply to a temporal one" {
    // N = -1 (spatial) times N = -1 gives N = +1 (temporal), so the set of
    // tachyonic elements is not closed under multiplication. This is the precise
    // version of "split-complex is the language of tachyons": it is the language
    // of the NULL CONE and the BOOST GROUP, not of tachyons.
    const a = Z.init(0, 1); // N(a) = −1
    const b = Z.init(0, 1);
    try std.testing.expect(Z.norm(a) < 0);
    const p = Z.mul(a, b);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), p.re, 1e-15);
    try std.testing.expect(Z.norm(p) > 0);
    // whereas N = 1 is closed: it is the boost group
    const g1 = Z.fromRapidity(0.7);
    const g2 = Z.fromRapidity(-1.3);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), Z.norm(Z.mul(g1, g2)), 1e-15);
}

test "REGRESSION: the inverse works at both ends of the range" {
    // Witnesses measured before the fix (probe against the real module):
    //   z=(1e-200,0) → returned error.NotInvertible, although 1/re = 1e200 exists
    //   z=(1e200,0)  → returned (0,-0), i.e. ZERO as the inverse, silently
    // Both cases must now give the correct inverse and z·z⁻¹ = 1.
    // Note: pairs with |re| = |im| are TRUE zero divisors (N = 0) and must raise
    // an error — hence the boundary cases are chosen off the null cone.
    const cases = [_]Z{
        Z.init(1e-200, 0),
        Z.init(1e200, 0),
        Z.init(1e-160, 0),
        Z.init(1e160, 0),
        Z.init(1e-200, 1e-201),
        Z.init(1e200, 1e199), // Z.init(1e200, 1e199), // under the naive version N overflowed to +inf
        Z.init(3, 1),
        Z.init(0, 1),
        Z.init(-2.5, 0.75),
    };
    for (cases) |z| {
        const iv = try Z.inv(z);
        try std.testing.expect(std.math.isFinite(iv.re) and std.math.isFinite(iv.im));
        const p = Z.mul(z, iv);
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), p.re, 1e-12);
        try std.testing.expectApproxEqAbs(@as(f64, 0.0), p.im, 1e-12);
    }

    // True zero divisors must still be rejected — and exactly.
    try std.testing.expectError(error.NotInvertible, Z.inv(Z.init(0, 0)));
    try std.testing.expectError(error.NotInvertible, Z.inv(Z.init(1, 1)));
    try std.testing.expectError(error.NotInvertible, Z.inv(Z.init(-3, 3)));
    try std.testing.expectError(error.NotInvertible, Z.inv(Z.init(1e-300, 1e-300)));
}
