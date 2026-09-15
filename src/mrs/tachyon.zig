//! MRS-LAB :: tachyon module
//!
//! Dispersion relation:  E² = p² + m²,  where m² ∈ R.
//!
//! A note on terminology, and a real correction to the original project spec:
//! for m² < 0 it is NOT enough to write "E² = p² + m² with m² < 0" and call it
//! done, because the equation has two qualitatively different regimes separated
//! by the threshold |p| = μ, where μ² = −m²:
//!
//!   1. |p| > μ  — E real, the mode propagates, and the group velocity
//!      v_g = |p| / E = |p| / √(p² − μ²) is ALWAYS greater than 1.
//!      That is superluminality without instability.
//!
//!   2. |p| < μ  — E² < 0, so E is imaginary. The mode does not propagate;
//!      it grows exponentially at rate γ = √(μ² − p²).
//!      That is the tachyonic (spinodal) instability.
//!
//!   3. |p| = μ  — the threshold; E = 0 and v_g → ∞. The singularity sits here.
//!
//! Link to the form module: for a four-momentum (E, p) we have
//! g(P,P) = E² − p² = m². A mode is a tachyon EXACTLY when its four-momentum
//! is in the `spatial` class with respect to the declared convention — so the
//! definition in `form` and the definition from the dispersion relation are
//! the same definition.
//!
//! Computational layer: dual numbers give the exact derivative dE/dp in a
//! single pass, with no finite-difference step. That is a measurable advantage
//! (see `scanPhaseAccuracy` and results/RESULTS.md, thesis T4).
//!
//! Scope note: this module is deliberately small. It contains what the engine
//! and the benchmarks actually use, and nothing else.

const std = @import("std");
const dual = @import("dual.zig");
const form = @import("form.zig");
const sig = @import("signature.zig");

pub const Tachyon = struct {
    /// μ² = −m² > 0
    mu2: f64,

    pub fn init(mu2: f64) error{NegativeMu2}!Tachyon {
        if (mu2 < 0.0) return error.NegativeMu2;
        return .{ .mu2 = mu2 };
    }

    pub fn mu(self: Tachyon) f64 {
        return @sqrt(self.mu2);
    }

    /// E² = p² + m² = p² − μ².
    pub fn energy2(self: Tachyon, p: f64) f64 {
        return p * p - self.mu2;
    }

    /// Does the mode propagate (is E real)?
    pub fn isPropagating(self: Tachyon, p: f64) bool {
        return self.energy2(p) > 0.0;
    }

    /// Is the mode unstable (E imaginary, exponential growth)?
    pub fn isUnstable(self: Tachyon, p: f64) bool {
        return self.energy2(p) < 0.0;
    }

    /// Instability threshold: |p| = μ.
    pub fn threshold(self: Tachyon) f64 {
        return self.mu();
    }

    /// E = √(p² − μ²). Error when the mode does not propagate.
    pub fn energy(self: Tachyon, p: f64) error{ImaginaryEnergy}!f64 {
        const e2 = self.energy2(p);
        if (e2 < 0.0) return error.ImaginaryEnergy;
        return @sqrt(e2);
    }

    /// Growth rate γ = √(μ² − p²) for |p| < μ; zero otherwise.
    pub fn growthRate(self: Tachyon, p: f64) f64 {
        const e2 = self.energy2(p);
        return if (e2 < 0.0) @sqrt(-e2) else 0.0;
    }

    /// Group velocity v_g = dE/dp = p/E. For |p| > μ it is always > 1 —
    /// superluminality here is generic, not an accident.
    pub fn groupVelocity(self: Tachyon, p: f64) error{ImaginaryEnergy}!f64 {
        return p / try self.energy(p);
    }

    /// Minkowski form of the four-momentum (E, p, 0, …, 0) in a given
    /// convention: g(P,P) = s_t·(E² − p²), where s_t = time_sign.
    /// The VALUE depends on the convention (−μ² for (+,−,−,−), +μ² for
    /// (−,+,+,+)), but the CLASS does not — see `momentumClass`.
    pub fn fourMomentumNorm2(self: Tachyon, E: f64, p: f64, ts: sig.TimeSign) f64 {
        _ = self;
        return ts.f() * (E * E - p * p);
    }

    /// Class of the four-momentum of a mode, inside the form module's classes.
    /// This ties the two modules into one definition instead of two
    /// independent ones.
    ///
    /// Note: on shell E² = p² − μ², so m² = E² − p² = −μ² < 0, i.e. the
    /// four-momentum is SPATIAL in both regimes (propagating and unstable).
    /// That is the definition of a tachyon and it does NOT depend on the sign
    /// convention, because the class depends on s_t·g and s_t·s_t·m² = m².
    pub fn momentumClass(self: Tachyon, p: f64) form.Class {
        _ = p;
        const m2 = -self.mu2; // E² − p² on shell
        if (m2 == 0.0) return .null_like;
        return if (m2 > 0.0) .temporal else .spatial;
    }
};

// ---------------------------------------------------------------------------
// Derivative of the dispersion relation: dual numbers
// ---------------------------------------------------------------------------

/// E(p) = √(p² − μ²) evaluated in dual numbers.
/// Returns (E, dE/dp) exactly — with no finite-difference step.
pub fn energyDual(p: dual.D, mu2: f64) dual.D {
    const p2 = dual.D.mul(p, p);
    const shifted = dual.D.sub(p2, dual.D.cst(mu2));
    return dual.D.sqrt(shifted);
}

/// dE/dp from dual algebra. Exact value: p/E.
pub fn energyDerivAD(p: f64, mu2: f64) f64 {
    return energyDual(dual.D.at(p), mu2).der;
}

/// dE/dp analytically — the reference formula.
pub fn energyDerivExact(p: f64, mu2: f64) f64 {
    return p / @sqrt(p * p - mu2);
}

// ---------------------------------------------------------------------------
// Composite function: energy with phase
// ---------------------------------------------------------------------------
//
// METHODOLOGICAL NOTE. For E(p) = √(p² − μ²) alone, the dual-algebra
// derivative comes out EXACTLY equal to the analytic formula p/E, because the
// chain rule reproduces that same formula. Such a test would be a tautology —
// it would check nothing beyond the correctness of the differentiation rules.
// That is why T4 measures on a function whose derivative is genuinely
// composite:
//
//     f(p) = E(p)·sin(p)
//     f'(p) = (p/E(p))·sin(p) + E(p)·cos(p)

pub fn phaseEnergyDual(p: dual.D, mu2: f64) dual.D {
    return dual.D.mul(energyDual(p, mu2), dual.D.sin(p));
}

pub fn phaseEnergyAt(p: f64, mu2: f64) f64 {
    return @sqrt(p * p - mu2) * @sin(p);
}

/// Derivative of f from dual algebra: one pass, no finite-difference step.
pub fn phaseEnergyDerivAD(p: f64, mu2: f64) f64 {
    return phaseEnergyDual(dual.D.at(p), mu2).der;
}

/// Derivative of f from central differences: two function calls per derivative.
pub fn phaseEnergyDerivFD(p: f64, mu2: f64, h: f64) f64 {
    return dual.centralDiff(phaseEnergyAt, p, mu2, h);
}

/// Analytic reference: f' = (p/E)·sin p + E·cos p.
pub fn phaseEnergyDerivExact(p: f64, mu2: f64) f64 {
    const E = @sqrt(p * p - mu2);
    return (p / E) * @sin(p) + E * @cos(p);
}

pub const PhaseAccuracy = struct {
    samples: usize,
    max_ad_error: f64,
    max_fd_error: f64,
    h: f64,
};

/// Derivative scan over a range of momenta. We compare the absolute error
/// against the analytic reference. This is thesis T4.
pub fn scanPhaseAccuracy(
    mu2: f64,
    p_from: f64,
    p_to: f64,
    samples: usize,
    h: f64,
) PhaseAccuracy {
    var max_ad: f64 = 0;
    var max_fd: f64 = 0;
    for (0..samples) |i| {
        const t = @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(samples - 1));
        const p = p_from + t * (p_to - p_from);
        const exact = phaseEnergyDerivExact(p, mu2);
        max_ad = @max(max_ad, @abs(phaseEnergyDerivAD(p, mu2) - exact));
        max_fd = @max(max_fd, @abs(phaseEnergyDerivFD(p, mu2, h) - exact));
    }
    return .{ .samples = samples, .max_ad_error = max_ad, .max_fd_error = max_fd, .h = h };
}

// ---------------------------------------------------------------------------
// Finding the superluminal momentum: Newton with dual derivatives
// ---------------------------------------------------------------------------

pub const NewtonResult = struct {
    p: f64,
    iterations: usize,
    converged: bool,
};

/// v_g(p) = p/E(p) in dual numbers.
pub fn groupVelocityDual(p: dual.D, mu2: f64) dual.D {
    return dual.D.div(p, energyDual(p, mu2));
}

/// Analytic solution of v_g(p) = k: p* = k·μ / √(k² − 1).
/// It serves as a reference for the Newton iteration, not as part of it.
pub fn superluminalMomentum(mu2: f64, k: f64) f64 {
    const mu = @sqrt(mu2);
    return k * mu / @sqrt(k * k - 1.0);
}

/// v_g(p) = p/E(p) as a plain function — analytic reference.
pub fn groupVelocityAt(p: f64, mu2: f64) f64 {
    return p / @sqrt(p * p - mu2);
}

/// Newton with derivatives from dual algebra.
///
/// NUMERICAL NOTE. Start just above the threshold |p| = μ. Starting on the
/// other side of the solution (p > p*) throws Newton past the asymptote at
/// p = μ onto the negative branch: v_g is an odd function, so the equation also
/// has the root −p*. That is a property of the problem, not an implementation
/// defect, and the test below pins it down.
pub fn newtonAD(mu2: f64, k: f64, p0: f64, tol: f64, max_iter: usize) NewtonResult {
    var p = p0;
    for (0..max_iter) |i| {
        const gv = groupVelocityDual(dual.D.at(p), mu2);
        const f = gv.val - k;
        if (@abs(f) < tol) return .{ .p = p, .iterations = i, .converged = true };
        if (gv.der == 0.0) break;
        p -= f / gv.der;
    }
    return .{ .p = p, .iterations = max_iter, .converged = false };
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "three dispersion regimes and the threshold" {
    const t = try Tachyon.init(9.0); // mu = 3
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), t.threshold(), 1e-15);

    // propagation: |p| > mu
    try std.testing.expect(t.isPropagating(5.0));
    try std.testing.expect(!t.isPropagating(2.0));
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), try t.energy(5.0), 1e-15); // sqrt(25-9)

    // instability: |p| < mu
    try std.testing.expect(t.isUnstable(2.0));
    try std.testing.expectError(error.ImaginaryEnergy, t.energy(2.0));

    // growth rate: gamma = sqrt(9 - 4) = sqrt(5)
    try std.testing.expectApproxEqAbs(@sqrt(@as(f64, 5.0)), t.growthRate(2.0), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), t.growthRate(5.0), 1e-15);

    // threshold: E -> 0, growth rate -> 0, v_g -> infinity
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), t.energy2(3.0), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), t.growthRate(3.0), 1e-15);
}

test "superluminality is generic: v_g > 1 for every |p| > mu" {
    const t = try Tachyon.init(4.0); // mu = 2
    var prng = std.Random.DefaultPrng.init(77);
    const rnd = prng.random();
    for (0..2000) |_| {
        const p = 2.0 + rnd.float(f64) * 50.0;
        const vg = try t.groupVelocity(p);
        try std.testing.expect(vg > 1.0);
    }
}

test "tachyon four-momentum is spatial (class spatial)" {
    const t = try Tachyon.init(9.0);
    const s = sig.minkowski_3_1;
    const f = form.DiagonalForm{ .signature = s };

    // propagating mode: E = 4, p = 5
    const E = try t.energy(5.0);
    // (+,−,−,−): g = 16 - 25 = -9 = m^2
    try std.testing.expectApproxEqAbs(-t.mu2, t.fourMomentumNorm2(E, 5.0, .mostly_minus), 1e-12);
    // (−,+,+,+): g = -16 + 25 = +9 = +mu^2 — different value, same class
    try std.testing.expectApproxEqAbs(t.mu2, t.fourMomentumNorm2(E, 5.0, .mostly_plus), 1e-12);

    // the class is convention invariant: always spatial
    try std.testing.expectEqual(form.Class.spatial, t.momentumClass(5.0));
    const four = [_]f64{ E, 5.0, 0.0, 0.0 };
    try std.testing.expectEqual(form.Class.spatial, f.classify(&four));

    // the unstable mode is spatial too, but with imaginary E
    try std.testing.expectEqual(form.Class.spatial, t.momentumClass(2.0));
}

test "composite function: dual derivative matches analytic, FD lags" {
    const acc = scanPhaseAccuracy(9.0, 4.0, 100.0, 512, 1e-5);
    try std.testing.expect(acc.max_ad_error < 1e-11);
    try std.testing.expect(acc.max_fd_error > 1e-9);
    try std.testing.expect(acc.max_ad_error * 100.0 < acc.max_fd_error);
}

test "dual derivative of E reproduces the analytic formula bit for bit" {
    // This is an exactness claim, and the exact comparison IS the assertion.
    // For this expression the chain rule composes the very same floating point
    // operations as the closed form p/E, so the two are identical bit for bit.
    // The value of MRS here is not a "better result than the analytic formula"
    // (impossible) but that the formula need not be derived or typed by hand.
    const mu2 = 9.0;
    for ([_]f64{ 4.0, 5.0, 10.0, 50.0, 1e6 }) |p| {
        try std.testing.expectEqual(energyDerivExact(p, mu2), energyDerivAD(p, mu2));
    }
}

test "Newton with dual derivative hits the analytic solution" {
    const mu2 = 9.0;
    const k = 3.0;
    const p_star = superluminalMomentum(mu2, k);
    try std.testing.expectApproxEqAbs(9.0 / @sqrt(@as(f64, 8.0)), p_star, 1e-12);

    // start just above the threshold; convergence is quadratic, so a few steps
    const p0 = @sqrt(mu2) * 1.01;
    const res = newtonAD(mu2, k, p0, 1e-12, 50);
    try std.testing.expect(res.converged);
    try std.testing.expect(res.iterations <= 20);
    try std.testing.expectApproxEqAbs(p_star, res.p, 1e-11);
    try std.testing.expectApproxEqAbs(k, groupVelocityAt(res.p, mu2), 1e-10);

    // Newton started on the wrong side of the asymptote runs onto the negative
    // branch — a property of v_g(p) = p/E(p), which is odd and has a pole at
    // p = mu
    const bad = newtonAD(mu2, k, p_star * 1.5, 1e-12, 50);
    try std.testing.expect(!bad.converged or bad.p < 0.0);
}
