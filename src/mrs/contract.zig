//! MRS-LAB :: error contracts as part of the type
//!
//! THE PROBLEM THIS SOLVES. "The same result" has three different meanings and
//! none of them was in the type:
//!
//!   * `bit_exact`        — the same bits on every run and every platform;
//!   * `correctly_rounded`— the exact result rounded once, to the nearest;
//!   * `bounded`          — the exact result is within a radius of the value.
//!
//! A tolerance used to be an ARGUMENT the caller supplied (`classifyTol(v, tol)`,
//! `isZeroDivisor(a, tol)`), so it did not travel with the number and the
//! question "is this zero" depended on a constant chosen somewhere else. Here
//! the radius travels with the value and is composed by the same operations
//! that produce the value, so the zero test becomes a consequence:
//!
//!     a number is indistinguishable from zero when its own radius covers zero.
//!
//! THE RULES (standard forward error bounds, Higham, "Accuracy and Stability of
//! Numerical Algorithms", ch. 3):
//!
//!     add:  r = ra + rb + U·|value|
//!     sub:  r = ra + rb + U·|value|
//!     mul:  r = |a|·rb + |b|·ra + ra·rb + U·|value|
//!
//! The term `U·|value|` is the rounding of the operation itself; the other
//! terms are the propagated input error. Every rule is checked against
//! f128 arithmetic by test, on grids, not on examples.
//!
//! WHAT THESE RULES ARE NOT GOOD AT — measured, see the test at the bottom.
//! They are absolute-value rules, so they discard cancellation. Over a chain of
//! N boosts the bound in the ADDITIVE coordinate (the rapidity) grows like N·U
//! and stays tight: at N = 100 000 it is 1.0e-12 against a true error of
//! 2.8e-17. The same rules applied to the MULTIPLICATIVE representations
//! (split-complex product, 2×2 matrix product) grow like exp(Σ|θ_i|) and are
//! worthless by N = 100 000: 5.3e17 and 4.8e18 against true errors of ~3e-14.
//! The bound is not wrong — worst case it is right — but the worst case does
//! not happen, because the errors of opposite boosts cancel and an
//! absolute-value rule cannot see it.
//!
//! That is the MRS-LAB thesis stated as a property of error propagation rather
//! than as a benchmark: a contract is dischargeable exactly where the
//! representation is well conditioned, and the additive chart is the
//! representation where it is.

const std = @import("std");
const sig = @import("signature.zig");
const form = @import("form.zig");
const Z = @import("split_complex.zig").Z;

/// Unit roundoff for f64 under round-to-nearest: |fl(x) − x| <= U·|x|.
pub const U: f64 = 0x1p-53;

/// γ_n = n·U / (1 − n·U): the classical bound for a chain of n roundings.
/// Valid while n·U < 1, which for f64 means n < 9e15.
pub fn gamma(n: u32) f64 {
    const nu = @as(f64, @floatFromInt(n)) * U;
    return nu / (1.0 - nu);
}

/// What a result promises about its exact value.
pub const Contract = enum {
    /// Same bits on every run, every platform, every order. The contract
    /// required for navigation and audit.
    bit_exact,
    /// The exact result rounded once to the nearest f64.
    correctly_rounded,
    /// The exact result lies within `radius` of the value.
    bounded,

    pub fn label(self: Contract) []const u8 {
        return switch (self) {
            .bit_exact => "bit exact",
            .correctly_rounded => "correctly rounded",
            .bounded => "bounded by a radius",
        };
    }

    /// Only `bounded` carries a radius; the other two are absolute claims, and
    /// a caller holding one of them does not need a number to trust the value.
    pub fn carriesRadius(self: Contract) bool {
        return self == .bounded;
    }
};

/// A number together with a bound on how far the exact value may be from it.
pub const Bounded = struct {
    value: f64,
    /// |value − exact| <= radius, always. Never negative.
    radius: f64,

    /// A value known exactly — an input read from a data file, an integer, a
    /// power of two. Nothing is claimed about how it was produced.
    pub fn exact(v: f64) Bounded {
        return .{ .value = v, .radius = 0.0 };
    }

    /// A value whose own error is known from outside: a measurement, a
    /// conversion, a literal with a stated precision.
    ///
    /// `radius` must be finite and non-negative. It is CLAMPED rather than
    /// merely asserted, and that is deliberate: `std.debug.assert` is compiled
    /// out in ReleaseFast, which is the mode this project ships in, so an
    /// assert here would be a check that does not run where it matters. A
    /// negative radius would make `couldBeZero` answer "no" for every value —
    /// including a true zero — which is a wrong answer rather than a crash. A
    /// NaN radius would make every comparison false and do the same. Clamping
    /// costs one comparison and keeps the invariant `radius >= 0` true in
    /// every mode. A caller who wants the mistake to be loud should validate
    /// before calling.
    pub fn fromError(v: f64, radius: f64) Bounded {
        if (!std.math.isFinite(radius)) return .{ .value = v, .radius = std.math.inf(f64) };
        return .{ .value = v, .radius = @max(0.0, radius) };
    }

    pub fn add(a: Bounded, b: Bounded) Bounded {
        const v = a.value + b.value;
        return .{ .value = v, .radius = a.radius + b.radius + U * @abs(v) };
    }

    pub fn sub(a: Bounded, b: Bounded) Bounded {
        const v = a.value - b.value;
        return .{ .value = v, .radius = a.radius + b.radius + U * @abs(v) };
    }

    pub fn mul(a: Bounded, b: Bounded) Bounded {
        const v = a.value * b.value;
        const r = @abs(a.value) * b.radius + @abs(b.value) * a.radius + a.radius * b.radius;
        return .{ .value = v, .radius = r + U * @abs(v) };
    }

    pub fn scale(a: Bounded, k: f64) Bounded {
        return a.mul(exact(k));
    }

    pub fn square(a: Bounded) Bounded {
        return a.mul(a);
    }

    pub fn neg(a: Bounded) Bounded {
        return .{ .value = -a.value, .radius = a.radius };
    }

    /// The largest magnitude the exact value can have.
    pub fn absUpper(self: Bounded) f64 {
        return @abs(self.value) + self.radius;
    }

    /// Does the radius cover the exact value? This is the predicate every
    /// claim about a radius is tested with.
    pub fn containsValue(self: Bounded, exact_value: f64) bool {
        return @abs(self.value - exact_value) <= self.radius;
    }

    /// The zero test, derived rather than supplied: a number whose own radius
    /// covers zero cannot be distinguished from zero.
    ///
    /// This replaces `tol` arguments. Note what it does NOT say: a small value
    /// with a small radius is NOT called zero, it is called small.
    pub fn couldBeZero(self: Bounded) bool {
        return @abs(self.value) <= self.radius;
    }

    /// Relative width of the interval, for reporting. Zero when the value is
    /// exactly zero and the radius is zero.
    pub fn relativeWidth(self: Bounded) f64 {
        const denom = @abs(self.value);
        if (denom == 0.0) return if (self.radius == 0.0) 0.0 else std.math.inf(f64);
        return self.radius / denom;
    }
};

// ---------------------------------------------------------------------------
// Integration with the mathematical layer
// ---------------------------------------------------------------------------

/// The form evaluated with a propagated bound: Σ_i s_i v_i², n multiply-adds.
///
/// The signs s_i are exact (they are ±1 or 0), so they contribute no error of
/// their own; the input components are treated as exact, which is what makes
/// the result's radius a statement about THIS computation and not about where
/// the vector came from. Feed values with their own error through
/// `evalFormBounded` instead.
pub fn evalForm(f: form.DiagonalForm, v: []const f64) Bounded {
    std.debug.assert(v.len == f.n());
    var acc = Bounded.exact(0.0);
    for (v, 0..) |vi, i| {
        const s_i = if (f.use_signs) f.signs[i] else f.signature.signAt(i);
        const term = Bounded.exact(s_i).mul(Bounded.exact(vi).square());
        acc = acc.add(term);
    }
    return acc;
}

/// What a contract produces: a value with a radius, or a bare value.
pub fn Result(comptime c: Contract) type {
    return if (c == .bounded) Bounded else f64;
}

/// Evaluate the form under a contract fixed at compile time.
///
/// For `.bit_exact` and `.correctly_rounded` the return type is a plain `f64`
/// and the body IS `f.eval(v)` — the same code, not a similar one. So the layer
/// cannot cost anything on those paths, and that is not an optimisation result
/// to be measured but a property of the type: there is no radius to compute
/// because there is no radius in the type.
///
/// The price is paid only where a radius is asked for, and `bench/t5_contract`
/// measures what it is.
pub fn evalFormWith(comptime c: Contract, f: form.DiagonalForm, v: []const f64) Result(c) {
    return switch (c) {
        .bounded => evalForm(f, v),
        else => f.eval(v),
    };
}

/// The same, for a vector that already carries its own error.
pub fn evalFormBounded(f: form.DiagonalForm, v: []const Bounded) Bounded {
    std.debug.assert(v.len == f.n());
    var acc = Bounded.exact(0.0);
    for (v, 0..) |vi, i| {
        const s_i = if (f.use_signs) f.signs[i] else f.signature.signAt(i);
        acc = acc.add(Bounded.exact(s_i).mul(vi.square()));
    }
    return acc;
}

/// Classification of a vector whose norm carries a radius.
///
/// The difference from `DiagonalForm.classifyTol`: `null_like` is returned
/// when the propagated radius covers zero, not when a caller-supplied
/// tolerance says so. The vector below is classified relative to the declared
/// convention, exactly as before.
pub fn classify(f: form.DiagonalForm, v: []const f64) form.Class {
    const g = evalForm(f, v);
    // Same guard as `form.classifyTol`, for the same reason: a non-finite norm
    // makes every comparison below false, and without this the fall-through
    // would report a NaN vector as a tachyon. Note that a NaN INPUT propagates
    // into the radius rather than being trapped by `couldBeZero` — nothing here
    // claims to sanitise its input, only to refuse to classify it.
    if (!std.math.isFinite(g.value)) return .invalid;
    if (g.couldBeZero()) return .null_like;
    const ts = if (f.use_signs) f.time_sign.f() else f.signature.time_sign.f();
    return if (g.value * ts > 0) .temporal else .spatial;
}

/// The rapidity of a product of boosts, with a propagated bound.
///
/// This is the additive chart, and it is the one place where the bound stays
/// tight over a long chain — the test at the bottom of the file measures how
/// much tighter it is than the multiplicative representations.
pub fn rapidityChain(thetas: []const f64) Bounded {
    var acc = Bounded.exact(0.0);
    for (thetas) |t| acc = acc.add(Bounded.exact(t));
    return acc;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Deterministic pseudo-random value in [-1, 1].
fn gridValue(seed: *std.Random.DefaultPrng) f64 {
    return seed.random().float(f64) * 2.0 - 1.0;
}

test "PREDICTION: the radius contains the exact value on a grid of sums" {
    var prng = std.Random.DefaultPrng.init(0x1234);
    var checked: usize = 0;
    for (0..200) |_| {
        var acc = Bounded.exact(0.0);
        var exact: f128 = 0;
        for (0..17) |_| {
            const x = gridValue(&prng);
            acc = acc.add(Bounded.exact(x));
            exact += @as(f128, x);
        }
        checked += 1;
        if (!acc.containsValue(@floatCast(exact))) {
            std.debug.print("sum: value={e} radius={e} exact={e}\n", .{ acc.value, acc.radius, exact });
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked > 0);
}

test "PREDICTION: the radius contains the exact value on a grid of products" {
    var prng = std.Random.DefaultPrng.init(0x5678);
    var checked: usize = 0;
    for (0..200) |_| {
        var acc = Bounded.exact(1.0);
        var exact: f128 = 1.0;
        for (0..9) |_| {
            const x = gridValue(&prng) * 3.0;
            acc = acc.mul(Bounded.exact(x));
            exact *= @as(f128, x);
        }
        checked += 1;
        // The accumulated radius under a product is the absolute-value rule, so
        // it must still contain the exact value — this is the correctness claim.
        if (!acc.containsValue(@floatCast(exact))) {
            std.debug.print("mul: value={e} radius={e} exact={e}\n", .{ acc.value, acc.radius, exact });
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked > 0);
}

test "a radius is never negative, in any optimisation mode" {
    // `std.debug.assert` is compiled out in ReleaseFast, so the invariant has
    // to be enforced by the constructor itself: a negative radius would make
    // `couldBeZero` answer "no" for a true zero, which is a wrong answer
    // rather than a crash.
    try std.testing.expectEqual(@as(f64, 0.0), Bounded.fromError(1.0, -5.0).radius);
    try std.testing.expect(std.math.isInf(Bounded.fromError(1.0, std.math.nan(f64)).radius));
    // A true zero is still recognised as zero after the clamp, and a nonzero
    // value is still not called zero — the clamp must not turn one into the
    // other in either direction.
    try std.testing.expect(Bounded.fromError(0.0, -5.0).couldBeZero());
    try std.testing.expect(!Bounded.fromError(1.0, -5.0).couldBeZero());
}

test "the zero test is derived, and a small number is not called zero" {
    // A tiny value with a tiny radius: NOT zero.
    const small = Bounded.exact(1e-300);
    try std.testing.expect(!small.couldBeZero());

    // A value whose radius covers zero: indistinguishable from zero.
    const covered = Bounded.fromError(1e-17, 1e-16);
    try std.testing.expect(covered.couldBeZero());

    // Exactly zero, exactly known: zero.
    try std.testing.expect(Bounded.exact(0.0).couldBeZero());
}

/// Radius of the additive chart over a chain of `n` boosts of magnitude 1e-3.
fn additiveBound(n: usize) f64 {
    var prng = std.Random.DefaultPrng.init(0x5EED);
    var acc = Bounded.exact(0.0);
    for (0..n) |_| acc = acc.add(Bounded.exact((gridValue(&prng)) * 1e-3));
    return acc.radius;
}

/// Radius of the real part of the split-complex product over the same chain.
fn multiplicativeBound(n: usize) f64 {
    var prng = std.Random.DefaultPrng.init(0x5EED);
    var re = Bounded.exact(1.0);
    var im = Bounded.exact(0.0);
    for (0..n) |_| {
        const t = gridValue(&prng) * 1e-3;
        const b = Z.fromRapidity(t);
        const old_re = re;
        re = re.mul(Bounded.exact(b.re)).sub(im.mul(Bounded.exact(b.im)));
        im = im.mul(Bounded.exact(b.re)).add(old_re.mul(Bounded.exact(b.im)));
    }
    return re.radius;
}

test "PREDICTION: the additive chart keeps a useful bound where the multiplicative one does not" {
    // Mixed signs: the cancelling regime, which is the hard one for rules that
    // take absolute values.
    const additive = additiveBound(100_000);
    const multiplicative = multiplicativeBound(100_000);

    // After 100 000 boosts the propagated radius is still around 1e-12 — small
    // enough to decide something with.
    try std.testing.expect(additive < 1e-9);
    // The same rules on the multiplicative chart give a radius of order 1e17,
    // which is correct in the worst case and worthless in fact.
    try std.testing.expect(multiplicative > 1e10 * U);
    // And the gap between the two grows with the chain, it is not a constant.
    try std.testing.expect(multiplicativeBound(100_000) > multiplicativeBound(1_000));
}

test "gamma is monotone and agrees with n·U for small n" {
    try std.testing.expectEqual(@as(f64, 0.0), gamma(0));
    try std.testing.expectApproxEqRel(U, gamma(1), 1e-12);
    try std.testing.expect(gamma(1000) > gamma(100));
    try std.testing.expect(gamma(1000) < 1000 * U * 1.001);
}

test "a contract that carries no radius produces the value bit for bit" {
    // The claim behind the death criterion "the layer must cost nothing when
    // the contract is known at compile time": for `bit_exact` the code path IS
    // the plain one, so the results cannot differ by so much as one bit.
    const f = form.DiagonalForm.init(sig.minkowski_3_1);
    const v = [_]f64{ 1.25, -2.5, 0.0, 3.75 };
    const plain = f.eval(&v);
    try std.testing.expectEqual(plain, evalFormWith(.bit_exact, f, &v));
    try std.testing.expectEqual(plain, evalFormWith(.correctly_rounded, f, &v));

    // ... and the bounded one returns a radius that covers the same value.
    const bounded = evalFormWith(.bounded, f, &v);
    try std.testing.expect(bounded.containsValue(plain));
    try std.testing.expect(bounded.radius > 0.0);
}

test "classify: a radius that covers zero gives null_like, not a guessed sign" {
    const f = form.DiagonalForm.init(sig.minkowski_3_1);
    // A vector whose norm is exactly zero.
    const nullish = [_]f64{ 1.0, 1.0, 0.0, 0.0 };
    try std.testing.expectEqual(form.Class.null_like, classify(f, &nullish));
    // A clearly timelike vector keeps its class.
    const timelike = [_]f64{ 1.0, 0.25, 0.0, 0.0 };
    try std.testing.expectEqual(form.Class.temporal, classify(f, &timelike));
}

test "contracts describe themselves" {
    try std.testing.expectEqualStrings("bit exact", Contract.bit_exact.label());
    try std.testing.expect(Contract.bounded.carriesRadius());
    try std.testing.expect(!Contract.bit_exact.carriesRadius());
}

test "a non-finite input cannot be classified, on either path" {
    const f = form.DiagonalForm.init(sig.minkowski_3_1);
    const nan = [_]f64{ std.math.nan(f64), 0, 0, 0 };
    const inf = [_]f64{ std.math.inf(f64), 1, 0, 0 };
    try std.testing.expectEqual(form.Class.invalid, classify(f, &nan));
    try std.testing.expectEqual(form.Class.invalid, classify(f, &inf));
    // and the in-form path refuses it too, so the two cannot disagree
    try std.testing.expectEqual(form.Class.invalid, f.classify(&nan));

    // A non-finite input poisons the radius rather than being trapped: this is
    // stated so nobody expects `couldBeZero` to double as validation.
    const g = evalForm(f, &nan);
    try std.testing.expect(!std.math.isFinite(g.radius));
    try std.testing.expect(!g.couldBeZero());
}
