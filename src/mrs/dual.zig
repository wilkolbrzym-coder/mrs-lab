//! MRS-LAB :: dual numbers (ε² = 0)
//!
//! An element a + ε·b with ε² = 0. This is the smallest algebra in which a
//! derivative can be written down: ε is infinitesimal and ε² vanishes.
//!
//! The practical consequence MRS-LAB exploits: **the exact derivative in a
//! single pass**. If f: R → R is extended to f: D → D (keeping the same
//! operations), then f(a + ε) = f(a) + ε·f'(a). That gives the derivative to
//! within rounding error, with no step size and no subtraction of nearby
//! numbers — so no O(h²) truncation error and no catastrophic cancellation.
//!
//! This is the same mechanism commercial libraries call forward-mode automatic
//! differentiation. MRS-LAB does not reinvent it; it makes it an elementary
//! operation of the description language.

const std = @import("std");

pub const D = struct {
    val: f64 = 0,
    der: f64 = 0,

    /// Independent variable: (x, 1).
    pub fn at(x: f64) D {
        return .{ .val = x, .der = 1.0 };
    }

    /// Constant: (c, 0).
    pub fn cst(c: f64) D {
        return .{ .val = c, .der = 0.0 };
    }

    pub fn add(a: D, b: D) D {
        return .{ .val = a.val + b.val, .der = a.der + b.der };
    }

    pub fn sub(a: D, b: D) D {
        return .{ .val = a.val - b.val, .der = a.der - b.der };
    }

    /// (a + εa')(b + εb') = ab + ε(ab' + a'b)
    pub fn mul(a: D, b: D) D {
        return .{ .val = a.val * b.val, .der = a.val * b.der + a.der * b.val };
    }

    pub fn div(a: D, b: D) D {
        return .{
            .val = a.val / b.val,
            .der = (a.der * b.val - a.val * b.der) / (b.val * b.val),
        };
    }

    pub fn scale(a: D, k: f64) D {
        return .{ .val = a.val * k, .der = a.der * k };
    }

    pub fn neg(a: D) D {
        return .{ .val = -a.val, .der = -a.der };
    }

    /// √a = √val + ε·a'/(2√val)
    pub fn sqrt(a: D) D {
        const r = @sqrt(a.val);
        return .{ .val = r, .der = a.der / (2.0 * r) };
    }

    /// sin/cos by the chain rule — needed for wave operators.
    pub fn sin(a: D) D {
        return .{ .val = @sin(a.val), .der = a.der * @cos(a.val) };
    }

    pub fn cos(a: D) D {
        return .{ .val = @cos(a.val), .der = -a.der * @sin(a.val) };
    }

    /// exp: (exp val, a'·exp val)
    pub fn exp(a: D) D {
        const e = std.math.exp(a.val);
        return .{ .val = e, .der = a.der * e };
    }
};

/// Derivative of a one-variable function given by a function pointer.
pub fn derivOf(f: *const fn (D) D, x: f64) f64 {
    return f(D.at(x)).der;
}

/// Derivative at x for a function that also needs a constant μ.
pub fn derivOfParam(f: *const fn (D, f64) D, x: f64, param: f64) f64 {
    return f(D.at(x), param).der;
}

/// Central difference as the reference baseline — this is what "textbook"
/// mathematics does when it has no dual algebra.
pub fn centralDiff(f: *const fn (f64, f64) f64, x: f64, param: f64, h: f64) f64 {
    return (f(x + h, param) - f(x - h, param)) / (2.0 * h);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

fn squareFn(x: D) D {
    return D.mul(x, x);
}

test "eps squared is zero" {
    const e = D.at(0.0); // (0,1) = eps
    const e2 = D.mul(e, e);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), e2.val, 0.0);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), e2.der, 0.0);
}

test "exact derivative in a single pass" {
    // d/dx x² = 2x
    try std.testing.expectApproxEqAbs(@as(f64, 6.0), derivOf(squareFn, 3.0), 1e-15);
    try std.testing.expectApproxEqAbs(@as(f64, -1.0), derivOf(squareFn, -0.5), 1e-15);
}

test "dual algebra versus floating point arithmetic" {
    // (a + εa')·(b + εb'): eps² = 0 makes the result an exact algebraic
    // identity, not an approximation.
    const a = D{ .val = 1e8, .der = 1.0 };
    const b = D{ .val = 1e-8, .der = 1.0 };
    const p = D.mul(a, b);
    try std.testing.expectApproxEqRel(@as(f64, 1.0), p.val, 1e-15);
    try std.testing.expectApproxEqRel(a.val + b.val, p.der, 1e-15);
}
