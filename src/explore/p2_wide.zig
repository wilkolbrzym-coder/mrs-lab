//! MRS-LAB :: P2a BEYOND the range where the rule was observed
//!
//! THE PREDICTION, DECLARED BEFORE THE MEASUREMENT — and this file is written
//! so that the declaration comes first, in source order and in the test order.
//!
//! The table reports the rule **P2a holds ⇔ p + q <= 2**, verified for every
//! signature with `p+q+r <= 5` (`zig build explore`). That is an observation
//! INSIDE a range. Nothing about it was ever checked beyond `n = 5`, and the
//! engine refused there for an implementation reason (32-bit masks), not
//! because the mathematics stops.
//!
//! So the rule makes a sharp prediction at `n = 6`, which is declared here:
//!
//!     (0,0,6)   p+q = 0   predicted: holds   — six nilpotent generators are
//!                                             invisible to the scalar norm
//!     (0,1,5)   p+q = 1   predicted: holds
//!     (0,2,4)   p+q = 2   predicted: holds
//!     (0,3,3)   p+q = 3   predicted: FAILS
//!     (0,6,0)   p+q = 6   predicted: FAILS
//!
//! Both outcomes are results, and neither is a failure of the experiment:
//!
//!   * if every line holds, the rule survived a factor-16 enlargement of the
//!     space in which it was found, and is no longer a fit to 55 points;
//!   * if any line fails, we have located the BOUNDARY of our own rule — new
//!     knowledge, of exactly the kind a laboratory is for.
//!
//! WHY THIS IS NOT A SECOND ENGINE. It is the same arithmetic and the same
//! determining set D = {0} ∪ {e_i} ∪ {e_i+e_j}; what differs is only the width
//! of a blade mask (64 instead of 32) and the fact that D is generated on the
//! fly instead of materialised, because at n = 6 the dense grid would be
//! 2081 × 64 coefficients. The conjugate signs come from `exact.zig`, the same
//! function the in-range path uses, and a cross-check test requires the two to
//! agree on all 55 signatures with p+q+r <= 5 — where both can run — before any
//! n = 6 answer is allowed to mean anything.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const exact = @import("exact.zig");
const nrm = @import("norm_mult.zig");
const sigs = @import("signatures.zig");

/// Largest number of blades this path handles: 2^6 = 64.
pub const MAX_M: usize = 64;
/// |D| at MAX_M: 1 + 64 + 64·63/2.
pub const MAX_GRID_WIDE: usize = 1 + MAX_M + (MAX_M * (MAX_M - 1)) / 2;
/// D has at most two blades; a product of two such vectors has at most four.
const MAX_TERMS: usize = 4;

/// A multivector as a list of nonzero terms, with 64-bit-wide masks.
pub const Terms = struct {
    mask: [MAX_TERMS]u8 = undefined,
    val: [MAX_TERMS]i64 = undefined,
    len: usize = 0,

    pub fn eql(a: Terms, b: Terms) bool {
        if (a.len != b.len) return false;
        for (0..a.len) |i| {
            if (a.mask[i] != b.mask[i] or a.val[i] != b.val[i]) return false;
        }
        return true;
    }
};

/// `e_i · e_j` for every pair, once per algebra.
const Table = struct {
    mask: [MAX_M * MAX_M]u8 = undefined,
    sign: [MAX_M * MAX_M]i8 = undefined,
    m: usize = 0,

    /// Refuses an algebra wider than MAX_M rather than writing past the table.
    /// The same defect the 0.1.5 review found in `exact.ProductTable.init`, and
    /// the same fix: a guard that runs in EVERY mode, not a `std.debug.assert`
    /// that ReleaseFast compiles out.
    fn init(alg: cl.Algebra) error{TooManyBlades}!Table {
        const m = alg.basisCount();
        if (m > MAX_M) return error.TooManyBlades;
        var t = Table{ .m = m };
        for (0..m) |i| {
            for (0..m) |j| {
                const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
                t.mask[i * m + j] = @intCast(bp.mask);
                t.sign[i * m + j] = bp.sign;
            }
        }
        return t;
    }
};

inline fn conjSign(mask: u8) i64 {
    return @as(i64, exact.cliffordConjSign(@intCast(mask)));
}

/// Product of two sparse factors, merging terms that land on the same blade.
fn mulTerms(t: *const Table, a: Terms, b: Terms, out: *Terms) void {
    out.len = 0;
    for (0..a.len) |i| {
        for (0..b.len) |j| {
            const idx = @as(usize, a.mask[i]) * t.m + b.mask[j];
            const s = t.sign[idx];
            if (s == 0) continue;
            const mk = t.mask[idx];
            const term = a.val[i] * b.val[j] * @as(i64, s);
            var merged = false;
            for (0..out.len) |k| {
                if (out.mask[k] == mk) {
                    out.val[k] += term;
                    merged = true;
                    break;
                }
            }
            if (!merged) {
                out.mask[out.len] = mk;
                out.val[out.len] = term;
                out.len += 1;
            }
        }
    }
}

/// Scalar part of `a·conj(a)`, computed the way `exact.normScalarTerms` does:
/// conjugate the second factor, multiply, then apply the sign.
fn scalarNorm(t: *const Table, a: Terms) i64 {
    var out: i64 = 0;
    for (0..a.len) |i| {
        for (0..a.len) |j| {
            const idx = @as(usize, a.mask[i]) * t.m + a.mask[j];
            const s = t.sign[idx];
            if (s == 0) continue;
            if (t.mask[idx] != 0) continue; // the scalar part only
            out += a.val[i] * (a.val[j] * conjSign(a.mask[j])) * @as(i64, s);
        }
    }
    return out;
}

pub const Decision = struct {
    holds: bool = false,
    grid_size: usize = 0,
    /// The first counterexample, when one exists: x and y as blade terms.
    x: Terms = .{},
    y: Terms = .{},
    has_witness: bool = false,
};

/// The k-th element of D = {0} ∪ {e_i} ∪ {e_i+e_j}, in the SAME order the
/// in-range `buildGrid` uses, generated instead of stored.
fn gridAt(m: usize, idx: usize, out: *Terms) void {
    if (idx == 0) {
        out.len = 0;
        return;
    }
    var k = idx - 1;
    if (k < m) {
        out.mask[0] = @intCast(k);
        out.val[0] = 1;
        out.len = 1;
        return;
    }
    k -= m;
    var i: usize = 0;
    while (i < m) : (i += 1) {
        const cnt = m - i - 1;
        if (k < cnt) {
            out.mask[0] = @intCast(i);
            out.val[0] = 1;
            out.mask[1] = @intCast(i + 1 + k);
            out.val[1] = 1;
            out.len = 2;
            return;
        }
        k -= cnt;
    }
    unreachable; // idx < gridSize(m)
}

pub fn gridSize(m: usize) usize {
    return 1 + m + (m * (m - 1)) / 2;
}

/// Decides P2a — `N_sc(x·y) = N_sc(x)·N_sc(y)` for all real x, y — by the same
/// determining-set argument the in-range engine uses, with 64-bit-wide masks.
pub fn decideScalarMultiplicative(alg: cl.Algebra) error{TooManyBlades}!Decision {
    const t = try Table.init(alg);
    const m = alg.basisCount();
    const k = gridSize(m);

    // The same three wastes T6 removed from the in-range path are removed here,
    // because a claim about the exact mode that rests on a slower version of it
    // would be inconsistent with itself: D is generated ONCE and its norms are
    // computed ONCE, instead of k times each inside the inner loop.
    var d: [MAX_GRID_WIDE]Terms = undefined;
    var norms: [MAX_GRID_WIDE]i64 = undefined;
    for (0..k) |i| {
        gridAt(m, i, &d[i]);
        norms[i] = scalarNorm(&t, d[i]);
    }

    var p: Terms = undefined;
    var i: usize = 0;
    while (i < k) : (i += 1) {
        var j: usize = 0;
        while (j < k) : (j += 1) {
            mulTerms(&t, d[i], d[j], &p);
            // Checked multiply, like the in-range path: the exact layer's rule
            // is that a coefficient beyond i64 panics rather than wrapping, and
            // a raw `*` here would quietly break that rule.
            const expected = std.math.mul(i64, norms[i], norms[j]) catch
                @panic("MRS p2_wide: overflow in exact arithmetic");
            if (scalarNorm(&t, p) != expected) {
                return .{
                    .holds = false,
                    .grid_size = k,
                    .x = d[i],
                    .y = d[j],
                    .has_witness = true,
                };
            }
        }
    }
    return .{ .holds = true, .grid_size = k };
}

/// Terms for a triple, for reporting.
pub fn tripleTerms(m: usize, idx: usize) Terms {
    var out: Terms = .{};
    gridAt(m, idx, &out);
    return out;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "the wide path is the same decision at n <= 4, where both engines run" {
    // Without this, an n = 6 answer from the wide path would be a claim about
    // an untested program. It must reproduce the in-range verdict exactly —
    // including on the signatures where the verdict is NO, because agreeing
    // only on "yes" would be agreement on nothing.
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, 4);
    var grid: [nrm.MAX_GRID]exact.IntVec = undefined;
    var holds: usize = 0;
    var fails: usize = 0;
    for (0..n) |t| {
        const alg = try (sigs.SigBuf.build(buf[t], .mostly_minus)).algebra();
        const wide = try decideScalarMultiplicative(alg);
        const in_range = nrm.checkScalarMultiplicativeBrute(alg, &grid);

        if (wide.holds != in_range.holds) {
            std.debug.print("wide vs in-range disagree on ({d},{d},{d}): wide={} in-range={}\n", .{
                buf[t].p, buf[t].q, buf[t].r, wide.holds, in_range.holds,
            });
            return error.TestUnexpectedResult;
        }
        if (wide.holds) holds += 1 else fails += 1;
    }
    try std.testing.expectEqual(@as(usize, 34), n);
    try std.testing.expect(holds > 0);
    try std.testing.expect(fails > 0);
}

test "PREDICTION DECLARED IN THE HEADER: the rule p+q <= 2 at n = 6" {
    // The five lines of the declaration, in the order it was written. Only the
    // outcome is reported; a failure here is a RESULT, not a broken test, so
    // the message says which one it is.
    //
    // ReleaseFast only: the determining set at n = 6 has |D| = 2081 elements and
    // |D|² = 4 329 961 pairs, which is 5 seconds here and would be minutes in
    // Debug. Debug still exercises this whole path through the cross-check test
    // above, on every signature with n <= 4.
    if (@import("builtin").mode != .ReleaseFast) return error.SkipZigTest;
    const cases = [_]struct { p: u5, q: u5, r: u5, predicted: bool }{
        .{ .p = 0, .q = 0, .r = 6, .predicted = true },
        .{ .p = 0, .q = 1, .r = 5, .predicted = true },
        .{ .p = 0, .q = 2, .r = 4, .predicted = true },
        .{ .p = 0, .q = 3, .r = 3, .predicted = false },
        .{ .p = 0, .q = 6, .r = 0, .predicted = false },
    };
    for (cases) |c| {
        const alg = try (sigs.SigBuf.build(.{ .p = c.p, .q = c.q, .r = c.r }, .mostly_minus)).algebra();
        const d = try decideScalarMultiplicative(alg);
        if (d.holds != c.predicted) {
            std.debug.print(
                "RULE BOUNDARY FOUND: ({d},{d},{d}) p+q={d} predicted holds={} but the " ++
                    "determining set says holds={s}\n",
                .{ c.p, c.q, c.r, c.p + c.q, c.predicted, if (d.holds) "true" else "FALSE" },
            );
            return error.TestUnexpectedResult;
        }
    }
}
