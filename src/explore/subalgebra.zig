//! MRS-LAB :: P4 — subalgebras, ideals, centre
//!
//! QUESTION. Which subspaces spanned by blades are closed under the product,
//! which are unital, which are commutative, which are two-sided ideals of the
//!
//! WHY THIS IS EXACT. Blades are linearly independent, and the product of two
//! blades is — up to sign — a single blade
//! (e_A·e_B = ±e_{A△B}) or zero. Hence a subspace spanned by a set of blades B
//! is closed under the product **if and only if** for every pair of blades in B
//! their product (or zero) lies in the span of B. The check therefore reduces
//! to a mask membership test — no floating point arithmetic and no tolerance.
//! zmiennoprzecinkowej i bez tolerancji.
//!
//! MEASURED RESULTS (verified below, exhaustively for n <= 4):
//!   * the blades containing a degenerate generator span a proper, nonzero,
//!     two-sided ideal I (the radical ideal);
//!   * I is nilpotent with index EXACTLY r + 1. In particular I² = 0 only for
//!     r = 1 — for r >= 2 distinct degenerate generators multiply to something
//!     nonzero, e.g. ζ₁·ζ₂ ≠ 0.
//!   * for r = 0 there are no proper nonzero blade-spanned ideals.
//!
//! SCOPE — read carefully, this is a boundary and not a detail. We enumerate
//! only subspaces SPANNED BY BLADES. That is a sublattice of all subalgebras,
//! and not all subalgebras. Consequence: in split algebras (e.g.
//! Cl(1,0) ≅ R⊕R) there exist proper ideals spanned by idempotents (1 ± ω)/2
//! which this engine **cannot see**, because idempotents are not blades. The
//! blades. The engine does not hide this — it reports the result with its scope.
//!
//! Exhaustive for n <= 4 (2^16 = 65 536 subsets). For n = 5 there are 2^32
//! subsets, so enumeration refuses — explicitly, not silently.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const exact = @import("exact.zig");
const sigs = @import("signatures.zig");

/// Largest number of blades at which subset enumeration is feasible.
pub const MAX_EXHAUSTIVE_BASIS: usize = 16;

pub const EnumerateError = error{TooManyBlades};

pub const Info = struct {
    /// Subset of blades; bit i means the blade with mask i.
    set: u32 = 0,
    dim: u5 = 0,
    unital: bool = false,
    commutative: bool = false,
    /// Two-sided ideal of the WHOLE algebra.
    ideal: bool = false,
    proper: bool = false,
    /// I·I = 0 (the product vanishes entirely) — a property of a nilpotent ideal.
    square_zero: bool = false,

    pub fn writeTo(self: Info, alg: cl.Algebra, w: anytype) !void {
        try w.print("dim={d} {{", .{self.dim});
        var first = true;
        const m = alg.basisCount();
        for (0..m) |i| {
            if (self.set & (@as(u32, 1) << @intCast(i)) == 0) continue;
            if (!first) try w.writeAll(", ");
            try alg.writeBlade(w, @intCast(i));
            first = false;
        }
        if (first) try w.writeAll("0");
        try w.writeAll("}");
    }
};

/// Is the span of a set of blades closed under the product?
pub fn isClosed(alg: cl.Algebra, set: u32) bool {
    const m = alg.basisCount();
    for (0..m) |i| {
        if (set & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (set & (@as(u32, 1) << @intCast(j)) == 0) continue;
            const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue; //             if (bp.sign == 0) continue; // zero belongs to every subspace
            if (set & (@as(u32, 1) << @intCast(bp.mask)) == 0) return false;
        }
    }
    return true;
}

/// Is the subalgebra commutative (its blades commute pairwise)?
pub fn isCommutative(alg: cl.Algebra, set: u32) bool {
    const m = alg.basisCount();
    for (0..m) |i| {
        if (set & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (set & (@as(u32, 1) << @intCast(j)) == 0) continue;
            const ab = cl.bladeMul(alg, @intCast(i), @intCast(j));
            const ba = cl.bladeMul(alg, @intCast(j), @intCast(i));
            if (ab.mask != ba.mask or ab.sign != ba.sign) return false;
        }
    }
    return true;
}

/// Is the span of the set a two-sided ideal of the whole algebra?
pub fn isTwoSidedIdeal(alg: cl.Algebra, set: u32) bool {
    if (set == 0) return true; //     if (set == 0) return true; // the zero ideal
    const m = alg.basisCount();
    for (0..m) |b| {
        if (set & (@as(u32, 1) << @intCast(b)) == 0) continue;
        for (0..m) |a| { //         for (0..m) |a| { // over all blades of the whole algebra
            const ab = cl.bladeMul(alg, @intCast(a), @intCast(b));
            if (ab.sign != 0 and set & (@as(u32, 1) << @intCast(ab.mask)) == 0) return false;
            const ba = cl.bladeMul(alg, @intCast(b), @intCast(a));
            if (ba.sign != 0 and set & (@as(u32, 1) << @intCast(ba.mask)) == 0) return false;
        }
    }
    return true;
}

/// Is I·I = 0 (checked on blades — sufficient, because the product is bilinear).
pub fn isSquareZero(alg: cl.Algebra, set: u32) bool {
    const m = alg.basisCount();
    for (0..m) |i| {
        if (set & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (set & (@as(u32, 1) << @intCast(j)) == 0) continue;
            if (cl.bladeMul(alg, @intCast(i), @intCast(j)).sign != 0) return false;
        }
    }
    return true;
}

/// Collects information about a set of blades.
pub fn info(alg: cl.Algebra, set: u32) Info {
    const m = alg.basisCount();
    return .{
        .set = set,
        .dim = @intCast(@popCount(set)),
        .unital = (set & 1) != 0,
        .commutative = if (isClosed(alg, set)) isCommutative(alg, set) else false,
        .ideal = isTwoSidedIdeal(alg, set),
        .proper = @popCount(set) < m,
        .square_zero = isSquareZero(alg, set),
    };
}

/// Exhaustive enumeration of closed subsets of blades.
/// Writes into `out` and returns the number found. `out` must be large
/// enough; the remainder is silently dropped, so use `countClosed` for
/// diagnostics.
pub fn enumerateClosed(alloc: std.mem.Allocator, alg: cl.Algebra) ![]Info {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;

    var list = std.ArrayList(Info).empty;
    errdefer list.deinit(alloc);

    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 0;
    while (s < total) : (s += 1) {
        const set: u32 = @intCast(s);
        if (!isClosed(alg, set)) continue;
        try list.append(alloc, info(alg, set));
    }
    return list.toOwnedSlice(alloc);
}

/// Number of closed subsets (without allocating the result).
pub fn countClosed(alg: cl.Algebra) EnumerateError!usize {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;
    var count: usize = 0;
    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 0;
    while (s < total) : (s += 1) {
        if (isClosed(alg, @intCast(s))) count += 1;
    }
    return count;
}

/// Number of proper, nonzero, two-sided blade-spanned ideals.
pub fn countProperIdeals(alg: cl.Algebra) EnumerateError!usize {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;
    var count: usize = 0;
    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 1; //     var s: u64 = 1; // skip the empty set (the zero ideal)
    while (s < total) : (s += 1) {
        const set: u32 = @intCast(s);
        const dim: u5 = @intCast(@popCount(set));
        if (dim == 0 or dim == m) continue; // trywialne
        if (isTwoSidedIdeal(alg, set)) count += 1;
    }
    return count;
}

// ---------------------------------------------------------------------------
// Ideal powers and nilpotency
// ---------------------------------------------------------------------------

/// Product of two subspaces spanned by blades.
///
/// STRUCTURAL FACT this rests on: the product of two blades is, up to sign, a
/// single blade. Therefore the product of blade-spanned subspaces is AGAIN
/// blade-spanned, and no linear algebra over a field is needed — a set of masks
/// is enough.
pub fn mulSets(alg: cl.Algebra, a: u32, b: u32) u32 {
    const m = alg.basisCount();
    var out: u32 = 0;
    for (0..m) |i| {
        if (a & (@as(u32, 1) << @intCast(i)) == 0) continue;
        for (0..m) |j| {
            if (b & (@as(u32, 1) << @intCast(j)) == 0) continue;
            const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue;
            out |= (@as(u32, 1) << @intCast(bp.mask));
        }
    }
    return out;
}

/// Nilpotency index: the smallest k >= 1 with I^k = {0}.
/// Returns null when I^k ≠ 0 for `max_steps` steps.
pub fn nilpotencyIndex(alg: cl.Algebra, set: u32, max_steps: u32) ?u32 {
    if (set == 0) return 1;
    var cur = set;
    var k: u32 = 1;
    while (k <= max_steps) : (k += 1) {
        cur = mulSets(alg, cur, set);
        if (cur == 0) return k + 1;
    }
    return null;
}

// ---------------------------------------------------------------------------
// Tests — they check PREDICTIONS, not merely that the code runs
// ---------------------------------------------------------------------------

test "even blades form a unital subalgebra of dimension 2^(n-1)" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        if (buf[i].n() > 4) continue;
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const ev = exact.evenBlades(alg);
        try std.testing.expect(isClosed(alg, ev));
        try std.testing.expect((ev & 1) != 0); // unitarna: zawiera skalar 1
        try std.testing.expectEqual(@as(u5, @intCast(alg.basisCount() / 2)), @as(u5, @intCast(@popCount(ev))));
    }
}

test "the centre is always a subalgebra" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        if (buf[i].n() > 4) continue;
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const ctr = exact.centerBasis(alg);
        try std.testing.expect(isClosed(alg, ctr));
        try std.testing.expect((ctr & 1) != 0);
    }
}

test "PREDICTION: the radical ideal is nilpotent with index EXACTLY r+1" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (!t.isDegenerate() or t.n() > 4) continue;
        checked += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const I = exact.bladesWithDegenerate(alg);

        const closed = isClosed(alg, I);
        const ideal = isTwoSidedIdeal(alg, I);
        const proper = @as(usize, @popCount(I)) < alg.basisCount();
        const idx = nilpotencyIndex(alg, I, 8);

        if (!closed or !ideal or !proper or idx == null or idx.? != @as(u32, t.r) + 1) {
            std.debug.print(
                "radical ideal: t=({d},{d},{d}) dim={d}/{d} closed={} ideal={} proper={} idx={?}\n",
                .{ t.p, t.q, t.r, @popCount(I), alg.basisCount(), closed, ideal, proper, idx },
            );
            return error.TestUnexpectedResult;
        }

        // I² = 0 holds EXACTLY for r = 1 — for r >= 2 distinct degenerate
        // generators multiply to something nonzero (e.g. ζ1·ζ2 ≠ 0).
        try std.testing.expectEqual(t.r == 1, isSquareZero(alg, I));

        const cnt = try countProperIdeals(alg);
        if (cnt < 1) {
            std.debug.print("no proper ideal: t=({d},{d},{d})\n", .{ t.p, t.q, t.r });
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked > 0);
}

test "PREDICTION: r = 0 gives no proper blade-spanned ideals" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.isDegenerate() or t.n() > 4) continue;
        checked += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        // Note: this is a theorem about the BLADE scope. Split cases
        // (e.g. Cl(1,0) ≅ R⊕R) have ideals spanned by idempotents, which this engine
        // cannot see — and that is a deliberate limitation.
        try std.testing.expectEqual(@as(usize, 0), try countProperIdeals(alg));
    }
    try std.testing.expect(checked > 0);
}

test "the trivial sets are always closed: {0} and the whole algebra" {
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus)).algebra();
    try std.testing.expect(isClosed(alg, 0));
    const full: u32 = @intCast((@as(u64, 1) << @intCast(alg.basisCount())) - 1);
    try std.testing.expect(isClosed(alg, full));
    try std.testing.expect(isTwoSidedIdeal(alg, full));
    // the whole algebra is not proper
    try std.testing.expect(!info(alg, full).proper);
}

test "enumeration refuses for n = 5 instead of silently running for half an hour" {
    const alg = try (sigs.SigBuf.build(.{ .p = 0, .q = 5 }, .mostly_minus)).algebra();
    try std.testing.expectEqual(@as(usize, 32), alg.basisCount());
    try std.testing.expectError(error.TooManyBlades, countClosed(alg));
    try std.testing.expectError(error.TooManyBlades, countProperIdeals(alg));
}

test "enumeration is deterministic and agrees with the counter" {
    const alloc = std.testing.allocator;
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_minus)).algebra();
    const a = try enumerateClosed(alloc, alg);
    defer alloc.free(a);
    const b = try enumerateClosed(alloc, alg);
    defer alloc.free(b);
    try std.testing.expectEqual(a.len, b.len);
    try std.testing.expectEqual(try countClosed(alg), a.len);
    for (a, b) |x, y| {
        try std.testing.expectEqual(x.set, y.set);
        try std.testing.expectEqual(x.dim, y.dim);
        try std.testing.expectEqual(x.ideal, y.ideal);
    }
    // results are independent of the sign convention (an algebraic property)
    const alg_plus = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_plus)).algebra();
    try std.testing.expectEqual(try countClosed(alg), try countClosed(alg_plus));
}
