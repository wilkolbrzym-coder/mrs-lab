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
//! engine does not hide this — it reports the result together with its scope.
//!
//! TWO ENUMERATION METHODS, and why the second one is not an approximation.
//!
//! `countsBrute` scans all 2^m blade subsets (m = 2^n). Exact, but its cost
//! grows as 2^(2^n): n = 4 is 65 536 subsets, n = 5 is 2^32 = 4.3e9.
//!
//! `countsByClosure` enumerates the FIXED POINTS of the closure operator
//! `S ↦ smallest closed superset of S`, by Ganter's next closure. The closed
//! sets form a closure system — the intersection of two closed sets is closed,
//! because a product of two elements of an intersection lies in both — so the
//! fixed points are exactly the closed sets, each visited exactly once. The
//! cost is proportional to the NUMBER OF CLOSED SETS, not to the number of
//! subsets: at n = 5 that is 375 for a non-degenerate algebra and 31 242 668
//! for `(0,0,5)`, both measured, against the 4.3e9 subsets a scan would touch.
//!
//! In a non-degenerate algebra the closed sets are exactly the GF(2)-linear
//! subspaces of the blade masks — the product of two blades is ± a single
//! blade `e_{A△B}`, so closure under the product IS closure under symmetric
//! difference — and their number is `1 + Σ_k [n choose k]_2`, the Gaussian
//! binomial sum: 68 at n = 4 and 375 at n = 5, which is what the table
//! reports. A degenerate generator breaks that equivalence (a shared
//! degenerate index makes the product vanish, so fewer products constrain the
//! set) and the count is far larger.
//!
//! The two methods agree on every signature with n <= 4 — all 34 of them,
//! asserted by test — which is what licenses using the second at n = 5.
//!
//! Enumeration stops at n = 5. Above that the masks no longer fit in a `u32`
//! and the degenerate counts are already in the tens of millions; the refusal
//! is explicit, not silent.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const exact = @import("exact.zig");
const sigs = @import("signatures.zig");

/// Largest number of blades at which the subset SCAN is feasible.
pub const MAX_EXHAUSTIVE_BASIS: usize = 16;

/// Largest number of blades the closure walk handles. The blade masks are held
/// in a `u32` throughout this module and m = 2^n, so 32 blades is n = 5.
pub const MAX_CLOSURE_BASIS: usize = 32;

pub const EnumerateError = error{TooManyBlades};

/// What P4 counts for one algebra: the closed blade-spanned subspaces, and the
/// PROPER two-sided ideals among them.
///
/// Every blade-spanned two-sided ideal is closed — for a, b in an ideal the
/// product a·b lies in it — so filtering the closed sets cannot lose an ideal.
/// That is what lets the ideal count ride on the closure walk instead of
/// scanning all 2^m subsets a second time.
pub const Counts = struct {
    closed: usize = 0,
    proper_ideals: usize = 0,
};

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
            if (bp.sign == 0) continue; // zero belongs to every subspace
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
    if (set == 0) return true; // the zero ideal
    const m = alg.basisCount();
    for (0..m) |b| {
        if (set & (@as(u32, 1) << @intCast(b)) == 0) continue;
        for (0..m) |a| { // over all blades of the whole algebra
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
///
/// This is the only entry point still limited to MAX_EXHAUSTIVE_BASIS, and the
/// reason is memory rather than time: it MATERIALISES one `Info` per closed
/// subspace, and `(0,0,5)` has 31 242 668 of them — around half a gigabyte.
/// Counting does not need the list, which is why `countClosed` reaches n = 5
/// and this does not. If you need the n = 5 list, walk it yourself from
/// `countsByClosure` and consume each set instead of collecting it.
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

/// Counts by scanning every subset. The ORACLE: it assumes nothing about the
/// structure, so it is the reference the closure walk is tested against.
pub fn countsBrute(alg: cl.Algebra) EnumerateError!Counts {
    const m = alg.basisCount();
    if (m > MAX_EXHAUSTIVE_BASIS) return error.TooManyBlades;
    var c = Counts{};
    const total: u64 = @as(u64, 1) << @intCast(m);
    var s: u64 = 0;
    while (s < total) : (s += 1) {
        const set: u32 = @intCast(s);
        if (isClosed(alg, set)) c.closed += 1;
        const dim: u5 = @intCast(@popCount(set));
        if (dim == 0 or dim == m) continue; // the trivial ideals
        if (isTwoSidedIdeal(alg, set)) c.proper_ideals += 1;
    }
    return c;
}

/// Working tables for the closure walk, sized for MAX_CLOSURE_BASIS blades so
/// they live on the stack and the walk allocates nothing.
const ClosureTable = struct {
    const ZERO_PRODUCT: u8 = 0xFF;

    /// `prod[i*m + j]`: the mask of `e_i·e_j`, or ZERO_PRODUCT when the
    /// product vanishes (a shared degenerate generator).
    prod: [MAX_CLOSURE_BASIS * MAX_CLOSURE_BASIS]u8 = undefined,
    /// `left[b]`: union of the masks of `a·b` over every blade `a` of the
    /// algebra; `right[b]` the same for `b·a`.
    left: [MAX_CLOSURE_BASIS]u32 = undefined,
    right: [MAX_CLOSURE_BASIS]u32 = undefined,
    m: usize = 0,

    fn init(alg: cl.Algebra) ClosureTable {
        const m = alg.basisCount();
        var t = ClosureTable{ .m = m };
        for (0..m) |i| {
            t.left[i] = 0;
            t.right[i] = 0;
            for (0..m) |j| {
                const bp = cl.bladeMul(alg, @intCast(i), @intCast(j));
                t.prod[i * m + j] = if (bp.sign == 0) ZERO_PRODUCT else @intCast(bp.mask);
            }
        }
        for (0..m) |b| {
            for (0..m) |a| {
                const ab = cl.bladeMul(alg, @intCast(a), @intCast(b));
                if (ab.sign != 0) t.left[b] |= @as(u32, 1) << @intCast(ab.mask);
                const ba = cl.bladeMul(alg, @intCast(b), @intCast(a));
                if (ba.sign != 0) t.right[b] |= @as(u32, 1) << @intCast(ba.mask);
            }
        }
        return t;
    }

    /// The smallest closed superset of `set`.
    ///
    /// The loop terminates after at most `m` rounds: a round either changes
    /// nothing (done) or adds at least one blade, and there are `m` of them.
    /// The `m` is therefore a proof of termination and not a budget — but the
    /// loop is written so that the guard is visible rather than implied.
    fn closure(self: *const ClosureTable, set: u32) u32 {
        var cur = set;
        var round: usize = 0;
        while (round <= self.m) : (round += 1) {
            var out = cur;
            var ii = cur;
            while (ii != 0) {
                const i: usize = @ctz(ii);
                ii &= ii - 1;
                var jj = cur;
                while (jj != 0) {
                    const j: usize = @ctz(jj);
                    jj &= jj - 1;
                    const pm = self.prod[i * self.m + j];
                    if (pm != ZERO_PRODUCT) out |= @as(u32, 1) << @intCast(pm);
                }
            }
            if (out == cur) return cur;
            cur = out;
        }
        return cur;
    }

    /// Is `set` a two-sided ideal? Requires `set` to be closed, which the walk
    /// guarantees for every set it emits.
    fn isIdeal(self: *const ClosureTable, set: u32) bool {
        var needed: u32 = 0;
        var ii = set;
        while (ii != 0) {
            const b: usize = @ctz(ii);
            ii &= ii - 1;
            needed |= self.left[b] | self.right[b];
        }
        return (needed & ~set) == 0;
    }
};

/// Counts by walking the fixed points of the closure operator (Ganter's next
/// closure) instead of scanning subsets. Same answer as `countsBrute` wherever
/// the two can both run; see the module header for why that is exact.
pub fn countsByClosure(alg: cl.Algebra) EnumerateError!Counts {
    const m = alg.basisCount();
    if (m > MAX_CLOSURE_BASIS) return error.TooManyBlades;
    const t = ClosureTable.init(alg);
    // The whole algebra is closed, so it is always the last set emitted.
    const full: u32 = if (m == MAX_CLOSURE_BASIS) ~@as(u32, 0) else (@as(u32, 1) << @intCast(m)) - 1;

    var c = Counts{};
    var a = t.closure(0); // the empty set is closed
    while (true) {
        c.closed += 1;
        const dim = @popCount(a);
        if (dim != 0 and dim != m and t.isIdeal(a)) c.proper_ideals += 1;
        if (a == full) return c;

        // Next closure: the largest i outside `a` whose closure adds nothing
        // below i. A closed set has exactly one successor, so the first such i
        // found scanning downwards is it.
        var next: ?u32 = null;
        var i: usize = m;
        while (i > 0) {
            i -= 1;
            const bit = @as(u32, 1) << @intCast(i);
            if (a & bit != 0) continue;
            const low: u32 = if (i == 0) 0 else (@as(u32, 1) << @intCast(i)) - 1;
            const b = t.closure((a & low) | bit);
            if ((b & ~a) & low == 0) {
                next = b;
                break;
            }
        }
        a = next orelse return c; // unreachable: `full` is closed
    }
}

/// Counts for one algebra, by whichever method reaches it.
pub fn counts(alg: cl.Algebra) EnumerateError!Counts {
    const m = alg.basisCount();
    if (m <= MAX_EXHAUSTIVE_BASIS) return countsBrute(alg);
    return countsByClosure(alg);
}

/// Number of closed subsets (without allocating the result).
pub fn countClosed(alg: cl.Algebra) EnumerateError!usize {
    return (try counts(alg)).closed;
}

/// Number of proper, nonzero, two-sided blade-spanned ideals.
pub fn countProperIdeals(alg: cl.Algebra) EnumerateError!usize {
    return (try counts(alg)).proper_ideals;
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
        try std.testing.expect((ev & 1) != 0); // unital: contains the scalar 1
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

test "enumeration refuses above n = 5 instead of silently running for hours" {
    // (0,6) has 64 blades, so a blade mask no longer fits in a u32.
    const alg6 = try (sigs.SigBuf.build(.{ .p = 0, .q = 6 }, .mostly_minus)).algebra();
    try std.testing.expectEqual(@as(usize, 64), alg6.basisCount());
    try std.testing.expectError(error.TooManyBlades, counts(alg6));
    try std.testing.expectError(error.TooManyBlades, countClosed(alg6));
    try std.testing.expectError(error.TooManyBlades, countProperIdeals(alg6));

    // (0,5) is INSIDE the walk's range and outside the scan's, which is the
    // whole point of having two methods.
    const alg5 = try (sigs.SigBuf.build(.{ .p = 0, .q = 5 }, .mostly_minus)).algebra();
    try std.testing.expectEqual(@as(usize, 32), alg5.basisCount());
    try std.testing.expectError(error.TooManyBlades, countsBrute(alg5));
    _ = try countsByClosure(alg5);
}

test "PREDICTION: the closure walk reproduces the subset scan on every n <= 4" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, 4);
    try std.testing.expectEqual(@as(usize, 34), n);
    for (0..n) |i| {
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const brute = try countsBrute(alg);
        const walk = try countsByClosure(alg);
        if (brute.closed != walk.closed or brute.proper_ideals != walk.proper_ideals) {
            std.debug.print("closure walk disagrees on ({d},{d},{d}): walk {d}/{d}, scan {d}/{d}\n", .{
                buf[i].p,            buf[i].q,
                buf[i].r,            walk.closed,
                walk.proper_ideals,  brute.closed,
                brute.proper_ideals,
            });
            return error.TestUnexpectedResult;
        }
    }
}

test "PREDICTION: the closure walk carries P4 past n = 4" {
    // r = 0: the closed sets are the blade-mask subspaces, 1 + Σ [n choose k]_2.
    const r0 = [_]usize{ 0, 3, 6, 17, 68, 375 };
    for (1..6) |dim| {
        const alg = try (sigs.SigBuf.build(.{ .q = @intCast(dim) }, .mostly_minus)).algebra();
        try std.testing.expectEqual(r0[dim], (try countsByClosure(alg)).closed);
        // ... and none of them is a proper two-sided ideal
        try std.testing.expectEqual(@as(usize, 0), (try countsByClosure(alg)).proper_ideals);
    }
    // The degenerate counts at n = 5, one per r. Every row of the committed
    // table is checked here, but the four large ones only in ReleaseFast:
    // together they cost about 50 seconds in Debug, which is the whole cost of
    // the test suite, and Debug would be exercising the same code path the two
    // small cases above already cover. ReleaseFast is the mode the engine is
    // built and measured in, so that is where the numbers of the report are
    // pinned.
    const deg = [_]usize{ 135534, 733827, 3033464, 10716233, 31242668 };
    const cheap_until = if (@import("builtin").mode == .ReleaseFast) deg.len else 2;
    for (deg[0..cheap_until], 1..) |expected, r| {
        const alg = try (sigs.SigBuf.build(.{ .p = 0, .q = @intCast(5 - r), .r = @intCast(r) }, .mostly_minus)).algebra();
        try std.testing.expectEqual(expected, (try countsByClosure(alg)).closed);
    }
    // (2,3,0) — the signature the audit asked about first, and the reason this
    // method exists at all.
    const alg231 = try (sigs.SigBuf.build(.{ .p = 2, .q = 3 }, .mostly_minus)).algebra();
    try std.testing.expectEqual(@as(usize, 375), (try countsByClosure(alg231)).closed);

    // The audit's measured claim, now guarded: the blade-scoped lattice is
    // blind to the time/space split, so the count depends on n and r and not
    // on how the remaining dimensions are divided between time and space.
    const splits = [_]sigs.Triple{
        .{ .p = 1, .q = 4 }, .{ .p = 2, .q = 3 },
        .{ .p = 3, .q = 2 }, .{ .p = 4, .q = 1 },
        .{ .p = 5, .q = 0 },
    };
    for (splits) |t| {
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        try std.testing.expectEqual(@as(usize, 375), (try countsByClosure(alg)).closed);
    }
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
