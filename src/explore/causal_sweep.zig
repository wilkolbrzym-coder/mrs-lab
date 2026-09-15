//! MRS-LAB :: the causal verdict swept over EVERY signature in range
//!
//! WHY THIS FILE EXISTS. The causality tests in `src/mrs/causal.zig` ran on a
//! handful of canonical signatures, and every one of them had at least one
//! spatial dimension. Signatures with two or three time dimensions and no space
//! — (2,0,0), (3,0,0) — were therefore never exercised, and the engine reported
//! them as PARTIAL ORDERS, contradicting its own Theorem 5.2 (⪯ is a partial
//! order ⟺ p = 1 and r = 0). Two defects hid behind that gap:
//!
//!   1. the antisymmetry witness for p >= 2 insisted on a spatial dimension,
//!      although d = e_t1 already gives u ⪯ u+d ⪯ u with d ≠ 0;
//!   2. for q = 0 the vector generator has no spatial dimension to work with, so
//!      it produced nothing — and an EMPTY search was reported as
//!      "transitive / convex".
//!
//! Both are fixed, and this file is the guard: it sweeps the whole range of
//! triples instead of a list of examples, so either one fails loudly. It also
//! records two structural facts about multi-time signatures that were measured
//! here for the first time: the strict-future sector cannot loop, and the
//! Lorentzian order survives inside a two-time signature as the restriction to a
//! (1,q) subspace.

const std = @import("std");
const mrs = @import("mrs");
const sig = mrs.signature;
const form = mrs.form;
const causal = mrs.causal;
const exact = @import("exact.zig");
const sigs = @import("signatures.zig");
const sub = @import("subalgebra.zig");

const Triple = sigs.Triple;
const MAX_DIM = causal.MAX_DIM;

/// LIFETIME. `Order` stores the signature as a slice into `SigBuf.roles`, so the
/// buffer must outlive the order. These helpers therefore take a pointer to a
/// buffer owned by the CALLER; building the buffer inside the helper and
/// returning the order leaves a dangling slice — that is undefined behaviour,
/// and it shows up as "switch on corrupt value" far away from the cause.
fn classifyTriple(sb: *const sigs.SigBuf, arrow: usize, trials: usize, rnd: std.Random) !causal.Causality {
    const o = try causal.Order.init(sb.signature(), arrow);
    return o.diagnose(trials, rnd);
}

fn orderOf(sb: *const sigs.SigBuf, arrow: usize) !causal.Order {
    return causal.Order.init(sb.signature(), arrow);
}

// ---------------------------------------------------------------------------
// The sweep
// ---------------------------------------------------------------------------

test "SWEEP: ⪯ is a partial order exactly for Lorentzian signatures, for every p >= 1" {
    var prng = std.Random.DefaultPrng.init(20260915);
    const rnd = prng.random();
    var buf: [64]Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);

    var checked: usize = 0;
    var multi_time: usize = 0;
    var no_space: usize = 0;
    var arrows: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.p == 0) continue; // no temporal dimension: no time arrow exists
        checked += 1;
        if (t.p >= 2) multi_time += 1;
        if (t.q == 0) no_space += 1;

        var sb = sigs.SigBuf.build(t, .mostly_minus);
        // Every temporal dimension is a legal arrow (`SigBuf` puts them first), and
        // the verdict must not depend on which one is chosen: the arrow changes the
        // relation, not whether the relation is an order.
        var arrow: usize = 0;
        while (arrow < t.p) : (arrow += 1) {
            const diag = try classifyTriple(&sb, arrow, 2000, rnd);
            if (diag.isPartialOrder() != t.isLorentzian()) {
                std.debug.print(
                    "SWEEP broken: ({d},{d},{d}) arrow={d} partial_order={} lorentzian={} method={s}\n",
                    .{ t.p, t.q, t.r, arrow, diag.isPartialOrder(), t.isLorentzian(), @tagName(diag.method) },
                );
                return error.TestUnexpectedResult;
            }
            // A sampled verdict must never rest on an empty search.
            if (diag.method == .sampling) try std.testing.expect(diag.trials_checked > 0);
            arrows += 1;
        }
    }
    try std.testing.expect(checked >= 35);
    try std.testing.expect(multi_time >= 10);
    try std.testing.expect(no_space >= 4);
    try std.testing.expect(arrows >= 45);
}

test "transitivity is decided where the engine can decide it, sampled only where it must" {
    var prng = std.Random.DefaultPrng.init(7);
    const rnd = prng.random();
    var buf: [64]Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);

    var half_space: usize = 0;
    var witness: usize = 0;
    var sampled: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.p == 0) continue;
        var sb = sigs.SigBuf.build(t, .mostly_minus);
        var arrow: usize = 0;
        while (arrow < t.p) : (arrow += 1) {
            const diag = try classifyTriple(&sb, arrow, 2000, rnd);

            if (t.q == 0) {
                // No spatial dimension: no nonzero vector is spatial, so the cone is
                // the closed half-space {v : v_arrow >= 0} and transitivity holds.
                try std.testing.expectEqual(causal.Method.half_space_proof, diag.method);
                try std.testing.expect(diag.transitive);
                try std.testing.expectEqual(@as(usize, 0), diag.trials_checked);
                half_space += 1;
            } else if (t.p >= 2) {
                // canonicalWitness: u ⪯ w ⪯ v with u ⋠ v, no sampling involved.
                try std.testing.expectEqual(causal.Method.constructive_witness, diag.method);
                try std.testing.expect(!diag.transitive);
                witness += 1;
            } else {
                try std.testing.expectEqual(causal.Method.sampling, diag.method);
                try std.testing.expect(diag.trials_checked > 0);
                try std.testing.expect(diag.transitive);
                sampled += 1;
            }
            try std.testing.expect(diag.isProved() or diag.trials_checked > 0);
        }
    }
    try std.testing.expect(half_space >= 4);
    try std.testing.expect(witness >= 6);
    try std.testing.expect(sampled >= 8);
}

test "q = 0: the generator produces nothing, the search is empty, and no claim is made" {
    var prng = std.Random.DefaultPrng.init(11);
    const rnd = prng.random();
    const t = sig.Role.temporal;
    const two_times: sig.Signature = .{ .roles = &.{ t, t } };
    const o = try causal.Order.init(two_times, 0);

    var v: [MAX_DIM]f64 = undefined;
    try std.testing.expect(!o.randomConeVec(&v, rnd, .timelike_future));
    try std.testing.expect(!o.randomConeVec(&v, rnd, .null_future));

    const probe = causal.probeTransitivity(o, 1000, rnd);
    try std.testing.expectEqual(@as(usize, 0), probe.checked);
    try std.testing.expect(!probe.convex); // an empty search is NOT convexity
    try std.testing.expect(!o.coneIsConvex(1000, rnd));

    // ... and yet the diagnosis still answers, because it does not rely on the
    // generator here: it uses the half-space argument.
    const diag = try o.diagnose(1000, rnd);
    try std.testing.expect(diag.transitive);
    try std.testing.expectEqual(causal.Method.half_space_proof, diag.method);
}

// ---------------------------------------------------------------------------
// Multi-time structure: loops, the arrow, and the Lorentzian core
// ---------------------------------------------------------------------------

test "multi-time: distinct points in both directions, and a 3-cycle" {
    var sb = sigs.SigBuf.build(.{ .p = 2, .q = 1 }, .mostly_minus);
    const o = try orderOf(&sb, 0);
    const zero = [_]f64{ 0, 0, 0 };
    // d = e_t1 + 0.5·e_s: arrow component 0, g = s_t·(1 − 0.25) ≠ 0, so d and −d
    // are both future directed.
    const d = [_]f64{ 0, 1, 0.5 };
    const two_d = [_]f64{ 0, 2, 1 };

    try std.testing.expect(o.leq(&zero, &d));
    try std.testing.expect(o.leq(&d, &zero)); // 2-cycle through distinct points
    try std.testing.expect(o.leq(&d, &two_d));
    try std.testing.expect(o.leq(&two_d, &zero)); // 0 ⪯ d ⪯ 2d ⪯ 0, all distinct
    try std.testing.expect(form.DiagonalForm.classify(o.f, &d) == .temporal);
    try std.testing.expect(form.DiagonalForm.classify(o.f, &two_d) == .temporal);
}

test "multi-time: the arrow coordinate forbids closed loops, addition still leaves the cone" {
    var sb = sigs.SigBuf.build(.{ .p = 2, .q = 1 }, .mostly_minus);
    const o = try orderOf(&sb, 0);
    // Two STRICTLY future steps (arrow component > 0) sum to a SPATIAL vector.
    const d1 = [_]f64{ 0.1, 1.0, 1.0 };
    const d2 = [_]f64{ 0.1, -1.0, 1.0 };
    const sum = [_]f64{ 0.2, 0.0, 2.0 };
    try std.testing.expect(o.inFutureCone(&d1));
    try std.testing.expect(o.inFutureCone(&d2));
    try std.testing.expect(o.f.classify(&d1) == .temporal);
    try std.testing.expect(o.f.classify(&d2) == .temporal);
    try std.testing.expect(d1[0] > 0.0 and d2[0] > 0.0);
    try std.testing.expect(!o.inFutureCone(&sum));

    // Consequence: a closed loop of STRICT steps is impossible, because the arrow
    // coordinate strictly increases and cannot come back to its starting value.
    // Random walk confirmation: the arrow coordinate is monotone along the walk.
    var prng = std.Random.DefaultPrng.init(31337);
    const rnd = prng.random();
    var steps: [causal.MAX_DIM]f64 = undefined;
    for (0..2000) |_| {
        var pos = [_]f64{ 0, 0, 0 };
        for (0..6) |_| {
            if (!o.randomConeVecAny(&steps, rnd)) break;
            if (steps[0] <= 0.0) continue; // keep only strict steps
            const before = pos[0];
            for (0..3) |k| pos[k] += steps[k];
            try std.testing.expect(pos[0] > before);
        }
    }
}

test "multi-time: the arrow is extra structure, and it changes the relation" {
    var sb = sigs.SigBuf.build(.{ .p = 2, .q = 3 }, .mostly_minus);
    const o0 = try orderOf(&sb, 0);
    const o1 = try orderOf(&sb, 1);
    var prng = std.Random.DefaultPrng.init(99);
    const rnd = prng.random();

    var differ: usize = 0;
    var witness: [5]f64 = undefined;
    for (0..20000) |_| {
        var x: [5]f64 = undefined;
        for (0..5) |i| x[i] = (rnd.float(f64) - 0.5) * 4.0;
        if (o0.inFutureCone(&x) != o1.inFutureCone(&x)) {
            if (differ == 0) witness = x;
            differ += 1;
        }
    }
    try std.testing.expect(differ > 0);
    // The witness is a vector both orders classify as in the cone / outside it,
    // in opposite ways: the form alone cannot pick one of the two arrows.
    try std.testing.expect(o0.inFutureCone(&witness) != o1.inFutureCone(&witness));
}

test "a Lorentzian order survives inside every two-time signature, as a (1,q) restriction" {
    // Ambient (2,3) against the Lorentzian (1,3) on the subspace spanned by one
    // temporal and all spatial dimensions (the second temporal component is 0).
    var sb_a = sigs.SigBuf.build(.{ .p = 2, .q = 3 }, .mostly_minus);
    var sb_s = sigs.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus);
    const o_a = try orderOf(&sb_a, 0);
    const o_s = try orderOf(&sb_s, 0);
    var prng = std.Random.DefaultPrng.init(4242);
    const rnd = prng.random();

    var pairs: usize = 0;
    for (0..20000) |_| {
        var ua: [5]f64 = undefined;
        var va: [5]f64 = undefined;
        for (0..5) |i| {
            ua[i] = rnd.float(f64) - 0.5;
            va[i] = rnd.float(f64) - 0.5;
        }
        ua[1] = 0.0;
        va[1] = 0.0; // stay inside the subspace
        const su = [4]f64{ ua[0], ua[2], ua[3], ua[4] };
        const sv = [4]f64{ va[0], va[2], va[3], va[4] };
        try std.testing.expectEqual(o_s.leq(&su, &sv), o_a.leq(ua[0..5], va[0..5]));
        pairs += 1;
    }
    try std.testing.expect(pairs > 15000);

    // Transitivity inside the subspace, measured with the AMBIENT relation:
    // 0 ⪯ w ⪯ v for cone vectors generated inside the subspace.
    var checked: usize = 0;
    for (0..20000) |_| {
        var d1: [MAX_DIM]f64 = undefined;
        var d2: [MAX_DIM]f64 = undefined;
        if (!o_s.randomConeVec(&d1, rnd, .timelike_future)) continue;
        if (!o_s.randomConeVec(&d2, rnd, .timelike_future)) continue;
        const zero = [5]f64{ 0, 0, 0, 0, 0 };
        const wv = [5]f64{ d1[0], 0, d1[1], d1[2], d1[3] };
        const vv = [5]f64{ d1[0] + d2[0], 0, d1[1] + d2[1], d1[2] + d2[2], d1[3] + d2[3] };
        try std.testing.expect(o_a.leq(&zero, &wv));
        try std.testing.expect(o_a.leq(&wv, &vv));
        try std.testing.expect(o_a.leq(&zero, &vv)); // transitivity holds here
        checked += 1;
    }
    try std.testing.expect(checked > 15000);
}

// ---------------------------------------------------------------------------
// The algebra side: the P4 lattice is blind to the time/space split
// ---------------------------------------------------------------------------

const CensusKey = struct {
    n: usize,
    r: u5,
    closed: usize,
    ideals: usize,
    centre: u5,
    first: Triple,
};

test "P4 counts depend only on n and r, never on how many dimensions are time" {
    var buf: [64]Triple = undefined;
    const n = sigs.enumerateTriples(&buf, 4);

    var keys: [16]CensusKey = undefined;
    var used: usize = 0;

    for (0..n) |i| {
        const t = buf[i];
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const closed = try sub.countClosed(alg);
        const ideals = try sub.countProperIdeals(alg);
        const centre = exact.dimOfSet(exact.centerBasis(alg));

        var found = false;
        for (0..used) |k| {
            if (keys[k].n != t.n() or keys[k].r != t.r) continue;
            found = true;
            if (keys[k].closed != closed or keys[k].ideals != ideals or keys[k].centre != centre) {
                std.debug.print(
                    "SPLIT MATTERS: ({d},{d},{d}) closed={d}/{d} ideals={d}/{d} centre={d}/{d} vs first ({d},{d},{d})\n",
                    .{ t.p, t.q, t.r, closed, keys[k].closed, ideals, keys[k].ideals, centre, keys[k].centre, keys[k].first.p, keys[k].first.q, keys[k].first.r },
                );
                return error.TestUnexpectedResult;
            }
        }
        if (!found) {
            keys[used] = .{ .n = t.n(), .r = t.r, .closed = closed, .ideals = ideals, .centre = centre, .first = t };
            used += 1;
        }
    }
    try std.testing.expect(used >= 8);
}
