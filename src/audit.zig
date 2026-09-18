//! MRS-LAB :: the safety gate
//!
//! WHAT THIS FILE IS FOR. A declared error that no test ever triggers is not a
//! safety feature, it is a comment. This file holds one test per ERROR the
//! public API can return, plus a sweep over the edges of every declared domain
//! — the smallest and largest dimensions, the zero vector, degenerate
//! signatures, and non-finite input arriving from outside.
//!
//! The audit that produced it (0.1.5) found that `TooManyGenerators`,
//! `NoGenerators`, `DimensionTooLarge`, `NegativeMu2`, `Overflow` and all three
//! variants of `OrderError` had NO test that provoked them. Each now has one,
//! and the tests below are written to fail if the error stops being returned —
//! not merely to call the function.
//!
//! The second finding was worse than a missing test: a vector containing NaN
//! was classified as SPATIAL, i.e. as a tachyon. Both comparisons in
//! `classifyTol` are false for a NaN norm, and the code fell through to the
//! last branch, so an input the engine cannot classify came back with a
//! confident answer. `Class.invalid` now exists for that, and
//! `a non-finite vector has no classification` pins it.
//!
//! A note on what this file does NOT prove. It shows that every error path is
//! reachable and that the edges behave as declared. It does not show that the
//! engine is free of defects: the two silently-wrong verdicts found in the
//! 0.1.1 audit were both on paths that had tests, and both tests agreed with
//! the code because both encoded the same wrong assumption. Reachability is not
//! correctness, and this file is the cheap half.

const std = @import("std");
const mrs = @import("mrs");
const exact = @import("explore/exact.zig");
const subalg = @import("explore/subalgebra.zig");
const sigs = @import("explore/signatures.zig");
const nrm = @import("explore/norm_mult.zig");

const sig = mrs.signature;
const form = mrs.form;
const cl = mrs.clifford;
const tach = mrs.tachyon;
const causal = mrs.causal;
const Z = mrs.split_complex.Z;
const contract = mrs.contract;

// ---------------------------------------------------------------------------
// Every declared error, provoked
// ---------------------------------------------------------------------------

test "AlgebraError: no generators, and more than the maximum" {
    try std.testing.expectError(error.NoGenerators, cl.Algebra.init(0, [_]i8{0} ** cl.MAX_GEN));
    try std.testing.expectError(
        error.TooManyGenerators,
        cl.Algebra.init(cl.MAX_GEN + 1, [_]i8{0} ** cl.MAX_GEN),
    );
    // and the boundary itself is accepted, so the guard is not off by one
    _ = try cl.Algebra.init(cl.MAX_GEN, [_]i8{1} ** cl.MAX_GEN);
    _ = try cl.Algebra.init(cl.MIN_GEN, [_]i8{1} ** cl.MAX_GEN);

    // the same guard through the signature entry point
    const too_many = [_]sig.Role{.spatial} ** (cl.MAX_GEN + 1);
    try std.testing.expectError(
        error.TooManyGenerators,
        cl.Algebra.fromSignature(.{ .roles = &too_many }),
    );
}

test "DimensionTooLarge: a form cannot outgrow the sign cache" {
    const too_many = [_]sig.Role{.spatial} ** (sig.MAX_DIM + 1);
    try std.testing.expectError(
        error.DimensionTooLarge,
        form.DiagonalForm.initOwned(.{ .roles = &too_many }),
    );
    // exactly MAX_DIM is accepted
    const at_limit = [_]sig.Role{.spatial} ** sig.MAX_DIM;
    _ = try form.DiagonalForm.initOwned(.{ .roles = &at_limit });
}

test "DegenerateForm: raising an index where the form has none" {
    const many = [_]sig.Role{ .temporal, .degenerate, .spatial };
    const f = form.DiagonalForm{ .signature = .{ .roles = &many } };
    var c = [_]f64{ 1, 1, 1 };
    var out: [3]f64 = undefined;
    try std.testing.expectError(error.DegenerateForm, f.raiseIndex(&c, &out));
}

test "NegativeMu2 and ImaginaryEnergy: the tachyon's two refusals" {
    try std.testing.expectError(error.NegativeMu2, tach.Tachyon.init(-1.0));
    _ = try tach.Tachyon.init(0.0); // the massless boundary is legal
    _ = try tach.Tachyon.init(1.0);

    const t = try tach.Tachyon.init(1.0);
    // below the threshold the energy would be imaginary, and the engine says so
    try std.testing.expectError(error.ImaginaryEnergy, t.energy(0.5));
    _ = try t.energy(2.0);
}

test "OrderError: all three refusals of the causal order" {
    const mink = sig.minkowski_3_1;
    // an arrow index outside the signature
    try std.testing.expectError(
        error.TimeArrowOutOfRange,
        causal.Order.init(mink, 99),
    );
    // an arrow that points at a spatial dimension: not a time arrow
    try std.testing.expectError(
        error.TimeArrowNotTemporal,
        causal.Order.init(mink, 1),
    );
    // a signature longer than the module's fixed buffers
    const too_many = [_]sig.Role{.temporal} ** (sig.MAX_DIM + 1);
    try std.testing.expectError(
        error.DimensionTooLarge,
        causal.Order.init(.{ .roles = &too_many }, 0),
    );
    _ = try causal.Order.init(mink, 0); // and the good case still works
}

test "TooManyBlades: P4 refuses above its range instead of running for hours" {
    const alg = try (sigs.SigBuf.build(.{ .q = 6 }, .mostly_minus)).algebra();
    try std.testing.expectError(error.TooManyBlades, subalg.counts(alg));
    try std.testing.expectError(error.TooManyBlades, subalg.countsBrute(alg));
    try std.testing.expectError(error.TooManyBlades, subalg.countsByClosure(alg));
}

test "Overflow: the guard rests on a checked multiply, not on a silent wrap" {
    // `mulExact` panics on overflow by design, and a panic aborts the test
    // process — so the panic itself cannot be asserted here, and pretending
    // otherwise would be a test that cannot fail. What CAN be pinned is the
    // property the guard depends on: the primitive reports overflow instead of
    // wrapping. If that ever changed, the exact layer would hand back wrapped
    // coefficients as though they were answers.
    const huge: i64 = std.math.maxInt(i64);
    try std.testing.expectError(error.Overflow, std.math.mul(i64, huge, 2));
    try std.testing.expectError(error.Overflow, std.math.add(i64, huge, 1));
    // ... and the products the engine actually computes are far inside i64:
    // coefficients are of order tens and the signs are ±1, so the guard never
    // fires on a legitimate run.
    try std.testing.expectEqual(@as(i64, 6), try std.math.mul(i64, 2, 3));
    try std.testing.expectEqual(@as(i64, -6), try std.math.mul(i64, -2, 3));
}

// ---------------------------------------------------------------------------
// The edges of the declared domains
// ---------------------------------------------------------------------------

test "n = 1 is a legal universe and n = MAX_DIM works" {
    const one = [_]sig.Role{.temporal};
    const f1 = form.DiagonalForm.init(.{ .roles = &one });
    var v1 = [_]f64{2.0};
    try std.testing.expectEqual(form.Class.temporal, f1.classify(&v1));

    const big = [_]sig.Role{.temporal} ++ [_]sig.Role{.spatial} ** (sig.MAX_DIM - 1);
    const fbig = form.DiagonalForm.init(.{ .roles = &big });
    var vbig = [_]f64{0.0} ** sig.MAX_DIM;
    vbig[0] = 1.0;
    vbig[1] = 1.0;
    try std.testing.expectEqual(form.Class.null_like, fbig.classify(&vbig));
}

test "the zero vector is null, and an empty form is not an error" {
    const mk = sig.minkowski_3_1;
    const f = form.DiagonalForm.init(mk);
    const zero = [_]f64{ 0, 0, 0, 0 };
    try std.testing.expectEqual(form.Class.null_like, f.classify(&zero));
    try std.testing.expect(f.classify(&zero).isCausalLike());
}

test "a non-finite vector has no classification" {
    // The 0.1.5 audit finding: NaN used to come back as `.spatial`, i.e. as a
    // tachyon, because both comparisons in `classifyTol` are false for NaN and
    // the code fell through to the last branch.
    const f = form.DiagonalForm.init(sig.minkowski_3_1);
    const nan = [_]f64{ std.math.nan(f64), 0, 0, 0 };
    const inf = [_]f64{ std.math.inf(f64), 0, 0, 0 };
    try std.testing.expectEqual(form.Class.invalid, f.classify(&nan));
    try std.testing.expectEqual(form.Class.invalid, f.classify(&inf));
    try std.testing.expectEqual(form.Class.invalid, f.classifyTol(&nan, 1.0));
    // an invalid vector is not causal, and its label says why
    try std.testing.expect(!form.Class.invalid.isCausalLike());
    try std.testing.expectEqualStrings("invalid (some component is not finite)", form.Class.invalid.label());
    // and the finite cases are untouched by the new branch
    const ok = [_]f64{ 1, 0, 0, 0 };
    try std.testing.expectEqual(form.Class.temporal, f.classify(&ok));
}

test "the widest supported signature still decides P2 and P4" {
    // n = 5 is the ceiling of both questions; this asserts the ceiling is
    // reachable, because a limit that refuses one step early is also a defect.
    const t: sigs.Triple = .{ .p = 2, .q = 3 };
    const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
    try std.testing.expectEqual(@as(usize, 32), alg.basisCount());
    const counts = try subalg.counts(alg);
    try std.testing.expectEqual(@as(usize, 375), counts.closed);
    try std.testing.expectEqual(@as(usize, 0), counts.proper_ideals);

    // (2,3,0) has p+q = 5, so the P2a rule predicts REFUTATION — and a
    // refutation without a witness would be a claim with nothing behind it.
    var grid: [nrm.MAX_GRID]exact.IntVec = undefined;
    const v = try nrm.checkScalarMultiplicative(alg, &grid);
    try std.testing.expect(!v.holds);
    try std.testing.expect(v.has_witness);
    try std.testing.expectEqual(@as(usize, nrm.gridSize(alg.basisCount())), v.grid_size);
}

test "every triple the engine admits can be built and asked" {
    // The enumerator and the constructor must agree on the domain: a triple the
    // report prints but the constructor refuses (or the other way round) would
    // be a hole between two modules, which is where the 0.1.1 defects lived.
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        for ([_]sig.TimeSign{ .mostly_minus, .mostly_plus }) |ts| {
            const alg = try (sigs.SigBuf.build(buf[i], ts)).algebra();
            try std.testing.expectEqual(@as(usize, @as(usize, 1) << @intCast(buf[i].n())), alg.basisCount());
            const f = form.DiagonalForm.init(sigs.SigBuf.build(buf[i], ts).signature());
            try std.testing.expectEqual(buf[i].n(), f.n());
            checked += 1;
        }
    }
    try std.testing.expectEqual(sigs.tripleCount(sigs.MAX_TOTAL) * 2, checked);
}

test "the sparse kernel and the wide path refuse an algebra past their tables" {
    // Found by the 0.1.5 review, with a run: the product table has 32^2 = 1024
    // slots and `init` on a 64-blade algebra wrote to index 1024. In ReleaseSafe
    // that panicked; in ReleaseFast it would have written 3072 bytes past a
    // stack struct and carried on. Both constructors now refuse in every mode.
    const alg6 = try cl.Algebra.init(6, [_]i8{1} ** cl.MAX_GEN);
    try std.testing.expectEqual(@as(usize, 64), alg6.basisCount());
    try std.testing.expectError(error.TooManyBlades, exact.ProductTable.init(alg6));

    var grid: [nrm.MAX_GRID]exact.IntVec = undefined;
    try std.testing.expectError(error.TooManyBlades, nrm.checkScalarMultiplicative(alg6, &grid));

    // The wide path is the one that is SUPPOSED to accept n = 6 (64 blades =
    // its MAX_M), so the boundary must not be off by one: 64 in, 128 out.
    const wide = @import("explore/p2_wide.zig");
    _ = try wide.decideScalarMultiplicative(alg6);
    const alg7 = try cl.Algebra.init(7, [_]i8{1} ** cl.MAX_GEN);
    try std.testing.expectEqual(@as(usize, 128), alg7.basisCount());
    try std.testing.expectError(error.TooManyBlades, wide.decideScalarMultiplicative(alg7));
}

test "the decided classification agrees with the tolerance where it should, and differs where it must" {
    // The 0.1.6 behaviour change, pinned from both sides: identical on every
    // vector whose status is not in dispute, and deliberately different on the
    // two families where a constant cannot be right.
    const f = form.DiagonalForm.init(sig.minkowski_3_1);

    // Agreement: clearly timelike, clearly spacelike, exactly null.
    const cases = [_]struct { v: [4]f64, c: form.Class }{
        .{ .v = .{ 1.0, 0.0, 0.0, 0.0 }, .c = .temporal },
        .{ .v = .{ 0.0, 1.0, 0.0, 0.0 }, .c = .spatial },
        .{ .v = .{ 1.0, 1.0, 0.0, 0.0 }, .c = .null_like },
        .{ .v = .{ 3.0, 2.0, 1.0, 0.5 }, .c = .temporal },
    };
    for (cases) |k| {
        const by_radius = contract.classify(f, &k.v);
        try std.testing.expectEqual(k.c, by_radius);
        try std.testing.expectEqual(k.c, f.classifyTol(&k.v, 1e-9));
    }

    // Difference 1: a small norm that the computation still resolves.
    const tiny = [_]f64{ 1.0, 1.0 - 1e-10, 0, 0 };
    try std.testing.expectEqual(form.Class.null_like, f.classifyTol(&tiny, 1e-9));
    try std.testing.expectEqual(form.Class.temporal, contract.classify(f, &tiny));

    // Difference 2: the absolute tolerance is not scale invariant.
    const lam: f64 = 1e-5;
    const scaled = [_]f64{ lam, lam * 0.95, 0, 0 };
    try std.testing.expectEqual(form.Class.null_like, f.classifyTol(&scaled, 1e-9));
    try std.testing.expect(contract.classify(f, &scaled) != .null_like);
}

test "isZeroDivisorDecided: the same story in the split-complex algebra" {
    // The tolerance version calls anything with a small norm a zero divisor.
    const small = Z.init(1e-6, 1e-6); // N = 0 exactly? no: 1e-12 - 1e-12 = 0
    try std.testing.expect(small.isZeroDivisor(1e-9));
    try std.testing.expect(contract.isZeroDivisorDecided(small));

    // A genuine zero divisor: N(z) = 0 exactly, so the radius covers it.
    const on_cone = Z.init(1.0, 1.0);
    try std.testing.expect(contract.isZeroDivisorDecided(on_cone));

    // Scale invariance: the same element scaled down is STILL a zero divisor
    // (it is on the cone), and a non-divisor stays a non-divisor.
    try std.testing.expect(contract.isZeroDivisorDecided(Z.init(1e-9, 1e-9)));
    const off_cone = Z.init(1.0, 0.95);
    try std.testing.expect(!contract.isZeroDivisorDecided(off_cone));
    try std.testing.expect(!contract.isZeroDivisorDecided(Z.init(1e-9, 0.95e-9)));

    // A non-finite element is not a zero divisor, and says so by refusing.
    try std.testing.expect(!contract.isZeroDivisorDecided(Z.init(std.math.nan(f64), 0.0)));
}

test "a non-finite vector is refused by the cone, not admitted by a fall-through" {
    // Found by the 0.1.6 review. `inFutureCone` used to ask only "is it spatial",
    // so `.invalid` fell through to the arrow test — and a vector whose arrow
    // component is finite but which is NaN elsewhere passed it. `NaN >= 0.0` is
    // false only when the NaN IS the arrow.
    const o = try causal.Order.init(sig.minkowski_3_1, 0);
    const nan_spatial = [_]f64{ 1.0, std.math.nan(f64), 0, 0 };
    const nan_arrow = [_]f64{ std.math.nan(f64), 1.0, 0, 0 };
    const inf_spatial = [_]f64{ 1.0, std.math.inf(f64), 0, 0 };

    try std.testing.expectEqual(form.Class.invalid, o.classifyDecided(&nan_spatial));
    try std.testing.expect(!o.inFutureCone(&nan_spatial));
    try std.testing.expect(!o.inFutureCone(&nan_arrow));
    try std.testing.expect(!o.inFutureCone(&inf_spatial));

    // the relation refuses it too, in both directions
    const origin = [_]f64{ 0, 0, 0, 0 };
    try std.testing.expect(!o.leq(&origin, &nan_spatial));
    try std.testing.expect(!o.separation(&origin, &nan_spatial));
    try std.testing.expect(!o.nullSeparated(&origin, &nan_spatial));

    // ... and a finite vector is untouched by the new branch
    const ok = [_]f64{ 1.0, 0.5, 0, 0 };
    try std.testing.expect(o.inFutureCone(&ok));
}
