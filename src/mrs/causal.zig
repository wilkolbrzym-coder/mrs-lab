//! MRS-LAB :: causality module
//!
//! Two relations that must be told apart — and conflating them was the most
//! serious error in the original project specification:
//!
//!   * SEPARATION (unoriented).  u ⊑ v  ⟺  g(v−u, v−u) ≤ 0.
//!     "The difference is temporal or null". Symmetric on null separations,
//!     therefore NOT antisymmetric, therefore not an order.
//!
//!   * CAUSAL ORDER (oriented). u ⪯ v ⟺ (v−u) lies in the closed future
//!     cone, i.e. the class of the difference is `temporal` or `null_like`
//!     AND the component along the chosen time arrow is >= 0.
//!
//! Orientation requires extra structure: a choice of time arrow. For p = 1 the
//! form determines it up to a sign; for p >= 2 it does not determine it at all.
//! That is not an implementation detail, it is the content of Theorem 5.2.
//!
//! ---------------------------------------------------------------------------
//! MAIN RESULT OF THE LABORATORY
//! ---------------------------------------------------------------------------
//! Theorem 5.1 (transitivity requires a single time).
//!   For r = 0 and q >= 1: ⪯ is transitive  ⟺  p = 1.
//!   (<=) The future cone is convex, and convexity of the cone gives
//!        transitivity: a, b in the cone imply a + b in the cone.
//!   (=>) For p >= 2 there is a witness — see `canonicalWitness`.
//!
//! Theorem 5.2 (⪯ is a partial order exactly for Lorentzian signatures).
//!   Reflexivity: always.
//!   Transitivity: ⟺ p = 1 (with r = 0).
//!   Antisymmetry: ⟺ r = 0 and p <= 1.
//!       Proof (<=) for p = 1, r = 0: if d and −d are both in the cone, the
//!       arrow component forces d_arrow = 0, and then g(d) >= 0 forces all
//!       spatial components to vanish, so d = 0.
//!       Proof (=>) for p >= 2: a vector with zero arrow component, a nonzero
//!       second temporal component and a sufficiently small spatial part lies
//!       in the cone together with its opposite.
//!       For r > 0: every radical vector is such a witness.
//!   Conclusion: ⪯ is a partial order ⟺ p = 1 and r = 0 ⟺ `isLorentzian()`.
//!   That means the `isLorentzian` predicate is not decoration — it is the
//!   necessary and sufficient condition for causality to be an order.
//!   causality is an order.
//!
//! Theorem 5.3 (r > 0 destroys causality).
//!   If r > 0 there is a radical vector d ≠ 0 with g(d,d) = 0 and zero arrow
//!   component, so u ⪯ u+d ⪯ u with d ≠ 0. The relation stops telling points
//!   apart. Causality requires r = 0.

const std = @import("std");
const sig = @import("signature.zig");
const form = @import("form.zig");
const contract = @import("contract.zig");

const Signature = sig.Signature;
const DiagonalForm = form.DiagonalForm;

pub const MAX_DIM: usize = 8;

pub const OrderError = error{
    TimeArrowOutOfRange,
    TimeArrowNotTemporal,
    DimensionTooLarge,
};

/// A witness of a property violation. The vectors live on the stack, because the
/// causality module is used inside simulation hot loops.
pub const Witness = struct {
    u: [MAX_DIM]f64 = [_]f64{0} ** MAX_DIM,
    w: [MAX_DIM]f64 = [_]f64{0} ** MAX_DIM,
    v: [MAX_DIM]f64 = [_]f64{0} ** MAX_DIM,
    len: usize = 0,

    pub fn uSlice(self: *const Witness) []const f64 {
        return self.u[0..self.len];
    }
    pub fn wSlice(self: *const Witness) []const f64 {
        return self.w[0..self.len];
    }
    pub fn vSlice(self: *const Witness) []const f64 {
        return self.v[0..self.len];
    }

    pub fn writeTo(self: *const Witness, w: anytype) !void {
        try w.writeAll("u=[");
        try writeVec(w, self.uSlice());
        try w.writeAll("]  w=[");
        try writeVec(w, self.wSlice());
        try w.writeAll("]  v=[");
        try writeVec(w, self.vSlice());
        try w.writeAll("]");
    }
};

fn writeVec(w: anytype, v: []const f64) !void {
    for (v, 0..) |x, i| {
        if (i != 0) try w.writeAll(", ");
        try w.print("{d:.3}", .{x});
    }
}

pub const ConeMode = enum {
    /// Class `temporal`, future directed.
    timelike_future,
    /// Class `null_like`, future directed.
    null_future,
};

/// How a verdict was reached. Scope is part of every answer, and for
/// transitivity the scope is the method: a proof is not a sample.
pub const Method = enum {
    /// Sampling only; `trials_checked` says how many trials produced a vector.
    sampling,
    /// `canonicalWitness` gives a concrete u ⪯ w ⪯ v with u ⋠ v.
    constructive_witness,
    /// q = 0: the cone is the closed half-space {v : v_arrow >= 0}.
    half_space_proof,
};

/// What goes wrong in a given signature.
pub const Causality = struct {
    reflexive: bool,
    transitive: bool,
    antisymmetric: bool,
    /// How `transitive` was obtained.
    method: Method = .sampling,
    /// Successful trials behind a sampled verdict — 0 when nothing was sampled.
    /// A sampled verdict with `trials_checked == 0` is not evidence, and
    /// `probeTransitivity` refuses to call it convex.
    trials_checked: usize = 0,

    pub fn isPartialOrder(self: Causality) bool {
        return self.reflexive and self.transitive and self.antisymmetric;
    }

    /// True when the verdict does not rest on sampling at all.
    pub fn isProved(self: Causality) bool {
        return self.method != .sampling;
    }

    pub fn label(self: Causality) []const u8 {
        if (self.isPartialOrder()) return "partial order";
        if (self.transitive) return "preorder (no antisymmetry)";
        return "neither transitive nor antisymmetric";
    }
};

/// The oriented causal order.
///
/// LIFETIME CONTRACT. `Order` stores the signature as a slice (`Signature.roles`),
/// so whatever buffer those roles point into — in the explore layer a `SigBuf` —
/// MUST outlive the order. Building that buffer inside a helper that RETURNS the
/// order leaves a dangling slice: undefined behaviour, whose symptom is a
/// "switch on corrupt value" panic far away from the cause.
pub const Order = struct {
    f: DiagonalForm,
    /// Roles copied BY VALUE, and a self-contained form built by `initOwned`.
    /// The reason is a lifetime hazard that was real: `Signature` borrows a
    /// `[]const Role` slice, so an `Order` built from a temporary buffer and
    /// returned from a helper used to dangle — it compiled and then failed at
    /// runtime with a corrupt-value switch. Copying removes the hazard at the
    /// type level instead of documenting it.
    roles: [MAX_DIM]sig.Role = undefined,
    len: u5 = 0,
    time_sign: sig.TimeSign = .mostly_minus,
    // REMOVED in 0.1.6: the field `tol: f64 = 1e-9`. It is worth saying why it
    // is gone rather than deprecated. Every relation on this type now decides
    // with the propagated radius, so nothing read the field any more — and a
    // public knob that silently does nothing is worse than no knob at all.
    // Callers who want a tolerance notion have `DiagonalForm.classifyTol(v, tol)`
    // and `split_complex.isZeroDivisor(a, tol)`; they are named as alternatives
    // on purpose, because they answer a different question.
    /// Index of the dimension chosen as the time arrow. For p >= 2 this is
    /// EXTRA structure that the form itself does not determine.
    time_arrow: usize,

    pub fn init(s: Signature, time_arrow: usize) OrderError!Order {
        if (s.n() > MAX_DIM) return error.DimensionTooLarge;
        if (time_arrow >= s.n()) return error.TimeArrowOutOfRange;
        if (s.roles[time_arrow] != .temporal) return error.TimeArrowNotTemporal;
        var o = Order{
            .f = try DiagonalForm.initOwned(s),
            .len = @intCast(s.n()),
            .time_sign = s.time_sign,
            .time_arrow = time_arrow,
        };
        for (0..s.n()) |i| o.roles[i] = s.roles[i];
        return o;
    }

    pub fn dim(self: Order) usize {
        return self.f.n();
    }

    /// A view into this Order's own copy of the roles, so it is valid for as
    /// long as the Order is. Takes a pointer, not a value: returning a slice
    /// into a by-value parameter would dangle.
    pub fn signature(self: *const Order) Signature {
        return .{ .roles = self.roles[0..self.len], .time_sign = self.time_sign };
    }

    pub fn norm2(self: Order, v: []const f64) f64 {
        return self.f.eval(v);
    }

    /// Closed future cone. The sign convention enters EXCLUSIVELY through
    /// `classify`, which is why this works for (+,−,−,−) as well,
    ///     and for (−,+,+,+).
    pub fn inFutureCone(self: Order, v: []const f64) bool {
        const cls = self.classifyDecided(v);
        // `.invalid` is REFUSED, not admitted. A vector with a non-finite
        // component has no norm, so no statement about the cone is true of it,
        // and the arrow test below cannot catch it on its own: a vector with a
        // finite arrow component and a NaN elsewhere used to pass, because
        // `NaN >= 0.0` is false only when the NaN is the arrow itself.
        if (cls == .spatial or cls == .invalid) return false;
        return v[self.time_arrow] >= 0.0;
    }

    /// Classification with the bound PROPAGATED from the computation, instead
    /// of a tolerance supplied from outside.
    ///
    /// Why this replaced `classifyTol(v, self.tol)` in the predicates above:
    /// the tolerance answered "is |g| smaller than 1e-9", which calls an
    /// ordinary vector with a small norm null. `1 - (1 - 1e-10)²` is about
    /// 2e-10, computed to a radius of order 1e-16 — the sign is decided and
    /// the vector is timelike, and only a constant said otherwise. Here
    /// `.null_like` means "the computation cannot separate the norm from zero",
    /// which is a statement about the computation and not about a threshold.
    ///
    /// A vector that is not finite has no classification: `Class.invalid`.
    pub fn classifyDecided(self: Order, v: []const f64) form.Class {
        return contract.classify(self.f, v);
    }

    /// Oriented relation: u ⪯ v.
    pub fn leq(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.inFutureCone(d[0..u.len]);
    }

    /// Unoriented separation: g(v−u, v−u) ≤ 0.
    /// This is NOT the same relation as ⪯ — it is symmetric on the null
    /// cone, so it is not an order.
    pub fn separation(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        const cls = self.classifyDecided(d[0..u.len]);
        return cls == .temporal or cls == .null_like;
    }

    /// Null separation: g(v−u, v−u) = 0 — an equivalence relation whose
    /// classes are the light rays.
    pub fn nullSeparated(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.classifyDecided(d[0..u.len]) == .null_like;
    }

    // -- cone vector generators ----------------------------------------------

    /// Random vector in the future cone. The construction is uniform across
    /// sign conventions:
    /// 1. draw the components outside the time arrow,
    /// 2. scale the spatial group so that its contribution to g is
    /// Sum_temporal_rest x² + 1 — then the contribution of all
    /// dimensions outside the arrow is exactly −s_arrow,
    /// 3. g = s_arrow(t² − 1), so t = 1 gives a null vector and t = 2
    /// timelike vector — in BOTH conventions, because "temporal" means
    /// "the temporal component dominates", not "g has a given sign".
    pub fn randomConeVec(
        self: Order,
        out: *[MAX_DIM]f64,
        rnd: std.Random,
        mode: ConeMode,
    ) bool {
        const s = self.signature();
        const nn = s.n();
        const s_arrow = s.signAt(self.time_arrow);

        for (0..nn) |i| out[i] = (rnd.float(f64) - 0.5) * 2.0;
        out[self.time_arrow] = 0;

        var q_space: f64 = 0;
        var q_time_rest: f64 = 0;
        for (0..nn) |i| {
            if (i == self.time_arrow) continue;
            switch (s.roles[i]) {
                .spatial => q_space += out[i] * out[i],
                .temporal => q_time_rest += out[i] * out[i],
                .degenerate => {},
            }
        }

        if (q_space <= 1e-12) return false; // no spatial dimension = no cone
        const f = @sqrt((q_time_rest + 1.0) / q_space);
        for (0..nn) |i| {
            if (s.roles[i] == .spatial) out[i] *= f;
        }

        var other: f64 = 0;
        for (0..nn) |i| {
            if (i == self.time_arrow) continue;
            other += s.signAt(i) * out[i] * out[i];
        }
        if (@abs(other + s_arrow) > 1e-9) return false; // construction failed

        out[self.time_arrow] = switch (mode) {
            .null_future => 1.0,
            .timelike_future => 2.0,
        };

        // Also decided rather than tolerance-checked: the generator must be able
        // to certify that the vector it built really is null, not merely that
        // the norm is under a threshold. If the construction fails, the caller
        // sees fewer successful trials — which the verdict now reports — rather
        // than a vector that only looked null.
        const cls = self.classifyDecided(out[0..nn]);
        return switch (mode) {
            .null_future => cls == .null_like,
            .timelike_future => cls == .temporal,
        };
    }

    /// Random vector in the cone with a random mode (timelike or null).
    /// Counterexamples to transitivity for p >= 2 require NULL vectors
    /// (proof: the sum of two timelike vectors with the same arrow sense
    /// stays in the cone), so the generator must be able to produce both.
    pub fn randomConeVecAny(self: Order, out: *[MAX_DIM]f64, rnd: std.Random) bool {
        const mode: ConeMode = if (rnd.float(f64) < 0.5) .timelike_future else .null_future;
        return self.randomConeVec(out, rnd, mode);
    }

    // -- witness search ------------------------------------------------------

    /// Transitivity violation ⪯: u ⪯ w, w ⪯ v, but u ⋠ v.
    pub fn searchTransitivityViolation(
        self: Order,
        trials: usize,
        rnd: std.Random,
    ) ?Witness {
        const nn = self.dim();
        var u: [MAX_DIM]f64 = undefined;
        var d1: [MAX_DIM]f64 = undefined;
        var d2: [MAX_DIM]f64 = undefined;
        var total: [MAX_DIM]f64 = undefined;

        for (0..trials) |_| {
            for (0..nn) |i| u[i] = (rnd.float(f64) - 0.5) * 4.0;
            if (!self.randomConeVecAny(&d1, rnd)) continue;
            if (!self.randomConeVecAny(&d2, rnd)) continue;
            for (0..nn) |i| total[i] = d1[i] + d2[i];
            if (!self.inFutureCone(total[0..nn])) {
                var wit = Witness{ .len = nn };
                @memcpy(wit.u[0..nn], u[0..nn]);
                for (0..nn) |i| {
                    wit.w[i] = u[i] + d1[i];
                    wit.v[i] = u[i] + d1[i] + d2[i];
                }
                return wit;
            }
        }
        return null;
    }

    /// Is the future cone convex (sampling)?
    /// Convexity ⟺ transitivity of ⪯.
    pub fn coneIsConvex(self: Order, trials: usize, rnd: std.Random) bool {
        const nn = self.dim();
        var a: [MAX_DIM]f64 = undefined;
        var b: [MAX_DIM]f64 = undefined;
        var sum: [MAX_DIM]f64 = undefined;
        var checked: usize = 0;
        for (0..trials) |_| {
            if (!self.randomConeVecAny(&a, rnd)) continue;
            if (!self.randomConeVecAny(&b, rnd)) continue;
            for (0..nn) |i| sum[i] = a[i] + b[i];
            checked += 1;
            if (!self.inFutureCone(sum[0..nn])) return false;
        }
        return checked > 0;
    }

    /// Constructive transitivity witness for p >= 2.
    /// Dimensions: t0 and t1 temporal, s spatial.
    ///   a   = e_t0 + e_s      → g = s_t − s_s = 0         (null)
    ///   v−a = e_t1 + e_s      → g = 0                      (null)
    ///   v   = e_t0 + e_t1 + 2·e_s → g = 2·s_t − 4·s_s < 0  (spatial)
    /// The inequality "2·s_t < 4·s_s" holds in both conventions,
    /// because s_s = −s_t.
    pub fn canonicalWitness(s: Signature) OrderError!?Witness {
        if (s.p() < 2) return null;
        var idx_t: [2]usize = undefined;
        var nt: usize = 0;
        var idx_s: usize = std.math.maxInt(usize);
        for (0..s.n()) |i| {
            switch (s.roles[i]) {
                .temporal => {
                    if (nt < 2) {
                        idx_t[nt] = i;
                        nt += 1;
                    }
                },
                .spatial => if (idx_s == std.math.maxInt(usize)) {
                    idx_s = i;
                },
                .degenerate => {},
            }
        }
        if (nt < 2 or idx_s == std.math.maxInt(usize)) return null;

        var wit = Witness{ .len = s.n() };
        wit.w[idx_t[0]] = 1.0;
        wit.w[idx_s] = 1.0;
        wit.v[idx_t[0]] = 1.0;
        wit.v[idx_t[1]] = 1.0;
        wit.v[idx_s] = 2.0;
        return wit;
    }

    /// Witness of missing antisymmetry of ⪯: u ≠ v with u ⪯ v and v ⪯ u.
    /// It exists exactly when r > 0 or p >= 2 (Theorem 5.2).
    /// Returns `null` for Lorentzian signatures — and that is consistent
    /// with the theorem, not a coincidence.
    pub fn searchAntisymmetryViolation(s: Signature) OrderError!?Witness {
        const nn = s.n();
        var idx_t: [2]usize = undefined;
        var nt: usize = 0;
        var idx_s: usize = std.math.maxInt(usize);
        var idx_rad: usize = std.math.maxInt(usize);
        for (0..nn) |i| {
            switch (s.roles[i]) {
                .temporal => {
                    if (nt < 2) {
                        idx_t[nt] = i;
                        nt += 1;
                    }
                },
                .spatial => if (idx_s == std.math.maxInt(usize)) {
                    idx_s = i;
                },
                .degenerate => if (idx_rad == std.math.maxInt(usize)) {
                    idx_rad = i;
                },
            }
        }

        var wit = Witness{ .len = nn };
        if (idx_rad != std.math.maxInt(usize)) {
            // radical vector: nonzero, zero norm, zero arrow component
            wit.v[idx_rad] = 1.0;
            return wit;
        }
        if (nt >= 2) {
            // d = e_t1 has arrow component 0 and g(d) = s_t ≠ 0, so d and −d both
            // lie in the closed future cone: u ⪯ u+d ⪯ u with d ≠ 0. No spatial
            // dimension is needed for this witness — requiring one used to hide
            // the violation, and (2,0,0), (3,0,0) were reported as partial orders
            // against Theorem 5.2.
            wit.v[idx_t[1]] = 1.0;
            if (idx_s != std.math.maxInt(usize)) wit.v[idx_s] = 0.5;
            return wit;
        }
        return null;
    }

    /// Full diagnosis: is ⪯ a partial order, and what exactly fails.
    ///
    /// Transitivity is DECIDED wherever the engine can decide it, and the method
    /// travels with the verdict (see `Causality.method`):
    ///   * q = 0 — with no spatial dimension no nonzero vector is spatial, so the
    ///     future cone is the closed half-space {v : v_arrow >= 0}, which is closed
    ///     under addition: transitivity HOLDS, by argument, not by sampling;
    ///   * p >= 2 — `canonicalWitness` exhibits u ⪯ w ⪯ v with u ⋠ v, so
    ///     transitivity FAILS constructively;
    ///   * p = 1 — sampling, with the number of successful trials in
    ///     `trials_checked`, so that an empty search cannot pass as a verdict.
    pub fn diagnose(self: Order, trials: usize, rnd: std.Random) OrderError!Causality {
        const s = self.signature();
        var c = Causality{
            .reflexive = true,
            .transitive = true,
            .antisymmetric = (try Order.searchAntisymmetryViolation(s)) == null,
        };
        if (s.q() == 0) {
            c.method = .half_space_proof;
        } else if (s.p() >= 2) {
            c.transitive = false;
            c.method = .constructive_witness;
        } else {
            const probe = probeTransitivity(self, trials, rnd);
            c.transitive = probe.convex;
            c.trials_checked = probe.checked;
        }
        return c;
    }
};

/// Experimental verification of Theorem 5.1 (sampling, not a decision).
pub const Transitivity = struct {
    violations: usize,
    trials: usize,
    checked: usize,
    /// True only when at least one trial produced a pair of cone vectors AND no
    /// violation was found. An empty search is not evidence: for q = 0 the vector
    /// generator has no spatial dimension to work with and reports `checked = 0`,
    /// which must never be read as "convex".
    convex: bool,
};

pub fn probeTransitivity(self: Order, trials: usize, rnd: std.Random) Transitivity {
    const nn = self.dim();
    var d1: [MAX_DIM]f64 = undefined;
    var d2: [MAX_DIM]f64 = undefined;
    var total: [MAX_DIM]f64 = undefined;
    var violations: usize = 0;
    var checked: usize = 0;
    for (0..trials) |_| {
        if (!self.randomConeVecAny(&d1, rnd)) continue;
        if (!self.randomConeVecAny(&d2, rnd)) continue;
        for (0..nn) |i| total[i] = d1[i] + d2[i];
        checked += 1;
        if (!self.inFutureCone(total[0..nn])) violations += 1;
    }
    return .{
        .violations = violations,
        .trials = trials,
        .checked = checked,
        .convex = checked > 0 and violations == 0,
    };
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "Minkowski 3+1: ⪯ is a partial order" {
    var prng = std.Random.DefaultPrng.init(2024);
    const rnd = prng.random();
    const o = try Order.init(sig.minkowski_3_1, 0);

    const zero = [_]f64{ 0, 0, 0, 0 };
    try std.testing.expect(o.leq(&zero, &zero)); // reflexivity

    const res = probeTransitivity(o, 20_000, rnd);
    try std.testing.expectEqual(@as(usize, 0), res.violations);
    try std.testing.expect(res.checked > 15_000);
    try std.testing.expect(res.convex);

    // antisymmetry: no witness, and that is the content of the theorem
    try std.testing.expectEqual(
        @as(?Witness, null),
        try Order.searchAntisymmetryViolation(sig.minkowski_3_1),
    );
    try std.testing.expectEqual(
        @as(?Witness, null),
        try Order.canonicalWitness(sig.minkowski_3_1),
    );

    const diag = try o.diagnose(20_000, rnd);
    try std.testing.expect(diag.isPartialOrder());
}

test "Theorem 5.1: p = 2 breaks transitivity" {
    const o = try Order.init(sig.two_times_2_1, 0);

    const wit = (try Order.canonicalWitness(sig.two_times_2_1)).?;
    const u = wit.uSlice();
    const wv = wit.wSlice();
    const v = wit.vSlice();

    try std.testing.expect(o.leq(u, wv)); // 0 ⪯ a
    try std.testing.expect(o.leq(wv, v)); // a ⪯ v
    try std.testing.expect(!o.leq(u, v)); // 0 ⋠ v  ← naruszenie

    // concrete numbers: a null, v−a null, v spatial
    try std.testing.expectEqual(form.Class.null_like, o.f.classify(wv));
    var diff: [3]f64 = undefined;
    for (0..3) |i| diff[i] = v[i] - wv[i];
    try std.testing.expectEqual(form.Class.null_like, o.f.classify(&diff));
    try std.testing.expectEqual(form.Class.spatial, o.f.classify(v));

    // independent confirmation: random trials also find a violation
    var prng = std.Random.DefaultPrng.init(31);
    try std.testing.expect(o.searchTransitivityViolation(5000, prng.random()) != null);
}

test "Theorem 5.2: no antisymmetry for p >= 2 and for r > 0" {
    // p = 2, r = 0
    {
        const o = try Order.init(sig.two_times_2_1, 0);
        const wit = (try Order.searchAntisymmetryViolation(sig.two_times_2_1)).?;
        const u = wit.uSlice();
        const v = wit.vSlice();
        try std.testing.expect(o.leq(u, v));
        try std.testing.expect(o.leq(v, u));
        try std.testing.expect(v[1] != 0.0); // u ≠ v
    }
    // r = 1
    {
        const o = try Order.init(sig.degenerate_2_1_1, 0);
        const wit = (try Order.searchAntisymmetryViolation(sig.degenerate_2_1_1)).?;
        const u = wit.uSlice();
        const v = wit.vSlice();
        try std.testing.expect(o.leq(u, v));
        try std.testing.expect(o.leq(v, u));
        try std.testing.expect(v[3] != 0.0); // radical vector, u ≠ v
    }
    // p = 1, r = 0: no witness — a partial order
    {
        try std.testing.expectEqual(
            @as(?Witness, null),
            try Order.searchAntisymmetryViolation(sig.minkowski_3_1),
        );
    }
    // Euclidean (p = 0): no cone, no time arrow
    try std.testing.expectError(
        error.TimeArrowNotTemporal,
        Order.init(sig.euclidean_4, 0),
    );
}

test "Theorem 5.2 — full characterisation: order ⟺ Lorentzian" {
    var prng = std.Random.DefaultPrng.init(4242);
    // Signatures with no spatial dimension are in this list on purpose: the
    // p >= 2 witness does not need one, and omitting them hid a wrong verdict
    // (see the test below).
    const t = sig.Role.temporal;
    const d_ = sig.Role.degenerate;
    const two_times_no_space: sig.Signature = .{ .roles = &.{ t, t } };
    const three_times_no_space: sig.Signature = .{ .roles = &.{ t, t, t } };
    const time_only: sig.Signature = .{ .roles = &.{t} };
    const time_and_radical: sig.Signature = .{ .roles = &.{ t, d_ } };
    const cases = [_]sig.Signature{
        sig.minkowski_3_1,
        sig.minkowski_3_1_flipped,
        sig.minkowski_1_1,
        sig.two_times_2_1,
        sig.degenerate_2_1_1,
        two_times_no_space,
        three_times_no_space,
        time_only,
        time_and_radical,
    };
    for (cases) |s| {
        const o = try Order.init(s, 0);
        const diag = try o.diagnose(4000, prng.random());
        try std.testing.expectEqual(s.isLorentzian(), diag.isPartialOrder());
    }
}

test "q = 0 is decided, never sampled: the cone is a half-space" {
    var prng = std.Random.DefaultPrng.init(77);
    const rnd = prng.random();
    const t = sig.Role.temporal;
    const two_times_no_space: sig.Signature = .{ .roles = &.{ t, t } };
    const o = try Order.init(two_times_no_space, 0);

    // The vector generator has no spatial dimension to work with, so it produces
    // nothing. Before the fix `convex` was then reported as true, i.e. an empty
    // search passed as evidence.
    var v: [MAX_DIM]f64 = undefined;
    try std.testing.expect(!o.randomConeVec(&v, rnd, .timelike_future));
    const probe = probeTransitivity(o, 1000, rnd);
    try std.testing.expectEqual(@as(usize, 0), probe.checked);
    try std.testing.expect(!probe.convex);

    // The diagnosis still answers, and says how: by the half-space argument.
    const diag = try o.diagnose(1000, rnd);
    try std.testing.expect(diag.transitive);
    try std.testing.expectEqual(Method.half_space_proof, diag.method);
    try std.testing.expect(diag.isProved());
    try std.testing.expect(!diag.isPartialOrder()); // 0 ⪯ e_t1 ⪯ 0 with e_t1 ≠ 0

    // Independent confirmation of the half-space reading: every vector with a
    // nonnegative arrow component is in the cone, and sums keep that property.
    var a: [MAX_DIM]f64 = undefined;
    var b: [MAX_DIM]f64 = undefined;
    for (0..2000) |_| {
        for (0..2) |i| {
            a[i] = rnd.float(f64) * 2.0 - 0.5; // arrow component >= -0.5
            b[i] = rnd.float(f64) * 2.0 - 0.5;
        }
        if (a[0] < 0.0 or b[0] < 0.0) continue;
        try std.testing.expect(o.inFutureCone(a[0..2]));
        try std.testing.expect(o.inFutureCone(b[0..2]));
        var sum: [2]f64 = undefined;
        for (0..2) |i| sum[i] = a[i] + b[i];
        try std.testing.expect(o.inFutureCone(&sum));
    }
}

test "the antisymmetry witness for p >= 2 does not need a spatial dimension" {
    const t = sig.Role.temporal;
    const two_times_no_space: sig.Signature = .{ .roles = &.{ t, t } };
    const wit = (try Order.searchAntisymmetryViolation(two_times_no_space)).?;
    const o = try Order.init(two_times_no_space, 0);
    const u = wit.uSlice();
    const v = wit.vSlice();
    try std.testing.expect(o.leq(u, v));
    try std.testing.expect(o.leq(v, u));
    try std.testing.expect(v[1] != 0.0); // u ≠ v
    // The witness is timelike, not null: g(e_t1) = s_t.
    try std.testing.expect(o.f.classify(v) == .temporal);
}

test "null separation is symmetric, so it is not an order" {
    const o = try Order.init(sig.minkowski_1_1, 0);
    const u = [_]f64{ 0, 0 };
    const v = [_]f64{ 1, 1 }; // null: g = 1 − 1 = 0

    // unoriented relation: both directions, u ≠ v → no antisymmetry
    try std.testing.expect(o.separation(&u, &v));
    try std.testing.expect(o.separation(&v, &u));
    try std.testing.expect(o.nullSeparated(&u, &v));
    try std.testing.expect(u[0] != v[0]);

    // oriented relation: v ⪯ u does NOT hold (the time arrow decreases)
    try std.testing.expect(o.leq(&u, &v));
    try std.testing.expect(!o.leq(&v, &u));

    // on the quotient by null separation the order is strict
    const x = [_]f64{ 2, 0 };
    try std.testing.expect(o.leq(&u, &x));
    try std.testing.expect(!o.leq(&x, &u));
}

test "the null vector generator really produces g = 0" {
    var prng = std.Random.DefaultPrng.init(4242);
    const rnd = prng.random();
    const o = try Order.init(sig.minkowski_3_1, 0);
    var v: [MAX_DIM]f64 = undefined;
    var hits: usize = 0;
    for (0..1000) |_| {
        if (o.randomConeVec(&v, rnd, .null_future)) {
            try std.testing.expectApproxEqAbs(@as(f64, 0.0), o.norm2(v[0..4]), 1e-9);
            hits += 1;
        }
    }
    try std.testing.expect(hits > 900);
}

test "the timelike vector generator works in both conventions" {
    var prng = std.Random.DefaultPrng.init(555);
    const rnd = prng.random();
    const cases = [_]sig.Signature{
        sig.minkowski_3_1,
        sig.minkowski_3_1_flipped,
        sig.two_times_2_1,
    };
    for (cases) |s| {
        const o = try Order.init(s, 0);
        var v: [MAX_DIM]f64 = undefined;
        var hits: usize = 0;
        for (0..500) |_| {
            if (o.randomConeVec(&v, rnd, .timelike_future)) {
                try std.testing.expect(o.inFutureCone(v[0..s.n()]));
                try std.testing.expectEqual(form.Class.temporal, o.f.classify(v[0..s.n()]));
                hits += 1;
            }
        }
        try std.testing.expect(hits > 450);
    }
}

test "Theorem 5.3: r > 0 degenerates causality" {
    const o = try Order.init(sig.degenerate_2_1_1, 0);

    const rad = [_]f64{ 0, 0, 0, 1 }; // radical vector
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), o.norm2(&rad), 1e-15);
    try std.testing.expect(o.nullSeparated(&rad, &[_]f64{ 0, 0, 0, 0 }));

    // causality stops telling points apart: 0 ⪯ rad ⪯ 0 with rad ≠ 0
    const zero = [_]f64{ 0, 0, 0, 0 };
    try std.testing.expect(o.leq(&zero, &rad));
    try std.testing.expect(o.leq(&rad, &zero));

    var out: [4]f64 = undefined;
    try std.testing.expectError(error.DegenerateForm, o.f.raiseIndex(&rad, &out));
}

test "the time arrow must be a temporal dimension" {
    try std.testing.expectError(
        error.TimeArrowNotTemporal,
        Order.init(sig.minkowski_3_1, 1),
    );
    try std.testing.expectError(
        error.TimeArrowOutOfRange,
        Order.init(sig.minkowski_3_1, 9),
    );
}

test "an Order survives the frame it was built from" {
    const makeOrderInLocalFrame = struct {
        fn call() OrderError!Order {
            // The signature borrows THIS local array. Order must copy the roles.
            var roles = [_]sig.Role{ .temporal, .spatial, .spatial, .spatial };
            const s = Signature{ .roles = &roles, .time_sign = .mostly_minus };
            return Order.init(s, 0);
        }
    }.call;

    var o = try makeOrderInLocalFrame();
    try std.testing.expectEqual(@as(usize, 4), o.dim());
    try std.testing.expect(o.signature().isLorentzian());
    try std.testing.expectEqual(@as(usize, 1), o.signature().p());

    const v = [_]f64{ 2, 0, 0, 0 };
    try std.testing.expectEqual(form.Class.temporal, o.f.classify(&v));
    try std.testing.expect(o.inFutureCone(&v));

    // the whole causal layer must work on a copied signature
    var prng = std.Random.DefaultPrng.init(11);
    const res = probeTransitivity(o, 2000, prng.random());
    try std.testing.expectEqual(@as(usize, 0), res.violations);
    try std.testing.expect(res.checked > 1500);
}

test "COUNTEREXAMPLE: the tolerance called it null, the computation decides it" {
    // This is the behaviour change of 0.1.6, stated as the smallest case that
    // shows it. v = (1, 1−1e-10, 0, 0) has g = 2e-10 − 1e-20 ≈ 2e-10, computed
    // to a radius of order 1e-16.
    const o = try Order.init(sig.minkowski_3_1, 0);
    var v = [_]f64{ 1.0, 1.0 - 1e-10, 0, 0 };

    // The old rule: |g| <= 1e-9, so the vector is reported as NULL.
    try std.testing.expectEqual(form.Class.null_like, o.f.classifyTol(&v, 1e-9));
    // The computation: 2e-10 against a radius of 2e-16 — decided, timelike.
    try std.testing.expectEqual(form.Class.temporal, o.classifyDecided(&v));
    // ... and the relation follows the sharper classification, not the constant:
    // a null vector is in its own future cone, and so is this one, but they are
    // now distinguishable by `nullSeparated`.
    try std.testing.expect(o.inFutureCone(&v));
    try std.testing.expect(!o.nullSeparated(&[_]f64{ 0, 0, 0, 0 }, &v));

    // A genuinely null vector is still null: the radius covers a computed zero.
    const nullish = [_]f64{ 1.0, 1.0, 0, 0 };
    try std.testing.expectEqual(form.Class.null_like, o.classifyDecided(&nullish));
    try std.testing.expect(o.nullSeparated(&[_]f64{ 0, 0, 0, 0 }, &nullish));
}
