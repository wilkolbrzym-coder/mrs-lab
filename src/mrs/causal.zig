//! MRS-LAB :: causality module
//!
//! Two relations that must be told apart — and conflating them was the most
//! serious error in the original project specification:
//!
//!   * SEPARACJA (niezorientowana).  u ⊑ v  ⟺  g(v−u, v−u) ≤ 0.
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
//! Twierdzenia 5.2.
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
//! lorentzowskich).
//!   Reflexivity: always.
//!   Transitivity: ⟺ p = 1 (with r = 0).
//!   Antysymetria: ⟺ r = 0 i p ≤ 1.
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
    ///     /// Class `temporal`, future directed.
    timelike_future,
    ///     /// Class `null_like`, future directed.
    null_future,
};

/// What goes wrong in a given signature.
pub const Causality = struct {
    reflexive: bool,
    transitive: bool,
    antisymmetric: bool,

    pub fn isPartialOrder(self: Causality) bool {
        return self.reflexive and self.transitive and self.antisymmetric;
    }

    pub fn label(self: Causality) []const u8 {
        if (self.isPartialOrder()) return "partial order";
        if (self.transitive) return "preorder (no antisymmetry)";
        return "neither transitive nor antisymmetric";
    }
};

pub const Order = struct {
    f: DiagonalForm,
    /// Tolerancja rozpoznania wektora zerowego. W arytmetyce zmiennoprzecinkowej
    ///     /// g(v,v) = 0 is realised as |g| ~ eps, so "null" MUST be a notion with
    ///     /// a tolerance. In MRS-LAB the tolerance is an explicit parameter of the
    ///     /// type, not a magic constant scattered through the code.
    tol: f64 = 1e-9,
    ///     /// Index of the dimension chosen as the time arrow. For p >= 2 this is
    ///     /// EXTRA structure that the form itself does not determine.
    time_arrow: usize,

    pub fn init(s: Signature, time_arrow: usize) OrderError!Order {
        if (s.n() > MAX_DIM) return error.DimensionTooLarge;
        if (time_arrow >= s.n()) return error.TimeArrowOutOfRange;
        if (s.roles[time_arrow] != .temporal) return error.TimeArrowNotTemporal;
        return .{ .f = DiagonalForm.init(s), .time_arrow = time_arrow };
    }

    pub fn dim(self: Order) usize {
        return self.f.n();
    }

    pub fn signature(self: Order) Signature {
        return self.f.signature;
    }

    pub fn norm2(self: Order, v: []const f64) f64 {
        return self.f.eval(v);
    }

    ///     /// Closed future cone. The sign convention enters EXCLUSIVELY through
    ///     /// `classify`, which is why this works for (+,−,−,−) as well,
    ///     and for (−,+,+,+).
    pub fn inFutureCone(self: Order, v: []const f64) bool {
        if (self.f.classifyTol(v, self.tol) == .spatial) return false;
        return v[self.time_arrow] >= 0.0;
    }

    /// Oriented relation: u ⪯ v.
    pub fn leq(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.inFutureCone(d[0..u.len]);
    }

    /// Separacja niezorientowana: g(v−u, v−u) ≤ 0.
    ///     /// This is NOT the same relation as ⪯ — it is symmetric on the null
    ///     /// cone, so it is not an order.
    pub fn separation(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.f.classifyTol(d[0..u.len], self.tol) != .spatial;
    }

    ///     /// Null separation: g(v−u, v−u) = 0 — an equivalence relation whose
    ///     /// classes are the light rays.
    pub fn nullSeparated(self: Order, u: []const f64, v: []const f64) bool {
        var d: [MAX_DIM]f64 = undefined;
        for (0..u.len) |i| d[i] = v[i] - u[i];
        return self.f.classifyTol(d[0..u.len], self.tol) == .null_like;
    }

    // // -- cone vector generators ----------------------------------------------

    ///     /// Random vector in the future cone. The construction is uniform across
    /// konwencji znaku:
    ///     ///   1. draw the components outside the time arrow,
    ///     ///   2. scale the spatial group so that its contribution to g is
    ///     ///      Sum_temporal_rest x² + 1 — then the contribution of all
    ///     ///      dimensions outside the arrow is exactly −s_arrow,
    ///     ///   3. g = s_arrow(t² − 1), so t = 1 gives a null vector and t = 2
    ///     ///      timelike vector — in BOTH conventions, because "temporal" means
    ///     ///      "the temporal component dominates", not "g has a given sign".
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

        if (q_space <= 1e-12) return false; //         if (q_space <= 1e-12) return false; // no spatial dimension = no cone
        const f = @sqrt((q_time_rest + 1.0) / q_space);
        for (0..nn) |i| {
            if (s.roles[i] == .spatial) out[i] *= f;
        }

        var other: f64 = 0;
        for (0..nn) |i| {
            if (i == self.time_arrow) continue;
            other += s.signAt(i) * out[i] * out[i];
        }
        if (@abs(other + s_arrow) > 1e-9) return false; //         if (@abs(other + s_arrow) > 1e-9) return false; // construction failed

        out[self.time_arrow] = switch (mode) {
            .null_future => 1.0,
            .timelike_future => 2.0,
        };

        const cls = self.f.classifyTol(out[0..nn], self.tol);
        return switch (mode) {
            .null_future => cls == .null_like,
            .timelike_future => cls == .temporal,
        };
    }

    ///     /// Random vector in the cone with a random mode (timelike or null).
    ///     /// Counterexamples to transitivity for p >= 2 require NULL vectors
    ///     /// (proof: the sum of two timelike vectors with the same arrow sense
    ///     /// stays in the cone), so the generator must be able to produce both.
    pub fn randomConeVecAny(self: Order, out: *[MAX_DIM]f64, rnd: std.Random) bool {
        const mode: ConeMode = if (rnd.float(f64) < 0.5) .timelike_future else .null_future;
        return self.randomConeVec(out, rnd, mode);
    }

    // // -- witness search ------------------------------------------------------

    ///     /// Transitivity violation ⪯: u ⪯ w, w ⪯ v, but u ⋠ v.
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

    ///     /// Is the future cone convex (sampling)?
    ///     /// Convexity ⟺ transitivity of ⪯.
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

    ///     /// Constructive transitivity witness for p >= 2.
    /// Wymiary: t0 i t1 czasowe, s przestrzenny.
    ///   a   = e_t0 + e_s      → g = s_t − s_s = 0         (zerowy)
    ///   v−a = e_t1 + e_s      → g = 0                      (zerowy)
    ///   v   = e_t0 + e_t1 + 2·e_s → g = 2·s_t − 4·s_s < 0  (przestrzenny)
    ///     /// The inequality "2·s_t < 4·s_s" holds in both conventions,
    /// bo s_s = −s_t.
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

    ///     /// Witness of missing antisymmetry of ⪯: u ≠ v with u ⪯ v and v ⪯ u.
    ///     /// It exists exactly when r > 0 or p >= 2 (Theorem 5.2).
    ///     /// Returns `null` for Lorentzian signatures — and that is consistent
    ///     /// with the theorem, not a coincidence.
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
        if (nt >= 2 and idx_s != std.math.maxInt(usize)) {
            // d: arrow 0, second temporal component 1, space 0.5
            // → g = s_t·(1 − 0.25), temporal class in both conventions
            wit.v[idx_t[1]] = 1.0;
            wit.v[idx_s] = 0.5;
            return wit;
        }
        return null;
    }

    ///     /// Full diagnosis: is ⪯ a partial order, and what exactly fails.
    /// zawodzi. Weryfikacja empiryczna Twierdzenia 5.2.
    pub fn diagnose(self: Order, trials: usize, rnd: std.Random) OrderError!Causality {
        const trans_viol = self.searchTransitivityViolation(trials, rnd) != null;
        const anti_viol = (try Order.searchAntisymmetryViolation(self.signature())) != null;
        return .{
            .reflexive = true,
            .transitive = !trans_viol,
            .antisymmetric = !anti_viol,
        };
    }
};

/// Eksperymentalna weryfikacja Twierdzenia 5.1.
pub const Transitivity = struct {
    violations: usize,
    trials: usize,
    checked: usize,
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
        .convex = violations == 0,
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
    try std.testing.expect(o.leq(&zero, &zero)); //     try std.testing.expect(o.leq(&zero, &zero)); // reflexivity

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

    // konkretne liczby: a zerowy, v−a zerowy, v przestrzenny
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
        try std.testing.expect(v[3] != 0.0); //     try std.testing.expect(v[3] != 0.0); // radical vector, u ≠ v
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
    const cases = [_]sig.Signature{
        sig.minkowski_3_1,
        sig.minkowski_3_1_flipped,
        sig.minkowski_1_1,
        sig.two_times_2_1,
        sig.degenerate_2_1_1,
    };
    for (cases) |s| {
        const o = try Order.init(s, 0);
        const diag = try o.diagnose(4000, prng.random());
        try std.testing.expectEqual(s.isLorentzian(), diag.isPartialOrder());
    }
}

test "null separation is symmetric, so it is not an order" {
    const o = try Order.init(sig.minkowski_1_1, 0);
    const u = [_]f64{ 0, 0 };
    const v = [_]f64{ 1, 1 }; // zerowy: g = 1 − 1 = 0

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

    const rad = [_]f64{ 0, 0, 0, 1 }; //     const rad = [_]f64{ 0, 0, 0, 1 }; // radical vector
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
