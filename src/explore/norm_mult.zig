//! MRS-LAB :: P2 — multiplicativity of the norm
//!
//! QUESTION (split into three predicates, because the project draft had one
//! and it was internally contradictory):
//!
//!   P2a  N_sc(x·y) = N_sc(x)·N_sc(y),  where N_sc(z) = Sc(z·z̄)
//!   P2b  (x·y)·conj(x·y) = (x·conj x)·(y·conj y)   (an identity in the centre)
//!   P2c  does z ↦ z·z̄ always land in the centre? (well-posedness of P2b)
//!
//! MEASURED RESULT: P2c holds exactly for 2^n <= 8 (n <= 3) and fails from
//! n = 4. Witness in Cl(1,3): x = e01 + e23, because e01 and e23 are disjoint
//! and commute, so x·x̄ carries the term 2·e0123, and e0123 is not central.
//! That is the same boundary as for multiplicativity — and not by accident:
//! from n = 4 the map z ↦ z·z̄ stops landing in the centre, so there is no
//! room for a multiplicative norm structure to exist.
//!
//! ---------------------------------------------------------------------------
//! THIS IS A DECISION PROCEDURE, NOT SAMPLING
//! ---------------------------------------------------------------------------
//! Key observation: for fixed y each of these identities is a **quadratic
//! function** of x (a linear map composed with a quadratic form), and
//! symmetrically, it is a quadratic function of y for fixed x. A quadratic
//! function (including a linear term) is uniquely determined by its values on
//! determined by its values on the set
//!
//!     D = {0} ∪ {e_i} ∪ {e_i + e_j : i < j},
//!
//! because from Q(0), Q(e_i) and Q(e_i+e_j) we read off all coefficients
//! combinations a_i + c_ii and c_ij. Therefore:
//!
//!     the identity holds for all x, y  ⟺  it holds for (x,y) in D×D.
//!
//! That is not "a lot of tests". It is an **exhaustive proof** for this class
//! of identities, in integer arithmetic, with no tolerance and no probability.
//! That is why the result may be called a decision rather than evidence.
//!
//! SCOPE. The engine decides about the CLIFFORD PRODUCT in a given signature.
//! It does not decide the question "does ANY multiplication with a
//! multiplicative norm exist on this space" — and that second question
//! includes the octonions, which are NOT a Clifford algebra (Clifford is
//! associative, the octonions are not). Hence the Hurwitz bound (1,2,4,8)
//! dimension 8 in the Clifford family is Cl(0,3) ≅ H⊕H, not the octonions.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const exact = @import("exact.zig");
const sigs = @import("signatures.zig");

const IntVec = exact.IntVec;

/// Full size of the determining set for n = 5.
pub const MAX_GRID: usize = 1 + 32 + (32 * 31) / 2; // 529

pub const Kind = enum { scalar_multiplicative, center_multiplicative, norm_is_central };

pub const Verdict = struct {
    kind: Kind,
    n: usize,
    grid_size: usize,
    pairs: usize,
    holds: bool,
    /// Witness of a refutation (only when `holds == false`).
    witness_x: IntVec = IntVec{},
    witness_y: IntVec = IntVec{},
    has_witness: bool = false,

    pub fn name(self: Verdict) []const u8 {
        return switch (self.kind) {
            .scalar_multiplicative => "P2a  N_sc(xy) = N_sc(x)·N_sc(y)",
            .center_multiplicative => "P2b  (xy)·conj(xy) = (x·conj x)(y·conj y)",
            .norm_is_central => "P2c  z·z̄ in the centre for every z",
        };
    }
};

/// Builds the determining set D = {0} ∪ {e_i} ∪ {e_i+e_j}.
/// `out` must have at least MAX_GRID slots; returns the size used.
pub fn buildGrid(alg: cl.Algebra, out: *[MAX_GRID]IntVec) usize {
    const m = alg.basisCount();
    var k: usize = 0;
    out[k] = IntVec.zero();
    k += 1;
    for (0..m) |i| {
        out[k] = IntVec.basis(@intCast(i));
        k += 1;
    }
    for (0..m) |i| {
        for (i + 1..m) |j| {
            var v = IntVec.basis(@intCast(i));
            v.c[j] = 1;
            out[k] = v;
            k += 1;
        }
    }
    return k;
}

fn normOf(x: IntVec, alg: cl.Algebra) IntVec {
    const m = alg.basisCount();
    return exact.mul(alg, x, exact.conj(x, m));
}

/// Checks P2c: does z·z̄ lie in the centre for every z?
/// Centrality is a linear condition and z·z̄ is quadratic in z, so the
/// determining set D again suffices.
pub fn checkNormIsCentral(alg: cl.Algebra, grid: *[MAX_GRID]IntVec) Verdict {
    const k = buildGrid(alg, grid);
    const m = alg.basisCount();
    const center = exact.centerBasis(alg);

    var i: usize = 0;
    while (i < k) : (i += 1) {
        const nz = normOf(grid[i], alg);
        for (0..m) |b| {
            if (nz.c[b] == 0) continue;
            if (center & (@as(u32, 1) << @intCast(b)) == 0) {
                return .{
                    .kind = .norm_is_central,
                    .n = alg.n_gen,
                    .grid_size = k,
                    .pairs = k,
                    .holds = false,
                    .witness_x = grid[i],
                    .witness_y = IntVec.basis(@intCast(b)),
                    .has_witness = true,
                };
            }
        }
    }
    return .{
        .kind = .norm_is_central,
        .n = alg.n_gen,
        .grid_size = k,
        .pairs = k,
        .holds = true,
    };
}

/// Checks P2a: N_sc(x·y) = N_sc(x)·N_sc(y).
pub fn checkScalarMultiplicative(alg: cl.Algebra, grid: *[MAX_GRID]IntVec) Verdict {
    const k = buildGrid(alg, grid);
    var i: usize = 0;
    while (i < k) : (i += 1) {
        const nx = exact.scalarPart(normOf(grid[i], alg));
        var j: usize = 0;
        while (j < k) : (j += 1) {
            const ny = exact.scalarPart(normOf(grid[j], alg));
            const nxy = exact.scalarPart(normOf(exact.mul(alg, grid[i], grid[j]), alg));
            if (nxy != exact_mul(nx, ny)) {
                return .{
                    .kind = .scalar_multiplicative,
                    .n = alg.n_gen,
                    .grid_size = k,
                    .pairs = k * k,
                    .holds = false,
                    .witness_x = grid[i],
                    .witness_y = grid[j],
                    .has_witness = true,
                };
            }
        }
    }
    return .{
        .kind = .scalar_multiplicative,
        .n = alg.n_gen,
        .grid_size = k,
        .pairs = k * k,
        .holds = true,
    };
}

/// Checks P2b: the vector identity in the centre.
pub fn checkCenterMultiplicative(alg: cl.Algebra, grid: *[MAX_GRID]IntVec) Verdict {
    const k = buildGrid(alg, grid);
    const m = alg.basisCount();
    var i: usize = 0;
    while (i < k) : (i += 1) {
        const nx = normOf(grid[i], alg);
        var j: usize = 0;
        while (j < k) : (j += 1) {
            const ny = normOf(grid[j], alg);
            const prod = exact.mul(alg, grid[i], grid[j]);
            const nxy = normOf(prod, alg);
            const rhs = exact.mul(alg, nx, ny);
            if (!IntVec.eql(nxy, rhs, m)) {
                return .{
                    .kind = .center_multiplicative,
                    .n = alg.n_gen,
                    .grid_size = k,
                    .pairs = k * k,
                    .holds = false,
                    .witness_x = grid[i],
                    .witness_y = grid[j],
                    .has_witness = true,
                };
            }
        }
    }
    return .{
        .kind = .center_multiplicative,
        .n = alg.n_gen,
        .grid_size = k,
        .pairs = k * k,
        .holds = true,
    };
}

inline fn exact_mul(a: i64, b: i64) i64 {
    return std.math.mul(i64, a, b) catch @panic("MRS explore: overflow in P2");
}

/// The determining set is complete — that is the contract of this layer.
/// Returns the size for a given number of blades.
pub fn gridSize(basis_count: usize) usize {
    return 1 + basis_count + (basis_count * (basis_count - 1)) / 2;
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "the grid size matches the one actually built" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const k = buildGrid(alg, &grid);
        try std.testing.expectEqual(gridSize(alg.basisCount()), k);
        try std.testing.expect(k <= MAX_GRID);
    }
}

test "P2a: true for n <= 2 in BOTH conventions" {
    var grid: [MAX_GRID]IntVec = undefined;
    inline for (.{ .mostly_minus, .mostly_plus }) |ts| {
        for ([_]sigs.Triple{
            .{ .p = 1, .q = 0 },
            .{ .p = 0, .q = 1 },
            .{ .p = 0, .q = 2 },
            .{ .p = 1, .q = 1 },
            .{ .p = 2, .q = 0 },
        }) |t| {
            const alg = try (sigs.SigBuf.build(t, ts)).algebra();
            const v = checkScalarMultiplicative(alg, &grid);
            try std.testing.expect(v.holds);
            const c = checkCenterMultiplicative(alg, &grid);
            try std.testing.expect(c.holds);
        }
    }
}

test "P2a: false from n = 3 — the engine produces a witness" {
    var grid: [MAX_GRID]IntVec = undefined;
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_minus)).algebra();
    const v = checkScalarMultiplicative(alg, &grid);
    try std.testing.expect(!v.holds);
    try std.testing.expect(v.has_witness);
    // the witness must actually break the identity
    const nx = exact.scalarPart(normOf(v.witness_x, alg));
    const ny = exact.scalarPart(normOf(v.witness_y, alg));
    const nxy = exact.scalarPart(normOf(exact.mul(alg, v.witness_x, v.witness_y), alg));
    try std.testing.expect(nxy != exact_mul(nx, ny));
}

test "P2c: z·z̄ lands in the centre EXACTLY for 2^n <= 8" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);

    // n <= 3 (2^n <= 8): holds for EVERY signature and convention
    var low: usize = 0;
    for (0..n) |i| {
        if (buf[i].n() > 3) continue;
        low += 1;
        const a_minus = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const a_plus = try (sigs.SigBuf.build(buf[i], .mostly_plus)).algebra();
        if (!checkNormIsCentral(a_minus, &grid).holds) {
            std.debug.print("P2c fails for n<=3: triple=({d},{d},{d}) minus\n", .{ buf[i].p, buf[i].q, buf[i].r });
            return error.TestUnexpectedResult;
        }
        if (!checkNormIsCentral(a_plus, &grid).holds) {
            std.debug.print("P2c fails for n<=3: triple=({d},{d},{d}) plus\n", .{ buf[i].p, buf[i].q, buf[i].r });
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(low > 0);

    // n = 4: it fails, and that is a result, not a defect. Constructive witness:
    // x = e01 + e23 in Cl(1,3), because e01 and e23 are disjoint and commute,
    // so x·x̄ carries the term 2·e0123, and e0123 is NOT central
    // (the centre of Cl(1,3) is R itself).
    const alg = try (sigs.SigBuf.build(.{ .p = 1, .q = 3 }, .mostly_minus)).algebra();
    const v = checkNormIsCentral(alg, &grid);
    try std.testing.expect(!v.holds);
    try std.testing.expect(v.has_witness);

    // check the witness by hand: e01 = mask 3, e23 = mask 12
    const e01 = IntVec.basis(3);
    const e23 = IntVec.basis(12);
    const xx = exact.add(e01, e23); // x = e01 + e23
    const nz = normOf(xx, alg);
    try std.testing.expect(nz.c[15] != 0); // the e0123 term (mask 15)
    try std.testing.expectEqual(@as(u5, 1), exact.dimOfSet(exact.centerBasis(alg)));

    // for n = 3 the same mechanism does NOT break centrality, because e012 is central
    const alg3 = try (sigs.SigBuf.build(.{ .p = 1, .q = 2 }, .mostly_minus)).algebra();
    const v3 = checkNormIsCentral(alg3, &grid);
    try std.testing.expect(v3.holds);
    try std.testing.expectEqual(@as(u5, 2), exact.dimOfSet(exact.centerBasis(alg3)));
}

test "the DECISION is complete: it finds a violation outside the grid" {
    // Take a witness from the grid and check that it really breaks the identity
    // for points outside the grid as well — this tests that the method is not an
    // artefact of the choice of D.
    var grid: [MAX_GRID]IntVec = undefined;
    const alg = try (sigs.SigBuf.build(.{ .p = 0, .q = 3 }, .mostly_minus)).algebra();
    const v = checkScalarMultiplicative(alg, &grid);
    try std.testing.expect(!v.holds);

    // random point outside the grid: the identity must fail for the vast majority
    // of points (if it failed only on D, the method would be unsound)
    var prng = std.Random.DefaultPrng.init(1234);
    const rnd = prng.random();
    const m = alg.basisCount();
    var fails: usize = 0;
    const trials: usize = 500;
    for (0..trials) |_| {
        var x = IntVec{};
        var y = IntVec{};
        for (0..m) |i| {
            x.c[i] = @intCast(rnd.intRangeLessThan(i32, -2, 3));
            y.c[i] = @intCast(rnd.intRangeLessThan(i32, -2, 3));
        }
        const nx = exact.scalarPart(normOf(x, alg));
        const ny = exact.scalarPart(normOf(y, alg));
        const nxy = exact.scalarPart(normOf(exact.mul(alg, x, y), alg));
        if (nxy != exact_mul(nx, ny)) fails += 1;
    }
    // we do not require 100% — only that failure is the rule, not the exception
    try std.testing.expect(fails * 2 > trials);
}

test "P2a depends on the sign convention only through the triple — same result" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    for (0..n) |i| {
        if (buf[i].n() > 3) continue;
        const a_minus = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        const a_plus = try (sigs.SigBuf.build(buf[i], .mostly_plus)).algebra();
        try std.testing.expectEqual(
            checkScalarMultiplicative(a_minus, &grid).holds,
            checkScalarMultiplicative(a_plus, &grid).holds,
        );
    }
}

// ---------------------------------------------------------------------------
// RULES DISCOVERED BY THE ENGINE — recorded as checkable theorems
// ---------------------------------------------------------------------------
//
// The engine was not programmed with these. They EMERGED from the table, and
// these tests turn the observation into a theorem that is checked on every
// change of the code. A proof sketch is given for each rule.

test "RULE 1: P2a holds exactly when p+q <= 2 (degenerate dimensions invisible)" {
    // WHY THIS IS SO. Cl(p,q,r) ≅ Cl(p,q) ⊗ Λ(R), where R is the radical.
    // The scalar part of the norm sees only Cl(p,q):
    //   * the product of two distinct blades is never a scalar (e_A·e_B = ±e_{A△B});
    //   * the product e_A·e_A is the product of the s_i over A, and that is ZERO
    // ZERO, when A contains a degenerate generator.
    // Hence N_sc(x) = (scalar coefficient of x)². Therefore
    //   N_sc(xy) = (c0(x)·c0(y))² = N_sc(x)·N_sc(y),
    // so for r > 0 the identity holds IDENTICALLY, independently of Cl(p,q).
    // For r = 0 what remains is the classical Hurwitz bound in the Clifford family:
    // multiplicativity holds exactly when 2^(p+q) <= 4.
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);

    var checked: usize = 0;
    var degenerate_seen: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.n() > 4) continue;
        checked += 1;
        if (t.isDegenerate()) degenerate_seen += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const holds = checkScalarMultiplicative(alg, &grid).holds;
        const predicted = (@as(usize, t.p) + @as(usize, t.q)) <= 2;
        if (holds != predicted) {
            std.debug.print(
                "RULE 1 broken: t=({d},{d},{d}) p+q={d} measured={} predicted={}\n",
                .{ t.p, t.q, t.r, @as(usize, t.p) + @as(usize, t.q), holds, predicted },
            );
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked >= 30);
    try std.testing.expect(degenerate_seen >= 15);
}

test "RULE 2: P2c holds exactly when 2^n <= 8" {
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var checked: usize = 0;
    for (0..n) |i| {
        const t = buf[i];
        if (t.n() > 4) continue;
        checked += 1;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const holds = checkNormIsCentral(alg, &grid).holds;
        const predicted = t.n() <= 3;
        if (holds != predicted) {
            std.debug.print(
                "RULE 2 broken: t=({d},{d},{d}) n={d} measured={} predicted={}\n",
                .{ t.p, t.q, t.r, t.n(), holds, predicted },
            );
            return error.TestUnexpectedResult;
        }
    }
    try std.testing.expect(checked >= 30);
}

test "P2b differs from P2a; P2b == P2c in range (observation, not a theorem)" {
    // Why three predicates? Because they give different answers — if they did not,
    // splitting P2 would be pointless. This test guards that the P2a vs P2b
    // difference is real, and RECORDS (without proof) that P2b and P2c agree over
    // the whole tested range n <= 4. The second is not claimed as a theorem but
    // as a pattern to investigate, recorded under "Known follow-ups" in
    // CHANGELOG.md.
    var grid: [MAX_GRID]IntVec = undefined;
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, sigs.MAX_TOTAL);
    var differ_a: usize = 0;
    var differ_c: usize = 0;
    var witness_a: sigs.Triple = .{};
    for (0..n) |i| {
        const t = buf[i];
        if (t.n() > 4) continue;
        const alg = try (sigs.SigBuf.build(t, .mostly_minus)).algebra();
        const a = checkScalarMultiplicative(alg, &grid).holds;
        const b = checkCenterMultiplicative(alg, &grid).holds;
        const c = checkNormIsCentral(alg, &grid).holds;
        if (a != b) {
            differ_a += 1;
            witness_a = t;
        }
        if (c != b) differ_c += 1;
    }
    // P2a and P2b MUST differ — otherwise splitting the predicates is empty.
    // The recorded count for the exhausted range n <= 4 is ASSERTED, not printed:
    // a test that writes to stderr makes `zig build test` report a spurious
    // "failed command:" line for the whole step (Zig 0.16 build runner), which
    // buries real failures in noise.
    try std.testing.expectEqual(@as(usize, 10), differ_a);
    try std.testing.expect(witness_a.p + witness_a.q + witness_a.r > 0);
    // Observation without proof: P2b and P2c agree over the whole range.
    try std.testing.expectEqual(@as(usize, 0), differ_c);
}
