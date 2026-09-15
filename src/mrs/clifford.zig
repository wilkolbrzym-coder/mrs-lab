//! MRS-0 :: Algebra Clifforda Cl(p,q) w reprezentacji bitmaskowej
//!
//! A multivector is a vector of coefficients over the basis of blades.
//! A blade is a subset of generators encoded as a bitmask: bit i means e_i.
//! This is exactly the representation the signature forces: the sign
//! e_i² = s_i is data, not the result of a matrix multiplication.
//!
//! Reordering rule (canonical sign). blade(mask) is the product
//! e_{i1}·e_{i2}·…·e_{ik} for i1 < i2 < … < ik. The product blade(mask)·e_i
//! is computed by appending e_i on the right and moving it leftwards:
//!
//!   * if bit i is already in the mask — e_i meets itself, gives s_i
//!     and the bit clears;
//!   * every generator with index greater than i that e_i passes
//!     picks up a sign of −1.
//!
//! Cost of one blade multiplication: O(number of generators), no allocation.
//!
//! Asymptotics of a multivector: a dense multivector has 2^n coefficients, so
//! a naive product costs 4^n. The product of multivectors with k nonzero
//! coefficients costs k² (up to the cost of merging), and that is thesis T3,
//! measured in `results/RESULTS.md`.

const std = @import("std");
const sig = @import("signature.zig");

const Signature = sig.Signature;

/// Maximum number of generators. Bounded deterministically, because the dense
/// representation has 2^n coefficients (n = 20 is already 1 048 576).
pub const MAX_GEN: u5 = 20;
pub const MIN_GEN: u5 = 1;

pub const AlgebraError = error{
    TooManyGenerators,
    NoGenerators,
};

pub const Algebra = struct {
    n_gen: u5,
    /// squares[i] = e_i² ∈ {+1, −1, 0}
    squares: [MAX_GEN]i8,

    pub fn init(n_gen: u5, squares: [MAX_GEN]i8) AlgebraError!Algebra {
        if (n_gen < MIN_GEN) return error.NoGenerators;
        if (n_gen > MAX_GEN) return error.TooManyGenerators;
        return .{ .n_gen = n_gen, .squares = squares };
    }

    /// Builds the algebra Cl from a signature: e_i² = s_i. Degenerate dimensions
    /// (s_i = 0) give nilpotent generators, i.e. a degenerate algebra — Cl(p,q,r),
    /// which is the quotient of Cl(p,q) by a nilpotent ideal.
    pub fn fromSignature(s: Signature) AlgebraError!Algebra {
        var sq = [_]i8{0} ** MAX_GEN;
        const nn = s.n();
        if (nn > MAX_GEN) return error.TooManyGenerators;
        for (0..nn) |i| {
            sq[i] = @intFromFloat(s.signAt(i));
        }
        return init(@intCast(nn), sq);
    }

    pub fn basisCount(self: Algebra) usize {
        return @as(usize, 1) << self.n_gen;
    }

    pub fn grade(self: Algebra, mask: u32) u32 {
        _ = self;
        return @popCount(mask);
    }

    /// Reverse sign of a multivector: (−1)^{k(k−1)/2}.
    /// The arithmetic wraps so that k = 0 does not underflow.
    pub fn reverseSign(self: Algebra, mask: u32) i8 {
        const k: u32 = self.grade(mask);
        const half = (k *% (k -% 1)) / 2;
        return if (half % 2 == 0) 1 else -1;
    }

    pub fn writeBlade(self: Algebra, w: anytype, mask: u32) !void {
        if (mask == 0) {
            try w.writeAll("1");
            return;
        }
        var first = true;
        for (0..self.n_gen) |i| {
            if (mask & (@as(u32, 1) << @intCast(i)) != 0) {
                if (!first) try w.writeAll("*");
                try w.print("e{d}", .{i});
                first = false;
            }
        }
    }
};

/// Product of two blades: e_A · e_B = sign · e_mask.
pub const BladeProd = struct {
    sign: i8,
    mask: u32,
};

/// Blade multiplication. Complexity O(popcount(B)) — at most n steps in practice,
/// without allocation and without 2^n × 2^n tables.
pub fn bladeMul(alg: Algebra, a_mask: u32, b_mask: u32) BladeProd {
    var sign: i8 = 1;
    var mask = a_mask;
    var b = b_mask;
    while (b != 0) {
        const i: u5 = @intCast(@ctz(b));
        b &= b - 1; // iterate over the generators of B in increasing order
        // generators with index > i that e_i passes
        const above: u32 = mask >> (i + 1);
        if (@popCount(above) & 1 == 1) sign = -sign;
        const bit: u32 = @as(u32, 1) << i;
        if (mask & bit != 0) {
            sign *= alg.squares[i];
            mask ^= bit;
        } else {
            mask |= bit;
        }
    }
    return .{ .sign = sign, .mask = mask };
}

// ---------------------------------------------------------------------------
// Sparse representation
// ---------------------------------------------------------------------------

pub const Term = struct {
    mask: u32,
    coeff: f64,
};

fn termLessThan(_: void, a: Term, b: Term) bool {
    return a.mask < b.mask;
}

/// Sparse multivector: a sorted list of (mask, coefficient), masks unique.
pub const Sparse = struct {
    /// A view over the coefficients. A mutable buffer can be passed as `[]Term`,
    /// but the multivector itself does not need write access — hence `const`.
    terms: []const Term,

    pub fn nnz(self: Sparse) usize {
        return self.terms.len;
    }

    pub fn deinit(self: Sparse, alloc: std.mem.Allocator) void {
        alloc.free(self.terms);
    }

    /// Coefficient at a mask (0 when absent).
    pub fn get(self: Sparse, mask: u32) f64 {
        for (self.terms) |t| {
            if (t.mask == mask) return t.coeff;
        }
        return 0.0;
    }

    pub fn writeTo(self: Sparse, alg: Algebra, w: anytype) !void {
        if (self.terms.len == 0) {
            try w.writeAll("0");
            return;
        }
        for (self.terms, 0..) |t, idx| {
            if (idx != 0) {
                if (t.coeff < 0) try w.writeAll(" - ") else try w.writeAll(" + ");
            }
            try w.print("{d:.6}*", .{@abs(t.coeff)});
            try alg.writeBlade(w, t.mask);
        }
    }
};

/// Geometric product of sparse multivectors.
/// Cost: O(nnz(a)·nnz(b)) blade operations + O(m log m) merging.
pub fn mulSparse(
    alloc: std.mem.Allocator,
    alg: Algebra,
    a: Sparse,
    b: Sparse,
) !Sparse {
    var acc = std.AutoHashMap(u32, f64).init(alloc);
    defer acc.deinit();

    for (a.terms) |ta| {
        for (b.terms) |tb| {
            const bp = bladeMul(alg, ta.mask, tb.mask);
            if (bp.sign == 0) continue; // nilpotent generator
            const contrib = ta.coeff * tb.coeff * @as(f64, @floatFromInt(bp.sign));
            const gop = try acc.getOrPut(bp.mask);
            if (!gop.found_existing) gop.value_ptr.* = 0.0;
            gop.value_ptr.* += contrib;
        }
    }

    const total: usize = @intCast(acc.count());
    var out = try alloc.alloc(Term, total);
    var idx: usize = 0;
    var it = acc.iterator();
    while (it.next()) |e| {
        if (e.value_ptr.* != 0.0) {
            out[idx] = .{ .mask = e.key_ptr.*, .coeff = e.value_ptr.* };
            idx += 1;
        }
    }
    const trimmed = try alloc.realloc(out, idx);
    std.mem.sort(Term, trimmed, {}, termLessThan);
    return .{ .terms = trimmed };
}

/// Reverse of a multivector.
pub fn reverseSparse(alloc: std.mem.Allocator, alg: Algebra, a: Sparse) !Sparse {
    const out = try alloc.alloc(Term, a.terms.len);
    for (a.terms, 0..) |t, i| {
        out[i] = .{
            .mask = t.mask,
            .coeff = t.coeff * @as(f64, @floatFromInt(alg.reverseSign(t.mask))),
        };
    }
    return .{ .terms = out };
}

/// Sparse product WITHOUT ALLOCATION in the hot path. Requires `scratch` of
/// length at least nnz(a)·nnz(b). Assumes both inputs are sorted by mask
/// posortowane po masce (co gwarantuje `mulSparse`).
///
/// This is the "maximally optimised" variant and it is the one used in
/// benchmark T3 — comparing against the allocating version would be unfair to MRS.
/// against MRS-LAB.
pub fn mulSparseScratch(alg: Algebra, a: Sparse, b: Sparse, scratch: []Term) Sparse {
    const need = a.terms.len * b.terms.len;
    std.debug.assert(scratch.len >= need);

    var len: usize = 0;
    for (a.terms) |ta| {
        for (b.terms) |tb| {
            const bp = bladeMul(alg, ta.mask, tb.mask);
            if (bp.sign == 0) continue;
            scratch[len] = .{
                .mask = bp.mask,
                .coeff = ta.coeff * tb.coeff * @as(f64, @floatFromInt(bp.sign)),
            };
            len += 1;
        }
    }

    std.mem.sort(Term, scratch[0..len], {}, termLessThan);

    // merge adjacent masks in place
    var out_len: usize = 0;
    var i: usize = 0;
    while (i < len) {
        const mask = scratch[i].mask;
        var acc: f64 = 0;
        while (i < len and scratch[i].mask == mask) : (i += 1) {
            acc += scratch[i].coeff;
        }
        if (acc != 0.0) {
            scratch[out_len] = .{ .mask = mask, .coeff = acc };
            out_len += 1;
        }
    }
    return .{ .terms = scratch[0..out_len] };
}

// ---------------------------------------------------------------------------
// Dense representation (reference)
// ---------------------------------------------------------------------------

/// Product in the dense representation: a loop over all blade pairs.
/// Cost 4^n — the price a representation pays when it does not know sparsity.
pub fn mulDense(
    alloc: std.mem.Allocator,
    alg: Algebra,
    a: []const f64,
    b: []const f64,
) ![]f64 {
    const m = alg.basisCount();
    std.debug.assert(a.len == m and b.len == m);
    const out = try alloc.alloc(f64, m);
    @memset(out, 0.0);
    for (0..m) |i| {
        if (a[i] == 0.0) continue;
        for (0..m) |j| {
            const bp = bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue;
            out[bp.mask] += a[i] * b[j] * @as(f64, @floatFromInt(bp.sign));
        }
    }
    return out;
}

/// Dense product with explicit zero skipping — the "obvious optimisation"
/// that anyone adds once they notice zeros in the data. Still O(4^n) checks.
pub fn mulDenseSkippingZeros(
    alloc: std.mem.Allocator,
    alg: Algebra,
    a: []const f64,
    b: []const f64,
) ![]f64 {
    const m = alg.basisCount();
    std.debug.assert(a.len == m and b.len == m);
    const out = try alloc.alloc(f64, m);
    @memset(out, 0.0);
    for (0..m) |i| {
        if (a[i] == 0.0) continue;
        for (0..m) |j| {
            if (b[j] == 0.0) continue;
            const bp = bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue;
            out[bp.mask] += a[i] * b[j] * @as(f64, @floatFromInt(bp.sign));
        }
    }
    return out;
}

/// Variants that allocate no result — a baseline must not be charged for an
/// allocation nobody would make in a hot loop.
pub fn mulDenseInto(alg: Algebra, a: []const f64, b: []const f64, out: []f64) void {
    const m = alg.basisCount();
    std.debug.assert(a.len == m and b.len == m and out.len == m);
    @memset(out, 0.0);
    for (0..m) |i| {
        if (a[i] == 0.0) continue;
        for (0..m) |j| {
            const bp = bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue;
            out[bp.mask] += a[i] * b[j] * @as(f64, @floatFromInt(bp.sign));
        }
    }
}

pub fn mulDenseSkippingZerosInto(
    alg: Algebra,
    a: []const f64,
    b: []const f64,
    out: []f64,
) void {
    const m = alg.basisCount();
    std.debug.assert(a.len == m and b.len == m and out.len == m);
    @memset(out, 0.0);
    for (0..m) |i| {
        if (a[i] == 0.0) continue;
        for (0..m) |j| {
            if (b[j] == 0.0) continue;
            const bp = bladeMul(alg, @intCast(i), @intCast(j));
            if (bp.sign == 0) continue;
            out[bp.mask] += a[i] * b[j] * @as(f64, @floatFromInt(bp.sign));
        }
    }
}

// ---------------------------------------------------------------------------
// Automatic strategy selection — the actual MRS-computational contribution
// ---------------------------------------------------------------------------

pub const Strategy = enum { sparse, dense };

/// Selects a strategy from the known number of nonzero coefficients and a
/// profitability threshold. Merging in the sparse path is logarithmic, so the
/// threshold is sharp: sparsity wins as long as k_a·k_b·C < 4^n.
/// k_a·k_b·C < 4^n.
pub fn chooseStrategy(
    ka: usize,
    kb: usize,
    basis_count: usize,
    merge_factor: f64,
) Strategy {
    const sparse_ops = @as(f64, @floatFromInt(ka)) *
        @as(f64, @floatFromInt(kb)) * merge_factor;
    const dense_ops = @as(f64, @floatFromInt(basis_count)) *
        @as(f64, @floatFromInt(basis_count));
    return if (sparse_ops < dense_ops) .sparse else .dense;
}

// ---------------------------------------------------------------------------
// Konwersje
// ---------------------------------------------------------------------------

pub fn sparseToDense(
    alloc: std.mem.Allocator,
    alg: Algebra,
    a: Sparse,
) ![]f64 {
    const m = alg.basisCount();
    const out = try alloc.alloc(f64, m);
    @memset(out, 0.0);
    for (a.terms) |t| out[t.mask] += t.coeff;
    return out;
}

/// A dense multivector as a list of nonzero coefficients.
pub fn denseToSparse(alloc: std.mem.Allocator, a: []const f64) !Sparse {
    var n: usize = 0;
    for (a) |x| {
        if (x != 0.0) n += 1;
    }
    const out = try alloc.alloc(Term, n);
    var idx: usize = 0;
    for (a, 0..) |x, i| {
        if (x != 0.0) {
            out[idx] = .{ .mask = @intCast(i), .coeff = x };
            idx += 1;
        }
    }
    return .{ .terms = out };
}

/// Dense table of signs and masks for the whole blade multiplication table.
/// Used to check that the bitmask representation realises the same algebra as
/// the matrix representation (translation proof).
pub fn buildBladeTable(
    alloc: std.mem.Allocator,
    alg: Algebra,
) !BladeTable {
    const m = alg.basisCount();
    const signs = try alloc.alloc(i8, m * m);
    const masks = try alloc.alloc(u32, m * m);
    for (0..m) |i| {
        for (0..m) |j| {
            const bp = bladeMul(alg, @intCast(i), @intCast(j));
            signs[i * m + j] = bp.sign;
            masks[i * m + j] = bp.mask;
        }
    }
    return .{ .signs = signs, .masks = masks };
}

pub const BladeTable = struct {
    signs: []i8,
    masks: []u32,

    pub fn deinit(self: BladeTable, alloc: std.mem.Allocator) void {
        alloc.free(self.signs);
        alloc.free(self.masks);
    }
};

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "defining relations of the Clifford algebra" {
    const alg = try Algebra.fromSignature(sig.minkowski_3_1);
    const n = alg.n_gen;

    // e_i e_j = −e_j e_i for i ≠ j
    for (0..n) |i| {
        for (0..n) |j| {
            if (i == j) continue;
            const p = bladeMul(alg, @as(u32, 1) << @intCast(i), @as(u32, 1) << @intCast(j));
            const q = bladeMul(alg, @as(u32, 1) << @intCast(j), @as(u32, 1) << @intCast(i));
            try std.testing.expectEqual(p.mask, q.mask);
            try std.testing.expectEqual(-p.sign, q.sign);
        }
    }

    // e_i² = s_i
    for (0..n) |i| {
        const p = bladeMul(alg, @as(u32, 1) << @intCast(i), @as(u32, 1) << @intCast(i));
        try std.testing.expectEqual(@as(u32, 0), p.mask);
        try std.testing.expectEqual(alg.squares[i], p.sign);
    }
}

test "reconstructing Cl(1,3) from Minkowski 3+1" {
    const alg = try Algebra.fromSignature(sig.minkowski_3_1);
    try std.testing.expectEqual(@as(i8, 1), alg.squares[0]); // temporal: e_0² = +1
    try std.testing.expectEqual(@as(i8, -1), alg.squares[1]); // space
    try std.testing.expectEqual(@as(i8, -1), alg.squares[2]);
    try std.testing.expectEqual(@as(i8, -1), alg.squares[3]);
    try std.testing.expectEqual(@as(usize, 16), alg.basisCount());

    // signature (-,+,+,+) gives a different e_0² — and that is the whole difference
    const alg_plus = try Algebra.fromSignature(sig.minkowski_3_1_flipped);
    try std.testing.expectEqual(@as(i8, -1), alg_plus.squares[0]);
}

test "associativity of blade multiplication (random triples)" {
    const alg = try Algebra.fromSignature(sig.two_times_2_1);
    var prng = std.Random.DefaultPrng.init(99);
    const rnd = prng.random();
    const m = alg.basisCount();
    for (0..5000) |_| {
        const a: u32 = rnd.intRangeLessThan(u32, 0, @intCast(m));
        const b: u32 = rnd.intRangeLessThan(u32, 0, @intCast(m));
        const c: u32 = rnd.intRangeLessThan(u32, 0, @intCast(m));
        const ab = bladeMul(alg, a, b);
        const bc = bladeMul(alg, b, c);
        const left = bladeMul(alg, ab.mask, c);
        const right = bladeMul(alg, a, bc.mask);
        try std.testing.expectEqual(left.mask, right.mask);
        // (ab)c = a(bc) as algebra elements: sign·e_mask
        try std.testing.expectEqual(ab.sign * left.sign, bc.sign * right.sign);
    }
}

test "nilpotent generator: a degenerate signature gives e_i² = 0" {
    const alg = try Algebra.fromSignature(sig.degenerate_2_1_1);
    try std.testing.expectEqual(@as(i8, 0), alg.squares[3]);
    const p = bladeMul(alg, 1 << 3, 1 << 3);
    try std.testing.expectEqual(@as(i8, 0), p.sign);
}

test "sparse and dense products give the same result" {
    const alloc = std.testing.allocator;
    const alg = try Algebra.fromSignature(sig.minkowski_3_1);

    // rotor: scalar + bivector e1*e2 (a typical element of the Spin group)
    const a = Sparse{ .terms = &.{
        .{ .mask = 0b0110, .coeff = 0.8 },
        .{ .mask = 0b0000, .coeff = 0.6 },
    } };
    const ad = try sparseToDense(alloc, alg, a);
    defer alloc.free(ad);

    const b = Sparse{
        .terms = &.{
            .{ .mask = 0b0011, .coeff = 1.5 }, // e0*e1
            .{ .mask = 0b0001, .coeff = -0.25 }, // e0
        },
    };
    const bd = try sparseToDense(alloc, alg, b);
    defer alloc.free(bd);

    const sp = try mulSparse(alloc, alg, a, b);
    defer sp.deinit(alloc);
    const dp = try mulDense(alloc, alg, ad, bd);
    defer alloc.free(dp);

    for (0..alg.basisCount()) |i| {
        try std.testing.expectApproxEqAbs(dp[i], sp.get(@intCast(i)), 1e-13);
    }
}

test "reverse: (e0*e1)~ = −e0*e1, a vector and a scalar unchanged" {
    const alg = try Algebra.fromSignature(sig.minkowski_3_1);
    try std.testing.expectEqual(@as(i8, 1), alg.reverseSign(0b0000)); // scalar
    try std.testing.expectEqual(@as(i8, 1), alg.reverseSign(0b0001)); // vector
    try std.testing.expectEqual(@as(i8, -1), alg.reverseSign(0b0011)); // biwektor
    try std.testing.expectEqual(@as(i8, -1), alg.reverseSign(0b0111)); // trivector
    try std.testing.expectEqual(@as(i8, 1), alg.reverseSign(0b1111)); // 4-vector
}

test "the allocation-free variant agrees with the allocating one" {
    const alloc = std.testing.allocator;
    const alg = try Algebra.fromSignature(sig.minkowski_3_1);

    const a = Sparse{ .terms = &.{
        .{ .mask = 0b0000, .coeff = 0.8 },
        .{ .mask = 0b0110, .coeff = 0.6 },
        .{ .mask = 0b1001, .coeff = -0.3 },
    } };
    const b = Sparse{ .terms = &.{
        .{ .mask = 0b0011, .coeff = 1.5 },
        .{ .mask = 0b0001, .coeff = -0.25 },
        .{ .mask = 0b0110, .coeff = 0.9 },
    } };

    const ref = try mulSparse(alloc, alg, a, b);
    defer ref.deinit(alloc);

    var scratch: [16]Term = undefined;
    const got = mulSparseScratch(alg, a, b, &scratch);

    try std.testing.expectEqual(ref.terms.len, got.terms.len);
    for (ref.terms, 0..) |t, i| {
        try std.testing.expectEqual(t.mask, got.terms[i].mask);
        try std.testing.expectApproxEqRel(t.coeff, got.terms[i].coeff, 1e-14);
    }
}

test "strategy choice: sparse wins, dense does not" {
    // n = 12 → 4096 coefficients, 4^n = 16.7M
    const basis = @as(usize, 1) << 12;
    try std.testing.expectEqual(Strategy.sparse, chooseStrategy(4, 4, basis, 1.0));
    try std.testing.expectEqual(Strategy.dense, chooseStrategy(basis, basis, basis, 1.0));
}

test "the blade table agrees with the direct product" {
    const alloc = std.testing.allocator;
    const alg = try Algebra.fromSignature(sig.minkowski_1_1);
    const table = try buildBladeTable(alloc, alg);
    defer table.deinit(alloc);
    const m = alg.basisCount();
    for (0..m) |i| {
        for (0..m) |j| {
            const bp = bladeMul(alg, @intCast(i), @intCast(j));
            try std.testing.expectEqual(bp.sign, table.signs[i * m + j]);
            try std.testing.expectEqual(bp.mask, table.masks[i * m + j]);
        }
    }
}
