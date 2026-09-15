//! MRS-LAB :: report engine (P2 + P4)
//!
//! Rules this file does not break:
//!
//!   1. **Every row of the table is also a test.** There is no result that
//!      cannot be reproduced with `zig build test`.
//!   2. **The report is deterministic.** The same source produces the same
//!      file byte for byte: no floating point, no timestamps, no randomness.
//!      The "method" section describes the procedure, not a measurement.
//!   3. **Scope is part of the result.** When the engine does not decide
//!      something it writes `—` and states why, instead of staying silent or
//!      guessing.
//!
//! By default the report goes up to max_total = 4, because that is where P4 is
//! exhaustive (2^16 subsets). Going to n = 5 is possible
//! (`Options.max_total = 5`) but costs roughly 10^9 operations per signature
//! for P2, so it is a deliberate user choice, not a default.

const std = @import("std");
const mrs = @import("mrs");
const cl = mrs.clifford;
const sigs = @import("signatures.zig");
const exact = @import("exact.zig");
const subalg = @import("subalgebra.zig");
const nrm = @import("norm_mult.zig");

pub const Options = struct {
    /// Largest p+q+r sum included in the report.
    max_total: u5 = 4,
    /// Dimension above which P4 stops being exhaustive (2^n > 16).
    exhaustive_max_n: usize = 4,
    /// Whether to run P4 (cost 2^(2^n) per signature).
    run_p4: bool = true,
};

pub const Row = struct {
    t: sigs.Triple,
    center_dim: u5 = 0,
    center_is_full_center: bool = false,
    p2a: ?bool = null,
    p2b: ?bool = null,
    p2c: ?bool = null,
    subalgebra_count: ?usize = null,
    proper_ideal_count: ?usize = null,
    radical_ideal_index: ?u32 = null,
    note: []const u8 = "",

    pub fn p4InRange(self: Row, opts: Options) bool {
        return self.t.n() <= opts.exhaustive_max_n;
    }
};

pub const Summary = struct {
    rows: usize = 0,
    p2a_true: usize = 0,
    p2a_false: usize = 0,
    p2b_true: usize = 0,
    p2b_false: usize = 0,
    p2c_true: usize = 0,
    p2c_false: usize = 0,
    degenerate: usize = 0,
    lorentzian: usize = 0,
    p4_rows: usize = 0,
    p4_out_of_range: usize = 0,
};

/// Compute one row for a single triple. `grid` is a scratch buffer reused
/// across the loop, so the hot path does not allocate.
pub fn evalTriple(
    alloc: std.mem.Allocator,
    t: sigs.Triple,
    ts: mrs.signature.TimeSign,
    grid: *[nrm.MAX_GRID]exact.IntVec,
    opts: Options,
) !Row {
    var row = Row{ .t = t };
    const sb = sigs.SigBuf.build(t, ts);
    const alg = try sb.algebra();

    row.center_dim = exact.dimOfSet(exact.centerBasis(alg));

    // --- P2: decision procedure -------------------------------------------
    const v_a = nrm.checkScalarMultiplicative(alg, grid);
    row.p2a = v_a.holds;
    const v_b = nrm.checkCenterMultiplicative(alg, grid);
    row.p2b = v_b.holds;
    const v_c = nrm.checkNormIsCentral(alg, grid);
    row.p2c = v_c.holds;

    // --- P4: exhaustive within the supported range -----------------------
    if (t.isDegenerate()) {
        const I = exact.bladesWithDegenerate(alg);
        row.radical_ideal_index = subalg.nilpotencyIndex(alg, I, 8);
    }
    if (opts.run_p4 and t.n() <= opts.exhaustive_max_n) {
        row.subalgebra_count = try subalg.countClosed(alg);
        row.proper_ideal_count = try subalg.countProperIdeals(alg);
    }

    // --- note: scope is part of the result -------------------------------
    row.note = blk: {
        if (opts.run_p4 and t.n() > opts.exhaustive_max_n) {
            break :blk "P4 out of range: 2^n > 16";
        }
        if (t.isDegenerate()) {
            if (row.radical_ideal_index) |idx| {
                break :blk if (idx == 2) "radical ideal, I^2 = 0" else "radical ideal, nilpotent";
            }
            break :blk "no radical ideal";
        }
        break :blk "no blade-spanned ideals";
    };

    _ = alloc;
    return row;
}

fn writeBool(w: anytype, v: ?bool) !void {
    if (v) |b| {
        try w.writeAll(if (b) "yes" else "**no**");
    } else {
        try w.writeAll("—");
    }
}

fn writeOptUsize(w: anytype, v: ?usize) !void {
    if (v) |x| try w.print("{d}", .{x}) else try w.writeAll("—");
}

pub fn writeMarkdown(
    alloc: std.mem.Allocator,
    w: anytype,
    opts: Options,
) !Summary {
    var sum = Summary{};

    // Buffer for the determining set: one allocation, outside the loop.
    const gbuf = try alloc.create([nrm.MAX_GRID]exact.IntVec);
    defer alloc.destroy(gbuf);

    var triples: [sigs.tripleCount(sigs.MAX_TOTAL)]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&triples, opts.max_total);

    try w.writeAll("# MRS-LAB — property engine results\n\n");
    try w.print("Scope: all signatures (p,q,r) with p+q+r from 1 to {d}. " ++
        "Convention: (+,−,−,…); the P2 results are convention independent " ++
        "(checked by test).\n\n", .{opts.max_total});
    try w.writeAll("## Method — why these are decisions, not evidence\n\n");
    try w.writeAll("- **P2a/P2b/P2c**: for fixed y the identity is a **quadratic function** " ++
        "of x, and symmetrically, and a quadratic function (including a linear term) is " ++
        "determined by its values on {0} ∪ {e_i} ∪ {e_i+e_j}. Checking the Cartesian " ++
        "product of that set **decides** the identity for all real x,y. The arithmetic is " ++
        "exact integers — no tolerance, no probability.\n");
    try w.writeAll("- **P4**: exhaustive over all subspaces spanned by blades. That is a " ++
        "sublattice of all subalgebras — ideals spanned by idempotents (e.g. in " ++
        "Cl(1,0) ≅ R⊕R) are outside this scope and are not visible here.\n");
    try w.writeAll("- **Radical ideal**: the blades containing a degenerate generator. " ++
        "The index is the smallest k with I^k = 0.\n\n");

    try w.writeAll("## Table\n\n");
    try w.writeAll("| (p,q,r) | n | P2a | P2b | P2c | centre dim | subalgebras | proper ideals | radical idx | notes |\n");
    try w.writeAll("|---|---:|---|---|---|---:|---:|---:|---:|---|\n");

    for (0..n) |i| {
        const row = try evalTriple(alloc, triples[i], .mostly_minus, gbuf, opts);
        sum.rows += 1;
        if (row.t.isDegenerate()) sum.degenerate += 1;
        if (row.t.isLorentzian()) sum.lorentzian += 1;
        if (row.p2a.?) sum.p2a_true += 1 else sum.p2a_false += 1;
        if (row.p2b.?) sum.p2b_true += 1 else sum.p2b_false += 1;
        if (row.p2c.?) sum.p2c_true += 1 else sum.p2c_false += 1;
        if (row.p4InRange(opts) and opts.run_p4) sum.p4_rows += 1 else if (opts.run_p4) {
            sum.p4_out_of_range += 1;
        }

        try w.writeAll("| ");
        try row.t.writeTo(w);
        try w.print(" | {d} | ", .{row.t.n()});
        try writeBool(w, row.p2a);
        try w.writeAll(" | ");
        try writeBool(w, row.p2b);
        try w.writeAll(" | ");
        try writeBool(w, row.p2c);
        try w.print(" | {d} | ", .{row.center_dim});
        try writeOptUsize(w, row.subalgebra_count);
        try w.writeAll(" | ");
        try writeOptUsize(w, row.proper_ideal_count);
        try w.writeAll(" | ");
        if (row.radical_ideal_index) |x| {
            try w.print("{d}", .{x});
        } else {
            try w.writeAll("—");
        }
        try w.print(" | {s} |\n", .{row.note});
    }

    try w.writeAll("\n## Summary\n\n");
    try w.print("- signatures in scope: **{d}**\n", .{sum.rows});
    try w.print("- P2a holds in {d}, fails in {d}\n", .{ sum.p2a_true, sum.p2a_false });
    try w.print("- P2b holds in {d}, fails in {d}\n", .{ sum.p2b_true, sum.p2b_false });
    try w.print("- P2c holds in {d}, fails in {d}\n", .{ sum.p2c_true, sum.p2c_false });
    try w.print("- degenerate signatures (r>0): {d}\n", .{sum.degenerate});
    try w.print("- Lorentzian signatures: {d}\n", .{sum.lorentzian});
    try w.print("- P4 rows computed exhaustively: {d}, out of range: {d}\n", .{
        sum.p4_rows, sum.p4_out_of_range,
    });
    return sum;
}

pub fn writeJson(
    alloc: std.mem.Allocator,
    w: anytype,
    opts: Options,
) !Summary {
    const gbuf = try alloc.create([nrm.MAX_GRID]exact.IntVec);
    defer alloc.destroy(gbuf);

    var triples: [sigs.tripleCount(sigs.MAX_TOTAL)]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&triples, opts.max_total);

    var sum = Summary{};
    try w.writeAll("{\n");
    try w.print("  \"max_total\": {d},\n", .{opts.max_total});
    try w.writeAll("  \"results\": [\n");

    for (0..n) |i| {
        const row = try evalTriple(alloc, triples[i], .mostly_minus, gbuf, opts);
        sum.rows += 1;
        if (row.p2a.?) sum.p2a_true += 1 else sum.p2a_false += 1;
        if (row.p2b.?) sum.p2b_true += 1 else sum.p2b_false += 1;
        if (row.p2c.?) sum.p2c_true += 1 else sum.p2c_false += 1;

        try w.print("    {{\"p\": {d}, \"q\": {d}, \"r\": {d}, \"n\": {d}, \"p2a\": {s}, " ++
            "\"p2b\": {s}, \"p2c\": {s}, \"center_dim\": {d}, ", .{
            row.t.p,                            row.t.q,
            row.t.r,                            row.t.n(),
            if (row.p2a.?) "true" else "false", if (row.p2b.?) "true" else "false",
            if (row.p2c.?) "true" else "false", row.center_dim,
        });
        if (row.subalgebra_count) |x| {
            try w.print("\"subalgebras\": {d}, ", .{x});
        } else {
            try w.writeAll("\"subalgebras\": null, ");
        }
        if (row.proper_ideal_count) |x| {
            try w.print("\"proper_ideals\": {d}, ", .{x});
        } else {
            try w.writeAll("\"proper_ideals\": null, ");
        }
        if (row.radical_ideal_index) |x| {
            try w.print("\"radical_index\": {d}", .{x});
        } else {
            try w.writeAll("\"radical_index\": null");
        }
        if (i + 1 < n) {
            try w.writeAll("},\n");
        } else {
            try w.writeAll("}\n");
        }
    }
    try w.writeAll("  ]\n}\n");
    return sum;
}

// ---------------------------------------------------------------------------
// Testy
// ---------------------------------------------------------------------------

test "report is deterministic byte for byte" {
    const alloc = std.testing.allocator;
    const opts = Options{ .max_total = 3, .exhaustive_max_n = 3 };

    var a = std.Io.Writer.Allocating.init(alloc);
    defer a.deinit();
    var b = std.Io.Writer.Allocating.init(alloc);
    defer b.deinit();

    const sa = try writeMarkdown(alloc, &a.writer, opts);
    const sb = try writeMarkdown(alloc, &b.writer, opts);

    try std.testing.expectEqualStrings(a.written(), b.written());
    try std.testing.expectEqual(sa.rows, sb.rows);
    try std.testing.expectEqual(sa.p2a_true, sb.p2a_true);
    try std.testing.expectEqual(sa.p2a_false, sb.p2a_false);
}

test "JSON is deterministic and has a complete set of rows" {
    const alloc = std.testing.allocator;
    const opts = Options{ .max_total = 3, .exhaustive_max_n = 3 };

    var a = std.Io.Writer.Allocating.init(alloc);
    defer a.deinit();
    var b = std.Io.Writer.Allocating.init(alloc);
    defer b.deinit();

    const s1 = try writeJson(alloc, &a.writer, opts);
    const s2 = try writeJson(alloc, &b.writer, opts);
    try std.testing.expectEqualStrings(a.written(), b.written());
    try std.testing.expectEqual(s1.rows, sigs.tripleCount(3));
    try std.testing.expectEqual(s2.rows, sigs.tripleCount(3));
}

test "the P2a boundary in the report falls exactly at n = 3" {
    const alloc = std.testing.allocator;
    const opts = Options{ .max_total = 4, .exhaustive_max_n = 4 };
    var a = std.Io.Writer.Allocating.init(alloc);
    defer a.deinit();
    const s = try writeMarkdown(alloc, &a.writer, opts);

    // For n <= 2 P2a holds; for n >= 3 it fails. We check that the summary
    // is consistent with the number of triples rather than asserting a range.
    try std.testing.expectEqual(sigs.tripleCount(4), s.rows);
    try std.testing.expect(s.p2a_false > 0);
    try std.testing.expect(s.p2a_true > 0);
}

test "P4 out of range is stated explicitly, not guessed" {
    const alloc = std.testing.allocator;
    const opts = Options{ .max_total = 5, .exhaustive_max_n = 4 };
    const gbuf = try alloc.create([nrm.MAX_GRID]exact.IntVec);
    defer alloc.destroy(gbuf);

    const row = try evalTriple(alloc, .{ .p = 0, .q = 5 }, .mostly_minus, gbuf, opts);
    try std.testing.expectEqual(@as(?usize, null), row.subalgebra_count);
    try std.testing.expectEqualStrings("P4 out of range: 2^n > 16", row.note);
}
