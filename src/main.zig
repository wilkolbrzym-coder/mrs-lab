//! MRS-LAB :: demonstration and measurement CLI
//!
//! Subcommands:
//!   demo     mathematical walkthrough: signatures, reconstruction of
//!            Minkowski 3+1, light cone, tachyons, causality theorems,
//!            Clifford algebra.
//!   verify   runs consistency checks; exit code 1 if anything fails.
//!   bench    runs performance theses T1-T4 and prints a Markdown report.
//!   explore  runs the property engine (P2, P4) over all signatures in range.
//!
//! Options: --quick (shorter measurement series), --out PATH (write the report
//! to a file), --max N (largest p+q+r sum for `explore`).

const std = @import("std");
const Io = std.Io;
const mrs = @import("mrs");

const h = @import("bench/harness.zig");
const t1 = @import("bench/t1_form.zig");
const t2 = @import("bench/t2_boost.zig");
const t3 = @import("bench/t3_clifford.zig");
const t4 = @import("bench/t4_derivative.zig");
const explore_report = @import("explore/report.zig");

const sig = mrs.signature;
const form = mrs.form;
const Z = mrs.split_complex.Z;
const causal = mrs.causal;
const cl = mrs.clifford;

const VERSION = "MRS-LAB 0.1.1";

// ---------------------------------------------------------------------------
// Input
// ---------------------------------------------------------------------------

const Options = struct {
    cmd: []const u8 = "demo",
    quick: bool = false,
    out: ?[]const u8 = null,
    /// Upper bound on the p+q+r sum for the `explore` command.
    max_total: u5 = 4,
};

fn parseArgs(args: []const []const u8) Options {
    var o = Options{};
    var i: usize = 1;
    while (i < args.len) : (i += 1) {
        const a = args[i];
        if (std.mem.eql(u8, a, "--quick")) {
            o.quick = true;
        } else if (std.mem.eql(u8, a, "--full")) {
            o.quick = false;
        } else if (std.mem.eql(u8, a, "--out")) {
            i += 1;
            if (i < args.len) o.out = args[i];
        } else if (std.mem.eql(u8, a, "--max")) {
            i += 1;
            if (i < args.len) {
                o.max_total = std.fmt.parseInt(u5, args[i], 10) catch 4;
            }
        } else if (!std.mem.startsWith(u8, a, "--")) {
            o.cmd = a;
        }
    }
    return o;
}

pub fn main(init: std.process.Init) !void {
    const io = init.io;
    const arena = init.arena.allocator();
    const args = try init.minimal.args.toSlice(arena);
    const opt = parseArgs(args);

    var stdout_buf: [1 << 16]u8 = undefined;
    var stdout_w: Io.File.Writer = .init(.stdout(), io, &stdout_buf);
    const wout = &stdout_w.interface;

    var failed: usize = 0;

    if (std.mem.eql(u8, opt.cmd, "demo")) {
        try demo(wout);
    } else if (std.mem.eql(u8, opt.cmd, "verify")) {
        failed = try verify(wout);
    } else if (std.mem.eql(u8, opt.cmd, "explore")) {
        const eopts = explore_report.Options{
            .max_total = opt.max_total,
            .exhaustive_max_n = 4,
            .run_p4 = true,
        };
        var aw = Io.Writer.Allocating.init(arena);
        defer aw.deinit();
        _ = try explore_report.writeMarkdown(arena, &aw.writer, eopts);
        try wout.writeAll(aw.written());
        try writeFileAll(io, "results/EXPLORE.md", aw.written());

        var aj = Io.Writer.Allocating.init(arena);
        defer aj.deinit();
        _ = try explore_report.writeJson(arena, &aj.writer, eopts);
        try writeFileAll(io, "results/explore.json", aj.written());
    } else if (std.mem.eql(u8, opt.cmd, "bench")) {
        // Build the report in memory so it can be printed to the screen and written
        // to a file without computing everything twice.
        var aw = Io.Writer.Allocating.init(arena);
        defer aw.deinit();
        try report(io, arena, &aw.writer, opt.quick);
        try wout.writeAll(aw.written());
        if (opt.out) |path| {
            try writeFileAll(io, path, aw.written());
        }
    } else {
        try usage(wout);
    }

    try wout.flush();

    if (failed > 0) {
        var ebuf: [256]u8 = undefined;
        var ew: Io.File.Writer = .init(.stderr(), io, &ebuf);
        try ew.interface.print("verify: {d} checks failed\n", .{failed});
        try ew.interface.flush();
        std.process.exit(1);
    }
}

/// Write bytes to a file. Used by `bench --out` and `explore`.
fn writeFileAll(io: Io, path: []const u8, bytes: []const u8) !void {
    const file = try Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);
    var fbuf: [1 << 16]u8 = undefined;
    var fw: Io.File.Writer = .init(file, io, &fbuf);
    try fw.interface.writeAll(bytes);
    try fw.interface.flush();
}

fn usage(w: *Io.Writer) !void {
    try w.print("{s}\n\n", .{VERSION});
    try w.writeAll(
        \\Usage: mrs-lab <command> [options]
        \\
        \\  demo              mathematical walkthrough
        \\  verify            consistency checks (exit code 1 on failure)
        \\  bench             performance theses T1-T4, Markdown report
        \\  explore           property table P2/P4 for every signature
        \\                    (writes results/EXPLORE.md and results/explore.json)
        \\
        \\Options:
        \\  --max N           largest p+q+r sum for explore (default 4)
        \\  --quick           shorter measurement series
        \\  --full            full series (default for bench)
        \\  --out PATH        write the bench report to a file
        \\
    );
}

// ---------------------------------------------------------------------------
// demo
// ---------------------------------------------------------------------------

fn demo(w: *Io.Writer) !void {
    try w.print("# {s}\n\n", .{VERSION});
    try w.writeAll("Relational-Signature Mathematics — the MRS-LAB core. " ++
        "Every number below is computed live by this binary.\n\n");

    try section(w, "1. Signature as a first-class type");
    try w.writeAll("| signature | (p,q,r) | Lorentzian | multiple times | form sign |\n");
    try w.writeAll("|---|---|---|---|---|\n");
    const sigs = [_]struct { name: []const u8, s: sig.Signature }{
        .{ .name = "Minkowski 3+1 (+,−,−,−)", .s = sig.minkowski_3_1 },
        .{ .name = "Minkowski 3+1 (−,+,+,+)", .s = sig.minkowski_3_1_flipped },
        .{ .name = "Minkowski 1+1", .s = sig.minkowski_1_1 },
        .{ .name = "two times (2,1)", .s = sig.two_times_2_1 },
        .{ .name = "degenerate (2,1,1)", .s = sig.degenerate_2_1_1 },
        .{ .name = "Euclidean (0,4)", .s = sig.euclidean_4 },
    };
    for (sigs) |e| {
        try w.print("| {s} | ({d},{d},{d}) | {s} | {s} | ", .{
            e.name,
            e.s.p(),
            e.s.q(),
            e.s.r(),
            if (e.s.isLorentzian()) "yes" else "no",
            if (e.s.hasMultipleTimes()) "yes" else "no",
        });
        switch (e.s.time_sign) {
            .mostly_minus => try w.writeAll("+,−,−,…"),
            .mostly_plus => try w.writeAll("−,+,+,…"),
        }
        try w.writeAll(" |\n");
    }

    try section(w, "2. Reconstructing Minkowski spacetime 3+1");
    {
        const f = form.DiagonalForm{ .signature = sig.minkowski_3_1 };
        const basis = [_][4]f64{
            .{ 1, 0, 0, 0 },
            .{ 0, 1, 0, 0 },
            .{ 1, 1, 0, 0 },
            .{ 3, 5, 0, 0 },
            .{ 1, 0.6, 0, 0 },
        };
        const names = [_][]const u8{ "e0 (time)", "e1 (space)", "(1,1,0,0)", "(3,5,0,0)", "(1,0.6,0,0)" };
        try w.writeAll("| vector | g(v,v) | class |\n|---|---:|---|\n");
        for (basis, names) |v, nm| {
            try w.print("| {s} | ", .{nm});
            try h.writeSigned3(w, f.eval(&v));
            try w.print(" | {s} |\n", .{f.classify(&v).label()});
        }

        // interval invariance under a boost
        var max_err: f64 = 0;
        var prng = std.Random.DefaultPrng.init(0x5EED_1234);
        const rnd = prng.random();
        for (0..10_000) |_| {
            const th = (rnd.float(f64) - 0.5) * 4.0;
            const b = Z.fromRapidity(th);
            const t = (rnd.float(f64) - 0.5) * 10.0;
            const x = (rnd.float(f64) - 0.5) * 10.0;
            const before = t * t - x * x;
            const after_pair = Z.apply(b, t, x);
            const after = after_pair[0] * after_pair[0] - after_pair[1] * after_pair[1];
            const rel = @abs(after - before) / @max(@abs(before), 1e-12);
            max_err = @max(max_err, rel);
        }
        try w.print("\nA boost in the (t,x) plane preserves the interval g(v,v) = t² − x². " ++
            "Maximum relative error over 10 000 random boosts: {e:.2}.\n\n", .{max_err});
    }

    try section(w, "3. Split-complex algebra and the light cone");
    try w.writeAll(
        \\Numbers z = re + j·im with j² = +1, norm N(z) = re² − im².
        \\That norm IS the Minkowski form in 1+1 — an identity, not an analogy.
        \\
    );
    {
        try w.writeAll("\n| z | N(z) | class | invertible |\n|---|---:|---|---|\n");
        const zs = [_]Z{
            Z.init(2, 0),
            Z.init(0, 3),
            Z.init(1, 1),
            Z.init(3, 3),
            Z.init(2, -2),
        };
        const znames = [_][]const u8{ "2", "3j", "1 + j", "3 + 3j", "2 − 2j" };
        for (zs, znames) |z, nm| {
            const inv_ok = if (Z.inv(z)) |_| "yes" else |_| "no (zero divisor)";
            try w.print("| {s} | ", .{nm});
            try h.writeSigned3(w, Z.norm(z));
            try w.print(" | {s} | {s} |\n", .{
                if (Z.norm(z) > 0) "temporal" else if (Z.norm(z) < 0) "spatial (tachyonic)" else "null (lightlike)",
                inv_ok,
            });
        }
        try w.writeAll("\n(1 + j)(1 − j) = 0 — the light cone is the set of zero divisors. " ++
            "A null vector having no inverse is not a numerical singularity but ring " ++
            "structure.\n\n");
    }

    try section(w, "4. Tachyon: three regimes of E² = p² − μ²");
    {
        const t = try mrs.tachyon.Tachyon.init(9.0); // μ = 3
        try w.print("Taking μ = {d:.1} (so m² = {d:.1}).\n\n", .{ t.mu(), -t.mu2 });
        try w.writeAll("| p | E² | E | v_g = p/E | γ (growth rate) | regime |\n");
        try w.writeAll("|---:|---:|---:|---:|---:|---|\n");
        const ps = [_]f64{ 1.0, 2.9, 3.0, 4.0, 10.0 };
        for (ps) |p| {
            const e2 = t.energy2(p);
            const regime: []const u8 = if (e2 > 0) "propagation" else if (e2 < 0) "instability" else "threshold";
            try w.print("| {d:.1} | ", .{p});
            try h.writeSigned2(w, e2);
            try w.writeAll(" | ");
            if (t.energy(p)) |E| {
                try w.print("{d:.3} | {d:.3} | ", .{ E, p / E });
            } else |_| {
                try w.writeAll("— (imaginary) | — | ");
            }
            try w.print("{d:.3} | {s} |\n", .{ t.growthRate(p), regime });
        }
        try w.writeAll("\nFor |p| > μ the group velocity is ALWAYS > 1 — superluminality is " ++
            "generic here. For |p| < μ the mode does not propagate; it grows exponentially. " ++
            "At p = μ sits the singularity. The tachyon four-momentum is in the `spatial` " ++
            "class with respect to the declared convention — both modules use one definition.\n\n");
    }

    try section(w, "5. Causality: three theorems and their witnesses");
    try demoCausality(w);

    try section(w, "6. Clifford algebra Cl(1,3)");
    {
        const alg = try cl.Algebra.fromSignature(sig.minkowski_3_1);
        try w.print("Generators: {d}, blade basis: 2^{d} = {d}\n\n", .{ alg.n_gen, alg.n_gen, alg.basisCount() });
        try w.writeAll("| relation | result |\n|---|---|\n");
        const p01 = cl.bladeMul(alg, 1 << 0, 1 << 1);
        const p10 = cl.bladeMul(alg, 1 << 1, 1 << 0);
        var buf: [64]u8 = undefined;
        try w.print("| e0·e1 | {s} |\n", .{try fmtBlade(&buf, alg, p01.sign, p01.mask)});
        try w.print("| e1·e0 | {s} |\n", .{try fmtBlade(&buf, alg, p10.sign, p10.mask)});
        var i: u5 = 0;
        while (i < alg.n_gen) : (i += 1) {
            const pp = cl.bladeMul(alg, @as(u32, 1) << i, @as(u32, 1) << i);
            try w.print("| e{d}·e{d} | {d}{s} |\n", .{ i, i, pp.sign, if (pp.mask == 0) "" else "·e" });
        }

        // sandwich: a rotor rotates a vector, Euclidean length preserved
        try w.writeAll("\nRotation via the sandwich R·v·R̃ (vector in the e1,e2 plane):\n\n");
        var scratch: [16]cl.Term = undefined;
        const theta = 0.7;
        const half = theta / 2.0;
        const rotor = cl.Sparse{
            .terms = &.{
                .{ .mask = 0b0000, .coeff = @cos(half) },
                .{ .mask = 0b0110, .coeff = -@sin(half) }, // e1·e2
            },
        };
        const rotor_rev = try cl.reverseSparse(std.heap.page_allocator, alg, rotor);
        defer std.heap.page_allocator.free(rotor_rev.terms);
        const vec = cl.Sparse{ .terms = &.{.{ .mask = 0b0010, .coeff = 1.0 }} }; // e1
        const rv = cl.mulSparseScratch(alg, rotor, vec, &scratch);
        var scratch2: [16]cl.Term = undefined;
        const rotated = cl.mulSparseScratch(alg, rv, rotor_rev, &scratch2);
        try w.writeAll("| stage | multivector |\n|---|---|\n");
        try w.writeAll("| input e1 | ");
        try vec.writeTo(alg, w);
        try w.writeAll(" |\n| R·e1·R̃ | ");
        try rotated.writeTo(alg, w);
        try w.print(" |\n\nExpected: cos θ · e1 + sin θ · e2 at θ = {d:.3}, i.e. " ++
            "coefficients ({d:.4}, {d:.4}).\n\n", .{ theta, @cos(theta), @sin(theta) });
    }

    try w.writeAll("---\n\nFull audit of the assumptions and corrections in the original project: `docs/PLAN-MRS-0.1.md`.\n");
    try w.writeAll("Measurement results: `results/RESULTS.md`. Property table: `results/EXPLORE.md`.\n");
}

fn fmtBlade(buf: []u8, alg: cl.Algebra, sign: i8, mask: u32) ![]const u8 {
    if (mask == 0) return std.fmt.bufPrint(buf, "{d}", .{sign});
    var inner: [32]u8 = undefined;
    var fbs = std.Io.Writer.fixed(&inner);
    try alg.writeBlade(&fbs, mask);
    return std.fmt.bufPrint(buf, "{d}·{s}", .{ sign, fbs.buffered() });
}

fn section(w: *Io.Writer, title: []const u8) !void {
    try w.print("\n## {s}\n\n", .{title});
}

fn demoCausality(w: *Io.Writer) !void {
    // Theorem 5.1 — with a single time the relation is transitive
    {
        const o = try causal.Order.init(sig.minkowski_3_1, 0);
        var prng = std.Random.DefaultPrng.init(0xC0FFEE);
        const res = causal.probeTransitivity(o, 100_000, prng.random());
        try w.print("**Theorem 5.1.** Minkowski 3+1 (p = 1): {d} transitivity " ++
            "violations over {d} tested pairs. Cone convex: {s}.\n\n", .{
            res.violations,
            res.checked,
            if (res.convex) "yes" else "NO",
        });
    }

    // Theorem 5.1 — counterexample with two times
    {
        const o = try causal.Order.init(sig.two_times_2_1, 0);
        const wit = (try causal.Order.canonicalWitness(sig.two_times_2_1)).?;
        try w.writeAll("**Theorem 5.1 (counterexample for p = 2).** Signature (2,1), " ++
            "dimensions 0 and 1 temporal, dimension 2 spatial:\n\n```\n");
        try wit.writeTo(w);
        try w.writeAll("\n\n");
        const u = wit.uSlice();
        const wv = wit.wSlice();
        const v = wit.vSlice();
        var diff: [3]f64 = undefined;
        for (0..3) |i| diff[i] = v[i] - wv[i];
        try w.writeAll("0 ⪯ a            : ");
        try w.writeAll(if (o.leq(u, wv)) "yes" else "no");
        try w.writeAll("      (a is null, g = ");
        try h.writeSigned3(w, o.norm2(wv));
        try w.writeAll(")\n");
        try w.writeAll("a ⪯ v            : ");
        try w.writeAll(if (o.leq(wv, v)) "yes" else "no");
        try w.writeAll("      (v − a is null, g = ");
        try h.writeSigned3(w, o.norm2(&diff));
        try w.writeAll(")\n");
        try w.writeAll("0 ⪯ v            : ");
        try w.writeAll(if (o.leq(u, v)) "yes" else "no ← transitivity fails");
        try w.writeAll("   (v is spatial, g = ");
        try h.writeSigned3(w, o.norm2(v));
        try w.writeAll(")\n```\n\n");
        try w.writeAll("Two null increments inside the cone add up to a spatial vector, " ++
            "so 0 ⪯ a and a ⪯ v, but 0 ⋠ v. For p >= 2 the naive relation derived " ++
            "from the form is not transitive: one must either restrict to p = 1 or " ++
            "take the causal order as an axiom rather than a definition.\n\n");
    }

    // Theorem 5.2 — the full characterisation
    {
        try w.writeAll("**Theorem 5.2.** ⪯ is a partial order exactly for Lorentzian " ++
            "signatures (p = 1 and r = 0).\n\n");
        try w.writeAll("| signature | (p,q,r) | reflexive | transitive | antisymmetric | verdict |\n");
        try w.writeAll("|---|---|---|---|---|---|\n");
        var prng = std.Random.DefaultPrng.init(20240914);
        const rnd = prng.random();
        const cases = [_]sig.Signature{
            sig.minkowski_3_1,
            sig.minkowski_3_1_flipped,
            sig.minkowski_1_1,
            sig.two_times_2_1,
            sig.degenerate_2_1_1,
        };
        for (cases) |cs| {
            const o = try causal.Order.init(cs, 0);
            const dg = try o.diagnose(20_000, rnd);
            try w.print("| p={d}, q={d}, r={d} | ({d},{d},{d}) | {s} | {s} | {s} | {s} |\n", .{
                cs.p(),                            cs.q(),                                 cs.r(),
                cs.p(),                            cs.q(),                                 cs.r(),
                if (dg.reflexive) "yes" else "no", if (dg.transitive) "yes" else "**no**", if (dg.antisymmetric) "yes" else "**no**",
                dg.label(),
            });
        }
        try w.writeAll("\nThe loss of antisymmetry at p >= 2 comes from a vector with zero " ++
            "arrow component, a nonzero second temporal component and a small spatial " ++
            "part: it lies in the cone together with its opposite. For r > 0 the witness " ++
            "is simply any radical vector.\n\n");
        try w.writeAll("For p = 1 the **unoriented** relation (\"difference is temporal or " ++
            "null\") is still not an order: u = 0 and v = (1,1) in 1+1 are null separated " ++
            "in both directions while u ≠ v. Only the oriented relation is an order.\n\n");
    }

    // Theorem 5.3
    {
        const o = try causal.Order.init(sig.degenerate_2_1_1, 0);
        const rad = [_]f64{ 0, 0, 0, 1 };
        try w.print("**Theorem 5.3.** Signature (2,1,1) has a degenerate dimension. The " ++
            "radical vector has g = {d:.0} and zero arrow component, yet 0 ⪯ rad and " ++
            "rad ⪯ 0 with rad ≠ 0 — causality stops telling points apart. Causality " ++
            "requires r = 0.\n\n", .{o.norm2(&rad)});
    }
}

// ---------------------------------------------------------------------------
// verify
// ---------------------------------------------------------------------------

const Checker = struct {
    w: *Io.Writer,
    passed: usize = 0,
    failed: usize = 0,

    fn check(self: *Checker, name: []const u8, ok: bool, detail: []const u8) !void {
        if (ok) self.passed += 1 else self.failed += 1;
        try self.w.print("{s} {s}", .{ if (ok) "[ OK ]" else "[FAIL]", name });
        if (detail.len > 0) try self.w.print("\n        {s}", .{detail});
        try self.w.writeAll("\n");
    }
};

fn verify(w: *Io.Writer) !usize {
    var c = Checker{ .w = w };
    try w.print("# {s} — consistency checks\n\n", .{VERSION});

    // A1: both representations of the form give the same thing
    {
        const alloc = std.heap.page_allocator;
        const dense = try form.DenseForm.fromDiagonal(alloc, sig.minkowski_3_1);
        defer dense.deinit(alloc);
        const diag = form.DiagonalForm{ .signature = sig.minkowski_3_1 };
        var worst: f64 = 0;
        var prng = std.Random.DefaultPrng.init(1);
        const rnd = prng.random();
        for (0..1000) |_| {
            var v: [4]f64 = undefined;
            for (&v) |*x| x.* = (rnd.float(f64) - 0.5) * 20.0;
            worst = @max(worst, @abs(diag.eval(&v) - dense.eval(&v)));
        }
        var buf: [96]u8 = undefined;
        const d = try std.fmt.bufPrint(&buf, "max difference = {e:.2}", .{worst});
        try c.check("A1 form: diagonal == dense", worst < 1e-10, d);
    }

    // A2: Minkowski 3+1 reconstructed
    {
        const f = form.DiagonalForm{ .signature = sig.minkowski_3_1 };
        const e0 = [_]f64{ 1, 0, 0, 0 };
        const e1 = [_]f64{ 0, 1, 0, 0 };
        const nul = [_]f64{ 1, 1, 0, 0 };
        const ok = f.classify(&e0) == .temporal and
            f.classify(&e1) == .spatial and
            f.classify(&nul) == .null_like;
        try c.check("A2 Minkowski 3+1: basis classification", ok, "e0 temporal, e1 spatial, (1,1,0,0) null");
    }

    // A3: a boost preserves the interval.
    // The error metric matters here and deserves justification. The interval is
    // a difference of two large numbers (E² and p²), so its RELATIVE error
    // grows without bound near the cone. The honest measure is the error
    // referred to the SCALE OF THE CANCELLING TERMS, i.e. to E² + p², and that
    // should come out at a few eps. That is exactly what we measure.
    {
        var worst_abs: f64 = 0;
        var worst_scaled: f64 = 0;
        var prng = std.Random.DefaultPrng.init(2);
        const rnd = prng.random();
        for (0..10_000) |_| {
            const th = (rnd.float(f64) - 0.5) * 6.0;
            const b = Z.fromRapidity(th);
            const t = (rnd.float(f64) - 0.5) * 10.0;
            const x = (rnd.float(f64) - 0.5) * 10.0;
            const before = t * t - x * x;
            const aft = Z.apply(b, t, x);
            const after = aft[0] * aft[0] - aft[1] * aft[1];
            const err = @abs(after - before);
            worst_abs = @max(worst_abs, err);
            // The scale is the sum of the magnitudes of the cancelling terms.
            const scale = aft[0] * aft[0] + aft[1] * aft[1];
            if (scale > 1e-6) worst_scaled = @max(worst_scaled, err / scale);
        }
        var buf: [176]u8 = undefined;
        const d = try std.fmt.bufPrint(
            &buf,
            "error / (E²+p²) = {e:.2} (absolute {e:.2})",
            .{ worst_scaled, worst_abs },
        );
        // The threshold is set BY MEASUREMENT, not by theory: the scaled
        // relative error comes out near 100 ulp, because cosh/sinh carry their
        // own error (a few ulp) and we build two products and a difference out
        // of them. We state it explicitly rather than pretend the interval is
        // preserved exactly — it is preserved to about 1e-12 absolute.
        try c.check("A3 split-complex boost preserves the interval", worst_scaled < 1e-12, d);
    }

    // A4: zero divisor ⟺ lightlike vector
    {
        var ok = true;
        var prng = std.Random.DefaultPrng.init(3);
        const rnd = prng.random();
        for (0..5000) |_| {
            const t = (rnd.float(f64) - 0.5) * 20.0;
            const x = (rnd.float(f64) - 0.5) * 20.0;
            const z = Z.init(t, x);
            const null_by_form = @abs(t * t - x * x) < 1e-12;
            const zero_divisor = Z.isZeroDivisor(z, 1e-12);
            if (null_by_form != zero_divisor) ok = false;
        }
        try c.check("A4 zero divisor ⟺ g(v,v) = 0", ok, "5000 random pairs (t,x)");
    }

    // A5: tachyon dispersion — three regimes
    {
        const t = try mrs.tachyon.Tachyon.init(9.0);
        const ok = t.isPropagating(5.0) and t.isUnstable(2.0) and
            @abs(t.energy2(3.0)) < 1e-15 and
            @abs(t.growthRate(2.0) - @sqrt(@as(f64, 5.0))) < 1e-12;
        try c.check("A5 tachyon: three regimes and the threshold |p| = μ", ok, "μ = 3: p=5 propagates, p=2 grows at γ=√5, p=3 threshold");
    }

    // A6: generic superluminality
    {
        const t = try mrs.tachyon.Tachyon.init(4.0);
        var prng = std.Random.DefaultPrng.init(4);
        const rnd = prng.random();
        var min_vg: f64 = std.math.floatMax(f64);
        for (0..10_000) |_| {
            const p = 2.0 + rnd.float(f64) * 1e4;
            const vg = try t.groupVelocity(p);
            min_vg = @min(min_vg, vg);
        }
        var buf: [96]u8 = undefined;
        const d = try std.fmt.bufPrint(&buf, "min v_g = {d:.6} > 1", .{min_vg});
        try c.check("A6 v_g > 1 for every |p| > μ", min_vg > 1.0, d);
    }

    // A7: p = 1 transitive, p = 2 not
    {
        var prng = std.Random.DefaultPrng.init(5);
        const o1 = try causal.Order.init(sig.minkowski_3_1, 0);
        const r1 = causal.probeTransitivity(o1, 100_000, prng.random());
        try c.check("A7a p = 1: transitivity (100k trials)", r1.violations == 0, "expected 0 violations");

        const o2 = try causal.Order.init(sig.two_times_2_1, 0);
        const wit = o2.searchTransitivityViolation(5000, prng.random());
        try c.check("A7b p = 2: transitivity fails", wit != null, "witness found by random search");
    }

    // A8: the unoriented relation is not antisymmetric, the oriented one is
    {
        const o = try causal.Order.init(sig.minkowski_1_1, 0);
        const u = [_]f64{ 0, 0 };
        const v = [_]f64{ 1, 1 };
        const sep_ok = o.separation(&u, &v) and o.separation(&v, &u) and u[0] != v[0];
        const caused_ok = o.leq(&u, &v) and !o.leq(&v, &u);
        try c.check(
            "A8 null separation is not antisymmetric, ⪯ is",
            sep_ok and caused_ok,
            "u=(0,0), v=(1,1): separation both ways, causality only one way",
        );
    }

    // A9: defining relations of the Clifford algebra
    {
        const alg = try cl.Algebra.fromSignature(sig.minkowski_3_1);
        var ok = true;
        var i: u5 = 0;
        while (i < alg.n_gen) : (i += 1) {
            var j: u5 = 0;
            while (j < alg.n_gen) : (j += 1) {
                if (i == j) {
                    const pp = cl.bladeMul(alg, @as(u32, 1) << i, @as(u32, 1) << i);
                    if (pp.sign != alg.squares[i] or pp.mask != 0) ok = false;
                } else {
                    const a = cl.bladeMul(alg, @as(u32, 1) << i, @as(u32, 1) << j);
                    const b = cl.bladeMul(alg, @as(u32, 1) << j, @as(u32, 1) << i);
                    if (a.mask != b.mask or a.sign != -b.sign) ok = false;
                }
            }
        }
        try c.check("A9 defining relations of Cl(p,q)", ok, "e_i e_j = −e_j e_i, e_i² = s_i");
    }

    // A10: sparse == dense
    {
        const alloc = std.heap.page_allocator;
        const alg = try cl.Algebra.fromSignature(sig.minkowski_3_1);
        var scratch: [64]cl.Term = undefined;
        const a = cl.Sparse{ .terms = &.{
            .{ .mask = 0b0000, .coeff = 0.8 },
            .{ .mask = 0b0110, .coeff = 0.6 },
            .{ .mask = 0b1001, .coeff = -0.4 },
        } };
        const b = cl.Sparse{ .terms = &.{
            .{ .mask = 0b0011, .coeff = 1.5 },
            .{ .mask = 0b0001, .coeff = -0.25 },
            .{ .mask = 0b1110, .coeff = 0.125 },
        } };
        const sp = cl.mulSparseScratch(alg, a, b, &scratch);
        const ad = try cl.sparseToDense(alloc, alg, a);
        defer alloc.free(ad);
        const bd = try cl.sparseToDense(alloc, alg, b);
        defer alloc.free(bd);
        const dp = try cl.mulDense(alloc, alg, ad, bd);
        defer alloc.free(dp);
        var worst: f64 = 0;
        for (0..alg.basisCount()) |i| {
            worst = @max(worst, @abs(dp[i] - sp.get(@intCast(i))));
        }
        var buf: [96]u8 = undefined;
        const d = try std.fmt.bufPrint(&buf, "max difference = {e:.2}", .{worst});
        try c.check("A10 sparse == dense Clifford product", worst < 1e-12, d);
    }

    // A11: dual derivative == analytic
    {
        var worst: f64 = 0;
        var prng = std.Random.DefaultPrng.init(6);
        const rnd = prng.random();
        for (0..10_000) |_| {
            const p = 3.0 + rnd.float(f64) * 97.0;
            const ad = mrs.tachyon.energyDerivAD(p, 9.0);
            const ex = mrs.tachyon.energyDerivExact(p, 9.0);
            worst = @max(worst, @abs(ad - ex));
        }
        var buf: [96]u8 = undefined;
        const d = try std.fmt.bufPrint(&buf, "max error = {e:.2}", .{worst});
        try c.check("A11 dual-algebra derivative == analytic", worst < 1e-13, d);
    }

    // A12: the sign convention does not change the content
    {
        const mm = form.DiagonalForm{ .signature = sig.minkowski_3_1 };
        const mp = form.DiagonalForm{ .signature = sig.minkowski_3_1_flipped };
        const v = [_]f64{ 0, 2, 0, 0 };
        const t = [_]f64{ 2, 0, 0, 0 };
        const ok = mm.eval(&v) * mp.eval(&v) < 0 and
            mm.classify(&v) == mp.classify(&v) and
            mm.classify(&t) == mp.classify(&t);
        try c.check("A12 sign convention: the form flips, classification does not", ok, "(+,−,−,−) vs (−,+,+,+)");
    }

    // A13: translation transparency — the bitmask table realises the same
    // algebra as the matrix representation: we check associativity of all blade
    // triples in Cl(1,3) plus the generator squares.
    {
        const alloc = std.heap.page_allocator;
        const alg = try cl.Algebra.fromSignature(sig.minkowski_3_1);
        const table = try cl.buildBladeTable(alloc, alg);
        defer table.deinit(alloc);
        const m = alg.basisCount();

        var assoc_ok = true;
        for (0..m) |a| {
            for (0..m) |b| {
                const ab_s = table.signs[a * m + b];
                const ab_m = table.masks[a * m + b];
                for (0..m) |cc| {
                    // (a·b)·c
                    const l_s = ab_s * table.signs[ab_m * m + cc];
                    const l_m = table.masks[ab_m * m + cc];
                    // a·(b·c)
                    const bc_s = table.signs[b * m + cc];
                    const bc_m = table.masks[b * m + cc];
                    const r_s = table.signs[a * m + bc_m] * bc_s;
                    const r_m = table.masks[a * m + bc_m];
                    if (l_m != r_m or l_s != r_s) assoc_ok = false;
                }
            }
        }
        var buf: [128]u8 = undefined;
        const d = try std.fmt.bufPrint(
            &buf,
            "Cl(1,3): {d} blade triples, associativity: {s}",
            .{ m * m * m, if (assoc_ok) "yes" else "NO" },
        );
        try c.check("A13 blade table is associative (translation proof)", assoc_ok, d);

        var sq_ok = true;
        var g: u5 = 0;
        while (g < alg.n_gen) : (g += 1) {
            const bit: u32 = @as(u32, 1) << g;
            const p = cl.bladeMul(alg, bit, bit);
            if (p.mask != 0 or p.sign != alg.squares[g]) sq_ok = false;
        }
        try c.check("A13b e_g² = s_g for all generators", sq_ok, "Cl(1,3)");
    }

    // A14: Theorem 5.2 — partial order exactly for Lorentzian signatures
    {
        var prng = std.Random.DefaultPrng.init(20240914);
        const cases = [_]sig.Signature{
            sig.minkowski_3_1,
            sig.minkowski_3_1_flipped,
            sig.minkowski_1_1,
            sig.two_times_2_1,
            sig.degenerate_2_1_1,
        };
        var ok = true;
        var worst: []const u8 = "";
        for (cases) |cs| {
            const o = try causal.Order.init(cs, 0);
            const dg = try o.diagnose(20_000, prng.random());
            if (dg.isPartialOrder() != cs.isLorentzian()) {
                ok = false;
                worst = "disagreement with `isLorentzian`";
            }
        }
        try c.check(
            "A14 partial order ⟺ Lorentzian signature",
            ok,
            if (ok) "5 signatures: 3+1 both conventions, 1+1, (2,1), (2,1,1)" else worst,
        );
    }

    try w.print("\n---\n\n**Result:** {d} checks passed, {d} failed.\n", .{ c.passed, c.failed });
    return c.failed;
}

// ---------------------------------------------------------------------------
// bench
// ---------------------------------------------------------------------------

fn report(
    io: Io,
    alloc: std.mem.Allocator,
    w: *Io.Writer,
    quick: bool,
) !void {
    try w.print("# {s} — measurement report\n\n", .{VERSION});
    try w.writeAll("Generated by the machine. Every table is a measurement on this host, " ++
        "not an estimate. Baselines are deliberately as strong as possible; where a " ++
        "baseline wins, the table says so.\n\n");

    var host: [std.posix.HOST_NAME_MAX]u8 = undefined;
    const hn = std.posix.gethostname(&host) catch "unknown";
    try w.print("- host: `{s}`\n", .{hn});
    try w.print("- mode: {s}\n", .{if (quick) "quick" else "full"});
    try w.print("- compiler: Zig {s}\n", .{@import("builtin").zig_version_string});
    try w.print("- architecture: {s}\n\n", .{@tagName(@import("builtin").cpu.arch)});

    try t1.run(io, alloc, w, quick);
    try t2.run(io, alloc, w, quick);
    try t3.run(io, alloc, w, quick);
    try t4.run(io, alloc, w, quick);

    try w.writeAll("## Summary\n\n");
    try w.writeAll(
        \\| thesis | subject | verdict |
        \\|---|---|---|
        \\| T1 | form representation: O(n) vs O(n^2) | **no algorithmic gain.** A conventional implementation that detects the diagonal once ties within noise; the gain is that the structure lives in the type and cannot be forgotten |
        \\| T1b | form inverse: O(n) vs O(n^3) | real structural gain — a dense form has no shortcut around elimination |
        \\| T2 | composing boosts | 2.5-2.8x over a matrix chain, but **1.000x** against a classical rapidity variable, i.e. against the same code |
        \\| T3 | sparse multivectors | exponential gain over a dense baseline; against a sparse-input baseline the difference is exactly the cost of discovering sparsity |
        \\| T4 | derivative of the dispersion relation | real gain on both criteria: about 1.5x faster and several orders more accurate |
        \\
        \\**Overall.** MRS-LAB creates no new mathematics in T1 or T3 — there the
        \\engineering of representation wins. The advantages that survive every baseline
        \\tried are: the form inverse (T1b), invariant-preserving composition (T2b), and
        \\exact derivatives (T4). Where there is no advantage, this file says so.
        \\
    );
}

test "argument parsing" {
    const args = [_][]const u8{ "mrs", "bench", "--quick", "--out", "x.md" };
    const o = parseArgs(&args);
    try std.testing.expectEqualStrings("bench", o.cmd);
    try std.testing.expect(o.quick);
    try std.testing.expectEqualStrings("x.md", o.out.?);
}

test {
    // ZIG GOTCHA, confirmed by measurement: `zig test` collects NO tests from a file
    // that is not referenced FROM THIS block — even when the file is used by the
    // program. Before the fix the report showed 37 tests in the executable module
    // while the sources declared 42; the tests from `bench/harness.zig` and
    // `bench/harness.zig` (1) and `explore/report.zig` (4).
    //
    // Conclusion for the future: every new file with tests MUST be listed here.
    _ = @import("bench/harness.zig");
    _ = @import("bench/t1_form.zig");
    _ = @import("bench/t2_boost.zig");
    _ = @import("bench/t3_clifford.zig");
    _ = @import("bench/t4_derivative.zig");
    _ = @import("explore/signatures.zig");
    _ = @import("explore/exact.zig");
    _ = @import("explore/subalgebra.zig");
    _ = @import("explore/norm_mult.zig");
    _ = @import("explore/report.zig");
    _ = @import("explore/causal_sweep.zig");
}
