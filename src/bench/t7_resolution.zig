//! MRS-LAB :: thesis T7 — how close to the null cone can the engine still decide?
//!
//! WHAT CHANGED IN 0.1.6. The causal layer used to ask `|g| <= 1e-9` to decide
//! whether a vector is null. That is a constant this project chose, and it has
//! two consequences that this thesis MEASURES rather than asserts:
//!
//!   1. RESOLUTION. A vector whose norm is genuinely 1e-10 is not near-null in
//!      any useful sense — it is a vector 1e-10 away from the cone, and the
//!      engine refused to say which side. Deciding from the propagated radius
//!      moves that boundary from ~5e-10 (half the tolerance) down to the point
//!      where f64 itself stops separating, here ~1e-16.
//!
//!   2. SCALE. An absolute tolerance is not a statement about geometry. Scale a
//!      vector down by 1e5 and its norm falls by 1e10, so a vector that is TEN
//!      PER CENT off the cone gets reported as null. The propagated radius
//!      scales with the computation, so the decision is scale invariant — which
//!      is what a statement about the cone must be.
//!
//! Both parts below are measurements on the same family, `v = λ·(1, 1−δ, 0, 0)`
//! in (+,−,−,−), whose exact norm is `λ²(2δ − δ²)`.

const std = @import("std");
const Io = std.Io;
const mrs = @import("mrs");
const h = @import("harness.zig");

const sig = mrs.signature;
const form = mrs.form;
const contract = mrs.contract;

/// The old rule, kept in the source so the comparison is against real code.
fn calledNullByTolerance(f: form.DiagonalForm, v: []const f64, tol: f64) bool {
    return f.classifyTol(v, tol) == .null_like;
}

pub fn run(
    io: Io,
    alloc: std.mem.Allocator,
    w: anytype,
    quick: bool,
) !void {
    _ = io;
    _ = alloc;
    const f = form.DiagonalForm.init(sig.minkowski_3_1);
    const tol: f64 = 1e-9; // the constant the causal layer used

    try w.writeAll("### T7 — the resolution of a classification\n\n");
    try w.writeAll("Family: `v = λ·(1, 1−δ, 0, 0)` in (+,−,−,−), exact norm `λ²(2δ − δ²)`. " ++
        "At `λ = 1` the second part holds and the first column is the distance from " ++
        "the cone in units of the vector's own square.\n\n");

    // --- Part 1: how small a norm is still decided -------------------------
    try w.writeAll("| λ | δ | exact `g` | radius | 0.1.5 rule (`\\|g\\| <= 1e-9`) | 0.1.6 rule |\n");
    try w.writeAll("|---:|---:|---:|---:|---|---|\n");
    const deltas = if (quick)
        &[_]f64{ 1e-2, 1e-8, 1e-10, 1e-16, 1e-18 }
    else
        &[_]f64{ 1e-2, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16, 1e-18 };
    for (deltas) |d| {
        const lam: f64 = 1.0;
        var v = [_]f64{ lam, lam * (1.0 - d), 0, 0 };
        const exact = lam * lam * (2.0 * d - d * d);
        const b = contract.evalForm(f, &v);
        const new_says_null = contract.classify(f, &v) == .null_like;
        try w.print("| {d:.0} | {e:.0} | {e:.2} | {e:.2} | {s} | **{s}** |\n", .{
            lam,
            d,
            exact,
            b.radius,
            if (calledNullByTolerance(f, &v, tol)) "null" else "decided",
            if (new_says_null) "null" else "decided",
        });
    }

    // --- Part 2: the absolute tolerance is not scale invariant -------------
    try w.writeAll("\nThe same family, scaled. `δ = 0.05` means the exact norm is " ++
        "9.75 % of the vector's own square — a vector that is nowhere near the cone.\n\n");
    try w.writeAll("| λ | δ | exact `g` | \\|g\\| / λ² | 0.1.5 rule | 0.1.6 rule |\n");
    try w.writeAll("|---:|---:|---:|---:|---|---|\n");
    const lambdas = [_]f64{ 1e6, 1.0, 1e-2, 1e-5 };
    for (lambdas) |lam| {
        const d: f64 = 0.05;
        var v = [_]f64{ lam, lam * (1.0 - d), 0, 0 };
        const exact = lam * lam * (2.0 * d - d * d);
        try w.print("| {e:.0} | {d:.2} | {e:.2} | {d:.3} | {s} | **{s}** |\n", .{
            lam,
            d,
            exact,
            @abs(exact) / (lam * lam),
            if (calledNullByTolerance(f, &v, tol)) "**null**" else "decided",
            if (contract.classify(f, &v) == .null_like) "null" else "decided",
        });
    }

    try w.writeAll("\n**Verdict.** The tolerance version fails twice, and the failures are visible " ++
        "in the tables above. It calls a vector null while the computation can still separate " ++
        "its norm from zero (a refusal that costs nothing), and it calls a vector null because " ++
        "it is SMALL, which has nothing to do with being ON THE CONE — scale `δ = 0.05` down by " ++
        "1e5 and a vector 9.75 % off the cone is reported as null. The radius version decides " ++
        "whenever f64 can, and it is scale invariant because the radius is built from the same " ++
        "products the value is.\n\n");
    try w.writeAll("The honest limit: the resolution is now set by f64 and not by us, so it is " ++
        "**worse** than the tolerance for a genuinely fuzzy input — a measurement with an error " ++
        "of 1e-8 has a decided classification of a norm it does not know to 1e-8. That is what " ++
        "`Bounded.fromError` is for: an input with a stated error carries it into the radius, " ++
        "and the classification goes back to being undecided, this time for a reason that is " ++
        "true.\n\n");
}
