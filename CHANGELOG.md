# Changelog

All notable changes to MRS-LAB are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Two rules from the README apply to every entry below: a scope limit is a scope
limit and not a negative result, and "not shown" never means "impossible".

## [Unreleased] — 0.1.1

Audit of the multi-time signatures (p >= 2), i.e. exactly the region where
causality is supposed to break: two and three time dimensions, with and without
spatial dimensions. Two defects were found there, both in the causal layer, one
of them silently wrong. Everything below is reproducible: `zig build test`
(285 tests, Debug + ReleaseSafe + ReleaseFast), `zig build verify` (16 checks),
`zig build explore` (results byte-identical to the committed ones).

### Fixed

- **A silently wrong verdict: signatures with two or three time dimensions and
  NO spatial dimension were reported as partial orders.** `(2,0,0)` and
  `(3,0,0)` came back with `isPartialOrder() == true`, contradicting Theorem 5.2
  of the module itself (⪯ is a partial order ⟺ p = 1 and r = 0). Cause: the
  antisymmetry witness for p >= 2 required a spatial dimension in order to build
  `d = e_t1 + ½·e_s`, although the purely temporal `d = e_t1` already has arrow
  component 0 and `g(d) = s_t ≠ 0`, so `d` and `−d` both lie in the closed future
  cone and `u ⪯ u+d ⪯ u` with `d ≠ 0`. The guard was a condition that the
  mathematics never needed. After the fix: `(2,0,0)` and `(3,0,0)` report
  `partial_order = false`, `antisymmetric = false`, `transitive = true` — the
  honest verdict "preorder".
  Guarded by `explore/causal_sweep.zig :: SWEEP: ⪯ is a partial order exactly for
  Lorentzian signatures, for every p >= 1`, which sweeps every signature with
  `p >= 1, p+q+r <= 5` instead of a list of examples, and by
  `mrs/causal.zig :: the antisymmetry witness for p >= 2 does not need a spatial
  dimension`. The old test list contained no signature without a spatial
  dimension, which is precisely why the defect survived it.

- **An empty search was reported as a verdict.** For `q = 0` the cone-vector
  generator has no spatial dimension to normalise against, so it produces
  nothing; `probeTransitivity` then returned `convex = (violations == 0)` = true
  with `checked == 0`, i.e. "no violation" from zero trials. The sibling
  function `coneIsConvex` already refused to claim convexity without a single
  successful trial — the two disagreed on the same question. `convex` now
  requires `checked > 0`.
  Guarded by `explore/causal_sweep.zig :: q = 0: the generator produces nothing,
  the search is empty, and no claim is made`.

### Fixed

- A self-contradicting sentence in the T2 report: the prose said "the ratio is 1"
  while the printed number was a noisy `1.036x`. The identity check between MRS-LAB
  and a classical rapidity variable is exact by construction, so the sentence now
  reads "the ratio is 1 up to timer noise" and matches whatever the timer prints.

### Changed

- **Transitivity is now decided where the engine can decide it, and the method
  travels with the verdict.** `Causality` gains `method`
  (`constructive_witness` | `half_space_proof` | `sampling`) and `trials_checked`,
  plus `isProved()`. `Order.diagnose` now decides instead of sampling in two
  cases it used to approximate:
  - `q = 0`: no nonzero vector is spatial, so the future cone is the closed
    half-space `{v : v_arrow >= 0}`, which is closed under addition —
    transitivity HOLDS, by argument, with zero trials;
  - `p >= 2` with a spatial dimension: `canonicalWitness` gives
    `u ⪯ w ⪯ v` with `u ⋠ v` — transitivity FAILS constructively, so the verdict
    no longer depends on a random draw finding a null-vector counterexample.
  A sampled verdict can no longer be vacuous: `method == .sampling` implies
  `trials_checked > 0`, which the sweep asserts for every signature in range.
  This is the README rule "proof or counterexample, never probably" applied to
  the one place where it was not yet true.

- `form.Class.label()` now returns English (`"temporal (timelike)"`,
  `"spatial (tachyonic)"`); it is public, printed by `demo` and `verify`, and
  used to leak Polish into the console of an English-only project.

### Fixed

- **Translation leftovers.** The English-only pass of commit b8af845 missed
  strings and comments that a line-oriented scan does not catch: debug output
  inside a test, a `@panic` message, test names, `Class.label()`, and 70 dodgy
  comment lines, 67 of which carried duplicated comment markers (`///  /// …`)
  left by the rewriting pass. All of them are translated; no Polish text remains
  in `src/` (checked against the system `pl_PL` dictionary, then read by hand).

- **A test that wrote to stderr made the whole build step look broken.** The
  P2b/P2c observation test printed its count with `std.debug.print` on the
  SUCCESS path, and in Zig 0.16 the build runner answers any stderr from a test
  process with `failed command: … --listen=-` lines attributed to the step —
  while the build still succeeds and all tests pass. Evidence is now carried by
  assertions (`expectEqual(10, differ_a)`), so `zig build test` prints a clean
  summary, and a real failure is no longer buried in noise. Diagnosis of the
  noise is reproducible: adding a `std.debug.print` to any passing test makes
  that step emit the same line.

### Added

- `src/explore/causal_sweep.zig` — eight tests that sweep the causal verdict over
  every signature with `p >= 1` instead of a sample of them. Beyond guarding the
  two defects above, it records three structural facts measured while auditing
  multi-time signatures:
  - **the blade-scoped P4 lattice is blind to the time/space split.** For every
    triple with `n <= 4`, the number of closed blade-spanned subalgebras, the
    number of proper two-sided ideals and the dimension of the centre depend only
    on `n` and `r`: all five signatures `(0,4,0) … (4,0,0)` give 68 closed
    subspaces, and `(2,2,0)` is indistinguishable from `(1,3,0)` down to the
    unital-commutative subalgebras. A two-time universe is therefore
    algebraically invisible to P4 — the time/space split shows up in the causal
    layer, not in the subalgebra layer;
  - **a Lorentzian order survives inside every two-time signature.** The
    restriction of the `(2,q)` relation to the subspace spanned by one temporal
    and all spatial dimensions agrees with the Lorentzian `(1,q)` order on every
    sampled pair, and transitivity holds inside that subspace, while the ambient
    relation is not transitive. There are `p` such hyperplanes — one per choice
    of the kept time — and they are the tractable answer to "where does a hard,
    deterministic order live inside a multi-time signature";
  - **strict future steps cannot close a loop, even though addition can leave the
    cone.** Two timelike future steps with a positive arrow component can sum to
    a SPATIAL vector (witness inside the test), so the relation is not
    transitive even in the strict sector; but the arrow coordinate strictly
    increases along a strict step, so no sequence of strict steps returns to its
    start. In a two-time signature the "time loop" that does exist is the
    2-cycle `u ⪯ u+d ⪯ u` with `d ≠ 0` and zero arrow component, not a
    chronology-violating path.
- Removed the lifetime hazard of `Order` at the type level instead of documenting
  it. `Signature` borrows a `[]const Role` slice, so an order built from a
  temporary buffer and returned from a helper used to dangle — it compiled and then
  failed at run time as `switch on corrupt value`. `Order` now copies the roles
  into a fixed array, and builds its form with the new `DiagonalForm.initOwned`,
  which copies the signs, the dimension and the convention and stores an EMPTY role
  slice; the order is therefore self-contained and safe to return from the frame
  that built it. `Order.signature()` takes a pointer rather than a value, because
  returning a slice into a by-value parameter would dangle in turn. Covered by two
  regression tests that build the object inside a local frame and use it after
  returning.

### Known follow-ups

- `results/RESULTS.md` is a committed benchmark snapshot measured on one machine.
  It has been regenerated from the current sources, because the report HEADINGS are
  part of the artifact and had to match the source; an earlier revision kept a stale
  snapshot on purpose, which left Polish headings in the published file and broke
  the README's claim that the reports are reproducible from the repository. Its
  timings remain a snapshot from that machine and are not a claim about any other.
- P4 above `p+q+r = 4` is still out of reach: `(2,3,0)` — the signature asked
  about first — has `n = 5`, so the engine refuses (`error.TooManyBlades`) by
  construction. The bounded substitute used in the audit (all blade subsets of
  dimension <= 4: 187 closed subspaces, 107 of them unital and commutative, no
  proper ideal) is a lower bound within a declared scope, not the full lattice.
- P1 remains unimplemented as a classification: "which subgroups `H ⊆ O(p,q)`
  admit an invariant pointed convex cone" still has no LP layer behind it.
