# Changelog

All notable changes to MRS-LAB are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Two rules from the README apply to every entry below: a scope limit is a scope
limit and not a negative result, and "not shown" never means "impossible".

## [0.1.3] — 2026-09-15

The tolerance stops being an argument the caller supplies and becomes part of
the value. Reproducible: `zig build test` (327 tests, Debug + ReleaseSafe +
ReleaseFast), `zig build verify` (18 checks), `zig build bench` (T1–T5 into
`results/RESULTS.md`).

### Added

- **`src/mrs/contract.zig` — the error bound travels with the number.** Until
  0.1.2 a classification needed a tolerance from outside: `classifyTol(v, tol)`
  in `form.zig`, `isZeroDivisor(a, tol)` in `split_complex.zig`, a constant in
  `causal.zig`. The number did not know how well it was known, so "is this
  vector null?" was answered by a constant chosen somewhere else. `Bounded`
  carries `value` and `radius` with `|value − exact| <= radius`, and the radius
  is composed by the same operations that produce the value:
  - `add`/`sub`: `ra + rb + U·|value|`;
  - `mul`: `|a|·rb + |b|·ra + ra·rb + U·|value|`;
  where `U = 2⁻⁵³` is the unit roundoff. Every rule is checked against `f128`
  arithmetic on grids, not on examples.
- **The zero test is now derived rather than supplied**: a number is
  indistinguishable from zero when its own radius covers zero. The other half
  of that rule matters as much — `1e-300` known exactly is **not** called zero.
  `verify` A16 guards both directions.
- **`Contract` — three meanings of "the same result", now distinguishable in
  the type**: `bit_exact`, `correctly_rounded`, `bounded`. The contract can be
  a compile-time parameter, and when it carries no radius the layer IS the old
  code path: `evalFormWith(.bit_exact, …)` returns a plain `f64` from the body
  `f.eval(v)`. Those contracts therefore cost nothing by construction rather
  than by an optimiser, which is the death criterion of this design discharged
  rather than measured. A test asserts the results are bit-identical.
- `contract.classify` — classification of a vector against the declared
  convention, with `null_like` returned when the radius covers zero.
  `classifyTol` stays for callers who want to supply a tolerance; the
  difference is documented at both.
- `src/bench/t5_contract.zig` — the price of the radius, measured.
- `demo` prints the contract section: `g((1,1,0,0)) = 0 ± 5.55e-16` → null by
  the radius, and `1e-300` → not null. `verify` gains A15 (the propagated
  radius contains the exact form value, 2000 evaluations over 4 signatures
  against `f128`) and A16 (the zero test).

### Measured — and this is the interesting part, because it killed the naive version

The bounds were propagated through a 100 000-step chain of boosts with mixed
signs BEFORE the module was written. The exact answer is known analytically (a
product of boosts has invariant 1), so no oracle was needed:

| representation | propagated radius | true error | useful? |
|---|---:|---:|---|
| additive chart (rapidity) | 1.0e-12 | 2.8e-17 | yes |
| split-complex product | 5.3e17 | 2.9e-14 | **no** |
| matrix product (determinant) | 4.8e18 | 3.4e-14 | **no** |

All three are CORRECT worst-case bounds. Only the first is useful. The rules
take absolute values, and an absolute-value rule cannot see the cancellation
that the additive coordinate makes explicit, so the bound grows like
`exp(Σ|θᵢ|)` in the multiplicative representations and like `N·U` in the
additive one. The conclusion is the thesis of this project reached from the
error-propagation side instead of from a stopwatch: **a contract is
dischargeable exactly where the representation is well conditioned.** The test
`the additive chart keeps a useful bound where the multiplicative one does not`
guards both halves of that claim.

Two prototype defects were found and fixed before the numbers above were
believed, and they are recorded because a wrong number is worse than no number:
the first version indexed the wrong radius array in the matrix chain, which
produced a spurious exponential blow-up, and the first version normalised the
random rapidities by dividing by their sum, which is ill-conditioned when that
sum is near zero and produced a second spurious blow-up.

### Changed

- `results/RESULTS.md` gains T5 and its summary table gains a row. T5 is not a
  speed comparison and is not presented as one: it prices a capability that
  earlier versions did not have (1.67×–2.10× where a radius is asked for).

### Known follow-ups

- **The contract is not yet wired into `causal.zig` or `split_complex.zig`.**
  Their tolerance arguments still exist and still behave as before; the new path
  is `contract.classify` and `contract.Bounded`. Converting them would change
  the verdicts those functions return for near-null vectors, which is a
  behaviour change that belongs in its own release with its own counterexamples.
- **The approximation half is deliberately absent.** Bounds that compose through
  a pipeline of approximate summaries (sketches, digests, confidence levels) are
  a different library with a different open question — how the confidence levels
  compose — and do not belong in a project about metric signatures.
- The radius rules are the simple forward ones. Sharper rules exist (running
  error bounds that keep the correlation between terms, as in compensated
  summation); the measured table above is what those would have to beat.

## [0.1.2] — 2026-09-15

P4 reaches dimension 5, by a change of enumeration method rather than by
waiting longer. Everything below is reproducible: `zig build test` (300 tests,
Debug + ReleaseSafe + ReleaseFast), `zig build verify` (16/16),
`zig build explore` (1 m 41 s on the machine that produced the committed table,
and machine dependent like every other timing here; `explore --max 4` gives the
previous 34-row table in under a second).

### Added

- **The closure walk: P4 now decides every signature the engine supports.**
  Until 0.1.1 the subalgebra question stopped at `p+q+r = 4`, because
  `countsBrute` scans all `2^(2^n)` blade subsets and `n = 5` means 4 294 967 296
  of them. `countsByClosure` walks the FIXED POINTS of the closure operator
  `S ↦ smallest closed superset of S` (Ganter's next closure) and visits each
  closed subspace exactly once, so its cost is the number of closed subspaces:
  at `n = 5` that is 375 for a non-degenerate algebra and 31 242 668 for
  `(0,0,5)`. `(2,3,0)` — the signature the 0.1.1 audit asked about first and had
  to leave unanswered — is now a decided row: 375 closed blade-spanned
  subalgebras, none of them a proper two-sided ideal.
- `subalgebra.counts` returns both counts in one pass, and `Counts` carries
  them. The report used to call `countClosed` and `countProperIdeals`
  separately, and at `n = 5` that would have doubled an already expensive walk.
  Keying on the fixed point is sound because every blade-spanned two-sided
  ideal is closed: for `a, b` in an ideal the product `a·b` lies in it.
- Three tests, one per claim:
  - `the closure walk reproduces the subset scan on every n <= 4` — all 34
    signatures, both counts, walk versus scan. This is the licence for using
    the walk where the scan cannot go.
  - `the closure walk carries P4 past n = 4` — the `r = 0` column
    `3, 6, 17, 68, 375` and the five degenerate counts at `n = 5`.
  - `a buffer that is too small is filled, not overrun` — pins the
    `enumerateTriples` contract and, with it, why the report clamps its range.

### Changed

- **The report's default range is the whole supported range.** `zig build
  explore` now writes the 55-row table for `p+q+r <= 5`. Reaching `n = 5` while
  the default stopped at 4 was a scope limit in the one place the project
  promises none. Cost, measured: 1 minute 41 seconds, almost all of it in the
  three most degenerate rows and doubled by the report being rendered twice
  (markdown and JSON).
- **`explore-full` is gone.** It existed only to reach `p+q+r = 5` while the
  default stopped at 4; with the default at 5 it would have been a second name
  for the same command. `zig build run -- explore --max 4` is the fast subset.
- The report's method section now describes both P4 methods and states the
  scope of the blade enumeration separately from the procedure, because the
  blade scope and the dimension limit are different limits and were sitting in
  one sentence.

### Fixed

- **Three Polish strings, one of them in the published report.** The
  English-only pass of b8af845 claimed no Polish text remained in `src/`; that
  was false. `src/mrs/signature.zig` (a comment), `src/bench/t4_derivative.zig`
  (a table row label, which reached `results/RESULTS.md` as
  `MRS: liczby dualne`) and the generated report itself. None of the three
  contains a diacritic, which is why both the `pl_PL` dictionary scan and the
  hand read walked past them. The 0.1.1 entry is corrected rather than left
  standing.
- **Eight references to a `docs/` directory that is not in the repository** —
  eight call sites, including one printed into the user-visible `demo` output
  (`docs/PLAN-MRS-0.1.md`). They now point at files that exist: `README.md`,
  `CHANGELOG.md`, `results/RESULTS.md`, and `mrs/causal.zig` for the theorem
  numbers.
- Three lines truncated by the translation pass, each of which had left its own
  fragment duplicated behind it (`The blades.`, `thesis T3 (docs/07).`, `a
  theorem: it is a pattern to investigate`).
- **A range wider than the engine supports used to truncate the table
  silently.** `enumerateTriples` stops at the end of the buffer it is given, so
  `explore --max 6` printed 55 rows under a heading claiming `p+q+r from 1 to
  6`. The range is now clamped to `MAX_TOTAL` before it is used, so the heading
  and the rows cannot disagree.

### Known follow-ups

- The report computes the sweep twice — once for the markdown table, once for
  the JSON. At range 4 that is invisible; at range 5 it is about half of the
  1.5 minutes. Rendering both from one pass is a contained refactor that was
  not worth taking on in the same change as the enumeration.
- `enumerateClosed` still refuses above `n = 4`, because it materialises one
  `Info` per closed subspace and `(0,0,5)` has 31 242 668 of them. Counting does
  not need the list; a caller who needs the `n = 5` list has to consume the sets
  as the walk produces them. The limitation is now stated at the function.
- The four largest `n = 5` assertions run only in ReleaseFast, because together
  they cost about 50 seconds in Debug — the entire cost of the suite — while
  walking exactly the code path the two small cases already cover.
- P1 remains unimplemented as a classification: "which subgroups `H ⊆ O(p,q)`
  admit an invariant pointed convex cone" still has no LP layer behind it.

## [0.1.1] — 2026-09-15

Audit of the multi-time signatures (p >= 2), i.e. exactly the region where
causality is supposed to break: two and three time dimensions, with and without
spatial dimensions. Two defects were found there, both in the causal layer, one
of them silently wrong. Everything below is reproducible: `zig build test`
(291 tests, Debug + ReleaseSafe + ReleaseFast), `zig build verify` (16 checks),
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

- A self-contradicting sentence in the T2 report: the prose said "the ratio is 1"
  while the printed number was a noisy `1.036x`. The identity check between MRS-LAB
  and a classical rapidity variable is exact by construction, so the sentence now
  reads "the ratio is 1 up to timer noise" and matches whatever the timer prints.

- **Translation leftovers.** The English-only pass of commit b8af845 missed
  strings and comments that a line-oriented scan does not catch: debug output
  inside a test, a `@panic` message, test names, `Class.label()`, and 70 dodgy
  comment lines, 67 of which carried duplicated comment markers (`///  /// …`)
  left by the rewriting pass. All of those are translated. **The claim that no
  Polish text remained in `src/` was too strong**: three strings carried no
  diacritics and survived both the dictionary check and the hand read. They were
  found and fixed in 0.1.2, which also corrects this entry.

- **A test that wrote to stderr made the whole build step look broken.** The
  P2b/P2c observation test printed its count with `std.debug.print` on the
  SUCCESS path, and in Zig 0.16 the build runner answers any stderr from a test
  process with `failed command: … --listen=-` lines attributed to the step —
  while the build still succeeds and all tests pass. Evidence is now carried by
  assertions (`expectEqual(10, differ_a)`), so `zig build test` prints a clean
  summary, and a real failure is no longer buried in noise. Diagnosis of the
  noise is reproducible: adding a `std.debug.print` to any passing test makes
  that step emit the same line.

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
