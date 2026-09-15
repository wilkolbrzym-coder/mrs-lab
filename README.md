# MRS-LAB

**A computational laboratory for testing hypotheses about metric signatures.**

MRS-LAB treats the metric signature — how many dimensions are time, how many are
space, how many carry no metric information at all — as a first-class type
parameter, and then answers structural questions about the resulting algebras and
orders **by machine**, returning either a proof or a counterexample.

It is a tool for testing hypotheses, not a theory. Every claim in the property
tables is backed by a test in this repository, and every table is reproducible
byte for byte.

**Status:** 0.1.0. The engine decides questions P2 and P4 below for every
signature with `p+q+r <= 4` in about 0.4 s. P1 is partially decided, P3 is
specified but not implemented.

---

## What you can test, and what you cannot

This is the boundary of the tool. Everything outside the left column will not
work, and the engine says so rather than guessing.

### Supported

| question | scope | method |
|---|---|---|
| **P2** norm multiplicativity, `N(xy) = N(x)N(y)` | every signature with `p+q+r <= 5` | **decision procedure**: exact integer arithmetic, no sampling, no tolerance (see below) |
| **P4** closed subalgebras, two-sided ideals, the centre, the nilpotency index of the radical ideal | every signature with `p+q+r <= 4` | **exhaustive**: all `2^(2^n) <= 65 536` blade subsets |
| **P1** is the causal order a partial order? | verified on five signatures | randomised witnesses plus a constructive counterexample |
| both sign conventions `(+,-,-,-)` and `(-,+,+,+)` | all of the above | algebraic results are convention invariant, and that invariance is itself tested |
| exact reproducibility | all reports | byte-identical output between runs |
| form evaluation, boost composition, derivatives | `n <= 1024` for the numeric layer | `f64`, with per-operation error contracts stated in the source |

### Not supported

| limitation | why |
|---|---|
| **P4 is limited to subspaces spanned by blades** | the enumeration covers subspaces spanned by blade subsets — a sublattice of all subalgebras. Split algebras such as `Cl(1,0) ≅ R⊕R` have proper ideals spanned by the idempotents `(1 ± ω)/2`, and the engine **cannot see them**, because idempotents are not blades. This is a scope limit, not a negative result. |
| **P4 above `p+q+r = 4`** | for `n = 5` there are `2^32` subsets. Enumeration refuses explicitly (`error.TooManyBlades`) instead of running silently for hours. |
| **the Hurwitz bound of 8 cannot be reproduced in full** | the engine handles Clifford algebras, which are associative. Dimensions 1, 2 and 4 are covered; dimension 8 in this family is `Cl(0,3) ≅ H⊕H`, **not** the octonions, which are not a Clifford algebra. The question "does any multiplication with a multiplicative norm exist" is therefore outside the tool. |
| **P1 as a classification** | "which subgroups `H ⊆ O(p,q)` admit an invariant pointed convex cone" is specified and computable by linear programming over the rationals, but **not implemented**. What is implemented is the partial-order test for a given signature. |
| **P3, conserved positive-definite energy** | specified, not implemented; it needs symbolic reasoning, not enumeration. |
| **physics** | no measurable predictions, none, at any point. |
| **signatures above `p+q+r = 5`** | the constants `MAX_TOTAL` and `MAX_DIM` would have to be raised, and the dense coefficient arrays grow as `2^n`. |
| **threads, GPU, other architectures** | nothing is thread-parallel; all measurements are x86_64 Linux. Timings are machine dependent and must not be compared across machines. |
| **priority and novelty** | the engine computes whether a statement holds. It does not check the literature, so no result here should be read as a claim of priority. |

### Why the P2 answers are decisions rather than evidence

For fixed `y` each norm identity is a **quadratic function** of `x` (a linear map
composed with a quadratic form), and symmetrically a quadratic function of `y`.
A quadratic function, including a linear term, is determined by its values on

```
D = {0} ∪ {e_i} ∪ {e_i + e_j : i < j}
```

so checking the identity on `D × D` **decides** it for all real `x, y`. The
arithmetic is exact integers. That is why these results are reported as decisions
and not as evidence, and it is the reason the word "probably" appears nowhere in
the output.

---

## Results at a glance

### What the engine found

Sweep over all 34 signatures with `p+q+r <= 4` (`zig build explore`):

| result | statement | status |
|---|---|---|
| **Causality** | `⪯` is a partial order **iff** the signature is Lorentzian (`p = 1`, `r = 0`) | verified on 5 signatures; for `p >= 2` the temporal set is connected, so "the future" does not exist intrinsically |
| **P2a rule** | norm multiplicativity holds **iff `p+q <= 2`** — degenerate dimensions are invisible to the scalar norm | clean and verifiable; literature priority not checked |
| **P2c rule** | `z·z̄` lands in the centre **iff `2^n <= 8`** | witness at `n = 4`: `e01 + e23` in `Cl(1,3)` yields `−2·e0123` |
| **Radical ideal** | the ideal spanned by blades containing a degenerate generator is nilpotent with index **exactly `r+1`**; `I² = 0` only when `r = 1` | not found in the literature while writing this |
| **P2b ≡ P2c** | the two agree across the whole tested range | recorded as an **unproven pattern**, not a theorem — an open question |
| subalgebra counts | e.g. 68 closed blade-spanned subalgebras for `Cl(1,3)`, 9238 for `(0,0,4)` | complete tables, rarely published |

### Honest performance verdict

Three performance claims of this project were killed by its own benchmarks. They
are reported in full in `results/RESULTS.md`.

| thesis | verdict |
|---|---|
| **T1** form evaluation, `O(n)` vs `O(n²)` | 1129× against a dense matrix, but **1.00× against a conventional implementation that detects the diagonal once**. No algorithmic gain: `O(n²)` versus `O(n)` is a difference of *representations*, and the fastest known classical method for a diagonal form is the same `O(n)` sum. |
| **T1b** form inverse, `O(n)` vs `O(n³)` | **432×** at `n = 128`. Real and structural: a dense form has no shortcut around elimination. |
| **T2** composing boosts | 2.6–2.8× over a matrix chain and 11.8–12.0× over a matrix+atanh chain, but **1.000× against a classical rapidity variable** — the best known classical method is literally the same code (an identity check, not a measurement). |
| **T2b** invariant drift | real quality difference: the additive chart holds the invariant at **2.2·10⁻¹⁶** regardless of chain length, while the matrix chain drifts diffusively (`6.2·10⁻¹⁵` at N=1000, `3.2·10⁻¹⁴` at N=100 000). |
| **T2c** conditioning | where MRS-LAB *loses*: for `\|θ\| >= 20` the pair `(cosh θ, sinh θ)` loses the invariant entirely. The rapidity must stay the primary coordinate. |
| **T3** sparse multivectors | **903×** over a dense product at `n = 12, k = 4`; **13×** over the same algorithm fed a dense array, which is exactly the cost of discovering sparsity. |
| **T3b** crossover | the cost constant derived from the data is `C = 1888`; the measured crossover sits between `k = 8` and `k = 32`, and the prediction agrees with the measurement in 8 of 8 points. |
| **T4** derivative | dual numbers give the derivative **exactly** — bit for bit identical to the analytic formula — against `2.9·10⁻⁸` for the best tuned central difference, and about 1.5× faster. |

Where MRS-LAB wins, the win survives every baseline tried: the form inverse
(T1b), invariant-preserving composition (T2b), exact derivatives (T4). Where it
does not win, this file says so.

---

## Requirements

**Zig 0.16.0.** Nothing else. No dependencies, no package manager, no network.

```bash
zig version        # must print 0.16.0
```

If Zig is not installed, download the 0.16.0 tarball for your platform from
ziglang.org/download, unpack it anywhere and add that directory to `PATH`.

---

## Quick start

```bash
zig build test          # all tests: Debug + ReleaseSafe + ReleaseFast
zig build explore       # property table P2/P4 for every signature, p+q+r <= 4
zig build explore-full  # same, up to p+q+r <= 5 (slower)
zig build bench         # performance theses T1-T4 -> results/RESULTS.md
zig build bench-check   # short benchmark run built ReleaseSafe (catches UB)
zig build demo          # mathematical walkthrough, printed to stdout
zig build verify        # consistency checks; exit code 1 on any failure
zig build run -- <cmd>  # any subcommand directly
```

Everything writes to `results/`:

| file | content |
|---|---|
| `results/EXPLORE.md` | the property table, human readable |
| `results/explore.json` | the same table, machine readable |
| `results/RESULTS.md` | benchmark report with the honest verdicts |

All three are deterministic: the same source produces the same bytes.

---

## The two relations, kept apart

This distinction drives several results above, and it is easy to get wrong:

- **separation** (unoriented): `u ⊑ v` iff `g(v−u, v−u) <= 0`. Symmetric on null
  separations, so it is *not* antisymmetric, so it is not an order.
- **causal order** (oriented): `u ⪯ v` iff `v−u` lies in the closed future cone
  with respect to a declared time arrow. The arrow is extra structure: for `p = 1`
  the form determines it up to a sign, for `p >= 2` it does not determine it at
  all.

---

## How to test your own hypothesis

Three levels, from "just run it" to "add a new question to the engine".

### 1. Look up the answer in the existing table

```bash
zig build explore
grep '(1,3)' results/EXPLORE.md
```

Every row answers: is the norm multiplicative (P2a/P2b/P2c), what is the centre
dimension, how many closed subalgebras, how many proper two-sided ideals spanned
by blades, and what is the nilpotency index of the radical ideal.

### 2. Write a property test in Zig

Hypotheses that are finite in scope are best written as a sweep over every
signature. `src/explore/` gives you the tools: `signatures.zig` enumerates
triples, `exact.zig` gives exact integer multivectors with no allocation,
`subalgebra.zig` answers closure and ideal questions, `norm_mult.zig` decides the
norm identities.

```zig
// src/explore/my_question.zig
const std = @import("std");
const mrs = @import("mrs");
const sigs = @import("signatures.zig");
const exact = @import("exact.zig");

/// HYPOTHESIS: with a degenerate dimension present, the centre is always
/// at least two-dimensional.
pub fn centreDimAtLeastTwo(alg: mrs.clifford.Algebra) bool {
    return exact.dimOfSet(exact.centerBasis(alg)) >= 2;
}

test "sweep the hypothesis over every signature with p+q+r <= 4" {
    var buf: [64]sigs.Triple = undefined;
    const n = sigs.enumerateTriples(&buf, 4);
    for (0..n) |i| {
        if (!buf[i].isDegenerate()) continue;
        const alg = try (sigs.SigBuf.build(buf[i], .mostly_minus)).algebra();
        if (!centreDimAtLeastTwo(alg)) {
            std.debug.print("REFUTED by ({d},{d},{d})\n", .{ buf[i].p, buf[i].q, buf[i].r });
            return error.HypothesisRefuted;
        }
    }
}
```

Then add the file to the test block in `src/main.zig` — **this is required**:
`zig test` collects no tests from a file that is not referenced from a test
block, and the failure mode is silent (the run reports fewer tests, it does not
complain). Run `zig build test`.

That particular hypothesis is in fact false: `(0,3,1)` has a one-dimensional
centre, because the volume element of the nondegenerate part does not commute
with the degenerate generator. The engine will tell you so, printing the
signature as the counterexample.

### 3. Add a new question to the engine

`src/explore/report.zig` has the extension point. `evalTriple` computes one row
per signature; add your predicate there and a column in `writeMarkdown`, and it
appears in `results/EXPLORE.md` and `results/explore.json` for every signature in
range. Requirements, taken from the method used by P2 and P4:

1. the predicate must be decidable within a declared range;
2. the scope must be printed next to the answer (never a bare "yes");
3. a failure must produce a witness — a signature, a vector, a matrix;
4. the row must be reproducible by `zig build test`.

Raising the range beyond `p+q+r = 4` means raising `MAX_TOTAL` (and `MAX_DIM`
for the cached sign vector); the dense coefficient arrays and the subset
enumeration grow as `2^n` and as `2^(2^n)` respectively, so `n = 5` for P4 is out
of reach by construction, not by tuning.

---

## Design rules

- One code path for every signature; the signature is a compile-time parameter.
- Conventions are always explicit, never defaulted.
- Exact arithmetic wherever the mathematics is exact.
- Proof or counterexample — never "probably".
- Scope is part of every answer, and is printed with it.
- Tests run in Debug **and** ReleaseSafe **and** ReleaseFast, because ReleaseFast
  disables the safety checks that turn undefined behaviour into a panic.
- Benchmarks are additionally smoke-run in a checked mode (`bench-check`),
  because the benchmark binary itself is built ReleaseFast.

## Repository layout

```
src/
  mrs/         signature-parameterized algebra and causality (base layer)
  explore/     the property engine: signatures, exact arithmetic, P2, P4, report
  bench/       benchmark harness and theses T1-T4
  main.zig     CLI: demo | verify | bench | explore
results/       generated tables and reports (committed on purpose)
```

## License

Apache License 2.0 — see [LICENSE](LICENSE).
