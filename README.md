# MRS-LAB

**A computational laboratory for testing hypotheses about metric signatures.**

MRS-LAB treats the metric signature — how many dimensions are time, how many are
space, how many carry no metric information at all — as a first-class type
parameter, and then answers structural questions about the resulting algebras and
orders **by machine**, returning either a proof or a counterexample.

It is a tool for testing hypotheses, not a theory. Every claim in the property
tables is backed by a test in this repository, and every table is reproducible
byte for byte.

**Status:** 0.1.0 — the engine answers P2 and P4 (below) for every signature with
`p+q+r <= 4` in about 0.4 s. P1 and P3 are specified but not fully decided.

---

## What it is / is not

| it is | it is not |
|---|---|
| a property engine that decides or refutes statements about `Cl(p,q,r)` | a physical theory — it makes no measurable predictions |
| a signature-parameterized algebra layer (forms, Clifford, split-complex, dual numbers) | a faster mathematics — see [honest verdict](#honest-performance-verdict) |
| a reproducible benchmark suite with deliberately strong baselines | a replacement for linear algebra, analysis or topology |
| a place to test *your* hypotheses with proof or counterexample | a source of new axioms — notation cannot create theorems |

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

This distinction is the reason for several results below, and it is easy to get
wrong:

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

### 2. Write a property test in Zig (the ordinary way)

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
with the degenerate generator. The engine will tell you so, with the signature
printed as the counterexample.

### 3. Add a new question to the engine

`src/explore/report.zig` has the extension point. `evalTriple` computes one row
per signature; add your predicate there and a column in `writeMarkdown`, and it
appears in `results/EXPLORE.md` and `results/explore.json` for every signature in
range. Requirements, taken from the method used by P2 and P4:

1. the predicate must be decidable within a declared range;
2. the scope must be printed next to the answer (never a bare "yes");
3. a failure must produce a witness — a signature, a vector, a matrix;
4. the row must be reproducible by `zig build test`.

### Why the P2 answers are decisions rather than evidence

For fixed `y` each norm identity is a **quadratic function** of `x`, and
symmetrically. A quadratic function, including a linear term, is determined by its
values on `{0} ∪ {e_i} ∪ {e_i + e_j}`. Checking the identity on the Cartesian
product of that set therefore **decides** it for all real `x, y`, and the
arithmetic is exact integers. No tolerance, no probability, no sampling. That is
why results are reported as decisions and not as evidence.

---

## What the engine found

Sweep over all 34 signatures with `p+q+r <= 4` (`zig build explore`):

| result | statement | notes |
|---|---|---|
| **Causality** | `⪯` is a partial order **iff** the signature is Lorentzian (`p = 1`, `r = 0`) | verified on 5 signatures; for `p >= 2` the temporal set is connected, so "the future" does not exist intrinsically |
| **P2a rule** | norm multiplicativity holds **iff `p+q <= 2`** — degenerate dimensions are invisible to the scalar norm | clean and verifiable; I have not checked the literature for priority |
| **P2c rule** | `z·z̄` lands in the centre **iff `2^n <= 8`** | witness at `n = 4`: `e01 + e23` in `Cl(1,3)` yields `−2·e0123` |
| **Radical ideal** | the ideal spanned by blades containing a degenerate generator is nilpotent with index **exactly `r+1`**; `I² = 0` only for `r = 1` | not found in the literature |
| **P2b ≡ P2c** | the two agree across the whole range | recorded as an **unproven pattern**, not a theorem — an open question |

---

## Honest performance verdict

Three performance claims of this project were killed by its own benchmarks. They
are reported in full in `results/RESULTS.md`.

| thesis | verdict |
|---|---|
| **T1** form evaluation, `O(n)` vs `O(n²)` | 1129× against a dense matrix, but **1.00× against a conventional implementation that detects the diagonal once**. No algorithmic gain: O(n²) versus O(n) is a difference of *representations*, and the fastest known classical method for a diagonal form is the same O(n) sum. |
| **T1b** form inverse, `O(n)` vs `O(n³)` | **432×** at n=128. Real and structural: a dense form has no shortcut around elimination. |
| **T2** composing boosts | 2.6–2.8× over a matrix chain and 11.8–12.0× over a matrix+atanh chain, but **1.000× against a classical rapidity variable** — the best known classical method is literally the same code. |
| **T2b** invariant drift | real quality difference: the additive chart holds the invariant at **2·10⁻¹⁶** regardless of chain length, while the matrix chain drifts diffusively. |
| **T2c** conditioning | where MRS-LAB *loses*: for `\|θ\| >= 20` the pair (cosh θ, sinh θ) loses the invariant entirely. The rapidity must stay the primary coordinate. |
| **T3** sparse multivectors | **903×** over a dense product; **13×** over the same algorithm fed a dense array, which is exactly the cost of discovering sparsity. |
| **T4** derivative | dual numbers give the derivative **exactly** (bit for bit identical to the analytic formula) against 2.9·10⁻⁸ for the best tuned central difference, and about 1.5× faster. |

Where MRS-LAB wins, the win survives every baseline tried: the form inverse
(T1b), invariant-preserving composition (T2b), exact derivatives (T4). Where it
does not win, this file says so.

---

## Design rules

- One code path for every signature; the signature is a compile-time parameter.
- Conventions are always explicit, never defaulted.
- Exact arithmetic wherever the mathematics is exact.
- Proof or counterexample — never "probably".
- Scope is part of every answer.
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
