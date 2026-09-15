# MRS-LAB — property engine results

Scope: all signatures (p,q,r) with p+q+r from 1 to 4. Convention: (+,−,−,…); the P2 results are convention independent (checked by test).

## Method — why these are decisions, not evidence

- **P2a/P2b/P2c**: for fixed y the identity is a **quadratic function** of x, and symmetrically, and a quadratic function (including a linear term) is determined by its values on {0} ∪ {e_i} ∪ {e_i+e_j}. Checking the Cartesian product of that set **decides** the identity for all real x,y. The arithmetic is exact integers — no tolerance, no probability.
- **P4**: exhaustive over all subspaces spanned by blades. That is a sublattice of all subalgebras — ideals spanned by idempotents (e.g. in Cl(1,0) ≅ R⊕R) are outside this scope and are not visible here.
- **Radical ideal**: the blades containing a degenerate generator. The index is the smallest k with I^k = 0.

## Table

| (p,q,r) | n | P2a | P2b | P2c | centre dim | subalgebras | proper ideals | radical idx | notes |
|---|---:|---|---|---|---:|---:|---:|---:|---|
| (0,0,1) | 1 | yes | yes | yes | 2 | 4 | 1 | 2 | radical ideal, I^2 = 0 |
| (0,1,0) | 1 | yes | yes | yes | 2 | 3 | 0 | — | no blade-spanned ideals |
| (1,0,0) | 1 | yes | yes | yes | 2 | 3 | 0 | — | no blade-spanned ideals |
| (0,0,2) | 2 | yes | yes | yes | 2 | 14 | 4 | 3 | radical ideal, nilpotent |
| (0,1,1) | 2 | yes | yes | yes | 1 | 10 | 1 | 2 | radical ideal, I^2 = 0 |
| (0,2,0) | 2 | yes | yes | yes | 1 | 6 | 0 | — | no blade-spanned ideals |
| (1,0,1) | 2 | yes | yes | yes | 1 | 10 | 1 | 2 | radical ideal, I^2 = 0 |
| (1,1,0) | 2 | yes | yes | yes | 1 | 6 | 0 | — | no blade-spanned ideals |
| (2,0,0) | 2 | yes | yes | yes | 1 | 6 | 0 | — | no blade-spanned ideals |
| (0,0,3) | 3 | yes | yes | yes | 5 | 136 | 18 | 4 | radical ideal, nilpotent |
| (0,1,2) | 3 | yes | yes | yes | 3 | 89 | 4 | 3 | radical ideal, nilpotent |
| (0,2,1) | 3 | yes | yes | yes | 2 | 46 | 1 | 2 | radical ideal, I^2 = 0 |
| (0,3,0) | 3 | **no** | yes | yes | 2 | 17 | 0 | — | no blade-spanned ideals |
| (1,0,2) | 3 | yes | yes | yes | 3 | 89 | 4 | 3 | radical ideal, nilpotent |
| (1,1,1) | 3 | yes | yes | yes | 2 | 46 | 1 | 2 | radical ideal, I^2 = 0 |
| (1,2,0) | 3 | **no** | yes | yes | 2 | 17 | 0 | — | no blade-spanned ideals |
| (2,0,1) | 3 | yes | yes | yes | 2 | 46 | 1 | 2 | radical ideal, I^2 = 0 |
| (2,1,0) | 3 | **no** | yes | yes | 2 | 17 | 0 | — | no blade-spanned ideals |
| (3,0,0) | 3 | **no** | yes | yes | 2 | 17 | 0 | — | no blade-spanned ideals |
| (0,0,4) | 4 | yes | **no** | **no** | 8 | 9238 | 166 | 5 | radical ideal, nilpotent |
| (0,1,3) | 4 | yes | **no** | **no** | 4 | 5004 | 18 | 4 | radical ideal, nilpotent |
| (0,2,2) | 4 | yes | **no** | **no** | 2 | 2220 | 4 | 3 | radical ideal, nilpotent |
| (0,3,1) | 4 | **no** | **no** | **no** | 1 | 654 | 1 | 2 | radical ideal, I^2 = 0 |
| (0,4,0) | 4 | **no** | **no** | **no** | 1 | 68 | 0 | — | no blade-spanned ideals |
| (1,0,3) | 4 | yes | **no** | **no** | 4 | 5004 | 18 | 4 | radical ideal, nilpotent |
| (1,1,2) | 4 | yes | **no** | **no** | 2 | 2220 | 4 | 3 | radical ideal, nilpotent |
| (1,2,1) | 4 | **no** | **no** | **no** | 1 | 654 | 1 | 2 | radical ideal, I^2 = 0 |
| (1,3,0) | 4 | **no** | **no** | **no** | 1 | 68 | 0 | — | no blade-spanned ideals |
| (2,0,2) | 4 | yes | **no** | **no** | 2 | 2220 | 4 | 3 | radical ideal, nilpotent |
| (2,1,1) | 4 | **no** | **no** | **no** | 1 | 654 | 1 | 2 | radical ideal, I^2 = 0 |
| (2,2,0) | 4 | **no** | **no** | **no** | 1 | 68 | 0 | — | no blade-spanned ideals |
| (3,0,1) | 4 | **no** | **no** | **no** | 1 | 654 | 1 | 2 | radical ideal, I^2 = 0 |
| (3,1,0) | 4 | **no** | **no** | **no** | 1 | 68 | 0 | — | no blade-spanned ideals |
| (4,0,0) | 4 | **no** | **no** | **no** | 1 | 68 | 0 | — | no blade-spanned ideals |

## Summary

- signatures in scope: **34**
- P2a holds in 21, fails in 13
- P2b holds in 19, fails in 15
- P2c holds in 19, fails in 15
- degenerate signatures (r>0): 20
- Lorentzian signatures: 4
- P4 rows computed exhaustively: 34, out of range: 0
