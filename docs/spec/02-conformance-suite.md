# Conformance suite

Language-neutral oracles so a second implementation (R, Python, Julia, Rust, C++) can run the same tests. This suite is the first executable deliverable. No `gpa()` implementation ships in Milestone 0.

Schema: [inst/conformance/schema.json](../../inst/conformance/schema.json).
Fixtures: [inst/conformance/fixtures/](../../inst/conformance/fixtures/).
Matrices: inline JSON for tiny cases; little-endian column-major `float64` binaries under [inst/conformance/matrices/](../../inst/conformance/matrices/) when noted.

---

## Fixture contract

Each fixture is one JSON object:

- `id`, `title`, `kind`, `tags`
- `problem`: transform group, scale, translation, loss, views, optional masks / weights / incidence
- `expected`: numerical status, optimality status, objective (and tolerance), rank / uniqueness, optional `R` / `s` / `t` / consensus **up to gauge**, optional certificate result, or an expected error

Views are row-major in the JSON text only as nested arrays `[[row1], [row2], ...]`. That is display order. Semantically they are \(n\times d\) matrices under the right-action convention.

Comparisons:

- scalars: absolute + relative tolerances
- rotations: Frobenius distance after the declared gauge, or `polar` / `signed-permutation` exactness
- consensus: Procrustes distance modulo the stored gauge group
- objectives: kernel laws from [00-mathematics.md](00-mathematics.md) §4 must hold to ~machine precision on pairwise exact cases

---

## Kinds

| `kind` | What is checked |
|---|---|
| `pairwise` | Closed-form kernel against an oracle |
| `moments` | Dense / sparse / block moment identity |
| `law` | Algebraic invariance (gauge, order, composition) |
| `error` | Invalid input must fail with a declared code |
| `gower_history` | Historical / Gower-style iteration identities |
| `identifiability` | Overlap graph and rank diagnostics |
| `certificate` | Dual test must pass or fail as labelled |

---

## Required coverage (Milestone 0 authored; Milestone 1 executes)

1. Pairwise \(d=1\), where \(O(1)\) and \(SO(1)\) diverge.
2. Rank-zero and rank-deficient \(C\); repeated singular values.
3. Reflections; huge and tiny scales; zero configuration weight.
4. One observed landmark; disconnected overlap graph.
5. Nonfinite values (must error).
6. Dense / sparse / block equivalence of weighted centered cross-products.
7. Gauge invariance of the objective.
8. Gower 1975 residual identity and algorithm-history laws (Psychometrika 40:33–51).

Additional families: proper-rotation determinant correction; similarity residual formulas \(F^\star=b-\gamma^2/a\) and \(a+b-2\gamma\); signed-permutation assignment; row-mask consensus mean; energy decomposition identity.

---

## Gower 1975 historical case

Gower, J. C. (1975). Generalized procrustes analysis. *Psychometrika* 40:33–51. Local copy: [docs/references/Gower_1975_generalized_procrustes.pdf](../references/Gower_1975_generalized_procrustes.pdf).

The paper supplies:

- the residual \(S_r=\sum_i\sum_j \Delta^2(P_j^{(i)},G_j)\);
- identity \(\sum_{i<j}\Delta^2(P^{(i)},P^{(j)})=K\sum_i\Delta^2(P^{(i)},G)\) (his (2), \(K=m\));
- translation to a common centroid (placeable at the origin);
- SVD / Eckart–Young rotation, with determinant correction when reflection is excluded;
- a global energy constraint to prevent collapse;
- a worked iteration whose successive residual sums of squares he published so others could check implementations;
- an energy table he did **not** equip with sampling degrees of freedom.

Milestone 0 ships:

- `gower-identity-pairwise-clusters` — identity (2) on a small complete example;
- `gower-1975-algorithm-laws` — monotone residual decrease, common-centroid law, energy constraint, and “monotone \(\neq\) global” documentation;
- `gower-1975-published-history` — Table 2 carcass scores (9 entities, 7 dimensions, 3 judges) and Table 5 successive \(S_r\). Table 1 in the paper is the ANOVA layout, not the coordinates. The published path is an oracle for Gower's 1975 rotation-then-scale schedule; `gpa()` need not reproduce that history.

Do not invent published digits.

---

## External package oracles

Language-neutral fixtures remain the primary suite. R, Python, and Julia packages may be used as **optional second implementations** of a problem that this specification already names. They are not a second guarantee matrix and they do not author new fixture digits.

A foreign call is an oracle only when every item below matches the compiled problem:

1. Transform group: \(O(d)\) versus \(SO(d)\) versus similarity (isotropic scale, translation on or off).
2. Right action \(T(X)=sXR+\mathbf{1}t\). Left-action or “rotate the second argument” APIs must be adapted at the call site, never by silently transposing stored views.
3. Centering and scale convention (raw, centered, unit centroid size / preshape, Gower energy). A matching residual after a global scale is a convention mismatch, not a failed test.
4. Weights and missingness. Complete rows only, unless the other API implements the same mask. Implicit zeros or `NA`-as-missing are never an oracle for cell masks.
5. The objective. Compare \(F\) and aligned configurations **modulo the stored gauge**. Do not require matching \(R_i\) in a raw basis, matching iteration counts, or matching status strings.

If any item fails, skip the comparison. “Both ran GPA” is not an oracle.

Optional adapters live in Suggests / extra dependencies and skip when the package is absent. An adapter must not copy a foreign status string onto a `gpa_fit`.

### Already in this repository

| Oracle | Role |
|---|---|
| NumPy SVD in `tools/generate_conformance.py` | Independent polar factor for pairwise exact fixtures |
| Gower 1975 Tables 2 and 5 | Published carcass scores and \(S_r\) path; history for *his* schedule, not ours |

### Pairwise \(O(d)\) / \(SO(d)\) / similarity

These check the closed-form kernel ([01-solver-guarantee-matrix.md](01-solver-guarantee-matrix.md)). Claim remains `exact_closed_form` from our polar, not from the foreign function.

| Source | Problem it actually solves | Adapter notes |
|---|---|---|
| SciPy `linalg.orthogonal_procrustes(A, B)` | \(\|AR-B\|_F\), \(R\in O(d)\) | Same right action. No translation, no scale, no \(SO(d)\). |
| `qc-procrustes` (`theochem/procrustes`) `orthogonal` / `rotational` | \(\|AQ-B\|_F\); rotational is \(SO(d)\); optional translate / scale / weights | Closest Python twin. Two-sided, permutation, and softassign methods are different problems. |
| R `MCMCpack::procrustes` | \(sXR+\mathbf{1}t^\top\approx X^\star\) | Same algebra. Defaults leave translation and dilation off. |
| R `vegan::procrustes` | Rotates **Y toward X** from SVD of \(X^\top Y\) | Flip the argument order. Always centers. `scale=TRUE` is dilation. `symmetric=TRUE` unit-Frobenius-normalizes both matrices — different \(F\). |
| R `shapes::procOPA` | Ordinary pairwise Procrustes | Map `reflect` onto `O` versus `SO`. |
| Julia `LinearAlgebra.svd` | Polar factor | Same role as NumPy. Many Julia APIs store points as columns (\(d\times n\)). |
| Julia `CoordinateTransformations.kabsch` | Weighted rigid or similarity | Always \(SO(d)\) (determinant correction). Points are columns. |
| Julia `Kabsch.jl` / BioStructures `superimpose!` | \(SO(3)\) RMSD | Only \(d=3\). Confirm layout before comparing. |

SciPy `spatial.procrustes` unit-trace-standardizes **both** matrices and reports disparity. Use it only as a smoke test that disparity is ~0 on a known similarity pair. It is not an oracle for \(F^\star=b-\gamma^2/a\).

### Multiway Gower / complete \(L_2\)

These are consensus oracles for complete, dense, isotropic GPA. Compare \(\sum_i\|T_i(X_i)-M\|_F^2\) after aligning one consensus to the other with pairwise \(O(d)\). Do not require matching iteration histories. The published Gower 1975 path remains the history oracle; `gpa()` need not reproduce it.

| Source | Usable when | Not an oracle for |
|---|---|---|
| R `shapes::procGPA` | Complete array; `scale=TRUE` is unit centroid size / preshape; `reflect=FALSE` is \(SO(d)\). Array layout is \(k\times m\times n\) (landmarks × dimensions × specimens). | Relative warps (`alpha ≠ 0`), their PCA, sliding. |
| R `geomorph::gpagen` | Same, with **no** `curves` / `surfaces` | Semilandmark sliding (charter non-goal). `Proj=TRUE` is a post-fit tangent projection. |
| R `Morpho::procSym` | Same family. Default `reflect=TRUE` is \(O(d)\); `CSinit=TRUE` rescales first. | Sliding, `pairedLM` symmetry, `orp` as if it were the GPA objective. |
| `qc-procrustes.generalized` | Orthogonal multiway \(\sum_{i<j}\|A_iT_i-A_jT_j\|_F^2\), right action | Translation / scale (their loop does not apply them). Our stationarity residual. Unequal \(n\) or masks. |

### Geometry, not solvers

`geomstats` Kendall / PreShape geometry may check `quotient_distance`, preshape projection, and gauge invariance of \(F\). It is not an oracle for `certified_global`.

### Explicit non-oracles

| Our method | Why a foreign package is not an oracle |
|---|---|
| Ling GPM and dual certificate | No other library ships Ling’s dual test. `manifoldalign` GPA is left-action GPM and is a caller, not a specification. |
| Cell masks / anisotropic \(\Sigma_d\) | Morphometric APIs drop landmarks (rows), not cells. Filling `NA` with 0 is forbidden here. |
| Huber / Tukey IRLS | No matching \(\rho\)-GPA to compare stationarity against. |
| LBW / TPS eigenproblem (Bai–Bartoli) | Morpho / geomorph TPS is a warp interpolant, not the constrained reference-space eigenproblem. |
| Signed permutation | Foreign permutation / QAP / softassign assigns rows or both sides, not the signed column permutation of \(C\). |
| Unknown correspondence, OT, ICP | Charter non-goals. |
| `vegan::protest` | Permutation test of ordinations, not a kernel law. |

### Execution order when adapters are added

1. SciPy `orthogonal_procrustes` and `qc-procrustes.rotational` against existing pairwise fixtures (second polar, plus a true \(SO(d)\) library).
2. `MCMCpack::procrustes` / `vegan::procrustes` on the same pairwise cases (R-only).
3. `shapes::procGPA` / `geomorph::gpagen` on a complete 2-D or 3-D landmark array versus `gpa()` with the matching group and scale convention — residual identity only.
4. Julia `svd` and `kabsch` for weighted pairwise \(SO(d)\) after transposing column-point storage.
5. GPM, IRLS, MM, and LBW stay on laws and this package’s residuals.

---

## Regenerating fixtures

```text
python3 tools/generate_conformance.py
```

The generator writes JSON only. It uses NumPy SVD as an independent polar-factor oracle for pairwise exact cases and hand-derived values for \(d=1\), identities, and error cases. External package oracles from the previous section are not fixture authors; they consume the same JSON or a live `gpa()` / `procrustes()` call.

The release campaign in `tests/testthat/test-conformance-campaign.R` executes remaining fixture kinds, lazy versus materialized stores, dense versus sparse moments, iterative status versus stationarity residual, and the frozen 1.0 export / S3 surface.
