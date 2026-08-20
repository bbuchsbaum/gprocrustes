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
- `gower-1975-published-history` — placeholder for the published successive \(S_r\) sequence, to be filled by transcription from Table 2 of the 1975 paper (coordinates in Table 1). Until that transcription is complete the fixture records the citation and the identities that any correct implementation must satisfy.

Do not invent published digits.

---

## Regenerating fixtures

```text
python3 tools/generate_conformance.py
```

The generator writes JSON only. It uses NumPy SVD as an independent polar-factor oracle for pairwise exact cases and hand-derived values for \(d=1\), identities, and error cases.
