# 1.0 API freeze

This is the public surface for version 1.0. Names fall out of the mathematics and [03-public-api.md](03-public-api.md). New families need a new constructor; they do not grow hidden arguments on `gpa()`.

The package has **one problem compiler** and **four engines**. It does not have one universal optimizer.

| Compiled family | Engine | Allowed claim |
|---|---|---|
| Pairwise \(O(d)\) / \(SO(d)\) / similarity | exact SVD / polar (`eigencore`) | `exact_closed_form` |
| Multiway orthogonal / similarity | GPM or Gower BCD | GPM: stationary, `certified_global` only after Ling’s dual test; Gower: `blockwise_stationary` |
| Weighted, robust, cell-masked | IRLS / MM | `blockwise_stationary` or `first_order_stationary` |
| Affine / LBW / TPS | variable projection + eigenproblem (`eigencore`) | `exact_closed_form` **only** for the constrained reference-space formulation |

Spectral work (pairwise polar, GPM certificates, LBW bottom eigenvectors, canonical axes) goes through **eigencore**. The package does not call base `svd()` / `eigen()` for those problems.

---

## Exported fitting

- `procrustes()`, `gpa()`, `proc_data()`, `proc_from_array()`
- `compile_proc_problem()`, `explain_solver()`
- `predict()`, `fitted()`, `residuals()`, `aligned()`, `block_apply()`

## Transformation constructors

- `proc_orthogonal()`, `proc_similarity()`, `proc_signed_permutation()`
- `proc_affine()`, `proc_lbw()`, `proc_tps()`

`O(d)` and `SO(d)` remain different typed groups. Affine / TPS are LBW constructors for \(\Phi\) and \(L\), not a second GPA world.

`reference_covariance` (\(\Lambda\) in \(M^\top M=\Lambda\)) is required to be explicit on LBW specs. It is never estimated invisibly.

## Metric, loss, gauge

- `proc_metric()` / `proc_weights()`, `proc_covariance()`
- `proc_squared_l2()`, `proc_huber()`, `proc_tukey()`
- `proc_gauge()`, `gpa_control()`

\(\Sigma_S\) is a fitting argument. \(\Sigma_M\) is an `infer()` argument.

## Accessors and geometry

- `consensus()`, `transformations()`, `coef()`
- `align_gauge()`, `canonicalize()`, `diagnose()`, `certify()`, `decompose()`
- `quotient_distance()`, `principal_angles()`, `subspace_distance()`
- `horizontal_component()`, `vertical_component()`, `tied_axis_blocks()`
- `tangent_coordinates()`, `tangent_project()`
- `consensus(fit, gauge = "native"|"canonical"|reference)`
- `apply_proc_transform()`, `inverse_proc_transform()`, `compose_proc_transform()`
- `datum_space_error()` — refused when no exact inverse exists

## Statistical layer

- `infer()`, `proc_shape_model()`
- `cross_validate()`, `tune_smoothness()`
- `tidy.gpa_fit()`, `glance.gpa_fit()`, `augment.gpa_fit()` — existing quantities only; broom registration on load

## Plots

- `plot_data()`, `autoplot.gpa_fit()`
- Types: overlay, residuals, support, convergence, influence, weights, uncertainty, decomposition, deformation, cv

## Not exported, not hidden inside `gpa()`

ICP, unknown correspondence, optimal transport, sliding semilandmarks, mesh repair, image registration, generic manifold or SDP modelling, PCA / CCA / NMF.

Semi-orthogonal / Stiefel alignment of unequal \(p_i\) is **post-1.0**.

---

## Guarantee vocabulary (do not collapse)

- Numerical status and optimality status stay separate.
- Ling’s dual certificate is not reused for \(SO(d)\), Huber / Tukey, cell masks, or LBW.
- LBW success is reported as a global eigen-solution **for the stated formulation**, never “GPA is globally solved.”
- Certificate unavailable never becomes “probably global.”
