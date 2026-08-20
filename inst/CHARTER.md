# gprocrustes charter

The goal is not “an R package containing many Procrustes functions.” It is a compact, mathematically explicit **alignment engine** whose results are fast, inspectable, numerically defensible, and honest about what has — and has not — been optimized.

The four source papers define four obligations:

- **Gower.** GPA must remain genuinely symmetric: all configurations move toward an estimated consensus. Translation, rotation, reflection, scaling, consensus, and residual variation are an interpretable decomposition, not merely an opaque optimizer.
- **Goodall.** Alignment is only the beginning. Weighting, covariance, robustness, shape-space geometry, estimation, and inference have to be distinguished cleanly. The metric used to superimpose observations need not be the statistical model metric used for inference.
- **Ling.** Generalized orthogonal Procrustes is nonconvex and generally NP-hard. An ordinary alternating solution must not be labelled globally optimal. Offer spectral initialization, the generalized power method, and — where possible — an optimality certificate.
- **Bai and Bartoli.** One disciplined deformable extension: linear-basis warps, especially affine and thin-plate splines, with regularization, incomplete shapes, cross-validation, and precise statements about the formulation for which the closed-form solution is global.

The existing R ecosystem (`vegan`, `shapes`, `geomorph`, `multiblock`) is not empty. Their interfaces are oriented toward ordination, landmark-array shape analysis, geometric morphometrics, or multiblock analysis. The opening is to build the first general, matrix-like, sparse-aware, missing-aware, solver-transparent GPA engine.

The mathematics is the primary specification: [docs/spec/00-mathematics.md](../docs/spec/00-mathematics.md).

---

## Eight governing rules

1. **Consensus first.** Multi-configuration fitting estimates a consensus. An arbitrary reference is allowed only as an explicitly requested anchored problem.
2. **Semantics before algorithms.** Transformation group, scaling convention, reflection policy, correspondence, weights, missingness, loss, regularization, and gauge are all represented explicitly.
3. **Guarantees are data.** Every fit reports whether its solution is exact, certified global, first-order stationary, merely monotone, or unconverged.
4. **Matrix-like, not array-bound.** The native representation is a named list of entity-by-dimension matrix-like objects — not a required \(p\times k\times n\) dense array.
5. **Missingness is native.** Missing landmarks and missing coordinates are represented by masks and incidence maps, not silently imputed before fitting.
6. **Sparse means computationally sparse.** Sparse inputs remain sparse through centering and sufficient-statistic calculations. The package never promises that a rotated output remains sparse.
7. **Statistics are separate from fitting.** The superimposition metric and model covariance are separate objects. Bootstrap, permutation, and analytic inference declare their assumptions.
8. **Plots are diagnostics.** Every principal plot should expose fit, uncertainty, missing support, overfitting, influence, or numerical behaviour — not merely produce an attractive overlay.

This should result in perhaps 20–25 important exported functions, not 150 loosely related utilities.

---

## Explicit non-goals

The package should not become: an ICP or unknown-correspondence library; an optimal-transport library; a complete geometric-morphometrics environment; a semilandmark digitization / sliding package; a mesh repair package; an image-registration application; a NIfTI or microscopy file-format package; a generic manifold-optimization or SDP modelling package; a PCA, CCA, NMF, or dimensionality-reduction package; a catalog of every deformation ever proposed.

Known correspondence is the core contract.

---

## Development sequence

- **Milestone 0.** Specification and conformance suite (this repository state).
- **Milestone 1.** Typed transforms, pairwise \(O(d)/SO(d)/\)similarity, sparse moments, fit object, lazy aligned views.
- **Milestone 2.** Symmetric Gower-style generalized similarity GPA.
- **Milestone 3.** Matrix-free GPM, dual certificate, overlap graph, backends.
- **Milestone 4.** Weights, robustness, cellwise missingness.
- **Milestone 5.** Statistical layer and adapters.
- **Milestone 6.** Affine / LBW / TPS module, then freeze the 1.0 API.
