# Public API and fit contract

About 20–25 exported functions. The API falls out of [00-mathematics.md](00-mathematics.md) and [00-conventions.md](00-conventions.md). No solvers are implemented in Milestone 0; this document freezes names and object shapes.

---

## Data

```r
dat <- proc_data(
  views = list(site_a = X_a, site_b = X_b),
  ids = list(site_a = ids_a, site_b = ids_b),
  observed = list(site_a = mask_a, site_b = mask_b)
)
```

Convenience: `gpa(list(X1, X2, X3), transform = "similarity")`.

Each view: rows = entities, columns = dimensions; local-to-global row map; optional observation mask; optional topology for plotting only; optional metadata.

---

## Constructors

Typed objects, not undocumented matrices or a single `weights` argument.

```r
proc_orthogonal(group = c("O", "SO"))
proc_similarity(group = c("SO", "O"), translation = TRUE, scaling = "isotropic")
proc_signed_permutation()
proc_affine(reference_covariance = 1)
proc_lbw(basis, penalty, smoothness = 0, reference_covariance = 1)
proc_tps(control_points = NULL, smoothness = 1, reference_covariance = 1)

proc_gauge(scale = c("gower", "preshape", "fixed_consensus", "none"),
           orientation = c("free", "principal"))

proc_metric(configuration = NULL, landmark = NULL, cell = NULL, precision = NULL)
proc_weights(...)  # alias that forwards to the four channels

proc_huber(level = "landmark", k = 1.345)
proc_tukey(level = "landmark", c = 4.685)
proc_squared_l2()

gpa_control(tolerance = 1e-8, max_iterations = 500L,
            keep_aligned = "lazy", certify = "auto")
```

String shortcuts (`transform = "similarity"`) expand to the corresponding constructor.

---

## Fitting

```r
fit <- gpa(
  dat,
  transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
  gauge = proc_gauge(scale = "gower", orientation = "free"),
  metric = proc_metric(configuration = configuration_weights, landmark = landmark_weights),
  loss = proc_huber(level = "landmark", k = 1.345),
  solver = "auto",
  control = gpa_control()
)
```

Anchored problems (one configuration treated as fixed) require an explicit argument such as `anchor = "site_a"`. The default is consensus-first.

`compile_proc_problem()` / `explain_solver()` expose the planner ([01-solver-guarantee-matrix.md](01-solver-guarantee-matrix.md)).

Advanced: `proc_from_statistics(moments, crossproducts, overlap_graph)` for streaming / map-reduce / privacy-preserving reduced statistics.

---

## Fit object (`gpa_fit`)

Must contain:

- problem specification
- data / correspondence summary
- consensus
- typed transformations
- gauge
- objective and energy decomposition
- convergence history
- rank and conditioning diagnostics
- missing-support / overlap-graph diagnostics
- final robust / precision weights
- solver and backend
- numerical timings
- numerical status and optimality status (separate)
- certificate, when available
- structured warnings
- call and package version

Must **not** automatically store every dense aligned matrix.

Print example:

```text
Generalized Procrustes fit
  Configurations: 48
  Consensus entities: 59,412
  Target dimension: 20
  Transformation: SO(20), no scale, implicit centering
  Missing observations: 7.8%
  Overlap graph: connected; minimum overlap rank 20
  Solver: spectral GPM
  Iterations: 14
  Numerical status: converged
  Optimality: certified global
  Dual minimum-eigenvalue bound: 2.4e-7
  Aligned data: lazy
```

---

## Generics

```r
print(fit)
summary(fit)

consensus(fit)
transformations(fit)
aligned(fit)                    # lazy proc_aligned_view
fitted(fit)
residuals(fit, level = "landmark")
coef(fit)
predict(fit, newdata = X_new)   # align to a fixed consensus

diagnose(fit)
certify(fit)
decompose(fit)
canonicalize(fit)
align_gauge(x, reference)

explain_solver(fit)             # also works on a compiled problem

tangent_coordinates(fit)
tangent_project(newdata, fit)   # labelled local

infer(fit, model = proc_shape_model(...), method = "bootstrap")
```

Broom:

```r
tidy(fit)     # transformation-level parameters
glance(fit)   # fit-level diagnostics
augment(fit)  # entity / configuration residuals and weights
```

---

## Transform algebra

```r
tr <- transformations(fit)[["site_a"]]
apply_proc_transform(tr, X_a)
inverse_proc_transform(tr)
compose_proc_transform(tr1, tr2)
```

`inverse_proc_transform()` on a noninvertible warp must say that no exact inverse exists. Do not construct a fiction.

`datum_space_error()` is a separate generic from the reference-space objective ([math §34](00-mathematics.md)).

---

## Lazy aligned views

`aligned(fit)` returns `proc_aligned_view` objects storing the original matrix-like object, centering, scale, rotation, translation, row map, and mask. They materialize only requested rows or blocks:

```r
aligned(fit)[["site_a"]][1:1000, ]
block_apply(aligned(fit)[["site_a"]], FUN)
```

A general rotation makes sparse columns dense. The package says so. It never promises that rotated output remains sparse.

---

## Plots and plot data

All 2D plots return ggplot2 objects. Core types diagnose something: `overlay`, `residuals`, `support`, `convergence`, `influence`, `weights`, `uncertainty`, `pairwise`, `spectrum`, `decomposition`, `deformation`, `cv`, `certificate`.

```r
autoplot(fit, type = "overlay")
plot_data(fit, type = "residuals")
```

Large data: aggregate or sample explicitly and report what was done.

---

## Component matching

```r
match_components(x, reference, transform = "signed_permutation")
```

Full orthogonal rotation of latent subspaces is mathematically legitimate and not always scientifically legitimate.

---

## Non-goals (not exported, not hidden inside `gpa()`)

ICP / unknown correspondence, optimal transport, sliding semilandmarks, mesh repair, image registration, NIfTI I/O, generic manifold or SDP modelling, PCA / CCA / NMF.
