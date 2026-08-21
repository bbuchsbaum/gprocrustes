# gprocrustes

[Charter](inst/CHARTER.md) · [Mathematics](docs/spec/00-mathematics.md) · [API freeze](docs/spec/04-api-freeze.md) · [Vignette](vignettes/gprocrustes.Rmd)

**gprocrustes** is an R package for pairwise and generalized Procrustes analysis that returns a shared consensus, the transform that took each configuration there, and a status saying what was actually solved.

Use it when the same entities — landmarks, stimuli, tasks, vertices — appear in more than one recording and the correspondence is known. The package will not invent a matching.

> **Status:** Source-only pre-release (`0.0.1.9000`). The 1.0 surface is frozen in the spec; it is not on CRAN yet.

## Quick start

```r
# install.packages("remotes")
# remotes::install_github("bbuchsbaum/gprocrustes")

library(gprocrustes)

X <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
rot90 <- matrix(c(0, 1, -1, 0), 2, 2, byrow = TRUE)
Y <- X %*% rot90

fit <- procrustes(X, Y, transform = proc_orthogonal("SO"))
fit
#> Pairwise Procrustes fit
#>   Group: SO
#>   Effective rank: 2 of 2
#>   Transform uniqueness: yes
#>   Objective uniqueness: yes
#>   Unidentified subspace dimension: 0
#>   Objective: 0
#>   Numerical status: exact
#>   Optimality: exact_closed_form
```

`O(d)` and `SO(d)` are different groups. This pair is an exact rotation, so the kernel reports `exact_closed_form` and the stored transform reproduces \(Y\).

Three recordings of the same four points, no noise:

```r
views <- list(
  A = X,
  B = Y,
  C = X %*% matrix(c(0, -1, 1, 0), 2, 2, byrow = TRUE)
)
g <- gpa(views, transform = proc_orthogonal("O"))
g$solver              # "gpm"
g$optimality_status   # "certified_global"
certify(g)$status     # "certified_global"
```

Unanchored multi-view problems estimate a consensus. A certificate is emitted only when Ling’s dual test succeeds. Ordinary descent is never relabelled “probably global.”

Rows are entities; columns are dimensions. The action is \(T(X)=sXR+\mathbf{1}t\). If your matrices are features-by-observations, transpose before you call the engine.

## Partial correspondence

Configurations need not share every entity. Name the rows; missing landmarks stay missing.

```r
dat <- proc_data(
  list(A = X[1:3, ], B = Y[c(1, 2, 4), ]),
  ids = list(A = c("nw", "ne", "se"), B = c("nw", "ne", "sw"))
)
partial <- gpa(dat, transform = proc_orthogonal("O"))
nrow(consensus(partial))   # 4 — the union of named entities
```

A missing **row** is not a missing **cell**. Cell masks are a different problem and get a different solver and a weaker claim.

## What it covers

- Pairwise \(O(d)\), \(SO(d)\), and similarity, including the exact residual formulas
- Symmetric consensus GPA (Gower) and matrix-free generalized power method
- An optimality certificate when the dual test has a validated lower bound; otherwise `not_certified` or `unavailable`
- Native row masks, cell masks, and sparse moments (sparse inputs stay sparse through centering)
- One deformable family: affine / thin-plate / linear-basis warps, with \(\Lambda\) stated explicitly
- Inference (`infer()`, `cross_validate()`) as a layer after fitting, not a second superimposition metric
- Diagnostics: overlay, residuals, support, influence, decomposition — not decorative plots

`vegan::procrustes` compares two ordinations. `shapes` and `geomorph` register landmark arrays for morphometrics. This package is the matrix-like engine underneath: typed groups, missingness, and claims you can read off the fit.

## Fit and boundaries

Good fit when correspondence is known, configurations are entity-by-dimension matrices (dense or sparse), and you need to know whether the number you got is exact, certified, or only stationary.

It does not do ICP, unknown matching, optimal transport, sliding semilandmarks, mesh repair, or image registration. Semi-orthogonal alignment of unequal column dimensions is post-1.0.

Install from GitHub. Spectral work goes through [eigencore](https://bbuchsbaum.github.io/eigencore/).

## Documentation

- [Aligning landmark configurations](vignettes/gprocrustes.Rmd) — the first guided workflow
- [Public API](docs/spec/03-public-api.md) and [1.0 freeze](docs/spec/04-api-freeze.md) — names and what may be claimed
- [Mathematics](docs/spec/00-mathematics.md) — primary specification (the R API falls out of it)
- [Conventions](docs/spec/00-conventions.md) — orientation, gauge, missingness, weights
- [Solver guarantee matrix](docs/spec/01-solver-guarantee-matrix.md) — the only allowed optimality claims
- [Conformance suite](docs/spec/02-conformance-suite.md) — language-neutral fixtures and external oracles
- [Source papers](docs/references/README.md) — Gower, Goodall, Ling, Bai–Bartoli

## License

MIT. Author: [Brad Buchsbaum](mailto:brad.buchsbaum@gmail.com).
