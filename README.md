# gprocrustes

**gprocrustes** is an R package that aligns two or more tables of corresponding points — landmarks, stimuli, vertices — into a shared shape. You get the aligned tables, the transform that produced each one, and a status that says whether the solution is exact or only a place the solver stopped.

Use it when you already know which row in one table is the same entity as which row in another. It will not guess a matching.

> **Status:** Pre-release, GitHub only (`0.0.1.9000`). Not on CRAN.

## Quick start

```r
# install.packages("remotes")
# remotes::install_github("bbuchsbaum/gprocrustes")

library(gprocrustes)

# Four points of a unit square. Rows are points; columns are x, y.
X <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
Y <- X %*% matrix(c(0, 1, -1, 0), 2, 2, byrow = TRUE)   # same square, rotated 90°

fit <- procrustes(X, Y, transform = proc_orthogonal("SO"))
fit$objective            # 0
fit$optimality_status    # "exact_closed_form"
```

`SO` means rotation only (no reflection). The stored transform applied to `X` recovers `Y`.

Several recordings of the same points go through `gpa()`. It estimates one consensus instead of treating one recording as the truth:

```r
g <- gpa(list(A = X, B = Y), transform = proc_orthogonal("O"))
g$optimality_status      # "certified_global" on this noiseless pair
```

`certified_global` means a separate check confirmed a global orthogonal solution. Most real fits report only that they converged to a stationary point. The package will not upgrade that to “probably global.”

If your matrices are features × observations, transpose them first.

## Missing landmarks

Name the rows when the tables do not share every point. The consensus has one row per named entity.

```r
dat <- proc_data(
  list(A = X[1:3, ], B = Y[c(1, 2, 4), ]),
  ids = list(A = c("nw", "ne", "se"), B = c("nw", "ne", "sw"))
)
nrow(consensus(gpa(dat, transform = proc_orthogonal("O"))))   # 4
```

## What it covers

- Align two configurations (rotation, optional reflection, scale, translation)
- Align many configurations to one consensus
- Partial overlap: some points present in only some recordings
- A status on every fit: exact, certified, stationary, or unconverged
- Overlay and residual plots for the fit, not a presentation layer

`vegan` compares two ordinations. `shapes` and `geomorph` register landmark arrays for morphometrics. This package is for matrix-shaped data when you need the correspondence, the missingness, and the claim to stay explicit.

It does not search for unknown matches (ICP, transport) and it is not a full geometric-morphometrics environment.

## Documentation

- [Aligning landmark configurations](vignettes/gprocrustes.Rmd) — first workflow
- [Mathematics and guarantees](docs/spec/00-mathematics.md) — the specification the API is built from
- [1.0 API](docs/spec/04-api-freeze.md) — exported names and allowed claims

## License

MIT. [Brad Buchsbaum](mailto:brad.buchsbaum@gmail.com).
