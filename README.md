# gprocrustes

A compact, mathematically explicit **alignment engine** for generalized Procrustes analysis. Results are meant to be fast, inspectable, numerically defensible, and honest about what has — and has not — been optimized.

This is not “an R package containing many Procrustes functions.” It is a solver-transparent GPA engine: consensus-first, missing-aware, sparse-aware, and gauge-aware.

**Current:** Milestone 6 — four engines (pairwise polar, Gower/GPM, IRLS/MM, LBW eigen) behind one compiler. Spectral work uses **eigencore**. The 1.0 API is frozen in [docs/spec/04-api-freeze.md](docs/spec/04-api-freeze.md).

```r
library(gprocrustes)
X <- matrix(rnorm(40), 20, 2)
Y <- X %*% matrix(c(0, -1, 1, 0), 2, 2)
fit <- procrustes(X, Y, transform = proc_orthogonal("SO"))
fit$optimality_status   # "exact_closed_form"

views <- list(A = X, B = Y, C = X %*% matrix(c(0, 1, -1, 0), 2, 2))
g <- gpa(views, transform = proc_orthogonal("O"))
g$solver                # "gpm"
g$optimality_status     # "certified_global" only if Ling's dual test succeeds
certify(g)
```

## Specification

The mathematics is the primary contract. The R API should fall out of it.

1. [docs/spec/00-mathematics.md](docs/spec/00-mathematics.md) — 42-section mathematical specification
2. [docs/spec/00-conventions.md](docs/spec/00-conventions.md) — orientation, groups, gauge, missingness, weights
3. [docs/spec/01-solver-guarantee-matrix.md](docs/spec/01-solver-guarantee-matrix.md) — allowed optimality claims
4. [docs/spec/02-conformance-suite.md](docs/spec/02-conformance-suite.md) — oracle fixtures
5. [docs/spec/03-public-api.md](docs/spec/03-public-api.md) — intended public surface
6. [inst/CHARTER.md](inst/CHARTER.md) — product charter and non-goals
7. [docs/references/](docs/references/README.md) — source PDFs (Gower, Goodall, Ling, Bai–Bartoli, …)

Convention: rows are entities, columns are dimensions, right action \(T(X)=sXR+\mathbf{1}t\).

## What it will do (1.0)

- Pairwise \(O(d)\) / \(SO(d)\) / similarity with exact closed-form kernels
- Symmetric consensus GPA (Gower), not an implicit reference
- Matrix-free generalized power method and, when the dual test succeeds, an optimality certificate
- Native row- and cell-missingness (different mathematics, different claims)
- One deformable module: linear-basis warps, including affine and thin-plate splines, with an explicit reference covariance \(\Lambda\)
- Inference as a separate layer from fitting
- SVD / eigenproblems through [eigencore](https://bbuchsbaum.github.io/eigencore/)

## What it will not do

ICP, unknown correspondence, optimal transport, sliding semilandmarks, mesh repair, image registration, or a catalog of historical warps.

## License

MIT. Author: Brad Buchsbaum.
