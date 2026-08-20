# gprocrustes

A compact, mathematically explicit **alignment engine** for generalized Procrustes analysis. Results are meant to be fast, inspectable, numerically defensible, and honest about what has — and has not — been optimized.

This is not “an R package containing many Procrustes functions.” It is a solver-transparent GPA engine: consensus-first, missing-aware, sparse-aware, and gauge-aware.

**Milestone 0 (current):** specification and language-neutral conformance suite. There is no `gpa()` yet.

## Specification

The mathematics is the primary contract. The R API should fall out of it.

1. [docs/spec/00-mathematics.md](docs/spec/00-mathematics.md) — 42-section mathematical specification
2. [docs/spec/00-conventions.md](docs/spec/00-conventions.md) — orientation, groups, gauge, missingness, weights
3. [docs/spec/01-solver-guarantee-matrix.md](docs/spec/01-solver-guarantee-matrix.md) — allowed optimality claims
4. [docs/spec/02-conformance-suite.md](docs/spec/02-conformance-suite.md) — oracle fixtures
5. [docs/spec/03-public-api.md](docs/spec/03-public-api.md) — intended public surface
6. [inst/CHARTER.md](inst/CHARTER.md) — product charter and non-goals

Convention: rows are entities, columns are dimensions, right action \(T(X)=sXR+\mathbf{1}t\).

## What it will do (1.0)

- Pairwise \(O(d)\) / \(SO(d)\) / similarity with exact closed-form kernels
- Symmetric consensus GPA (Gower), not an implicit reference
- Matrix-free generalized power method and, when the dual test succeeds, an optimality certificate
- Native row- and cell-missingness (different mathematics, different claims)
- One deformable module: linear-basis warps, including affine and thin-plate splines
- Inference as a separate layer from fitting

## What it will not do

ICP, unknown correspondence, optimal transport, sliding semilandmarks, mesh repair, image registration, or a catalog of historical warps.

## License

MIT. Author: Brad Buchsbaum.
