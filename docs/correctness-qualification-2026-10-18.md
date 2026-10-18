# Correctness qualification, 2026-10-18

Scope: the eight high-priority groups in the review of commit
`07a58e4bd59651a832bc262fff83bbb72d9e80f2`. This is numerical regression
evidence for the specified behavior, not a proof for every possible dataset or
a claim of global optimality for iterative GPA.

| Review group | Implemented behavior | Independent evidence |
| --- | --- | --- |
| Signed permutations | Source-to-target permutation agrees with the scored assignment | Exhaustive permutation/sign enumeration and literal weighted residuals in 3D and 4D, including non-involutive cycles |
| Pairwise families | Affine/LBW/TPS requests reject instead of inheriting the polar solver's exactness claim | Public typed-error tests for all three unsupported families |
| Transform and preshape constraints | Translation-prohibited similarity uses raw moments; conflicting centering rejects; preshape requests normalize before rotation fitting and cannot force an incompatible engine | Independent weighted raw SVD, shift invariance, unit norm/zero-centroid and exact Gram-matrix laws, zero-size and forced-engine counterexamples |
| TPS penalty | Constrained radial basis and penalty share the affine-orthogonal coefficient subspace; 3D bending energy is nonnegative | Independent augmented TPS saddle systems with literal kernels in 2D/3D, off-training evaluation, PSD/nullspace tests, non-affine positive energy, high-smoothness affine limit, prediction/gauge closure |
| Entity bootstrap | Multiplicity weights refit transforms; each replicate is evaluated at the fixed original landmark identities | Noiseless zero-uncertainty law; independent multiplicity-weight refits and polar alignment; shuffled local ID order with nonuniform weights; duplicate-only configuration resamples |
| Permutation inference | Values move relative to fixed IDs; smaller fitted losses form the alternative | Named/unnamed identical seeded null distributions, nontrivial null residuals, perfect correspondence p-value `1/13` for 12 draws |
| Observed-data engine | True vector loss, cell precision and robust score determine all updates; stationarity includes every free constrained block | Scalar global SO(2) search after analytic elimination; analytic Huber consensus; literal-objective central differences across rotation, translation, scale and consensus variables for L2/Huber/Tukey; public constrained tangent checks; reflection and zero-weight edge cases |
| Anchored LBW | Reject the unsupported homogeneous eigenproblem instead of dropping the anchor's linear term | Public affine/LBW/TPS anchor rejection and compiler capability checks |

Tests live in `tests/testthat/test-review-{pairwise,gauge,tps,infer,mm}.R`.
The observed-data optimizer uses local Cayley charts and limited-memory BFGS
with analytic derivatives; it does not store a dense Hessian over landmarks.
Its successful first-order status requires the true residual threshold, not
objective stagnation.

The new regressions were also run against an archive of the original commit as
a negative control: all five files detected defects, with 106 failed assertions
and seven test errors at that milestone. The original implementation failed
independent mathematical checks, not just new status/shape expectations. Final
boundary checks additionally cover conflicting centering, unsupported anchored
bootstraps, and invalid/overflowing resample counts.

Final built-package validation:

- R 4.5.1, aarch64 macOS, eigencore 1.3.0, explicit `en_US.UTF-8` locale.
- `R CMD build --no-manual` succeeded, including vignette construction.
- `R CMD check --as-cran --no-manual`: **0 errors, 0 warnings, 1 note**.
- Installed-package tests: **1,179 passing assertions, 0 failures, 0 warnings,
  0 skips**. The one check note is new-submission/development-version metadata.
- The checked tarball SHA-256 is
  `a7959b006a4e0e27dd0d58e68c9fbb9547e1f1cd3ad131dcefcabdb7aa60955b`.
- Raw check, test, and negative-control receipts are in
  `/private/tmp/gprocrustes-fixes.SJNmsN/` for this session.

Unsupported pairwise/anchor combinations and bootstrap models now fail
explicitly. Bootstrap supports complete squared-L2 orthogonal and unscaled
similarity fits with optional configuration/landmark weights. Entity bootstrap
estimates conditional alignment uncertainty at fixed original landmarks; it is
not uncertainty over new landmark locations. Robust, cell-masked, precision,
scale-constrained, anchored, and LBW bootstrap requests remain unsupported.

The four lower-priority review findings remain separate issues: ordinary Gower
collective-scale translation invariance (`gpr-am7`), matching named configuration
weights (`gpr-t3s`), weighted/mapped cross-validation (`gpr-6ql`), and explicit
weighted decomposition semantics (`gpr-knk`). No hosted CI, PDF manual,
performance benchmark, or external scientific replication is claimed here.
