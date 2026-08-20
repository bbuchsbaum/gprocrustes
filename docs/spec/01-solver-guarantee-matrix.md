# Solver and guarantee matrix

This is the query planner’s source of truth. Each row states applicability, numerical method, convergence criterion, and the **only** allowed optimality claim. Derivations: [00-mathematics.md](00-mathematics.md) §§4–20, 28–37, 42.

`explain_solver(compiled)` must print the selected row, the reasons, and the execution plan. A fit must never transform “certificate unavailable” into “probably global.”

---

## Status vocabulary

Numerical status and optimality status are separate fields on every `gpa_fit` ([00-conventions.md](00-conventions.md) §11).

Iterative solvers always record, compactly:

- iteration, objective, relative change, consensus change (modulo gauge), stationarity residual, orthogonality residual, elapsed time, whether an accelerated step was accepted;
- an **exact** final objective recomputed from the fitted transforms (not the incremental surrogate).

Convergence is assessed **modulo gauge**. Tied subspaces use projectors / principal angles.

---

## Capability table

### Pairwise \(O(d)\), complete, squared \(L_2\)

- **Applies when:** two views; same \(d\); row or configuration weights only; no cell mask; loss is \(\|W^{1/2}(sXR+1t-Y)\|_F^2\) with \(s=1\) or omitted; \(R\in O(d)\).
- **Method:** weighted moments → thin SVD of \(C=X_c^\top WY_c\) → polar \(R^\star=UV^\top\).
- **Stopping:** closed form.
- **Claim:** `exact_closed_form`. Unique transform iff \(\operatorname{rank}(C)=d\).

### Pairwise \(SO(d)\), complete, squared \(L_2\)

- **Applies when:** as above, \(R\in SO(d)\).
- **Method:** same SVD, then \(R^\star=UDV^\top\) with last diagonal entry \(\det(UV^\top)\).
- **Stopping:** closed form.
- **Claim:** `exact_closed_form`. Uniqueness is not the \(O(d)\) criterion; the determinant constraint can identify one extra dimension. Compute it.

### Pairwise similarity, standard row weights

- **Applies when:** translation and/or isotropic \(s>0\); otherwise as pairwise orthogonal.
- **Method:** eliminate \(t^\star=\bar y-s\bar x R\); polar or proper polar; \(s^\star=\operatorname{tr}((R^\star)^\top C)/a\), projected onto \(s\ge 0\).
- **Stopping:** closed form.
- **Claim:** `exact_closed_form`. Residual \(F^\star=b-\gamma^2/a\) (with scale) or \(a+b-2\gamma\) (without).

### Signed permutation

- **Applies when:** \(R=PD\) with permutation \(P\) and sign diagonal \(D\).
- **Method:** linear assignment on \(|C_{jk}|\), then \(D_{kk}=\operatorname{sign}(C_{P(k),k})\).
- **Stopping:** closed form.
- **Claim:** `exact_closed_form`.

### Generalized orthogonal, complete \(L_2\) (GOPP)

- **Applies when:** \(K\ge 3\) (or any multi-view orthogonal); complete rows after implicit centering; squared Frobenius; no coordinate-specific precision; \(R_i\in O(d)\).
- **Method:** spectral initialization (leading \(d\)-eigenspace of the implicit operator \(C\), never forming \(C=B^\top B\)) + generalized power method: \(Z=\sum_j X_jR_j\), \(B_i=X_i^\top Z\), \(R_i\leftarrow\operatorname{polar}(B_i)\). Optional dual certificate after a fixed point.
- **Stopping:** gauge-invariant \(d_F(S^{(t+1)},S^{(t)})\); relative objective change; \(r_{\mathrm{stat}}=\max_i\|\operatorname{skew}\{R_i^\top(CS)_i\}\|_F\); \(r_{\mathrm{orth}}=\max_i\|R_i^\top R_i-I\|_F\).
- **Claim:** `first_order_stationary` or `blockwise_stationary` by default. `certified_global` **only** if \(CS=\Lambda S\) and a lower bound on \(\lambda_{\min}(\Lambda-C)\) is nonnegative (Ling). Uniqueness modulo \(O(d)\) only if \(\lambda_{d+1}(\Lambda-C)>0\). High-SNR linear convergence is a theorem about Ling’s model, not a universal guarantee.

### Generalized similarity (Gower BCD)

- **Applies when:** multi-view similarity; Gower / preshape / fixed-consensus scale convention declared; row-complete or row-masked \(L_2\).
- **Method:** safeguarded block coordinate descent. Each transform update is exact pairwise Procrustes against the current consensus; consensus is the weighted mean. Inits: medoid, sequential, spectral; optional deterministic multistart. Safeguarded Anderson acceleration with rollback if the true objective rises.
- **Stopping:** monotone \(F^{(t+1)}\le F^{(t)}\) after rollback; gauge-invariant consensus / transform change; blockwise stationarity residual.
- **Claim:** `blockwise_stationary` if converged. **No general global claim.** Converged \(\neq\) globally optimal.

### Row-masked or robust GPA

- **Applies when:** whole-landmark masks and/or landmark-vector \(\rho\); not arbitrary cell masks.
- **Method:** majorization / IRLS. Inner updates are exact weighted pairwise or consensus formulas ([math §§18–19, 23](00-mathematics.md)).
- **Stopping:** monotone surrogate and nonincrease of the true objective; IRLS weight change; stationarity of the weighted problem.
- **Claim:** `blockwise_stationary` (monotone surrogate descent). Stronger claims depend on convexity of \(\rho\). Tukey is nonconvex.

### Cell mask or anisotropic coordinate metric

- **Applies when:** \(\Omega_i\) differs by coordinate, or \(\Sigma_d\not\propto I\).
- **Method:** MM (fill unobserved target coordinates with current fitted values; solve the complete surrogate; accept only a true-objective decrease) or product-manifold Riemannian optimization.
- **Stopping:** true observed-data objective decrease; first-order Riemannian residual.
- **Claim:** `first_order_stationary`. Never claim a closed-form SVD solution.

### Affine / LBW, Bai–Bartoli reference-space constraints

- **Applies when:** \(T_i(X)=\Phi_i B_i\); quadratic \(\mu_i\operatorname{tr}(B_i^\top L_i B_i)\); constraints \(\mathbf{1}^\top M=0\), \(M^\top M=\Lambda\); free-translation check \(\Phi_i a_i=\mathbf{1}\), \(L_i a_i=0\) passes; reference-space loss.
- **Method:** variable projection. Eliminate \(B_i^\star=(\Phi_i^\top\Phi_i+\mu_i L_i)^\dagger\Phi_i^\top M\); form implicit \(P=\sum_i(I-H_i)\); bottom eigenvectors of \(P\) (or \(P_\nu=P+\nu\mathbf{1}\mathbf{1}^\top\)). Never form explicit inverses.
- **Stopping:** eigen-gap / residual of the constrained eigenproblem.
- **Claim:** `exact_closed_form` **only** for the stated formulation after assumptions are checked. Call this a global eigen-solution for that formulation, not “GPA is globally solved.”

### LBW outside those assumptions

- **Applies when:** datum-space loss, failed free-translation check, or other broken hypotheses.
- **Method:** alternating or manifold optimization.
- **Stopping:** first-order residual of the declared objective.
- **Claim:** `first_order_stationary`.

### Small orthogonal GPA requiring proof

- **Applies when:** user requests `certify="sdp"` or the problem is tiny and a dual bound from GPM is inconclusive.
- **Method:** SDP relaxation or low-rank Burer–Monteiro.
- **Stopping:** solver tolerance plus tightness check.
- **Claim:** `certified_global` if the relaxation / certificate is tight; otherwise `converged_unverified`.

### Semi-orthogonal / Stiefel (unequal \(p_i\))

- **Applies when:** \(X_i\in\mathbb R^{n\times p_i}\), \(Q_i\in\operatorname{St}(p_i,q)\).
- **Method:** Riemannian optimization on \(\prod_i\operatorname{St}(p_i,q)\). Whitened inputs may reuse the orthogonal trace-sum path.
- **Stopping:** Riemannian gradient residual.
- **Claim:** `first_order_stationary`.

---

## `explain_solver()` shape

```text
Selected solver: generalized power method
Reasons:
  transformation is orthogonal
  loss is squared L2
  configurations are complete after implicit centering
  no coordinate-specific precision
  problem exceeds dense block-Gram threshold

Execution plan:
  spectral initialization via implicit operator
  matrix-free GPM
  blockwise q x q polar projections
  dual certificate after convergence

Allowed optimality claim:
  first-order stationary; certified global only if dual test succeeds
```

---

## Certificate report

On success:

```text
Optimality: certified global
Dual stationarity residual: ...
Lower bound on lambda_min(Lambda - C): ...
Nullity: ..., expected ...
Uniqueness modulo global O(d): supported | not supported
```

On failure:

```text
Optimality: not certified
Reason: estimated smallest dual eigenvalue = ...
Numerical convergence: yes | no
Interpretation: stationary candidate only
```

The \(O(d)\) dual certificate is not reused for \(SO(d)\), robust losses, or cellwise missingness.
