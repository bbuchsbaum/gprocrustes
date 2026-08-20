# Conventions and mathematical contract

This document is the implementer-facing contract. Derivations live in [00-mathematics.md](00-mathematics.md). If the two disagree, the mathematics document wins.

Package metadata (frozen here, not rediscovered later):

- Name: `gprocrustes`
- Public API: S3
- License: MIT
- Hard dependency later: `Matrix`
- Native representation: named list of entity-by-dimension matrix-like objects, not a required \(p\times k\times n\) array

---

## 1. Matrix orientation and action

Rows are entities (landmarks, vertices, stimuli, samples, genes, time points). Columns are dimensions (coordinates, features, components).

The similarity action is **right multiplication**:

\[
T(X)=sXR+\mathbf{1}t,
\qquad t\in\mathbb R^{1\times d}.
\]

Left-action callers (`Q %*% X`, features \(\times\) observations) transpose at an adapter boundary. The engine never silently flips orientation.

Gower’s SVD convention in the 1975 paper writes \(X_1^\top X_2=U^\top FV\) and \(H=V^\top U\) for fitting \(X_2H\) to \(X_1\). In this package’s polar notation that is \(C=X_c^\top WY_c=U\Sigma V^\top\) and \(R^\star=UV^\top\), so \(XR\) fits \(X\) to \(Y\).

---

## 2. Correspondence

Safe rules:

1. Identical unique row names are used automatically.
2. Differing or partial row names require an explicit map \(\pi_i\) / incidence \(P_i\).
3. Unnamed matrices of equal size may use positional correspondence. The fit **records** that assumption.
4. Duplicated entity identifiers are an error unless an aggregation rule is supplied.

`NA` is not the mathematical definition of missingness. Missingness is a mask or incidence map.

---

## 3. Transformation groups

These are different typed objects, not a Boolean `reflection=` switch:

| Family | Constraint | Status |
|---|---|---|
| Orthogonal \(O(d)\) | \(R^\top R=I\) | Essential |
| Proper rotation \(SO(d)\) | \(R^\top R=I\), \(\det R=1\) | Essential |
| Similarity | \(s>0\), translation, \(O(d)\) or \(SO(d)\) | Essential |
| Signed permutation | \(R=PD\), \(D=\operatorname{diag}(\pm1)\) | Small high-value utility |
| Semi-orthogonal / Stiefel | \(Q^\top Q=I_q\), possibly \(p_i\neq q\) | Post-1.0 acceptable |
| Affine | \(XA+\mathbf{1}t\) under reference constraints | Essential advanced |
| Linear-basis warp | \(\Phi B\) with quadratic penalty | Essential advanced |
| Thin-plate spline | The one built-in nonlinear warp, via LBW | Essential advanced |

A new family must supply: `identity`, `apply`, `fit_pair`, `compose`, `inverse` (or declare unavailable), `degrees_of_freedom`, `constraint_residual`, `regularization`, `capabilities`.

---

## 4. Scale conventions

Unconstrained \(\{s_i,M\}\) admits the collapse \(s_i=0\), \(M=0\). The package distinguishes three legitimate models ([math §9](00-mathematics.md)):

- **Preshape.** Center and set \(\|X_i\|_F=1\), then rotate only.
- **Gower collective scaling.** Estimate all \(s_i\) subject to \(\sum_i s_i^2\operatorname{tr}(X_i^\top X_i)=c\). After rotations, \(Gs=\lambda Ds\).
- **Fixed consensus scale.** e.g. \(\|M\|_F^2=1\), then estimate each \(s_i\).

These are not interchangeable.

---

## 5. Gauge

Estimation returns an **equivalence class**. Canonicalization returns a **presentation**.

For centered orthogonal GPA,

\[
(M,\{R_i\})\sim(MQ,\{R_iQ\}),\qquad Q\in O(d).
\]

Ling’s comparison metric is \(d_F(S,T)=\min_{Q\in O(d)}\|S-TQ\|_F\).

Tied or nearly tied singular values are compared via projectors, principal angles, subspace distances, and blockwise alignment inside the tied eigenspace — never by elementwise singular-vector signs.

An optimizer must not secretly impose a principal-axis orientation while solving. `canonicalize()` and `align_gauge()` are separate operations.

A fit stores a `gauge` object describing remaining freedoms (rotation, translation, scale, and disconnected-component gauges).

---

## 6. Missingness

Distinguish:

- a genuine numerical zero;
- a missing whole landmark (row mask);
- a missing coordinate within an observed landmark (cell mask \(\Omega\));
- a landmark absent because it is outside a view’s domain;
- a structurally masked region (e.g. cortical medial wall).

Sparse implicit zeros are **never** missing.

**Row-masked** problems preserve the orthogonal trace reduction and remain a quadratic synchronization problem ([math §18](00-mathematics.md)).

**Cell-masked** problems destroy that reduction ([math §20](00-mathematics.md)). They use MM or a product-manifold solver. The slower, weaker guarantee must be reported.

Default fitting is **direct observed-data**. Completion is an explicit model component, never invisible preprocessing.

---

## 7. Weights and metrics

Never a single `weights` argument. Channels:

- `configuration` — \(\alpha_i\), reliability of a whole view;
- `landmark` — diagonal \(W_i\), reliability of entities;
- `cell` — coordinate-specific reliability;
- `precision` — operator \(Q_i\) in \(\operatorname{vec}(E)^\top Q_i\operatorname{vec}(E)\).

Separately:

- \(\Sigma_S\) — superimposition / fitting metric;
- \(\Sigma_M\) — model covariance for inference.

They need not be equal ([math §22](00-mathematics.md)).

Row metric \(\otimes I_d\) keeps the SVD fast path. Anisotropic \(\Sigma_d\not\propto I\) requires a manifold / numerical solver.

---

## 8. Losses

Default squared \(L_2\): \(L(Y,Z)=\|Y-Z\|_F^2\).

Robust losses act on the **Euclidean norm of each landmark residual vector** \(r_{ij}=\|e_{ij}\|_2\), not on coordinates. Coordinatewise Huber changes under rotation.

Core losses:

- squared \(L_2\);
- Huber, landmark-vector;
- Tukey bisquare, marked nonconvex;
- trimmed configuration loss as an advanced option.

IRLS must expose final weights and effective sample sizes.

---

## 9. Master objective

\[
\min_{M,\theta}
\sum_i\alpha_i
L_i\!\bigl(T_{\theta_i}(X_i),P_iM\bigr)
+
\sum_i\lambda_i\Omega_i(\theta_i)
\]

subject to explicit gauge and transformation constraints.

Pairwise kernel laws ([math §4](00-mathematics.md)) — treat as machine-precision tests:

- \(t^\star=\bar y-s\bar x R\)
- \(C=X_c^\top WY_c\)
- \(R_O^\star=UV^\top\), \(\gamma_O=\|C\|_*\)
- \(R_{SO}^\star=UDV^\top\), \(\gamma_{SO}=\sum_{j<d}\sigma_j+\det(UV^\top)\sigma_d\)
- \(s^\star=\operatorname{tr}((R^\star)^\top C)/a\)
- \(F^\star=b-\gamma^2/a\) with scale; \(F^\star=a+b-2\gamma\) without
- Sparse identity: \(X_c^\top WY_c=X^\top WY-(X^\top w)(Y^\top w)^\top/(\mathbf{1}^\top w)\)

Consensus for fixed transforms: \(M^\star=A^{-1}\sum_i\alpha_i Y_i\).

Energy identity ([math §27](00-mathematics.md)):

\[
\sum_i\alpha_i\|Y_i\|_F^2
=
A\|M\|_F^2
+
\sum_i\alpha_i\|Y_i-M\|_F^2.
\]

Label this an energy decomposition, not ANOVA, unless a sampling model supplies degrees of freedom.

---

## 10. Identifiability

Before fitting, construct the configuration-overlap graph. Nodes are configurations. An edge exists when two configurations share enough informative entities. Edge metadata: overlap count and centered overlap rank.

Report: connected components; configurations disconnected from the main problem; rank-deficient overlaps; entities observed nowhere; entities observed only once; whether relative transforms are estimable; whether components have unrelated gauges.

Default for a disconnected problem: **error with a diagnostic**. An explicit option may fit each component separately.

Every pairwise fit reports effective rank \(r_\tau=\#\{j:\sigma_j>\tau\sigma_1\}\), transform uniqueness, objective uniqueness, and unidentified subspace dimension.

---

## 11. Status fields

Two statuses, never collapsed:

**Numerical status:** `converged` | `maximum_iterations` | `numerical_failure` | `invalid_problem`

**Optimality status:** `exact_closed_form` | `certified_global` | `first_order_stationary` | `blockwise_stationary` | `converged_unverified` | `not_converged`

Allowed claims are listed in [01-solver-guarantee-matrix.md](01-solver-guarantee-matrix.md).
