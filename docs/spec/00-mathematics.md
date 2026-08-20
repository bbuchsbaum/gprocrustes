# Primary mathematical specification

The mathematics is the **primary specification** of `gprocrustes`. The R API should fall out of it.

There is not one “Procrustes problem.” There is a family of related optimization problems distinguished by transformation group, scale convention, correspondence structure, metric, missingness, robustness, and whether the reference itself is constrained.

Gower’s original formulation already has this broader spirit: multiple configurations move symmetrically toward a consensus rather than declaring one configuration fixed, with translation, rotation/reflection, and scale handled explicitly. Goodall embeds this in statistical shape analysis, including weighted metrics and covariance models. Modern work exposes both the nonconvex optimization structure of orthogonal GPA and closed-form eigenproblems for certain deformable formulations.

Notation here uses \(K\) configurations (the charter’s \(m\)) and target dimension \(d\) (the charter’s \(q\)).

---

## 1. A single convention

For the library, use

\[
X_i\in\mathbb R^{n_i\times p_i},
\qquad i=1,\ldots,K,
\]

where rows are corresponding **entities** and columns are coordinates, features, or dimensions.

That convention works for landmarks, cortical vertices, cells, genes, stimuli, embedding observations, factor scores, and similar objects.

For the usual GPA problem,

\[
p_i=d,\qquad X_i\in\mathbb R^{n\times d}.
\]

The consensus is

\[
M\in\mathbb R^{n\times d}.
\]

Bai and Bartoli use transposed shapes \(D_i\in\mathbb R^{d\times m}\). Their formulas are transposed here into the row-observation convention.

For partial correspondence, introduce a row-selection / incidence matrix

\[
P_i\in\{0,1\}^{n_i\times n},
\]

so that \(P_iM\) selects from the global consensus the entities observed in configuration \(i\).

This is preferable to treating `NA` as the mathematical definition of missingness.

Left-action callers (features \(\times\) observations, \(QX\)) must transpose into this convention at the adapter boundary.

---

## 2. The master optimization problem

At the highest level the library solves

\[
\min_{M,\theta_1,\ldots,\theta_K}
\sum_{i=1}^{K}
\alpha_i
L_i\!\left(
T_{\theta_i}(X_i),
P_iM
\right)
+
\sum_{i=1}^{K}\lambda_i\Omega_i(\theta_i)
\]

subject to transformation-specific and gauge constraints.

Here \(\alpha_i\ge 0\) is a configuration weight, \(T_{\theta_i}\) is a transformation, \(L_i\) is the fitting loss or metric, and \(\Omega_i\) is an optional regularizer.

Bai and Bartoli start from almost exactly the unweighted missing-landmark version

\[
\min_{\{T_i\},S}
\sum_{i,j}
\gamma_{ij}
\|T_i(D_{ij})-S_j\|_2^2,
\]

with constraints needed to avoid degeneracy.

For the basic isotropic squared-error problem,

\[
L_i(Y,Z)=\|Y-Z\|_F^2.
\]

For row weights \(W_i=\operatorname{diag}(w_{i1},\ldots,w_{in_i})\),

\[
L_i(Y,Z)
=
\|W_i^{1/2}(Y-Z)\|_F^2.
\]

For a general statistical metric,

\[
L_i(Y,Z)
=
\operatorname{vec}(Y-Z)^\top
Q_i
\operatorname{vec}(Y-Z),
\qquad Q_i\succeq0.
\]

This master formulation is valuable because the package can compile a requested model into the strongest solver available for that exact mathematical special case.

---

## 3. Transformation groups

The fundamental similarity transformation is

\[
T_i(X)
=
s_iXR_i+\mathbf 1t_i,
\]

where \(t_i\in\mathbb R^{1\times d}\).

For reflection-allowing alignment,

\[
R_i\in O(d)
=
\{R:R^\top R=I\}.
\]

For genuine rotations,

\[
R_i\in SO(d)
=
\{R:R^\top R=I,\det R=1\}.
\]

And \(s_i>0\) for similarity transforms.

This distinction must be mathematically explicit. Goodall’s definition of shape excludes reflection and uses positive isotropic scale, whereas reflection shape quotients by the larger orthogonal group. He stresses that the \(O(d)\) versus \(SO(d)\) distinction is important in shape theory.

The package regards these as different transformation groups, not a Boolean post-processing option.

---

## 4. Pairwise Procrustes: the irreducible numerical kernel

Everything starts here.

Suppose

\[
X,Y\in\mathbb R^{n\times d}
\]

with row weights \(w_j\ge0\),

\[
W=\operatorname{diag}(w_1,\ldots,w_n),
\qquad
w_+=\mathbf1^\top W\mathbf1.
\]

Define weighted centroids

\[
\bar x
=
\frac{\mathbf1^\top WX}{w_+},
\qquad
\bar y
=
\frac{\mathbf1^\top WY}{w_+},
\]

and weighted-centered matrices

\[
X_c=X-\mathbf1\bar x,
\qquad
Y_c=Y-\mathbf1\bar y.
\]

We solve

\[
\min_{s,R,t}
\left\|
W^{1/2}
(sXR+\mathbf1t-Y)
\right\|_F^2.
\]

### Translation

For any \(s,R\), differentiating in \(t\) gives

\[
t^\star=\bar y-s\bar xR.
\]

So translation can be eliminated analytically. This is the modern weighted form of Gower’s observation that, at the optimum, the configurations share a centroid, which may without loss of generality be placed at the origin.

After translation elimination:

\[
\min_{s,R}
\|W^{1/2}(sX_cR-Y_c)\|_F^2.
\]

Expand:

\[
F(s,R)
=
s^2a+b-2s\,\operatorname{tr}(R^\top C),
\]

where

\[
a=\operatorname{tr}(X_c^\top WX_c),
\qquad
b=\operatorname{tr}(Y_c^\top WY_c),
\]

and

\[
C=X_c^\top WY_c.
\]

### Orthogonal solution

Take

\[
C=U\Sigma V^\top,
\qquad
\Sigma=\operatorname{diag}(\sigma_1,\ldots,\sigma_d).
\]

Then

\[
\max_{R\in O(d)}
\operatorname{tr}(R^\top C)
\]

has the exact solution

\[
R^\star=UV^\top.
\]

The optimum trace is

\[
\gamma_O
=
\sum_j\sigma_j
=
\|C\|_*.
\]

This is the polar factor of \(C\).

### Proper rotation

For \(SO(d)\),

\[
D=
\operatorname{diag}
\left(
1,\ldots,1,
\det(UV^\top)
\right),
\]

and

\[
R^\star=UDV^\top.
\]

The corresponding optimum is

\[
\gamma_{SO}
=
\sum_{j=1}^{d-1}\sigma_j
+
\det(UV^\top)\sigma_d.
\]

Gower explicitly describes the determinant correction required when reflection is excluded.

### Scale

Once \(R^\star\) is known,

\[
s^\star
=
\frac{\operatorname{tr}\left((R^\star)^\top C\right)}
     {\operatorname{tr}(X_c^\top WX_c)}
\]

for unconstrained positive one-sided scale, with projection onto \(s\ge0\) if necessary.

Thus the minimized one-sided similarity residual is

\[
F^\star
=
b-\frac{\gamma^2}{a}.
\]

Without scaling,

\[
F^\star=a+b-2\gamma.
\]

These formulas are **kernel-level laws** and must be verified to essentially machine precision.

---

## 5. Rank deficiency and uniqueness

This deserves explicit mathematics because it is routinely glossed over.

For \(C=U\Sigma V^\top\), if \(\operatorname{rank}(C)=d\), the orthogonal polar factor \(UV^\top\) is unique.

If \(r=\operatorname{rank}(C)<d\), only the action on the identified \(r\)-dimensional subspace is determined. There is an arbitrary orthogonal mapping between the null spaces.

The library must distinguish

\[
\text{objective uniquely minimized}
\]

from

\[
\text{transformation uniquely identified}.
\]

For \(SO(d)\), determinant constraints can sometimes resolve one dimension that would remain ambiguous under \(O(d)\), so the exact uniqueness criterion must be calculated rather than inferred from a crude landmark count.

Numerically, the relevant quantities are

\[
\sigma_1(C)\ge\cdots\ge\sigma_d(C),
\]

the effective rank

\[
r_\tau
=
\#\{j:\sigma_j>\tau\sigma_1\},
\]

and the spectral gaps.

That rank diagnostic must appear in every pairwise fit.

---

## 6. Sparse centering without making a sparse matrix dense

Although \(X_c=X-\mathbf1\bar x\) is generally dense even when \(X\) is sparse, its cross-product need never be formed explicitly:

\[
X_c^\top WY_c
=
X^\top WY
-
\frac{(X^\top w)(Y^\top w)^\top}{\mathbf1^\top w}.
\]

Likewise

\[
X_c^\top WX_c
=
X^\top WX
-
\frac{(X^\top w)(X^\top w)^\top}{\mathbf1^\top w}.
\]

The computational primitive is **weighted moments**, not `scale(X)`.

---

## 7. Classical GPA as estimation of a consensus

Let \(Y_i=T_i(X_i)\). The symmetric GPA objective is

\[
F(M,\Theta)
=
\sum_{i=1}^{K}
\alpha_i
\|Y_i-M\|_F^2.
\]

For fixed transformations,

\[
M^\star
=
\frac{1}{A}
\sum_i\alpha_iY_i,
\qquad
A=\sum_i\alpha_i.
\]

That is the mathematical reason the reference is a consensus rather than a privileged observation.

Gower’s original criterion is exactly this idea geometrically: corresponding points across configurations form clusters and are fitted to their centroid.

There is also the identity

\[
\sum_i\alpha_i\|Y_i-M\|_F^2
=
\frac{1}{A}
\sum_{i<j}
\alpha_i\alpha_j
\|Y_i-Y_j\|_F^2
\]

when \(M\) is the weighted mean.

Thus “fit everything to the consensus” and “minimize all pairwise disagreement” are two views of the same quadratic objective.

---

## 8. Why GPA has a gauge

Suppose the problem is centered orthogonal GPA.

If \(\{R_i,M\}\) is a solution, then for every \(Q\in O(d)\),

\[
R_i' = R_iQ,
\qquad
M'=MQ
\]

gives exactly the same objective. Therefore

\[
(M,\{R_i\})
\sim
(MQ,\{R_iQ\}).
\]

This is not numerical nonidentifiability. It is a **true symmetry of the mathematical problem**.

Ling defines a distance modulo the common right orthogonal action:

\[
d_F(S,T)
=
\min_{Q\in O(d)}
\|S-TQ\|_F.
\]

With translation included there is additionally a common translation gauge. Before scale is fixed there is a common scale gauge.

Gower notes that the entire fitted system may be rotated without altering the criterion, and proposes expressing the final result in the principal axes of the consensus only after estimation.

Therefore:

\[
\text{optimization and canonicalization are different operations}.
\]

An optimizer must not secretly force a principal-axis orientation while it is solving the problem.

---

## 9. Scaling and the collapse degeneracy

If all \(s_i\) and \(M\) are unconstrained,

\[
s_i=0,\qquad M=0
\]

is a trivial perfect solution.

Gower solves this by constraining the total transformed energy:

\[
\sum_i
s_i^2
\operatorname{tr}(X_i^\top X_i)
=
\text{constant}.
\]

There are three legitimate scale semantics the package must distinguish.

### Preshape GPA

Normalize each centered configuration to \(\|X_i\|_F=1\) and subsequently optimize rotations only.

### Gower collective scaling

Estimate all \(s_i\) simultaneously subject to a global energy constraint.

### Fixed-consensus-scale similarity GPA

Require, for example, \(\|M\|_F^2=1\) or another prescribed centroid size and estimate each \(s_i\).

These are not identical models.

A clean derivation of Gower / Ten Berge scaling follows after rotations are fixed. Let

\[
Z_i=X_iR_i,
\qquad
G_{ij}=\langle Z_i,Z_j\rangle_F,
\]

and

\[
D=\operatorname{diag}
(\|Z_1\|_F^2,\ldots,\|Z_K\|_F^2).
\]

With scales \(s=(s_1,\ldots,s_K)^\top\),

\[
M=\frac1K\sum_i s_iZ_i.
\]

Eliminating \(M\),

\[
F(s)
=
s^\top Ds
-
\frac1K s^\top Gs.
\]

Under \(s^\top Ds=c\), minimizing \(F\) is equivalent to

\[
\max_s s^\top Gs
\quad
\text{s.t.}\quad
s^\top Ds=c.
\]

Hence \(Gs=\lambda Ds\).

When each configuration has been standardized so \(D=I\), the optimal scale vector is the leading eigenvector of the cross-configuration Gram / correlation matrix. That is the eigenvector scaling step described by Goodall’s account of Gower / Ten Berge GPA.

---

## 10. The classical alternating GPA solver

For similarity GPA, the basic iteration is

\[
M^{(t)}
\longrightarrow
\{R_i^{(t+1)},s_i^{(t+1)}\}
\longrightarrow
M^{(t+1)}.
\]

Every transformation update is an exact pairwise Procrustes solution against the current consensus. Every consensus update is an exact least-squares mean. Consequently \(F^{(t+1)}\le F^{(t)}\).

Gower emphasizes this monotonic decrease but also correctly observes that monotonic decrease and a lower bound do **not** by themselves prove convergence to a global optimum.

That distinction must survive unchanged:

\[
\text{converged}
\neq
\text{globally optimal}.
\]

---

## 11. Orthogonal GPA as group synchronization

Assume complete, centered data, \(X_i\in\mathbb R^{n\times d}\), \(R_i\in O(d)\), and solve

\[
\min_{\{R_i\},M}
\sum_i
\|X_iR_i-M\|_F^2.
\]

Eliminate \(M=\frac1K\sum_iX_iR_i\). The nonconstant part becomes

\[
\max_{\{R_i\}}
\left\|
\sum_iX_iR_i
\right\|_F^2.
\]

Expanding,

\[
\max_{\{R_i\}}
\sum_{i,j}
\operatorname{tr}
\left(
R_i^\top
X_i^\top X_j
R_j
\right).
\]

Define \(C_{ij}=X_i^\top X_j\) and the block matrix \(C\in\mathbb R^{Kd\times Kd}\). Stack

\[
S=
\begin{bmatrix}
R_1\\
\vdots\\
R_K
\end{bmatrix}
\in\mathbb R^{Kd\times d}.
\]

Then

\[
\max_{S\in O(d)^K}
\langle C,SS^\top\rangle
=
\max_S\operatorname{tr}(S^\top CS).
\]

This is Ling’s GOPP formulation. The problem is NP-hard in general, even at \(d=1\).

**There is no general SVD formula for multi-configuration orthogonal GPA.**

---

## 12. Spectral initialization

Relax the blockwise constraints temporarily and solve

\[
\max_{U^\top U=I_d}
\operatorname{tr}(U^\top CU).
\]

The solution is the leading \(d\)-dimensional eigenspace of \(C\). Partition \(U\) into blocks \(U_i\in\mathbb R^{d\times d}\), then project each block to its polar factor:

\[
R_i^{(0)}
=
\operatorname{polar}(U_i).
\]

That is Ling’s spectral initialization. Different eigenvector bases spanning the same leading subspace are gauge-equivalent.

---

## 13. Generalized power method

Given \(S^{(t)}\), compute \(B=CS^{(t)}\). Partition \(B\) into blocks \(B_i\in\mathbb R^{d\times d}\), then

\[
R_i^{(t+1)}
=
\operatorname{polar}(B_i),
\qquad
S^{(t+1)}
=
\mathcal P_{O(d)^K}
(CS^{(t)}).
\]

If \(B_i=U_i\Sigma_iV_i^\top\), then \(R_i^{(t+1)}=U_iV_i^\top\).

An iteration is essentially one matrix multiplication followed by \(K\) small polar projections.

Under Ling’s high-SNR signal-plus-noise assumptions, spectral initialization plus GPM converges linearly to the global solution. That is a theorem about that model, not a universal guarantee for arbitrary GPA inputs.

---

## 14. Never form the giant matrix \(C\)

Write \(B=[X_1\;X_2\;\cdots\;X_K]\in\mathbb R^{n\times Kd}\). Then \(C=B^\top B\), but \(CS=B^\top(BS)\). Even better,

\[
BS
=
\sum_j X_jR_j,
\]

so block \(i\) of \(CS\) is

\[
[CS]_i
=
X_i^\top
\left(
\sum_jX_jR_j
\right).
\]

A GPM step is therefore \(Z=\sum_jX_jR_j\), followed by \(B_i=X_i^\top Z\) and \(K\) small SVDs.

Memory becomes approximately \(O(nd+Kd^2)\) rather than \(O(K^2d^2)\).

---

## 15. First-order stationarity on \(O(d)^K\)

A solver must return more than “relative objective changed by \(10^{-8}\).”

For \(f(S)=\operatorname{tr}(S^\top CS)\), the Euclidean gradient is \(\nabla f=2CS\).

For one block \(R_i\), the tangent space is \(\{R_i\Omega:\Omega^\top=-\Omega\}\). Projection of an arbitrary matrix \(G_i\) onto the tangent space is

\[
\Pi_{R_i}(G_i)
=
G_i-
R_i\operatorname{sym}(R_i^\top G_i).
\]

The Riemannian gradient block is

\[
\operatorname{grad}_if
=
2
\left[
(CS)_i-
R_i\operatorname{sym}
\{R_i^\top(CS)_i\}
\right].
\]

A stationary solution satisfies

\[
\operatorname{skew}
\{R_i^\top(CS)_i\}
=0
\]

for every \(i\).

Report

\[
r_{\mathrm{stat}}
=
\max_i
\left\|
\operatorname{skew}
\{R_i^\top(CS)_i\}
\right\|_F,
\qquad
r_{\mathrm{orth}}
=
\max_i
\|R_i^\top R_i-I\|_F.
\]

Those are more meaningful than iteration count.

---

## 16. The semidefinite relaxation

Set \(Z=SS^\top\). Then \(Z\succeq0\) and its diagonal blocks obey \(Z_{ii}=I_d\). Dropping the nonconvex rank condition gives

\[
\max_Z \langle C,Z\rangle
\quad
\text{s.t.}
\quad
Z\succeq0,\quad Z_{ii}=I_d.
\]

A full SDP solver is not the default engine. The dual mathematics supplies the certificate.

---

## 17. A genuine global-optimality certificate

Let a candidate \(S\in O(d)^K\) be given. Construct a symmetric block-diagonal matrix \(\Lambda=\operatorname{blkdiag}(\Lambda_1,\ldots,\Lambda_K)\) such that \(CS=\Lambda S\).

At a GPM fixed point one natural choice is

\[
\Lambda_i
=
\left(
[CS]_i[CS]_i^\top
\right)^{1/2}.
\]

Then define the dual slack \(L=\Lambda-C\).

Ling’s global-optimality theorem: \(CS=\Lambda S\) and \(\Lambda-C\succeq0\) is sufficient for \(SS^\top\) to be globally optimal for the SDR, and hence for the original orthogonal problem when the relaxation is tight. If \(\operatorname{rank}(\Lambda-C)=(K-1)d\), the global solution is unique modulo the common orthogonal gauge.

Numerically calculate

\[
r_{\mathrm{dual}}
=
\frac{\|CS-\Lambda S\|_F}
     {1+\|CS\|_F}
\]

and the lowest eigenvalues of \(\Lambda-C\). Because the known gauge contributes \(d\) zero eigenvalues, a useful strict uniqueness diagnostic is \(\lambda_{d+1}(\Lambda-C)>0\).

This gives the package the ability to say either **certified global** or **converged, but global optimality not established**.

Caveat: that certificate belongs to the \(O(d)\) trace problem. Do not silently reuse it for \(SO(d)\), arbitrary robust losses, or cellwise missingness.

---

## 18. Missing rows: exploitable structure

Let \(P_i\in\{0,1\}^{n_i\times n}\) select the rows observed by configuration \(i\), and let \(W_i\succeq0\) be diagonal row weights.

Consider centered orthogonal GPA:

\[
F=
\sum_i
\alpha_i
\|
W_i^{1/2}
(X_iR_i-P_iM)
\|_F^2.
\]

For fixed \(R_i\), define \(D=\sum_i\alpha_iP_i^\top W_iP_i\). For ordinary row masks, \(D\) is a diagonal matrix of total support at every global entity.

The normal equations for \(M\) are \(DM=\sum_i\alpha_iP_i^\top W_iX_iR_i\), hence

\[
M^\star
=
D^\dagger
\sum_i
\alpha_iP_i^\top W_iX_iR_i.
\]

Elementwise this is

\[
M_j
=
\frac{
\sum_{i:j\in i}
\alpha_iw_{ij}(X_iR_i)_j
}{
\sum_{i:j\in i}\alpha_iw_{ij}
}.
\]

Eliminate \(M\). Define \(B(S)=\sum_i\alpha_iP_i^\top W_iX_iR_i\). Then, apart from terms constant under orthogonal \(R_i\),

\[
F
=
\text{const}
-
\operatorname{tr}
\left[
B(S)^\top
D^\dagger
B(S)
\right].
\]

Therefore even **row-missing orthogonal GPA remains a quadratic synchronization problem**, with block couplings

\[
C_{ij}
=
\alpha_i\alpha_j
X_i^\top W_i
P_i
D^\dagger
P_j^\top W_j
X_j.
\]

That is the natural variable-elimination extension of the supplied papers. Missing rows need not automatically condemn the solver to imputation. And \(C\) need not be formed.

---

## 19. Translation plus partial overlap

Translation can also be eliminated exactly. For configuration \(i\),

\[
\min_{t_i}
\|
W_i^{1/2}
(X_iR_i+\mathbf1t_i-P_iM)
\|_F^2.
\]

Define \(c_i=\mathbf1^\top W_i\mathbf1\) and the weighted residual-centering matrix

\[
K_i
=
W_i
-
W_i\mathbf1
c_i^{-1}
\mathbf1^\top W_i.
\]

After minimizing over \(t_i\),

\[
F_i
=
\operatorname{tr}
\left[
(X_iR_i-P_iM)^\top
K_i
(X_iR_i-P_iM)
\right].
\]

Then \(D=\sum_i P_i^\top K_iP_i\), \(B=\sum_iP_i^\top K_iX_iR_i\), and \(M^\star=D^\dagger B\).

Here \(D\) is generally sparse rather than diagonal and has at least the expected global translation nullspace.

Translation is not “preprocessing.” It is a nuisance parameter that can often be **projected out analytically**.

---

## 20. Arbitrary cellwise missingness is genuinely harder

Suppose observation masks differ by coordinate: \(\Omega_i\in\{0,1\}^{n_i\times d}\). Then

\[
F_i
=
\|
\Omega_i\odot
(X_iR_i-P_iM)
\|_F^2.
\]

Now \(\|\Omega_i\odot X_iR_i\|_F^2\) depends on \(R_i\): rotation mixes observed and unobserved coordinate components. The elegant orthogonal trace reduction disappears.

\[
\text{missing rows} \neq \text{missing arbitrary cells}.
\]

One principled solver is majorization. Let \(E_i(R,M)=X_iR-P_iM\). At iteration \(t\), construct

\[
Q_i(R,M\mid R^{(t)},M^{(t)})
=
\|\Omega_i\odot E_i(R,M)\|_F^2
+
\|(1-\Omega_i)\odot
(E_i(R,M)-E_i^{(t)})\|_F^2.
\]

Then \(Q_i\ge F_i\) and \(Q_i(R^{(t)},M^{(t)}\mid\cdot)=F_i^{(t)}\). Minimizing the surrogate yields monotone descent in the true observed-data criterion.

Alternatively, optimize the true objective over the product manifold.

Pretending arbitrary cell missingness retains the ordinary SVD solution is not defensible.

---

## 21. Weighted Procrustes in Goodall’s sense

Goodall’s major statistical generalization is \(r=\operatorname{vec}(E)\), \(E=T(X)-Y\), with

\[
L(E)=r^\top\Sigma_S^{-1}r.
\]

This is weighted / generalized least squares.

A particularly important covariance structure is \(\Sigma=\Sigma_N\otimes\Sigma_d\), where \(\Sigma_N\) describes dependence between landmarks / entities and \(\Sigma_d\) dependence between coordinate dimensions.

If \(\Sigma_d=\sigma^2I_d\), then with \(Q^\top Q=\Sigma_N^{-1}\), the weighted problem reduces to ordinary Procrustes on \(QX\), \(QY\). The pairwise SVD remains exact.

For genuinely anisotropic coordinate covariance, \(\Sigma_d\not\propto I\), rotation and the coordinate metric no longer commute. In general there is no ordinary one-SVD solution.

Solver selection:

- row metric \(\otimes I_d\) \(\rightarrow\) SVD fast path;
- arbitrary coordinate metric \(\rightarrow\) manifold / numerical solver.

---

## 22. Fitting metric and statistical model metric are different objects

Let \(\Sigma_S\) denote the **superimposition metric**, and \(\Sigma_M\) the **model covariance** used for statistical inference. They need not be equal.

Goodall explicitly separates these quantities and warns that inference based uncritically on the fitting metric can be biased.

The library must allow

\[
\hat M
=
\arg\min_M
L_{\Sigma_S}(M),
\]

followed by statistical analysis under \(\operatorname{vec}(E_i)\sim N(0,\Sigma_M)\).

When \(\Sigma_S=\Sigma_M\), the fitting procedure can have maximum-likelihood / efficiency interpretations under the corresponding Gaussian model.

Do not have one ambiguous argument named `weights`.

---

## 23. Robust GPA

For geometric robustness, the residual should usually be robustified at the **landmark-vector level**.

Let \(e_{ij}=(T_i(X_i))_j-M_j\) and \(r_{ij}=\|e_{ij}\|_2\). Use

\[
F
=
\sum_{ij}
\alpha_i
\rho(r_{ij}).
\]

This is preferable to \(\sum_{ijk}\rho(e_{ijk})\), because a coordinatewise loss changes when the coordinate system rotates.

For a differentiable \(\rho\), IRLS uses \(q(r)=\psi(r)/r\) with \(\psi(r)=\rho'(r)\).

Huber:

\[
q(r)
=
\begin{cases}
1,&r\le c,\\[4pt]
c/r,&r>c.
\end{cases}
\]

Tukey bisquare (nonconvex):

\[
q(r)
=
\begin{cases}
\left[1-(r/c)^2\right]^2,&r<c,\\
0,&r\ge c.
\end{cases}
\]

Iteration \(t+1\) solves a weighted GPA with \(w_{ij}^{(t+1)}=q(r_{ij}^{(t)})\).

Goodall discusses robustness weights based on the lengths of landmark residual vectors and also notes the possibility of an entire configuration being an outlier. The library must expose those final \(w_{ij}\).

---

## 24. Shape space and preshape geometry

Given \(X\in\mathbb R^{n\times d}\), center it: \(X_c=HX\) with \(H=I-\frac1n\mathbf1\mathbf1^\top\). Its centroid size is \(r=\|X_c\|_F\). The preshape is

\[
Z=\frac{X_c}{\|X_c\|_F},
\qquad
\|Z\|_F=1.
\]

Shape identifies all rotations: \([Z]=\{ZR:R\in SO(d)\}\).

The full Procrustes chord distance is

\[
d_P^2([Z_1],[Z_2])
=
\min_{R\in SO(d)}
\|Z_1-Z_2R\|_F^2.
\]

For normalized preshapes,

\[
d_P^2
=
2-2
\max_{R\in SO(d)}
\operatorname{tr}(R^\top Z_2^\top Z_1).
\]

This is the quotient-space interpretation that makes alignment different from comparing matrix entries.

---

## 25. Tangent-space mathematics

At a centered consensus \(M\), infinitesimal nuisance transformations span the vertical / orbit directions.

- Translation: \(\mathbf1a^\top\)
- Scale: \(cM\)
- Rotation: \(MA\), \(A^\top=-A\)

The similarity orbit tangent space is

\[
\mathcal V_M
=
\{
\mathbf1a^\top+cM+MA:
a\in\mathbb R^d,\;
c\in\mathbb R,\;
A^\top=-A
\}.
\]

Under the Frobenius metric, a perturbation \(Z\) is horizontal if \(\mathbf1^\top Z=0\), \(\langle Z,M\rangle_F=0\), and \(Z^\top M\) is symmetric.

The resulting shape dimension is

\[
nd
-
d
-
1
-
\frac{d(d-1)}2
=
nd-\frac{d(d+1)}2-1.
\]

That agrees with Goodall’s shape-dimension count. Shape space is not generally affine, so tangent-space procedures are local approximations.

---

## 26. A general weighted tangent projection

Construct a matrix \(B_M\) whose columns vectorize the nuisance directions: \(\mathbf1e_1^\top,\ldots,\mathbf1e_d^\top\), \(M\), and \(MA_{jk}\) for a basis of skew-symmetric matrices.

If the model / fitting precision is \(Q\), the weighted projection onto the nuisance tangent is

\[
P_V
=
B_M
(B_M^\top QB_M)^\dagger
B_M^\top Q.
\]

The horizontal projection is \(P_H=I-P_V\). Residual shape variation is \(z_H=P_H\operatorname{vec}(E)\).

That is one engine for tangent coordinates, residual covariance, shape PCA, approximate Gaussian inference, and degrees-of-freedom calculations.

---

## 27. Gower’s energy decomposition

Once configurations are aligned, \(M=\frac1A\sum_i\alpha_iY_i\), we have

\[
\sum_i\alpha_i\|Y_i\|_F^2
=
A\|M\|_F^2
+
\sum_i\alpha_i\|Y_i-M\|_F^2.
\]

Thus total aligned energy = consensus energy + residual energy.

Residual energy by configuration: \(R_i=\alpha_i\|Y_i-M\|_F^2\).

By entity: \(R_j=\sum_i\alpha_i\|(Y_i)_j-M_j\|_2^2\).

Gower used this for an analysis-of-variance-style summary but did not attach sampling degrees of freedom. Call this an **energy decomposition** unless a statistical model has justified inferential ANOVA.

---

## 28. Affine GPA

An affine transform is \(T_i(X_i)=X_iA_i+\mathbf1t_i\). Use homogeneous coordinates \(\widetilde X_i=[X_i\;\mathbf1]\) and \(B_i=\begin{bmatrix}A_i\\t_i\end{bmatrix}\), so \(T_i(X_i)=\widetilde X_iB_i\).

The naive problem \(\min_{B_i,M}\sum_i\|\widetilde X_iB_i-M\|_F^2\) is degenerate because everything can collapse toward zero.

Bai and Bartoli resolve this using constraints on the reference shape:

\[
\mathbf1^\top M=0,
\qquad
M^\top M=\Lambda,
\]

where \(\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_d)\) specifies the covariance eigenvalues / axis scales of the reference.

Given \(M\), \(B_i^\star=\widetilde X_i^\dagger M\). Define the projector \(H_i=\widetilde X_i\widetilde X_i^\dagger\). Then

\[
\min_{B_i}
\|\widetilde X_iB_i-M\|_F^2
=
\|(I-H_i)M\|_F^2,
\]

so \(F(M)=\operatorname{tr}(M^\top PM)\) with \(P=\sum_i(I-H_i)\). Affine GPA has become an eigenvalue problem.

---

## 29. The Brockett / Stiefel eigenproblem

Minimize \(\operatorname{tr}(M^\top PM)\) subject to \(M^\top M=\Lambda\). Write \(M=Q\Lambda^{1/2}\) with \(Q^\top Q=I\). Then

\[
F(Q)
=
\operatorname{tr}
\left(
\Lambda^{1/2}Q^\top PQ\Lambda^{1/2}
\right).
\]

This is a Brockett cost on the Stiefel manifold.

If \(\lambda_1\ge\cdots\ge\lambda_d\), the global minimum pairs the **largest reference variance** with the **smallest eigenvalue** of \(P\). If \(Pu_j=\eta_ju_j\) with \(\eta_1\le\eta_2\le\cdots\), then

\[
M^\star
=
[u_1,\ldots,u_d]\Lambda^{1/2},
\]

after excluding the translation eigenvector where appropriate.

That is a different sort of “global solution” from classical rigid GPA: the transform family is richer, but the constrained reference-space formulation permits variable projection and a global eigenvalue solution.

---

## 30. Linear-basis warps

Let \(\phi_i(x)\in\mathbb R^{\ell_i}\) be a basis expansion. For all rows of \(X_i\), collect these into \(\Phi_i\in\mathbb R^{n_i\times\ell_i}\). A linear-basis warp is

\[
T_i(X_i)=\Phi_iB_i,
\qquad
B_i\in\mathbb R^{\ell_i\times d}.
\]

This includes affine maps when \(\Phi_i=[X_i\;\mathbf1]\).

Introduce a quadratic regularizer \(\Omega_i(B_i)=\operatorname{tr}(B_i^\top L_iB_i)\), \(L_i\succeq0\). The GPA problem becomes

\[
\min_{\{B_i\},M}
\sum_i
\|\Phi_iB_i-M\|_F^2
+
\sum_i
\mu_i
\operatorname{tr}(B_i^\top L_iB_i)
\]

with reference constraints.

For fixed \(M\),

\[
B_i^\star
=
(\Phi_i^\top\Phi_i+\mu_iL_i)^\dagger
\Phi_i^\top M.
\]

Define the penalized smoother \(H_i=\Phi_i(\Phi_i^\top\Phi_i+\mu_iL_i)^\dagger\Phi_i^\top\). After eliminating \(B_i\),

\[
F(M)
=
\operatorname{tr}
\left[
M^\top
\left(
\sum_i(I-H_i)
\right)
M
\right].
\]

The hard problem collapses to a constrained eigenproblem.

---

## 31. Thin-plate splines inside LBW

Schematically,

\[
f(x)
=
b+Ax+
\sum_{k=1}^c
a_k\phi(\|x-c_k\|),
\]

so \(f(x)=B^\top\phi_{\mathrm{TPS}}(x)\) is linear in its coefficients even though it is nonlinear in \(x\).

The bending energy \(J(f)=\int\|\nabla^2f(x)\|_F^2\,dx\) becomes \(J(f)=\operatorname{tr}(B^\top LB)=\|ZB\|_F^2\).

The affine part lies in the nullspace of the bending penalty: \(J(f)=0\) iff \(f\) is affine for the usual TPS construction.

The package does not need a bespoke TPS architecture. It needs a good LBW abstraction plus a TPS basis / penalty constructor.

---

## 32. The free-translation condition

For a basis warp, translation should be representable without penalty. There must exist \(a_i\) satisfying

\[
\Phi_i a_i=\mathbf1,
\qquad
L_i a_i=0.
\]

Under this condition, \(P\mathbf1=0\). The translation direction is a null direction of the reduced eigenproblem and can be removed.

Bai and Bartoli call this **free-translations** and show that affine and TPS warps satisfy the property. This is a runtime mathematical check for user-supplied warp bases.

---

## 33. Partial shapes plus LBWs

Let \(P_i\) select the observed consensus entities and let \(W_i\) be visibility / weights. Solve

\[
\min_{B_i,M}
\sum_i
\|
W_i^{1/2}
(\Phi_iB_i-P_iM)
\|_F^2
+
\mu_i
\operatorname{tr}(B_i^\top L_iB_i).
\]

For fixed \(M\), \(B_i^\star=A_i^\dagger\Phi_i^\top W_iP_iM\) where \(A_i=\Phi_i^\top W_i\Phi_i+\mu_iL_i\). Substitution gives \(F(M)=\operatorname{tr}(M^\top PM)\) with

\[
P
=
\sum_i
\left[
P_i^\top W_iP_i
-
P_i^\top W_i
\Phi_i
A_i^\dagger
\Phi_i^\top W_i
P_i
\right].
\]

Then impose \(M^\top M=\Lambda\) and centering either as \(\mathbf1^\top M=0\) or a penalty \(\nu\|\mathbf1^\top M\|_2^2\), which changes the eigenoperator to \(P_\nu=P+\nu\mathbf1\mathbf1^\top\).

The solution is again the appropriate bottom eigenvectors.

---

## 34. Reference-space versus datum-space objectives

The reference-space loss is \(L_r=\sum_i\|T_i(X_i)-P_iM\|^2\).

A datum-space / generative loss is \(L_d=\sum_i\|X_i-T_i^{-1}(P_iM)\|^2\).

For Euclidean isometries these coincide. For affine or nonlinear transforms they generally do not. TPS may not even possess a global exact inverse.

`inverse(transform)` and `datum_space_error()` must have strong mathematical contracts. Do not invent a plausible inverse.

---

## 35. Deformation overfitting and cross-validation

A richer warp can always improve its training objective. \(\mu\downarrow\Rightarrow L_r\downarrow\) does **not** imply a better warp.

Bai and Bartoli use held-out landmarks: fit the warp without a group \(g\), predict those landmarks, gauge-align the cross-validation consensus to the full consensus, and evaluate held-out discrepancy.

If \(H_g\) is the held-out entity set,

\[
CV(\mu)
=
\sum_{i,g}
\left\|
\widehat T_{i,-g}
(X_{i,H_g})
-
\widehat M_{-g,H_g}
\right\|^2
\]

after gauge alignment. That is the appropriate criterion for selecting TPS smoothness or basis complexity.

---

## 36. Different feature dimensions

Suppose \(X_i\in\mathbb R^{n\times p_i}\) and we want a common \(q\)-dimensional representation. Use \(Q_i\in\operatorname{St}(p_i,q)\) and solve

\[
\min_{M,\{Q_i\}}
\sum_i
\|X_iQ_i-M\|_F^2.
\]

This is generalized semi-orthogonal Procrustes.

If the inputs are whitened, \(X_i^\top X_i=c_iI\), the self-energy terms are constant and the problem resembles the orthogonal trace-sum problem.

For arbitrary \(X_i\), \(\operatorname{tr}(Q_i^\top X_i^\top X_iQ_i)\) depends on \(Q_i\), so the problem lives on \(\prod_i\operatorname{St}(p_i,q)\) and should be solved with Riemannian optimization rather than forced into an ordinary GPA formula.

---

## 37. Signed-permutation Procrustes

For PCA / ICA / NMF-like components, arbitrary rotations may destroy interpretation. Restrict \(R=PD\) where \(P\) is a permutation matrix and \(D=\operatorname{diag}(\pm1,\ldots,\pm1)\).

Given cross-covariance \(C\), \(\max_{P,D}\operatorname{tr}(R^\top C)\) reduces to an assignment problem with reward \(|C_{jk}|\). After solving the permutation, \(D_{kk}=\operatorname{sign}(C_{P(k),k})\).

The global solution is obtained by a Hungarian / linear-assignment algorithm.

---

## 38. Overlap and identifiability are graph-theoretic

With partial configurations, form an overlap graph \(G=(V,E)\) where \(V=\{1,\ldots,K\}\) and \(i\sim j\) when the configurations have sufficient shared information to constrain their relative transform.

Connectivity is necessary for one globally linked gauge. If the graph has \(c>1\) components, there are at least \(c\) independent orientation gauges.

Connectivity alone is not sufficient. Each overlap has a cross-covariance \(C_{ij}\) whose rank determines which dimensions of the relative transformation are actually identified.

Diagnostics must involve both graph connectivity and \(\operatorname{rank}(C_{ij})\). This is more reliable than “three points are enough.”

Default behaviour for a disconnected problem is an error with a diagnostic.

---

## 39. Statistical uncertainty must respect the gauge

Suppose bootstrap replicate \(b\) produces \(M^{(b)}\). It is meaningless to average these matrices directly because every bootstrap solution may have a different arbitrary gauge.

First solve \(Q_b=\arg\min_{Q\in G}\|M^{(b)}Q-M^{(0)}\|_F^2\), then replace \(M^{(b)}\leftarrow M^{(b)}Q_b\). Only then calculate coordinatewise uncertainty.

If the reference has repeated or nearly repeated eigenvalues, individual canonical axes are unstable. Summarize an eigenspace through its projector \(P=UU^\top\) or principal angles \(\theta_j=\arccos\sigma_j(U^\top V)\).

---

## 40. Large-data sketching (optional)

For extremely large \(n\), a subspace embedding \(S\in\mathbb R^{r\times n}\) may satisfy \((1-\varepsilon)\|Zv\|_2^2\le\|SZv\|_2^2\le(1+\varepsilon)\|Zv\|_2^2\) for vectors in the relevant joint column space, so \(X_i^\top X_j\approx(SX_i)^\top(SX_j)\).

The right use is: sketch \(\rightarrow\) spectral / GPA initialization \(\rightarrow\) exact full-data refinement. An approximate sketch is not normally the final answer.

---

## 41. Numerical linear algebra primitives

The implementation reduces to a compact set of primitives:

| Primitive | Mathematical operation |
|---|---|
| weighted moments | \(X^\top WY,\ X^\top w\) |
| polar factor | \(UV^\top\) from small SVD |
| proper polar | \(UDV^\top\) |
| orthogonal projection | \(I-A(A^\top A)^\dagger A^\top\) |
| ridge smoother | \(A(A^\top A+\lambda L)^\dagger A^\top\) |
| linear solve | \(Ax=b\), never explicit \(A^{-1}\) |
| partial eigensystem | extreme eigenvectors of symmetric operators |
| matrix-free multiplication | \(v\mapsto Av\) |
| pseudoinverse | rank-revealing SVD / QR |
| tangent projection | \(G-R\operatorname{sym}(R^\top G)\) |
| gauge alignment | another small Procrustes problem |
| PSD certificate | smallest eigenvalues of \(\Lambda-C\) |

Equations may contain \(A^{-1}\) for exposition. The implementation should almost never calculate an inverse.

---

## 42. The mathematical guarantee lattice

This classification is part of the formal specification, not documentation trivia.

| Mathematical case | Strongest legitimate claim |
|---|---|
| Pairwise \(O(d)\), squared loss | exact global |
| Pairwise \(SO(d)\), squared loss | exact global |
| Pairwise similarity, standard weights | exact global |
| Signed permutation | exact global |
| Classical multiway similarity GPA | monotone iterative / stationary candidate |
| Complete \(O(d)\) GOPP + successful dual test | certified global |
| Complete \(O(d)\) GOPP without certificate | stationary / probable candidate only if converged |
| Arbitrary cell-mask GPA | stationary / MM solution |
| Robust GPA | stationary IRLS / MM; stronger claims depend on loss |
| Affine reference-space GPA with Bai–Bartoli constraints | global eigen-solution for stated formulation |
| LBW / TPS satisfying their assumptions | global eigen-solution for stated formulation |
| Generic Stiefel / semi-orthogonal GPA | stationary manifold solution |

Never upgrade “certificate unavailable” to “probably global.”

---

## Four engines

Almost everything reduces to four mathematical engines:

1. **SVD / polar Procrustes** — exact pairwise and block updates.
2. **Orthogonal-group synchronization** — rigid multiway GPA, spectral initialization, GPM, and certification.
3. **Quadratic variable projection \(\rightarrow\) Stiefel eigenproblem** — affine / LBW / TPS GPA.
4. **Weighted / MM / Riemannian optimization** — cases that break those closed forms: arbitrary covariance, robust losses, arbitrary missing cells, heterogeneous dimensions.

Around those sits one shared geometry of gauge, quotient spaces, rank, missing support, weighting, and statistical covariance.
