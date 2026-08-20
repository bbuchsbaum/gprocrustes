#!/usr/bin/env python3
"""Generate language-neutral gprocrustes conformance fixtures.

Pairwise oracles use NumPy SVD as an independent polar-factor implementation
of docs/spec/00-mathematics.md §4. Historical Gower digits are not invented.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIXDIR = ROOT / "inst" / "conformance" / "fixtures"
MATDIR = ROOT / "inst" / "conformance" / "matrices"
RANK_TOL = 1e-10
OBJ_ATOL = 1e-12


def mat(a: np.ndarray) -> list:
    out = np.asarray(a, dtype=float).tolist()
    return _json_ready(out)


def _json_ready(obj):
    if isinstance(obj, float):
        if not np.isfinite(obj):
            return None
        return obj
    if isinstance(obj, list):
        return [_json_ready(x) for x in obj]
    return obj


def write_fixture(doc: dict) -> None:
    path = FIXDIR / f"{doc['id']}.json"
    path.write_text(json.dumps(doc, indent=2) + "\n")


def write_binary(name: str, A: np.ndarray) -> dict:
    MATDIR.mkdir(parents=True, exist_ok=True)
    path = MATDIR / name
    # column-major float64 little-endian (Fortran / R)
    path.write_bytes(np.asfortranarray(A, dtype="<f8").tobytes(order="F"))
    return {
        "path": f"matrices/{name}",
        "nrow": int(A.shape[0]),
        "ncol": int(A.shape[1]),
        "dtype": "float64_le",
        "order": "column-major",
    }


def weighted_center(X: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    w = np.asarray(w, dtype=float)
    wp = float(w.sum())
    xbar = (w @ X) / wp
    return X - xbar, xbar


def moments(
    X: np.ndarray, Y: np.ndarray, w: np.ndarray
) -> tuple[np.ndarray, float, float]:
    Xc, _ = weighted_center(X, w)
    Yc, _ = weighted_center(Y, w)
    W = np.diag(w)
    C = Xc.T @ W @ Yc
    a = float(np.trace(Xc.T @ W @ Xc))
    b = float(np.trace(Yc.T @ W @ Yc))
    return C, a, b


def sparse_centered_cross(X: np.ndarray, Y: np.ndarray, w: np.ndarray) -> np.ndarray:
    # math §6: never form Xc
    XtWY = X.T @ (w[:, None] * Y)
    Xtw = X.T @ w
    Ytw = Y.T @ w
    return XtWY - np.outer(Xtw, Ytw) / w.sum()


def polar(C: np.ndarray, group: str) -> tuple[np.ndarray, np.ndarray, float]:
    U, s, Vt = np.linalg.svd(C, full_matrices=True)
    # Pad to square polar if C is square; for rectangular, thin polar
    d = C.shape[1]
    R = U[:, :d] @ Vt[:d, :]
    if group == "SO":
        if np.linalg.det(R) < 0:
            U = U.copy()
            U[:, d - 1] *= -1
            R = U[:, :d] @ Vt[:d, :]
    gamma = float(np.trace(R.T @ C))
    return R, s, gamma


def effective_rank(s: np.ndarray, tau: float = RANK_TOL) -> int:
    if s.size == 0 or s[0] <= 0:
        return 0
    return int(np.sum(s > tau * s[0]))


def pairwise_expected(X, Y, w, group: str, scaling: str, translation: bool) -> dict:
    C, a, b = moments(X, Y, w)
    R, svals, gamma = polar(C, group)
    r = effective_rank(svals)
    d = X.shape[1]
    if scaling == "isotropic":
        s_star = 0.0 if a == 0 else max(gamma / a, 0.0)
        obj = b - (0.0 if a == 0 else (gamma * gamma) / a)
    else:
        s_star = 1.0
        obj = a + b - 2.0 * gamma
    _, xbar = weighted_center(X, w)
    _, ybar = weighted_center(Y, w)
    t = (ybar - s_star * (xbar @ R)) if translation else np.zeros(d)
    unique = r == d
    # SO can identify one extra direction vs O when rank is d-1 and last
    # singular value is positive after determinant correction — report raw rank.
    return {
        "numerical_status": "exact",
        "optimality_status": "exact_closed_form",
        "objective": float(obj),
        "objective_atol": OBJ_ATOL,
        "objective_rtol": 0.0,
        "scale": float(s_star),
        "translation": mat(t) if t.ndim else [float(t)],
        "rotation": mat(R),
        "rotation_compare": "polar",
        "effective_rank": r,
        "transform_unique": bool(unique),
        "objective_unique": True if (group == "O" and r >= 1) or r == d else True,
        "unidentified_subspace_dimension": int(max(d - r, 0)),
        "gamma": float(gamma),
        "a": float(a),
        "b": float(b),
    }


def view(name, X, ids=None, row_mask=None):
    v = {"name": name, "matrix": mat(X)}
    if ids is not None:
        v["ids"] = list(ids)
    if row_mask is not None:
        v["row_mask"] = [bool(x) for x in row_mask]
    return v


def pairwise_doc(
    fid,
    title,
    X,
    Y,
    group,
    scaling="none",
    translation=False,
    w=None,
    tags=None,
    notes=None,
):
    w = np.ones(X.shape[0]) if w is None else np.asarray(w, dtype=float)
    exp = pairwise_expected(X, Y, w, group, scaling, translation)
    doc = {
        "id": fid,
        "title": title,
        "kind": "pairwise",
        "tags": tags or [],
        "problem": {
            "transform": {
                "group": "similarity" if scaling != "none" or translation else group,
                "translation": translation,
                "scaling": scaling,
            },
            "loss": "squared_l2",
            "source": mat(X),
            "target": mat(Y),
            "row_weights": [float(x) for x in w],
            "views": [view("source", X), view("target", Y)],
        },
        "expected": exp,
    }
    if group in ("O", "SO") and doc["problem"]["transform"]["group"] == "similarity":
        doc["problem"]["transform"]["rotation_group"] = group
    elif group in ("O", "SO"):
        doc["problem"]["transform"]["group"] = group
    if notes:
        doc["notes"] = notes
    return doc


def main() -> None:
    FIXDIR.mkdir(parents=True, exist_ok=True)
    MATDIR.mkdir(parents=True, exist_ok=True)
    for old in FIXDIR.glob("*.json"):
        old.unlink()

    rng = np.random.default_rng(20260820)
    n_written = 0

    def emit(doc):
        nonlocal n_written
        write_fixture(doc)
        n_written += 1

    # --- d = 1: O(1) vs SO(1) ---
    X = np.array([[1.0], [2.0], [3.0], [4.0]])
    Y_same = X.copy()
    Y_flip = -X
    emit(
        pairwise_doc(
            "pw-o1-same",
            "Pairwise O(1), identical 1-d configurations",
            X,
            Y_same,
            "O",
            tags=["d1", "O"],
        )
    )
    emit(
        pairwise_doc(
            "pw-so1-same",
            "Pairwise SO(1), identical 1-d configurations",
            X,
            Y_same,
            "SO",
            tags=["d1", "SO"],
        )
    )
    emit(
        pairwise_doc(
            "pw-o1-flip",
            "Pairwise O(1) recovers reflection on the line",
            X,
            Y_flip,
            "O",
            tags=["d1", "O", "reflection"],
        )
    )
    emit(
        pairwise_doc(
            "pw-so1-flip",
            "Pairwise SO(1) cannot reflect; residual is 4||X||^2 after centering",
            X,
            Y_flip,
            "SO",
            tags=["d1", "SO", "reflection"],
            notes="SO(1)={+1}. O(1) chooses -1. This is the sharp d=1 divergence.",
        )
    )

    # --- exact rotation in 2d / 3d ---
    def rot2(th):
        c, s = np.cos(th), np.sin(th)
        return np.array([[c, -s], [s, c]])

    X2 = rng.normal(size=(8, 2))
    Rtrue = rot2(0.7)
    Y2 = X2 @ Rtrue
    emit(
        pairwise_doc(
            "pw-o2-exact-rotation",
            "Pairwise O(2) recovers a known rotation",
            X2,
            Y2,
            "O",
            tags=["O", "exact"],
        )
    )
    emit(
        pairwise_doc(
            "pw-so2-exact-rotation",
            "Pairwise SO(2) recovers a known rotation",
            X2,
            Y2,
            "SO",
            tags=["SO", "exact"],
        )
    )

    Strue = np.array([[1.0, 0.0], [0.0, -1.0]])  # reflection
    Y2r = X2 @ Strue
    emit(
        pairwise_doc(
            "pw-o2-reflection",
            "Pairwise O(2) allows a reflection",
            X2,
            Y2r,
            "O",
            tags=["O", "reflection"],
        )
    )
    emit(
        pairwise_doc(
            "pw-so2-reflection",
            "Pairwise SO(2) rejects a reflection",
            X2,
            Y2r,
            "SO",
            tags=["SO", "reflection"],
        )
    )

    X3 = rng.normal(size=(10, 3))
    # Householder reflection in 3d
    u = np.array([1.0, 0.2, -0.4])
    u = u / np.linalg.norm(u)
    Href = np.eye(3) - 2 * np.outer(u, u)
    Y3r = X3 @ Href
    emit(
        pairwise_doc(
            "pw-o3-reflection",
            "Pairwise O(3) allows a Householder reflection",
            X3,
            Y3r,
            "O",
            tags=["O", "reflection"],
        )
    )
    emit(
        pairwise_doc(
            "pw-so3-reflection",
            "Pairwise SO(3) applies determinant correction",
            X3,
            Y3r,
            "SO",
            tags=["SO", "reflection"],
        )
    )

    # proper 3d rotation
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    Y3 = X3 @ Q
    emit(
        pairwise_doc(
            "pw-so3-exact-rotation",
            "Pairwise SO(3) recovers a proper rotation",
            X3,
            Y3,
            "SO",
            tags=["SO", "exact"],
        )
    )

    # --- translation and scale ---
    ttrue = np.array([3.0, -2.0])
    Y2t = X2 @ Rtrue + ttrue
    emit(
        pairwise_doc(
            "pw-sim-translation",
            "Similarity with translation only",
            X2,
            Y2t,
            "SO",
            scaling="none",
            translation=True,
            tags=["similarity", "translation"],
        )
    )
    Y2s = 4.5 * (X2 @ Rtrue) + ttrue
    emit(
        pairwise_doc(
            "pw-sim-scale-translation",
            "Similarity with isotropic scale and translation",
            X2,
            Y2s,
            "SO",
            scaling="isotropic",
            translation=True,
            tags=["similarity", "scale", "translation"],
        )
    )
    # huge / tiny scales
    emit(
        pairwise_doc(
            "pw-sim-huge-scale",
            "Isotropic scale 1e6",
            X2,
            1e6 * (X2 @ Rtrue),
            "O",
            scaling="isotropic",
            translation=False,
            tags=["similarity", "scale", "extreme"],
        )
    )
    emit(
        pairwise_doc(
            "pw-sim-tiny-scale",
            "Isotropic scale 1e-6",
            X2,
            1e-6 * (X2 @ Rtrue),
            "O",
            scaling="isotropic",
            translation=False,
            tags=["similarity", "scale", "extreme"],
        )
    )

    # no-scale residual formula vs scale formula
    Xn = rng.normal(size=(6, 2))
    Yn = rng.normal(size=(6, 2))
    emit(
        pairwise_doc(
            "pw-o2-no-scale-formula",
            "O(2) residual a+b-2gamma",
            Xn,
            Yn,
            "O",
            tags=["O", "formula"],
        )
    )
    emit(
        pairwise_doc(
            "pw-sim-scale-formula",
            "Similarity residual b-gamma^2/a",
            Xn,
            Yn,
            "O",
            scaling="isotropic",
            tags=["similarity", "formula"],
        )
    )

    # --- weighted ---
    w = np.array([1.0, 1.0, 1.0, 0.25, 0.25, 4.0, 4.0, 1.0])
    emit(
        pairwise_doc(
            "pw-weighted-rows",
            "Row-weighted O(2) polar factor",
            X2,
            Y2,
            "O",
            w=w,
            tags=["weights"],
        )
    )
    w0 = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    emit(
        pairwise_doc(
            "pw-zero-row-weight",
            "Zero row weights drop those entities",
            X2,
            Y2,
            "O",
            w=w0,
            tags=["weights", "zero-weight"],
        )
    )

    # --- rank deficiency ---
    Xrk = np.array([[1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]])
    Yrk = np.array([[0.0, 1.0], [0.0, 2.0], [0.0, 3.0], [0.0, 4.0]])
    emit(
        pairwise_doc(
            "pw-rank1-crosscov",
            "Rank-1 cross-covariance in d=2: transform not unique",
            Xrk,
            Yrk,
            "O",
            tags=["rank", "nonunique"],
        )
    )
    Xz = np.zeros((5, 3))
    Yz = np.zeros((5, 3))
    emit(
        pairwise_doc(
            "pw-rank0-zeros",
            "Rank-zero (all zeros) configurations",
            Xz,
            Yz,
            "O",
            tags=["rank", "degenerate"],
        )
    )
    # repeated singular values: two equal columns after rotation
    Xrep = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.5, 0.5, 0.0],
        ]
    )
    Yrep = Xrep @ np.diag([1.0, 1.0, 0.0])
    emit(
        pairwise_doc(
            "pw-repeated-singular",
            "Nearly repeated / zero trailing singular value",
            Xrep,
            Yrep,
            "O",
            tags=["rank", "repeated-sigma"],
        )
    )

    # one landmark
    X1 = np.array([[1.2, -0.4]])
    Y1 = np.array([[0.3, 1.1]])
    emit(
        pairwise_doc(
            "pw-one-landmark",
            "Single observed landmark in d=2",
            X1,
            Y1,
            "O",
            translation=True,
            tags=["rank", "one-landmark"],
        )
    )

    # --- signed permutation ---
    Cperm = np.array([[0.1, 0.9, 0.0], [0.8, 0.05, 0.1], [0.0, 0.2, -0.7]])
    # assignment on |C|: col0->row1 (0.8), col1->row0 (0.9), col2->row2 (0.7)
    P = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    D = np.diag([1.0, 1.0, -1.0])
    Rsp = P @ D
    emit(
        {
            "id": "pw-signed-perm-3",
            "title": "Signed-permutation assignment on a 3x3 cross-covariance",
            "kind": "pairwise",
            "tags": ["signed_permutation"],
            "problem": {
                "transform": {
                    "group": "signed_permutation",
                    "translation": False,
                    "scaling": "none",
                },
                "loss": "squared_l2",
                "source": mat(np.eye(3)),
                "target": mat(Cperm.T),
                "notes_crosscovariance": mat(Cperm),
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "exact_closed_form",
                "rotation": mat(Rsp),
                "rotation_compare": "exact",
                "objective_unique": True,
                "transform_unique": True,
            },
            "notes": "Reward |C_jk|; D_kk = sign(C_{P(k),k}).",
        }
    )

    # identity / composition laws as pairwise exact identities
    emit(
        pairwise_doc(
            "pw-identity-map",
            "Identity configurations: R = I, objective 0",
            X2,
            X2,
            "SO",
            tags=["law", "identity"],
        )
    )

    # --- moments: dense vs sparse identity ---
    Xs = np.zeros((12, 4))
    Xs[[0, 2, 5, 7, 11], 0] = [1.0, -2.0, 0.5, 3.0, -1.0]
    Xs[[1, 4, 8], 2] = [4.0, -0.5, 2.5]
    Ys = np.zeros((12, 4))
    Ys[[0, 3, 5, 9], 1] = [2.0, 1.0, -3.0, 0.5]
    Ys[[2, 6, 10], 3] = [1.5, -1.0, 2.0]
    ws = np.ones(12)
    C_explicit, a_s, b_s = moments(Xs, Ys, ws)
    C_sparse = sparse_centered_cross(Xs, Ys, ws)
    emit(
        {
            "id": "moments-dense-sparse-agree",
            "title": "Centered cross-product: explicit Xc vs sufficient statistics",
            "kind": "moments",
            "tags": ["sparse", "moments"],
            "problem": {
                "transform": {"group": "none"},
                "source": mat(Xs),
                "target": mat(Ys),
                "row_weights": [float(x) for x in ws],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "crossproduct": mat(C_sparse),
                "a": float(a_s),
                "b": float(b_s),
                "objective_atol": 1e-14,
                "laws": [
                    "Xc.T W Yc == X.T W Y - (X.T w)(Y.T w).T / 1.T w",
                    "dense and sparse formulas agree",
                ],
            },
        }
    )
    assert np.allclose(C_explicit, C_sparse)

    # weighted sparse identity
    ww = np.linspace(0.5, 2.0, 12)
    C_e, a_w, b_w = moments(Xs, Ys, ww)
    C_s = sparse_centered_cross(Xs, Ys, ww)
    emit(
        {
            "id": "moments-weighted-sparse-agree",
            "title": "Weighted centered cross-product without forming Xc",
            "kind": "moments",
            "tags": ["sparse", "moments", "weights"],
            "problem": {
                "transform": {"group": "none"},
                "source": mat(Xs),
                "target": mat(Ys),
                "row_weights": [float(x) for x in ww],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "crossproduct": mat(C_s),
                "a": float(a_w),
                "b": float(b_w),
                "objective_atol": 1e-14,
                "laws": ["weighted sufficient-statistic identity"],
            },
        }
    )
    assert np.allclose(C_e, C_s)

    # binary-backed pair
    Xb = rng.normal(size=(20, 5))
    Qb, _ = np.linalg.qr(rng.normal(size=(5, 5)))
    Yb = Xb @ Qb
    emit(
        {
            "id": "pw-binary-o5",
            "title": "Pairwise O(5) with column-major float64 binaries",
            "kind": "pairwise",
            "tags": ["binary", "O"],
            "problem": {
                "transform": {"group": "O", "translation": False, "scaling": "none"},
                "loss": "squared_l2",
                "source": write_binary("pw-binary-o5-X.bin", Xb),
                "target": write_binary("pw-binary-o5-Y.bin", Yb),
                "row_weights": [1.0] * 20,
            },
            "expected": pairwise_expected(Xb, Yb, np.ones(20), "O", "none", False),
        }
    )

    # --- gauge invariance ---
    Qg = rot2(1.1)
    emit(
        {
            "id": "law-gauge-objective",
            "title": "Right-multiplying source and target by the same Q leaves the residual unchanged",
            "kind": "law",
            "tags": ["gauge"],
            "problem": {
                "transform": {"group": "O", "translation": False, "scaling": "none"},
                "loss": "squared_l2",
                "views": [
                    view("X", X2),
                    view("Y", Y2),
                    view("XQ", X2 @ Qg),
                    view("YQ", Y2 @ Qg),
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "F(X,Y)=F(XQ,YQ) for Q in O(d)",
                    "optimization and canonicalization are different operations",
                ],
                "objective": pairwise_expected(X2, Y2, np.ones(8), "O", "none", False)[
                    "objective"
                ],
                "objective_atol": OBJ_ATOL,
            },
        }
    )

    # configuration order invariance of pairwise residual
    emit(
        {
            "id": "law-entity-order",
            "title": "Permuting rows of X and Y together does not change the pairwise residual",
            "kind": "law",
            "tags": ["invariance"],
            "problem": {
                "transform": {"group": "O"},
                "loss": "squared_l2",
                "source": mat(X2),
                "target": mat(Y2),
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "exact_closed_form",
                "objective": pairwise_expected(X2, Y2, np.ones(8), "O", "none", False)[
                    "objective"
                ],
                "objective_atol": OBJ_ATOL,
                "laws": [
                    "row permutation of (X,Y) together is a gauge of the pairwise problem"
                ],
            },
        }
    )

    # energy identity for three aligned copies
    A1 = rng.normal(size=(5, 2))
    A2 = A1 + 0.1 * rng.normal(size=(5, 2))
    A3 = A1 + 0.2 * rng.normal(size=(5, 2))
    M = (A1 + A2 + A3) / 3.0
    total = float(
        np.linalg.norm(A1) ** 2 + np.linalg.norm(A2) ** 2 + np.linalg.norm(A3) ** 2
    )
    cons = 3.0 * float(np.linalg.norm(M) ** 2)
    resid = float(
        np.linalg.norm(A1 - M) ** 2
        + np.linalg.norm(A2 - M) ** 2
        + np.linalg.norm(A3 - M) ** 2
    )
    emit(
        {
            "id": "law-gower-energy",
            "title": "Gower energy identity: total = consensus + residual",
            "kind": "law",
            "tags": ["gower", "decomposition"],
            "problem": {
                "transform": {"group": "none"},
                "views": [view("A1", A1), view("A2", A2), view("A3", A3)],
                "configuration_weights": [1.0, 1.0, 1.0],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "a": total,
                "b": cons,
                "objective": resid,
                "objective_atol": 1e-12,
                "laws": ["sum_i ||Yi||^2 = K||M||^2 + sum_i ||Yi-M||^2"],
            },
        }
    )

    # pairwise disagreement identity
    pair = 0.0
    for U, V in ((A1, A2), (A1, A3), (A2, A3)):
        pair += float(np.linalg.norm(U - V) ** 2)
    emit(
        {
            "id": "law-consensus-vs-pairwise",
            "title": "Fit-to-consensus residual equals pairwise disagreement / K",
            "kind": "law",
            "tags": ["gower"],
            "problem": {
                "transform": {"group": "none"},
                "views": [view("A1", A1), view("A2", A2), view("A3", A3)],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "objective": resid,
                "gamma": pair / 3.0,
                "objective_atol": 1e-12,
                "laws": ["sum_i ||Yi-M||^2 = (1/K) sum_{i<j} ||Yi-Yj||^2"],
            },
        }
    )

    # --- Gower historical ---
    # Figure-1 geometry: m=3, n=4, p=2. Synthetic, not claimed as Table 1 digits.
    G1 = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    G2 = G1 @ rot2(0.3) + np.array([0.1, -0.2])
    G3 = G1 @ rot2(-0.4) * 1.2 + np.array([-0.05, 0.15])
    emit(
        {
            "id": "gower-identity-pairwise-clusters",
            "title": "Gower identity (2) on a 3-configuration 4-point example",
            "kind": "gower_history",
            "citation": "Gower, J. C. (1975). Psychometrika 40:33-51, equation (2).",
            "tags": ["gower", "historical"],
            "problem": {
                "transform": {
                    "group": "similarity",
                    "translation": True,
                    "scaling": "gower",
                },
                "views": [
                    view("X1", G1, ids=["a", "b", "c", "d"]),
                    view("X2", G2, ids=["a", "b", "c", "d"]),
                    view("X3", G3, ids=["a", "b", "c", "d"]),
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "sum_{i<j} ||Pi-Pj||^2 = K * sum_i ||Pi-G||^2 at the centroid G",
                    "after optimal translations the configurations share a centroid",
                ],
            },
            "notes": "Geometry matches Gower Figure 1 (m=3,n=4,p=2). Coordinates are not claimed as Table 1.",
        }
    )
    emit(
        {
            "id": "gower-1975-algorithm-laws",
            "title": "Gower 1975 algorithm laws (no invented table digits)",
            "kind": "gower_history",
            "citation": "Gower, J. C. (1975). Psychometrika 40:33-51.",
            "tags": ["gower", "historical"],
            "problem": {
                "transform": {
                    "group": "similarity",
                    "translation": True,
                    "scaling": "gower",
                },
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "each rotation/scale update is an exact pairwise Procrustes step against the current consensus",
                    "each consensus update is the arithmetic (or weighted) mean",
                    "the residual is monotone nonincreasing",
                    "a global energy constraint prevents the collapse s_i=0, M=0",
                    "monotone decrease does not prove a global optimum",
                    "the energy table is not ANOVA unless a sampling model is declared",
                ],
            },
        }
    )
    # Table 2: 9 carcasses (rows) by 7 characters (columns), three judges.
    # Table 5: successive Sr. Digits transcribed from Psychometrika 40:33-51;
    # Table 5 verified against a render of page 49. Gower's rotation/scale
    # schedule is not our BCD history; the published path is an oracle for
    # implementations that follow his steps, not a claim that gpa() matches it.
    judge1 = np.array(
        [
            [47, 44, 49, 38, 35, 40, 40],
            [72, 45, 41, 77, 72, 73, 35],
            [61, 49, 40, 58, 58, 62, 30],
            [66, 56, 45, 55, 53, 46, 30],
            [37, 72, 50, 27, 30, 33, 25],
            [76, 76, 53, 81, 79, 75, 45],
            [64, 59, 51, 72, 61, 66, 40],
            [21, 70, 43, 27, 22, 26, 20],
            [71, 70, 34, 72, 72, 71, 35],
        ],
        dtype=float,
    )
    judge2 = np.array(
        [
            [31, 39, 33, 29, 48, 38, 42],
            [30, 60, 36, 22, 36, 34, 39],
            [27, 55, 30, 18, 28, 22, 42],
            [48, 52, 53, 27, 21, 30, 31],
            [20, 55, 28, 22, 33, 27, 35],
            [21, 42, 31, 46, 76, 33, 42],
            [30, 52, 53, 35, 44, 30, 44],
            [5, 57, 53, 12, 13, 6, 31],
            [55, 63, 53, 77, 79, 57, 49],
        ],
        dtype=float,
    )
    judge3 = np.array(
        [
            [43, 46, 44, 22, 53, 44, 29],
            [53, 79, 75, 79, 73, 52, 27],
            [22, 85, 83, 19, 27, 17, 22],
            [28, 89, 78, 13, 29, 20, 24],
            [75, 86, 85, 34, 75, 55, 38],
            [53, 79, 82, 72, 78, 74, 38],
            [15, 85, 85, 46, 75, 52, 35],
            [5, 95, 95, 3, 20, 2, 24],
            [27, 78, 85, 89, 92, 81, 41],
        ],
        dtype=float,
    )
    emit(
        {
            "id": "gower-1975-published-history",
            "title": "Gower 1975 Table 2 carcass scores and Table 5 successive Sr",
            "kind": "gower_history",
            "citation": "Gower, J. C. (1975). Psychometrika 40:33-51, Tables 2 and 5.",
            "tags": ["gower", "historical", "transcribed"],
            "problem": {
                "transform": {
                    "group": "similarity",
                    "translation": True,
                    "scaling": "gower",
                },
                "views": [
                    view("judge_1", judge1, ids=list(range(1, 10))),
                    view("judge_2", judge2, ids=list(range(1, 10))),
                    view("judge_3", judge3, ids=list(range(1, 10))),
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "successive residual sums of squares must match Table 5 within published rounding when the 1975 rotation-then-scale schedule is followed",
                    "Gower stops when successive Sr differ by less than 0.0001 after the scaling step",
                    "Sr increased slightly after rotation steps 6 and 7 in the published path; Gower attributes that to 1975 numerical inaccuracy",
                ],
                "published_sr": {
                    "criterion": "Sr",
                    "initial": 0.661438,
                    "with_scaling": [
                        {
                            "iteration": 1,
                            "after_rotation": 0.657312,
                            "after_scaling": 0.616714,
                        },
                        {
                            "iteration": 2,
                            "after_rotation": 0.616620,
                            "after_scaling": 0.604215,
                        },
                        {
                            "iteration": 3,
                            "after_rotation": 0.604201,
                            "after_scaling": 0.600456,
                        },
                        {
                            "iteration": 4,
                            "after_rotation": 0.600452,
                            "after_scaling": 0.599322,
                        },
                        {
                            "iteration": 5,
                            "after_rotation": 0.599322,
                            "after_scaling": 0.598978,
                        },
                        {
                            "iteration": 6,
                            "after_rotation": 0.598980,
                            "after_scaling": 0.598875,
                        },
                        {
                            "iteration": 7,
                            "after_rotation": 0.598877,
                            "after_scaling": 0.598842,
                        },
                    ],
                    "without_scaling": [
                        {"iteration": 0, "after_rotation": 0.661438},
                        {"iteration": 1, "after_rotation": 0.657312},
                        {"iteration": 2, "after_rotation": 0.657158},
                        {"iteration": 3, "after_rotation": 0.657137},
                    ],
                },
            },
            "notes": (
                "Transcribed from Table 2 (9 carcasses x 7 characters, three judges) "
                "and Table 5 (successive Sr). Table 1 in the paper is the ANOVA layout, "
                "not the coordinates. gpa() uses consensus-first BCD and need not "
                "reproduce Gower's 1975 iteration path."
            ),
        }
    )

    # --- identifiability ---
    emit(
        {
            "id": "id-disconnected-overlap",
            "title": "Two configurations with no shared entities: disconnected overlap graph",
            "kind": "identifiability",
            "tags": ["overlap", "error-default"],
            "problem": {
                "transform": {"group": "O"},
                "views": [
                    view("A", rng.normal(size=(4, 2)), ids=["a", "b", "c", "d"]),
                    view("B", rng.normal(size=(4, 2)), ids=["e", "f", "g", "h"]),
                ],
            },
            "expected": {
                "numerical_status": "invalid_problem",
                "optimality_status": "not_applicable",
                "error_code": "disconnected_overlap_graph",
                "connected_components": 2,
                "laws": [
                    "default is an error; an explicit option may fit components separately"
                ],
            },
        }
    )
    emit(
        {
            "id": "id-connected-partial",
            "title": "Partial overlap that still connects three views",
            "kind": "identifiability",
            "tags": ["overlap"],
            "problem": {
                "transform": {"group": "O"},
                "views": [
                    view("A", rng.normal(size=(3, 2)), ids=["a", "b", "c"]),
                    view("B", rng.normal(size=(3, 2)), ids=["c", "d", "e"]),
                    view("C", rng.normal(size=(3, 2)), ids=["e", "f", "a"]),
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "connected_components": 1,
            },
        }
    )
    emit(
        {
            "id": "id-entity-observed-nowhere",
            "title": "Global id present in the map but observed by no view",
            "kind": "identifiability",
            "tags": ["overlap"],
            "problem": {
                "transform": {"group": "O"},
                "views": [
                    view("A", rng.normal(size=(2, 2)), ids=["a", "b"]),
                    view("B", rng.normal(size=(2, 2)), ids=["a", "b"]),
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "connected_components": 1,
                "laws": [
                    "entities observed nowhere must be reported; they do not enter M"
                ],
            },
            "notes": "A declared global entity 'z' with no observations is a diagnostic, not a row of M.",
        }
    )

    # row mask consensus mean
    obs = np.array(
        [
            [[1.0, 0.0], [2.0, 0.0], [np.nan, np.nan]],
            [[1.2, 0.1], [np.nan, np.nan], [3.0, 0.0]],
        ]
    )
    emit(
        {
            "id": "rowmask-consensus-mean",
            "title": "Row-masked consensus is the observed-data weighted mean",
            "kind": "law",
            "tags": ["missing", "row-mask"],
            "problem": {
                "transform": {"group": "O", "translation": False, "scaling": "none"},
                "views": [
                    view(
                        "A",
                        np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 0.0]]),
                        ids=["u", "v", "w"],
                        row_mask=[True, True, False],
                    ),
                    view(
                        "B",
                        np.array([[1.2, 0.1], [0.0, 0.0], [3.0, 0.0]]),
                        ids=["u", "v", "w"],
                        row_mask=[True, False, True],
                    ),
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "M_u = mean of observed transformed rows for u",
                    "M_v uses only view A",
                    "M_w uses only view B",
                    "implicit zeros are not missing",
                ],
            },
        }
    )

    # cell mask is harder
    emit(
        {
            "id": "cellmask-not-svd",
            "title": "Coordinate-wise masks destroy the ordinary SVD update",
            "kind": "law",
            "tags": ["missing", "cell-mask"],
            "problem": {
                "transform": {"group": "O"},
                "views": [
                    {
                        "name": "A",
                        "matrix": mat(np.array([[1.0, 2.0], [3.0, 4.0]])),
                        "cell_mask": mat(np.array([[1.0, 0.0], [1.0, 1.0]])),
                    }
                ],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "missing rows != missing arbitrary cells",
                    "ordinary polar factor is not an exact block update",
                    "allowed claim is first-order stationary or MM",
                ],
            },
        }
    )

    # --- errors ---
    emit(
        {
            "id": "err-nan",
            "title": "Nonfinite values are an error",
            "kind": "error",
            "tags": ["error", "nonfinite"],
            "problem": {
                "transform": {"group": "O"},
                "source": mat(np.array([[1.0, 0.0], [np.nan, 1.0]])),
                "target": mat(np.array([[0.0, 1.0], [1.0, 0.0]])),
            },
            "expected": {
                "numerical_status": "error",
                "optimality_status": "not_applicable",
                "error_code": "nonfinite_values",
            },
        }
    )
    emit(
        {
            "id": "err-inf",
            "title": "Infinite values are an error",
            "kind": "error",
            "tags": ["error", "nonfinite"],
            "problem": {
                "transform": {"group": "O"},
                "source": mat(np.array([[1.0, 0.0], [np.inf, 1.0]])),
                "target": mat(np.array([[0.0, 1.0], [1.0, 0.0]])),
            },
            "expected": {
                "numerical_status": "error",
                "optimality_status": "not_applicable",
                "error_code": "nonfinite_values",
            },
        }
    )
    emit(
        {
            "id": "err-duplicate-ids",
            "title": "Duplicate entity identifiers without an aggregation rule",
            "kind": "error",
            "tags": ["error", "correspondence"],
            "problem": {
                "transform": {"group": "O"},
                "views": [view("A", np.eye(2), ids=["a", "a"])],
            },
            "expected": {
                "numerical_status": "error",
                "optimality_status": "not_applicable",
                "error_code": "duplicate_entity_ids",
            },
        }
    )
    emit(
        {
            "id": "err-zero-config-weight-all",
            "title": "All configuration weights zero is invalid",
            "kind": "error",
            "tags": ["error", "weights"],
            "problem": {
                "transform": {"group": "O"},
                "views": [view("A", X2), view("B", Y2)],
                "configuration_weights": [0.0, 0.0],
            },
            "expected": {
                "numerical_status": "invalid_problem",
                "optimality_status": "not_applicable",
                "error_code": "zero_total_configuration_weight",
            },
        }
    )
    emit(
        {
            "id": "err-dimension-mismatch",
            "title": "Pairwise O(d) with unequal column counts is invalid without Stiefel",
            "kind": "error",
            "tags": ["error"],
            "problem": {
                "transform": {"group": "O"},
                "source": mat(rng.normal(size=(5, 2))),
                "target": mat(rng.normal(size=(5, 3))),
            },
            "expected": {
                "numerical_status": "invalid_problem",
                "optimality_status": "not_applicable",
                "error_code": "dimension_mismatch",
            },
        }
    )

    # extra pairwise variety for volume
    for i, th in enumerate((0.0, 0.2, 1.0, np.pi / 2, np.pi, 2.3)):
        Xi = rng.normal(size=(7, 2))
        Yi = Xi @ rot2(th) + 0.01 * rng.normal(size=Xi.shape)
        emit(
            pairwise_doc(
                f"pw-o2-angle-{i}",
                f"O(2) noisy rotation angle {th:.2f}",
                Xi,
                Yi,
                "O",
                tags=["O", "batch"],
            )
        )
        emit(
            pairwise_doc(
                f"pw-so2-angle-{i}",
                f"SO(2) noisy rotation angle {th:.2f}",
                Xi,
                Yi,
                "SO",
                tags=["SO", "batch"],
            )
        )

    for i, d in enumerate((2, 3, 4, 6)):
        Xi = rng.normal(size=(15, d))
        Qi, _ = np.linalg.qr(rng.normal(size=(d, d)))
        if np.linalg.det(Qi) < 0:
            Qi[:, 0] *= -1
        Yi = Xi @ Qi
        emit(
            pairwise_doc(
                f"pw-so-d{d}-clean",
                f"SO({d}) exact proper rotation",
                Xi,
                Yi,
                "SO",
                tags=["SO", "batch"],
            )
        )

    # q=1 additional weights
    X1w = rng.normal(size=(9, 1))
    emit(
        pairwise_doc(
            "pw-o1-weighted",
            "O(1) with unequal row weights",
            X1w,
            -3.0 * X1w,
            "O",
            w=np.linspace(0.2, 2.0, 9),
            tags=["d1", "weights"],
        )
    )
    emit(
        pairwise_doc(
            "pw-so1-weighted",
            "SO(1) with unequal row weights cannot flip",
            X1w,
            -3.0 * X1w,
            "SO",
            w=np.linspace(0.2, 2.0, 9),
            tags=["d1", "weights"],
        )
    )

    # certificate documentation fixture (no digits invented for a large GOPP)
    emit(
        {
            "id": "cert-od-unavailable-default",
            "title": "Certificate applies only to the O(d) GOPP trace problem",
            "kind": "certificate",
            "tags": ["certificate"],
            "problem": {
                "transform": {
                    "group": "SO",
                    "translation": True,
                    "scaling": "isotropic",
                },
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "certificate": "unavailable",
                "laws": [
                    "Ling dual certificate is for O(d) GOPP",
                    "do not reuse for SO(d), robust losses, or cell masks",
                    "certificate unavailable is not probably global",
                ],
            },
        }
    )
    emit(
        {
            "id": "cert-gpm-fixed-point-shape",
            "title": "Certificate report fields at a GPM fixed point",
            "kind": "certificate",
            "tags": ["certificate"],
            "problem": {"transform": {"group": "O"}},
            "expected": {
                "numerical_status": "converged",
                "optimality_status": "first_order_stationary",
                "certificate": "not_certified",
                "laws": [
                    "report r_dual = ||CS-Lambda S|| / (1+||CS||)",
                    "report a lower bound on lambda_min(Lambda-C)",
                    "uniqueness if lambda_{d+1}(Lambda-C)>0",
                    "without a nonnegative bound the status stays stationary",
                ],
            },
        }
    )

    # more laws
    emit(
        {
            "id": "law-orthogonal-preserves-norms",
            "title": "O(d) action preserves Frobenius norms and pairwise distances",
            "kind": "law",
            "tags": ["law", "O"],
            "problem": {"transform": {"group": "O"}, "source": mat(X2)},
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "||XR||_F = ||X||_F",
                    "row-wise Euclidean distances are preserved",
                ],
            },
        }
    )
    emit(
        {
            "id": "law-so-det-plus-one",
            "title": "Every SO(d) factor has determinant +1",
            "kind": "law",
            "tags": ["law", "SO"],
            "problem": {"transform": {"group": "SO"}},
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": ["det(R)=+1 for every proper-rotation fit"],
            },
        }
    )
    emit(
        {
            "id": "law-similarity-distance-ratios",
            "title": "Similarity preserves distance ratios",
            "kind": "law",
            "tags": ["law", "similarity"],
            "problem": {
                "transform": {
                    "group": "similarity",
                    "translation": True,
                    "scaling": "isotropic",
                }
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "d(T(x),T(y)) / d(T(u),T(v)) = d(x,y)/d(u,v) when denominators are nonzero"
                ],
            },
        }
    )
    emit(
        {
            "id": "law-huber-rotationally-invariant",
            "title": "Landmark-vector Huber is invariant to a common rotation of residuals",
            "kind": "law",
            "tags": ["law", "robust"],
            "problem": {"transform": {"group": "O"}, "loss": "huber"},
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "rho(||e_j||) is unchanged if every residual is right-multiplied by Q in O(d)",
                    "coordinatewise Huber is not rotationally invariant",
                ],
            },
        }
    )
    emit(
        {
            "id": "law-preshape-unit-frobenius",
            "title": "Preshape configurations have unit Frobenius norm after centering",
            "kind": "law",
            "tags": ["law", "preshape"],
            "problem": {"transform": {"group": "SO", "scaling": "preshape"}},
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "not_applicable",
                "laws": [
                    "||HX / ||HX||_F||_F = 1",
                    "subsequent optimization is rotation-only",
                ],
            },
        }
    )

    # additional pairwise with translation + weights
    emit(
        pairwise_doc(
            "pw-sim-weighted-translated",
            "Weighted similarity with translation",
            X2,
            2.2 * (X2 @ Rtrue) + np.array([5.0, -1.5]),
            "SO",
            scaling="isotropic",
            translation=True,
            w=w,
            tags=["similarity", "weights"],
        )
    )

    # zero configuration weight on one view (the other remains)
    emit(
        {
            "id": "pw-zero-one-config-weight",
            "title": "A single zero configuration weight drops that view",
            "kind": "pairwise",
            "tags": ["weights"],
            "problem": {
                "transform": {"group": "O"},
                "views": [view("keep", X2), view("drop", Y2)],
                "configuration_weights": [1.0, 0.0],
            },
            "expected": {
                "numerical_status": "exact",
                "optimality_status": "exact_closed_form",
                "laws": [
                    "zero configuration weight removes that view from the consensus"
                ],
            },
        }
    )

    manifest = {
        "n_fixtures": n_written,
        "schema": "inst/conformance/schema.json",
        "generator": "tools/generate_conformance.py",
        "seed": 20260820,
    }
    (ROOT / "inst" / "conformance" / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(f"wrote {n_written} fixtures")


if __name__ == "__main__":
    main()
