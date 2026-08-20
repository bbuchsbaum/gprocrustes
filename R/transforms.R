#' Orthogonal transformation family.
#'
#' @param group `"O"` (reflections allowed) or `"SO"` (proper rotations).
#' @return A `proc_transform_spec` object.
#' @export
proc_orthogonal <- function(group = c("O", "SO")) {
  group <- match.arg(group)
  structure(
    list(family = "orthogonal", group = group, translation = FALSE, scaling = "none"),
    class = c("proc_orthogonal", "proc_transform_spec")
  )
}

#' Similarity transformation family.
#'
#' @param group Orthogonal subgroup for the linear part.
#' @param translation Logical; estimate a translation.
#' @param scaling `"none"` or `"isotropic"`.
#' @return A `proc_transform_spec` object.
#' @export
proc_similarity <- function(group = c("SO", "O"),
                            translation = TRUE,
                            scaling = c("isotropic", "none")) {
  group <- match.arg(group)
  scaling <- match.arg(scaling)
  structure(
    list(
      family = "similarity",
      group = group,
      translation = isTRUE(translation),
      scaling = scaling
    ),
    class = c("proc_similarity", "proc_transform_spec")
  )
}

#' Signed-permutation transformation family.
#'
#' @return A `proc_transform_spec` object.
#' @export
proc_signed_permutation <- function() {
  structure(
    list(family = "signed_permutation", group = "signed_permutation",
         translation = FALSE, scaling = "none"),
    class = c("proc_signed_permutation", "proc_transform_spec")
  )
}

#' Affine maps under Bai-Bartoli reference-space constraints.
#'
#' Right action \eqn{T(X)=XA+\mathbf{1}t}, written as an LBW with a homogeneous
#' column. `reference_covariance` is \eqn{\Lambda} in \eqn{M^\top M=\Lambda}.
#' It is a modeling choice, never estimated invisibly.
#'
#' @param reference_covariance Scalar or length-\eqn{d} vector \eqn{\lambda}.
#' @export
proc_affine <- function(reference_covariance = 1) {
  structure(
    list(
      family = "affine",
      group = "affine",
      translation = TRUE,
      scaling = "none",
      smoothness = 0,
      reference_covariance = reference_covariance,
      warp = "affine"
    ),
    class = c("proc_affine", "proc_lbw_spec", "proc_transform_spec")
  )
}

#' Linear-basis warp: \eqn{T(X)=\Phi B} with quadratic penalty \eqn{\operatorname{tr}(B^\top LB)}.
#'
#' @param basis Function `function(X)` returning \eqn{\Phi}, or `"affine"`.
#' @param penalty Symmetric penalty \eqn{L}, a function of \eqn{\Phi}, or `NULL`.
#' @param smoothness Multiplier \eqn{\mu\ge 0}.
#' @param reference_covariance \eqn{\Lambda} in \eqn{M^\top M=\Lambda}.
#' @export
proc_lbw <- function(basis,
                     penalty = NULL,
                     smoothness = 0,
                     reference_covariance = 1) {
  structure(
    list(
      family = "lbw",
      group = "lbw",
      translation = TRUE,
      scaling = "none",
      basis = basis,
      penalty = penalty,
      smoothness = smoothness,
      reference_covariance = reference_covariance,
      warp = "lbw"
    ),
    class = c("proc_lbw", "proc_lbw_spec", "proc_transform_spec")
  )
}

#' Thin-plate spline as an LBW constructor.
#'
#' @param control_points Optional shared control-point matrix. `NULL` uses
#'   each view's observed landmarks.
#' @param smoothness Bending-energy weight \eqn{\mu}.
#' @param reference_covariance \eqn{\Lambda} in \eqn{M^\top M=\Lambda}.
#' @export
proc_tps <- function(control_points = NULL,
                     smoothness = 1,
                     reference_covariance = 1) {
  structure(
    list(
      family = "tps",
      group = "tps",
      translation = TRUE,
      scaling = "none",
      control_points = control_points,
      smoothness = smoothness,
      reference_covariance = reference_covariance,
      warp = "tps"
    ),
    class = c("proc_tps", "proc_lbw_spec", "proc_transform_spec")
  )
}

#' Expand a string shortcut into a transform specification.
#'
#' @param transform A `proc_transform_spec` or a string such as `"similarity"`.
#' @return A `proc_transform_spec`.
#' @export
as_proc_transform <- function(transform) {
  if (inherits(transform, "proc_transform_spec")) {
    return(transform)
  }
  if (is.list(transform) && !is.null(transform$group)) {
    grp <- transform$group
    trans <- isTRUE(transform$translation)
    scaling <- transform$scaling %||% "none"
    rot <- transform$rotation_group %||% if (grp %in% c("O", "SO")) grp else "O"
    if (identical(grp, "signed_permutation")) {
      return(proc_signed_permutation())
    }
    if (identical(grp, "none")) {
      return(structure(list(family = "none", group = "none",
                            translation = FALSE, scaling = "none"),
                       class = c("proc_none", "proc_transform_spec")))
    }
    if (identical(grp, "similarity") || trans || !identical(scaling, "none")) {
      g <- if (grp %in% c("O", "SO")) grp else rot
      if (!g %in% c("O", "SO")) g <- "O"
      sc <- if (scaling %in% c("isotropic", "none")) scaling else "none"
      return(proc_similarity(group = g, translation = trans, scaling = sc))
    }
    if (grp %in% c("O", "SO")) {
      return(proc_orthogonal(group = grp))
    }
  }
  if (!is.character(transform) || length(transform) != 1L) {
    .gproc_stop("invalid_problem", "Unknown transformation specification.")
  }
  switch(
    transform,
    orthogonal = proc_orthogonal("O"),
    O = proc_orthogonal("O"),
    SO = proc_orthogonal("SO"),
    rotation = proc_orthogonal("SO"),
    similarity = proc_similarity(),
    signed_permutation = proc_signed_permutation(),
    affine = proc_affine(),
    tps = proc_tps(),
    lbw = proc_lbw("affine"),
    .gproc_stop("invalid_problem", sprintf("Unknown transform shortcut '%s'.", transform))
  )
}

#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Identity fitted transform.
#'
#' @param spec A transform specification.
#' @param d Target dimension.
#' @export
proc_identity <- function(spec, d) {
  spec <- as_proc_transform(spec)
  .proc_fitted(
    spec = spec,
    R = diag(d),
    s = 1,
    t = rep(0, d)
  )
}

#' Apply a fitted transform: \eqn{T(X) = sXR + \mathbf{1}t}.
#'
#' @param tr A `proc_fitted_transform`.
#' @param X An \eqn{n \times d} matrix-like object.
#' @return The transformed matrix (dense if `R` mixes sparse columns).
#' @export
apply_proc_transform <- function(tr, X) {
  X <- .gproc_as_numeric_matrix(X, "X")
  if (is.function(tr$apply_warp)) {
    return(tr$apply_warp(X))
  }
  if (.gproc_ncol(X) != nrow(tr$R)) {
    .gproc_stop("dimension_mismatch", "X does not match the transform dimension.")
  }
  Y <- tr$s * (X %*% tr$R)
  if (any(tr$t != 0)) {
    Y <- sweep(as.matrix(Y), 2L, tr$t, "+")
  }
  Y
}

#' Inverse of a fitted similarity / orthogonal transform.
#'
#' @param tr A `proc_fitted_transform`.
#' @return A `proc_fitted_transform`.
#' @export
inverse_proc_transform <- function(tr) {
  if (!inherits(tr, "proc_fitted_transform")) {
    .gproc_stop("invalid_problem", "inverse_proc_transform() needs a fitted transform.")
  }
  if (identical(tr$spec$family, "affine") && is.null(tr$apply_warp)) {
    if (abs(det(tr$R)) <= 1e-12) {
      .gproc_stop("inverse_unavailable", "Affine linear part is singular; no exact inverse.")
    }
    Rinv <- solve(tr$R)
    tinv <- as.numeric(-tr$t %*% Rinv)
    return(.proc_fitted(spec = tr$spec, R = Rinv, s = 1, t = tinv))
  }
  if (!tr$spec$family %in% c("orthogonal", "similarity", "signed_permutation")) {
    .gproc_stop(
      "inverse_unavailable",
      "No exact global inverse exists for this transform family; do not invent one."
    )
  }
  if (!(tr$s > 0)) {
    .gproc_stop("inverse_unavailable", "Scale is zero; the transform is not invertible.")
  }
  Rinv <- t(tr$R)
  sinv <- 1 / tr$s
  tinv <- as.numeric(-sinv * tr$t %*% Rinv)
  .proc_fitted(spec = tr$spec, R = Rinv, s = sinv, t = tinv)
}

#' Compose two fitted transforms: `tr2` after `tr1`.
#'
#' @param tr1,tr2 Fitted transforms.
#' @export
compose_proc_transform <- function(tr1, tr2) {
  if (is.function(tr1$apply_warp) || is.function(tr2$apply_warp) ||
      tr1$spec$family %in% c("tps", "lbw") || tr2$spec$family %in% c("tps", "lbw")) {
    .gproc_stop(
      "inverse_unavailable",
      "Composition is not defined for a general linear-basis warp."
    )
  }
  R <- tr1$R %*% tr2$R
  s <- tr1$s * tr2$s
  t <- as.numeric(tr2$s * tr1$t %*% tr2$R + tr2$t)
  spec <- .gproc_composed_spec(tr1, tr2, R, s, t)
  .proc_fitted(spec = spec, R = R, s = s, t = t)
}

#' Spec of the composed action, not a copy of `tr2$spec`.
#'
#' @noRd
.gproc_composed_spec <- function(tr1, tr2, R, s, t) {
  g1 <- tr1$spec$group
  g2 <- tr2$spec$group
  detR <- if (is.matrix(R) && nrow(R) == ncol(R) && nrow(R) > 0L) det(R) else 1
  so <- identical(g1, "SO") && identical(g2, "SO") && is.finite(detR) && detR > 0
  group <- if (g1 %in% c("O", "SO") && g2 %in% c("O", "SO")) {
    if (so) "SO" else "O"
  } else {
    g2
  }
  if (is.finite(detR) && detR < 0 && identical(group, "SO")) {
    group <- "O"
  }
  has_s <- abs(s - 1) > 1e-12 ||
    identical(tr1$spec$scaling, "isotropic") ||
    identical(tr2$spec$scaling, "isotropic")
  has_t <- any(abs(t) > 1e-12) ||
    isTRUE(tr1$spec$translation) ||
    isTRUE(tr2$spec$translation)
  if (identical(tr1$spec$family, "signed_permutation") &&
      identical(tr2$spec$family, "signed_permutation") && !has_s && !has_t) {
    return(proc_signed_permutation())
  }
  if (has_s || has_t) {
    g <- if (group %in% c("O", "SO")) group else "O"
    return(proc_similarity(
      group = g,
      translation = has_t,
      scaling = if (has_s) "isotropic" else "none"
    ))
  }
  if (group %in% c("O", "SO")) {
    return(proc_orthogonal(group))
  }
  tr2$spec
}

#' @keywords internal
.proc_fitted <- function(spec, R, s, t, B = NULL, apply_warp = NULL, extra = NULL) {
  structure(
    list(
      spec = spec,
      R = as.matrix(R),
      s = as.numeric(s)[1L],
      t = as.numeric(t),
      B = B,
      apply_warp = apply_warp,
      extra = extra
    ),
    class = "proc_fitted_transform"
  )
}

#' @export
print.proc_transform_spec <- function(x, ...) {
  cat(sprintf(
    "Procrustes transform spec: %s (%s), translation=%s, scaling=%s\n",
    x$family, x$group, x$translation, x$scaling
  ))
  invisible(x)
}

#' @export
print.proc_fitted_transform <- function(x, ...) {
  cat(sprintf(
    "Fitted %s transform: d=%d, s=%.6g, translation=%s\n",
    x$spec$group, nrow(x$R), x$s, any(x$t != 0)
  ))
  invisible(x)
}

#' Degrees of freedom of a transform family in dimension `d`.
#'
#' @param spec Transform specification.
#' @param d Dimension.
#' @export
degrees_of_freedom <- function(spec, d) {
  spec <- as_proc_transform(spec)
  rot <- if (identical(spec$group, "SO")) d * (d - 1) / 2 else d * (d - 1) / 2
  # O(d) and SO(d) have the same dimension; O has two components.
  sc <- if (identical(spec$scaling, "isotropic")) 1 else 0
  tr <- if (isTRUE(spec$translation)) d else 0
  if (identical(spec$family, "signed_permutation")) {
    return(0)
  }
  if (identical(spec$family, "affine")) {
    return(as.integer(d * d + d))
  }
  if (spec$family %in% c("lbw", "tps")) {
    return(NA_integer_)
  }
  as.integer(rot + sc + tr)
}

#' Constraint residual \eqn{\|R^\top R - I\|_F} (and det for SO).
#'
#' @param tr Fitted transform.
#' @export
constraint_residual <- function(tr) {
  if (tr$spec$family %in% c("affine", "lbw", "tps")) {
    return(list(orthogonality = NA_real_, determinant = det(tr$R)))
  }
  R <- tr$R
  orth <- sqrt(sum((crossprod(R) - diag(nrow(R)))^2))
  det_err <- if (identical(tr$spec$group, "SO")) abs(det(R) - 1) else 0
  list(orthogonality = orth, determinant = det_err)
}
