#' Weighted column sums without forming a centered copy.
#'
#' Computes \eqn{X^\top w} for a vector of row weights.
#'
#' @param X Matrix-like object, \eqn{n \times d}.
#' @param w Nonnegative row weights of length \eqn{n}.
#' @return Numeric vector of length \eqn{d}.
#' @export
weighted_col_sums <- function(X, w) {
  X <- .gproc_as_numeric_matrix(X, "X")
  w <- .gproc_row_weights(w, .gproc_nrow(X))
  as.numeric(Matrix::crossprod(X, w))
}

#' Centered weighted cross-product via sufficient statistics.
#'
#' Implements the identity
#' \eqn{X_c^\top W Y_c = X^\top W Y - (X^\top w)(Y^\top w)^\top / (\mathbf{1}^\top w)}
#' so a sparse \eqn{X} is never explicitly centered.
#'
#' @param X,Y Matrix-like objects with the same number of rows.
#' @param w Nonnegative row weights. Recycled to 1 if omitted.
#' @return A \eqn{p_X \times p_Y} numeric matrix.
#' @export
centered_crossprod <- function(X, Y, w = NULL) {
  X <- .gproc_as_numeric_matrix(X, "X")
  Y <- .gproc_as_numeric_matrix(Y, "Y")
  if (.gproc_nrow(X) != .gproc_nrow(Y)) {
    .gproc_stop("dimension_mismatch", "X and Y must have the same number of rows.")
  }
  n <- .gproc_nrow(X)
  w <- .gproc_row_weights(w, n)
  wp <- sum(w)
  if (!(wp > 0)) {
    .gproc_stop("zero_total_configuration_weight", "Row weights must have positive sum.")
  }
  wy <- .gproc_row_scale(Y, w)
  raw <- as.matrix(Matrix::crossprod(X, wy))
  xtw <- as.numeric(Matrix::crossprod(X, w))
  ytw <- as.numeric(Matrix::crossprod(Y, w))
  raw - tcrossprod(xtw, ytw) / wp
}

#' Weighted second-moment scalar \eqn{\operatorname{tr}(X_c^\top W X_c)}.
#'
#' @param X Matrix-like object.
#' @param w Nonnegative row weights. Recycled to 1 if omitted.
#' @return A single nonnegative number.
#' @export
centered_trace <- function(X, w = NULL) {
  Cxx <- centered_crossprod(X, X, w)
  sum(diag(Cxx))
}

#' Weighted centroids without forming \eqn{X - 1\bar x}.
#'
#' @param X Matrix-like object.
#' @param w Nonnegative row weights. Recycled to 1 if omitted.
#' @return Numeric vector of column means.
#' @export
weighted_centroid <- function(X, w = NULL) {
  X <- .gproc_as_numeric_matrix(X, "X")
  w <- .gproc_row_weights(w, .gproc_nrow(X))
  wp <- sum(w)
  if (!(wp > 0)) {
    .gproc_stop("zero_total_configuration_weight", "Row weights must have positive sum.")
  }
  weighted_col_sums(X, w) / wp
}

#' Pairwise sufficient statistics used by the Procrustes kernel.
#'
#' @inheritParams centered_crossprod
#' @param center If `TRUE`, return translation-eliminated moments
#'   \eqn{X_c^\top W Y_c}. If `FALSE`, return the raw moments
#'   \eqn{X^\top W Y} used when translation is prohibited.
#' @return A list with `C`, `a`, `b`, `xbar`, `ybar`, `w_sum`, and `centered`.
#' @export
proc_moments <- function(X, Y, w = NULL, center = TRUE) {
  X <- .gproc_as_numeric_matrix(X, "X")
  Y <- .gproc_as_numeric_matrix(Y, "Y")
  .gproc_check_finite(X, "X")
  .gproc_check_finite(Y, "Y")
  w <- .gproc_row_weights(w, .gproc_nrow(X))
  xbar <- weighted_centroid(X, w)
  ybar <- weighted_centroid(Y, w)
  if (isTRUE(center)) {
    list(
      C = centered_crossprod(X, Y, w),
      a = centered_trace(X, w),
      b = centered_trace(Y, w),
      xbar = xbar,
      ybar = ybar,
      w_sum = sum(w),
      centered = TRUE
    )
  } else {
    list(
      C = .gproc_raw_crossprod(X, Y, w),
      a = .gproc_raw_trace(X, w),
      b = .gproc_raw_trace(Y, w),
      xbar = xbar,
      ybar = ybar,
      w_sum = sum(w),
      centered = FALSE
    )
  }
}

#' @noRd
.gproc_raw_crossprod <- function(X, Y, w) {
  as.matrix(Matrix::crossprod(X, .gproc_row_scale(Y, w)))
}

#' @noRd
.gproc_raw_trace <- function(X, w) {
  X <- as.matrix(X)
  sum(w * rowSums(X * X))
}

#' Weighted residual of a stored transform against a target.
#'
#' @noRd
.gproc_transform_residual <- function(X, Y, tr, w) {
  Yhat <- as.matrix(apply_proc_transform(tr, X))
  Ym <- as.matrix(Y)
  sum(w * rowSums((Yhat - Ym) * (Yhat - Ym)))
}

#' @keywords internal
.gproc_row_weights <- function(w, n) {
  if (is.null(w)) {
    return(rep(1, n))
  }
  w <- as.numeric(w)
  if (length(w) == 1L) {
    w <- rep(w, n)
  }
  if (length(w) != n) {
    .gproc_stop("invalid_problem", "Row weights must have length 1 or nrow(X).")
  }
  if (any(!is.finite(w)) || any(w < 0)) {
    .gproc_stop("invalid_problem", "Row weights must be finite and nonnegative.")
  }
  w
}

#' @keywords internal
.gproc_row_scale <- function(Y, w) {
  if (inherits(Y, "Matrix")) {
    return(Matrix::Diagonal(x = w) %*% Y)
  }
  Y * w
}
