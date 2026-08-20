#' Pairwise Procrustes kernel.
#'
#' Exact closed-form fit of `X` to `Y` under \(O(d)\), \(SO(d)\), or
#' similarity. Translation is eliminated analytically; centering uses
#' sufficient statistics and never forms \(X_c\) explicitly.
#'
#' @param X Source configuration, \(n \times d\).
#' @param Y Target configuration, \(n \times d\).
#' @param transform Transform specification or shortcut.
#' @param weights Row weights \(w\).
#' @param rank_tol Relative singular-value threshold \(\tau\).
#' @return A `proc_pair_fit`.
#' @export
procrustes <- function(X,
                       Y,
                       transform = proc_orthogonal("O"),
                       weights = NULL,
                       rank_tol = 1e-10) {
  spec <- as_proc_transform(transform)
  X <- .gproc_as_numeric_matrix(X, "X")
  Y <- .gproc_as_numeric_matrix(Y, "Y")
  .gproc_check_finite(X, "X")
  .gproc_check_finite(Y, "Y")
  if (.gproc_nrow(X) != .gproc_nrow(Y)) {
    .gproc_stop("dimension_mismatch", "X and Y must have the same number of rows.")
  }
  if (.gproc_ncol(X) != .gproc_ncol(Y)) {
    .gproc_stop(
      "dimension_mismatch",
      "Pairwise O(d)/SO(d)/similarity requires equal column counts."
    )
  }
  w <- .gproc_row_weights(weights, .gproc_nrow(X))
  if (identical(spec$family, "signed_permutation")) {
    return(.procrustes_signed_perm(X, Y, w, spec, rank_tol))
  }
  moms <- proc_moments(X, Y, w)
  polar <- .polar_factor(moms$C, group = spec$group, rank_tol = rank_tol)
  R <- polar$R
  gamma <- sum(R * moms$C)
  a <- moms$a
  b <- moms$b
  if (identical(spec$scaling, "isotropic")) {
    s <- if (a > 0) max(gamma / a, 0) else 0
  } else {
    s <- 1
  }
  t <- if (isTRUE(spec$translation)) {
    as.numeric(moms$ybar - s * (moms$xbar %*% R))
  } else {
    rep(0, length(moms$xbar))
  }
  objective <- .gproc_centered_residual(X, Y, R, s, w, moms$xbar, moms$ybar)
  tr <- .proc_fitted(spec, R, s, t)
  structure(
    list(
      transform = tr,
      moments = moms,
      gamma = gamma,
      objective = unname(objective),
      rank = polar$rank,
      singular_values = polar$singular_values,
      transform_unique = polar$transform_unique,
      objective_unique = TRUE,
      unidentified_subspace_dimension = polar$unidentified,
      numerical_status = "exact",
      optimality_status = "exact_closed_form",
      solver = "pairwise_polar"
    ),
    class = "proc_pair_fit"
  )
}

#' Polar / proper-polar factor of a cross-covariance.
#'
#' @param C Cross-covariance matrix.
#' @param group `"O"` or `"SO"`.
#' @param rank_tol Relative rank threshold.
#' @return List with `R`, singular values, rank, and uniqueness.
#' @export
polar_factor <- function(C, group = c("O", "SO"), rank_tol = 1e-10) {
  group <- match.arg(group)
  .polar_factor(as.matrix(C), group, rank_tol)
}

#' @keywords internal
.polar_factor <- function(C, group, rank_tol = 1e-10) {
  C <- as.matrix(C)
  storage.mode(C) <- "double"
  d <- ncol(C)
  if (nrow(C) != d) {
    .gproc_stop("dimension_mismatch", "Pairwise polar factor expects a square cross-covariance.")
  }
  if (d == 0L) {
    return(list(R = C, singular_values = numeric(), rank = 0L,
                transform_unique = TRUE, unidentified = 0L))
  }
  if (!all(is.finite(C))) {
    .gproc_stop("nonfinite_values", "Nonfinite values in the cross-covariance.")
  }
  sv <- .gproc_svd(C, nu = d, nv = d)
  sigma <- sv$d
  sigma1 <- if (length(sigma)) max(sigma[1L], 0) else 0
  rank <- if (sigma1 <= 0) 0L else as.integer(sum(sigma > rank_tol * sigma1))
  if (rank == 0L) {
    R <- diag(d)
  } else {
    R <- sv$u %*% t(sv$v)
    if (identical(group, "SO") && det(R) < 0) {
      sv$u[, d] <- -sv$u[, d]
      R <- sv$u %*% t(sv$v)
    }
  }
  unique <- rank == d
  list(
    R = R,
    singular_values = sigma,
    rank = rank,
    transform_unique = unique,
    unidentified = as.integer(max(d - rank, 0L))
  )
}

#' @keywords internal
.procrustes_signed_perm <- function(X, Y, w, spec, rank_tol) {
  # Component matching uses the raw cross-covariance, not the centered one.
  C <- as.matrix(Matrix::crossprod(X, .gproc_row_scale(Y, w)))
  moms <- proc_moments(X, Y, w)
  moms$C <- C
  d <- ncol(C)
  if (d > 8L) {
    .gproc_stop("invalid_problem", "Signed-permutation brute-force solver supports d <= 8.")
  }
  perms <- .gproc_permutations(d)
  best_score <- -Inf
  best_p <- seq_len(d)
  for (p in perms) {
    vals <- C[cbind(p, seq_len(d))]
    score <- sum(abs(vals))
    if (score > best_score) {
      best_score <- score
      best_p <- p
    }
  }
  P <- diag(d)[best_p, , drop = FALSE]
  signs <- sign(C[cbind(best_p, seq_len(d))])
  signs[signs == 0] <- 1
  D <- diag(signs, d, d)
  R <- P %*% D
  tr <- .proc_fitted(spec, R, 1, rep(0, d))
  structure(
    list(
      transform = tr,
      moments = moms,
      gamma = sum(R * C),
      objective = moms$a + moms$b - 2 * sum(R * C),
      rank = .polar_factor(C, "O", rank_tol)$rank,
      singular_values = .gproc_svd(C, nu = 0L, nv = 0L)$d,
      transform_unique = TRUE,
      objective_unique = TRUE,
      unidentified_subspace_dimension = 0L,
      numerical_status = "exact",
      optimality_status = "exact_closed_form",
      solver = "signed_permutation"
    ),
    class = "proc_pair_fit"
  )
}

#' @keywords internal
.gproc_centered_residual <- function(X, Y, R, s, w, xbar, ybar) {
  XR <- as.matrix(X %*% (s * R))
  Ym <- as.matrix(Y)
  shift <- as.numeric(s * as.numeric(xbar) %*% R - ybar)
  dlt <- sweep(XR - Ym, 2L, shift, "-")
  sum(w * rowSums(dlt * dlt))
}

#' @keywords internal
.gproc_permutations <- function(d) {
  if (d == 1L) {
    return(list(1L))
  }
  out <- list()
  rec <- function(head, rest) {
    if (!length(rest)) {
      out[[length(out) + 1L]] <<- head
      return(invisible())
    }
    for (i in seq_along(rest)) {
      rec(c(head, rest[i]), rest[-i])
    }
  }
  rec(integer(), seq_len(d))
  out
}

#' Align `x` to `reference` by a signed permutation of columns.
#'
#' @param x,reference Matrices with the same dimensions.
#' @param transform Ignored unless `"signed_permutation"`.
#' @export
match_components <- function(x, reference, transform = "signed_permutation") {
  if (!identical(as_proc_transform(transform)$family, "signed_permutation")) {
    .gproc_stop("invalid_problem", "match_components() currently supports signed_permutation only.")
  }
  procrustes(x, reference, transform = proc_signed_permutation())
}

#' @export
print.proc_pair_fit <- function(x, ...) {
  d <- nrow(x$transform$R)
  cat(sprintf(
    paste0(
      "Pairwise Procrustes fit\n",
      "  Group: %s\n",
      "  Effective rank: %d of %d\n",
      "  Transform uniqueness: %s\n",
      "  Objective uniqueness: %s\n",
      "  Unidentified subspace dimension: %d\n",
      "  Objective: %.6g\n",
      "  Numerical status: %s\n",
      "  Optimality: %s\n"
    ),
    x$transform$spec$group,
    x$rank, d,
    if (x$transform_unique) "yes" else "no",
    if (x$objective_unique) "yes" else "no",
    x$unidentified_subspace_dimension,
    x$objective,
    x$numerical_status,
    x$optimality_status
  ))
  invisible(x)
}
