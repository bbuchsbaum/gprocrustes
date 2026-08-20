#' Nuisance / orbit tangent basis at a consensus (math §§25–26).
#'
#' Columns are \(\operatorname{vec}(\mathbf{1}e_k^\top)\), \(\operatorname{vec}(M)\)
#' when scale is free, and \(\operatorname{vec}(MA)\) for a skew-symmetric basis.
#'
#' @noRd
.gproc_nuisance_basis <- function(M, spec) {
  n <- nrow(M)
  d <- ncol(M)
  cols <- list()
  if (isTRUE(spec$translation) || identical(spec$family, "orthogonal") ||
      identical(spec$family, "similarity")) {
    for (k in seq_len(d)) {
      A <- matrix(0, n, d)
      A[, k] <- 1
      cols[[length(cols) + 1L]] <- as.vector(A)
    }
  }
  if (identical(spec$scaling, "isotropic")) {
    cols[[length(cols) + 1L]] <- as.vector(M)
  }
  if (d >= 2L && spec$family %in% c("orthogonal", "similarity")) {
    for (i in seq_len(d - 1L)) {
      for (j in seq.int(i + 1L, d)) {
        A <- matrix(0, d, d)
        A[i, j] <- 1
        A[j, i] <- -1
        cols[[length(cols) + 1L]] <- as.vector(M %*% A)
      }
    }
  }
  if (!length(cols)) {
    return(matrix(0, n * d, 0))
  }
  do.call(cbind, cols)
}

#' Apply Kronecker precision \(Q=\Sigma_d^{-1}\otimes\Sigma_N^{-1}\) to \(\operatorname{vec}(Z)\).
#'
#' @noRd
.gproc_apply_vec_precision <- function(z, n, d, covar) {
  Z <- matrix(z, n, d)
  We <- diag(n)
  Wd <- diag(d)
  if (!is.null(covar) && !is.null(covar$entity)) {
    A <- (covar$entity + t(covar$entity)) / 2
    We <- if (isTRUE(covar$as_precision)) A else {
      ev <- .gproc_eigen_sym(A)
      ev$vectors %*% (ifelse(ev$values > 1e-12, 1 / ev$values, 0) * t(ev$vectors))
    }
  }
  if (!is.null(covar) && !is.null(covar$coordinate)) {
    A <- (covar$coordinate + t(covar$coordinate)) / 2
    Wd <- if (isTRUE(covar$as_precision)) A else {
      ev <- .gproc_eigen_sym(A)
      ev$vectors %*% (ifelse(ev$values > 1e-12, 1 / ev$values, 0) * t(ev$vectors))
    }
  }
  as.vector(We %*% Z %*% Wd)
}

#' Horizontal projection \(z_H=P_H\operatorname{vec}(Z)\) at \(M\).
#'
#' @noRd
.gproc_horizontal_project <- function(Z, M, spec, covar = NULL) {
  n <- nrow(M)
  d <- ncol(M)
  B <- .gproc_nuisance_basis(M, spec)
  z <- as.vector(Z)
  if (!ncol(B)) {
    return(Z)
  }
  QB <- apply(B, 2L, function(col) .gproc_apply_vec_precision(col, n, d, covar))
  if (!is.matrix(QB)) {
    QB <- matrix(QB, nrow = nrow(B), ncol = ncol(B))
  }
  BtQB <- crossprod(B, QB)
  Qz <- .gproc_apply_vec_precision(z, n, d, covar)
  sv <- .gproc_svd(BtQB)
  pos <- sv$d > 1e-10 * max(sv$d, 0)
  if (!any(pos)) {
    return(Z)
  }
  coef <- sv$v[, pos, drop = FALSE] %*% ((t(sv$u[, pos, drop = FALSE]) %*% crossprod(B, Qz)) / sv$d[pos])
  z_h <- z - as.vector(B %*% coef)
  matrix(z_h, n, d)
}

#' Tangent / horizontal coordinates of aligned residuals.
#'
#' Residual matrices are projected off the similarity (or fitted-group) orbit
#' at the consensus. The result is a local chart, not a global shape space.
#'
#' @param object A `gpa_fit`.
#' @param model Optional `proc_shape_model` whose covariance supplies \(Q\).
#'   Fitting-metric \(Q\) is used only if the model says so.
#' @param ... Unused.
#' @return A list of entity-by-dimension tangent residual matrices.
#' @export
tangent_coordinates <- function(object, ...) {
  UseMethod("tangent_coordinates")
}

#' @export
tangent_coordinates.gpa_fit <- function(object, model = NULL, ...) {
  M <- consensus(object)
  spec <- object$problem$transform
  covar <- .gproc_inference_covariance(object, model)
  Y <- .gproc_fitted_global(object)
  out <- lapply(Y, function(Yi) {
    E <- Yi - M
    E[!is.finite(E)] <- 0
    .gproc_horizontal_project(E, M, spec, covar)
  })
  attr(out, "nuisance_dimension") <- ncol(.gproc_nuisance_basis(M, spec))
  attr(out, "approximation") <- "local tangent chart at the fitted consensus"
  out
}

#' Project a residual or new configuration into the horizontal space.
#'
#' @param newdata A residual matrix, a configuration, or `proc_data` with one view.
#' @param fit A `gpa_fit` that supplies the consensus and gauge.
#' @param as_residual If `TRUE`, `newdata` is already \(Y-M\). If `FALSE`,
#'   it is aligned to the consensus first.
#' @export
tangent_project <- function(newdata, fit, as_residual = FALSE) {
  if (!inherits(fit, "gpa_fit")) {
    .gproc_stop("invalid_problem", "tangent_project() requires a gpa_fit.")
  }
  M <- consensus(fit)
  spec <- fit$problem$transform
  if (inherits(newdata, "proc_data")) {
    newdata <- newdata$views[[1L]]
  }
  Z <- as.matrix(newdata)
  if (!isTRUE(as_residual)) {
    if (nrow(Z) == nrow(M) && ncol(Z) == ncol(M)) {
      pair <- procrustes(Z, M, transform = spec)
      Z <- apply_proc_transform(pair$transform, Z) - M
    } else {
      .gproc_stop("invalid_problem", "newdata must match the consensus dimensions.")
    }
  }
  .gproc_horizontal_project(Z, M, spec, NULL)
}
