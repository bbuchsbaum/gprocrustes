#' Thin eigencore wrappers so SVD / eigen work never calls base `svd()` / `eigen()`.
#'
#' @noRd
.gproc_svd <- function(A, nu = min(nrow(A), ncol(A)), nv = min(nrow(A), ncol(A))) {
  A <- as.matrix(A)
  storage.mode(A) <- "double"
  m <- nrow(A)
  n <- ncol(A)
  nu <- as.integer(nu)[1L]
  nv <- as.integer(nv)[1L]
  if (m == 0L || n == 0L) {
    return(list(
      d = numeric(),
      u = matrix(0, m, max(nu, 0L)),
      v = matrix(0, n, max(nv, 0L))
    ))
  }
  rank <- max(1L, min(m, n, max(nu, nv, 1L)))
  vecs <- if (nu <= 0L && nv <= 0L) "none" else "both"
  fit <- eigencore::svd_partial(A, rank = rank, vectors = vecs, certify = FALSE)
  d <- as.numeric(eigencore::values(fit))
  u <- if (nu <= 0L) {
    matrix(0, m, 0L)
  } else {
    U <- eigencore::left_vectors(fit)
    if (is.null(U)) {
      matrix(0, m, 0L)
    } else {
      U[, seq_len(min(nu, ncol(U))), drop = FALSE]
    }
  }
  v <- if (nv <= 0L) {
    matrix(0, n, 0L)
  } else {
    V <- eigencore::right_vectors(fit)
    if (is.null(V)) {
      matrix(0, n, 0L)
    } else {
      V[, seq_len(min(nv, ncol(V))), drop = FALSE]
    }
  }
  list(d = d, u = u, v = v)
}

#' Symmetric eigen, eigenvalues in decreasing order (base `eigen` convention).
#'
#' @noRd
.gproc_eigen_sym <- function(A) {
  A <- as.matrix(A)
  storage.mode(A) <- "double"
  A <- (A + t(A)) / 2
  if (!nrow(A)) {
    return(list(values = numeric(), vectors = matrix(0, 0, 0)))
  }
  fit <- eigencore::eig_full(A, vectors = TRUE)
  vals <- as.numeric(eigencore::values(fit))
  vecs <- eigencore::vectors(fit)
  o <- order(vals, decreasing = TRUE)
  list(values = vals[o], vectors = vecs[, o, drop = FALSE])
}

#' Smallest `k` eigenpairs of a symmetric matrix, via eigencore.
#'
#' @noRd
.gproc_eigen_smallest <- function(A, k) {
  A <- as.matrix(A)
  storage.mode(A) <- "double"
  A <- (A + t(A)) / 2
  k <- min(as.integer(k)[1L], nrow(A))
  if (k < 1L) {
    return(list(values = numeric(), vectors = matrix(0, nrow(A), 0L)))
  }
  fit <- eigencore::eig_partial(
    A,
    k = k,
    target = eigencore::smallest(),
    vectors = TRUE,
    certify = FALSE
  )
  list(
    values = as.numeric(eigencore::values(fit)),
    vectors = eigencore::vectors(fit)
  )
}

#' PSD / symmetric solve without forming an explicit inverse.
#'
#' @noRd
.gproc_psd_solve <- function(A, B, tol = 1e-10) {
  B <- as.matrix(B)
  ev <- .gproc_eigen_sym(A)
  s1 <- if (length(ev$values)) max(ev$values[[1L]], 0) else 0
  keep <- ev$values > tol * max(s1, 1)
  if (!any(keep)) {
    return(matrix(0, nrow(A), ncol(B)))
  }
  U <- ev$vectors[, keep, drop = FALSE]
  U %*% (crossprod(U, B) / ev$values[keep])
}
