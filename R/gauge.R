#' Quotient / Procrustes distance after the optimal gauge alignment.
#'
#' For two consensuses this is \eqn{\min_Q\|XQ-Y\|_F}. For two rotation
#' stacks it is Ling's \eqn{d_F(S,T)=\min_Q\|S-TQ\|_F}.
#'
#' @param x,y Matrices with the same column dimension.
#' @param group `"O"` or `"SO"`.
#' @return Nonnegative quotient distance.
#' @export
#' @examples
#' X <- scale(matrix(rnorm(20), 10, 2), scale = FALSE)
#' Q <- matrix(c(0, -1, 1, 0), 2, 2)
#' quotient_distance(X, X %*% Q)
quotient_distance <- function(x, y, group = c("O", "SO")) {
  group <- match.arg(group)
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (ncol(x) != ncol(y) || nrow(x) != nrow(y)) {
    .gproc_stop("dimension_mismatch", "quotient_distance() needs matching dimensions.")
  }
  pair <- procrustes(x, y, transform = proc_orthogonal(group), center = FALSE)
  dlt <- x %*% pair$transform$R - y
  sqrt(sum(dlt * dlt))
}

#' Principal angles between two subspaces.
#'
#' @param U,V Matrices whose columns span the subspaces.
#' @return Numeric vector of principal angles in \eqn{[0, \pi/2]}.
#' @export
#' @examples
#' principal_angles(diag(2), matrix(c(0, 1, 1, 0), 2, 2))
principal_angles <- function(U, V) {
  U <- as.matrix(U)
  V <- as.matrix(V)
  if (!nrow(U) || !nrow(V)) {
    return(numeric())
  }
  QU <- qr.Q(qr(U))
  QV <- qr.Q(qr(V))
  k <- min(ncol(QU), ncol(QV))
  if (!k) {
    return(numeric())
  }
  sv <- .gproc_svd(crossprod(QU, QV), nu = 0L, nv = 0L)
  acos(pmin(1, pmax(0, sv$d[seq_len(k)])))
}

#' Grassmann distance \eqn{\bigl(\sum_j\sin^2\theta_j\bigr)^{1/2}}.
#'
#' @inheritParams principal_angles
#' @return Nonnegative Grassmann distance.
#' @export
#' @examples
#' subspace_distance(diag(2), matrix(c(0, 1, 1, 0), 2, 2))
subspace_distance <- function(U, V) {
  theta <- principal_angles(U, V)
  sqrt(sum(sin(theta)^2))
}

#' Horizontal part of a consensus perturbation (math sections 25-26).
#'
#' @param Z Perturbation matrix.
#' @param M Consensus.
#' @param spec Transform family that defines the orbit.
#' @return A matrix the same shape as `Z`.
#' @export
#' @examples
#' M <- scale(matrix(rnorm(20), 10, 2), scale = FALSE)
#' Z <- matrix(rnorm(20), 10, 2)
#' horizontal_component(Z, M)
horizontal_component <- function(Z, M, spec = proc_orthogonal("O")) {
  spec <- as_proc_transform(spec)
  .gproc_horizontal_project(as.matrix(Z), as.matrix(M), spec, NULL)
}

#' Vertical / orbit part of a consensus perturbation.
#'
#' @inheritParams horizontal_component
#' @return A matrix the same shape as `Z`.
#' @export
#' @examples
#' M <- scale(matrix(rnorm(20), 10, 2), scale = FALSE)
#' Z <- matrix(rnorm(20), 10, 2)
#' vertical_component(Z, M)
vertical_component <- function(Z, M, spec = proc_orthogonal("O")) {
  Z <- as.matrix(Z)
  Z - horizontal_component(Z, M, spec)
}

#' Eigenvalue clusters of \eqn{M^\top M} (tied canonical axes).
#'
#' @param M Consensus matrix.
#' @param rel_tol Relative gap below which consecutive eigenvalues are tied.
#' @return A list of column-index blocks.
#' @export
#' @examples
#' tied_axis_blocks(diag(c(3, 2, 2)))
tied_axis_blocks <- function(M, rel_tol = 1e-6) {
  M <- as.matrix(M)
  if (!ncol(M)) {
    return(list())
  }
  ev <- .gproc_eigen_sym(crossprod(M))
  lam <- ev$values
  blocks <- list()
  start <- 1L
  for (i in seq_len(length(lam) - 1L)) {
    scale <- max(abs(lam[[i]]), abs(lam[[i + 1L]]), 1e-15)
    if (abs(lam[[i]] - lam[[i + 1L]]) > rel_tol * scale) {
      blocks[[length(blocks) + 1L]] <- start:i
      start <- i + 1L
    }
  }
  blocks[[length(blocks) + 1L]] <- start:length(lam)
  blocks
}

#' @noRd
.gproc_canonical_consensus <- function(M, group = "O") {
  d <- ncol(M)
  sv <- .gproc_svd(M, nu = 0L, nv = d)
  Q <- sv$v
  if (ncol(Q) && det(Q) < 0 && identical(group, "SO")) {
    Q[, d] <- -Q[, d]
  }
  M %*% Q
}
