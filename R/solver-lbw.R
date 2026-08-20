#' Variable-projection LBW / affine / TPS GPA (Bai–Bartoli).
#'
#' Eliminates \(B_i\) and recovers \(M\) from the bottom eigenvectors of \(P\).
#' The only allowed closed-form claim is for this constrained reference-space
#' formulation after the free-translation check.
#'
#' @noRd
.lbw_gpa <- function(data, spec, gauge, metric, alphas, control, anchor = NULL) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  if (any(vapply(data$views, ncol, integer(1)) != d)) {
    .gproc_stop("dimension_mismatch", "LBW GPA requires a common target dimension.")
  }
  views <- .lbw_prepare_views(data, spec, metric)
  free <- vapply(views, `[[`, logical(1), "free_translation")
  if (!all(free)) {
    .gproc_stop(
      "invalid_problem",
      "Free-translation check failed: Phi a = 1 and L a = 0 must hold."
    )
  }
  P <- .lbw_form_P(views, alphas, n)
  ones <- rep(1, n)
  p1 <- as.numeric(P %*% ones)
  free_num <- sqrt(sum(p1^2)) <= 1e-8 * max(1, sqrt(sum(P^2)))
  ev <- .gproc_eigen_smallest(P, k = min(n, d + 1L))
  if (ncol(ev$vectors) < d + 1L) {
    .gproc_stop("invalid_problem", "LBW eigenproblem did not yield d shape directions.")
  }
  keep <- seq_len(d) + 1L
  U <- ev$vectors[, keep, drop = FALSE]
  eta <- ev$values[keep]
  Lambda <- .lbw_lambda(spec$reference_covariance, d)
  M <- U %*% diag(sqrt(Lambda), d, d)
  storage.mode(M) <- "double"
  transforms <- vector("list", length(nms))
  names(transforms) <- nms
  Y <- vector("list", length(nms))
  names(Y) <- nms
  pen <- 0
  for (nm in nms) {
    v <- views[[nm]]
    if (!is.null(anchor) && identical(nm, anchor)) {
      transforms[[nm]] <- proc_identity(spec, d)
      G <- matrix(NA_real_, n, d)
      G[v$map, ] <- apply_proc_transform(transforms[[nm]], v$X)
      Y[[nm]] <- G
      next
    }
    Mloc <- M[v$map, , drop = FALSE]
    B <- .gproc_psd_solve(v$A, crossprod(v$Phi, v$w * Mloc))
    pen <- pen + v$mu * sum(B * (v$L %*% B))
    tr <- .lbw_fitted(spec, B, v, d)
    transforms[[nm]] <- tr
    G <- matrix(NA_real_, n, d)
    G[v$map, ] <- v$Phi %*% B
    unobs <- v$map[!v$mask]
    if (length(unobs)) G[unobs, ] <- NA_real_
    Y[[nm]] <- G
  }
  obj_data <- 0
  for (nm in nms) {
    dlt <- Y[[nm]] - M
    dlt[!is.finite(dlt)] <- 0
    obj_data <- obj_data + alphas[[nm]] * sum(dlt^2)
  }
  gap <- ev$values[[d + 1L]] - ev$values[[d]]
  list(
    transforms = transforms,
    work_transforms = transforms,
    consensus = M,
    aligned = Y,
    objective = obj_data + pen,
    data_objective = obj_data,
    penalty = pen,
    history = data.frame(
      iteration = 1L,
      objective = obj_data + pen,
      relative_change = 0,
      eigen_gap = gap,
      elapsed_time = 0
    ),
    numerical_status = "converged",
    optimality_status = if (isTRUE(free_num)) "exact_closed_form" else "first_order_stationary",
    stationarity = sqrt(sum(p1^2)),
    free_translation = free_num,
    reference_covariance = Lambda,
    eigen_values = eta,
    eigen_gap = gap,
    folding = .lbw_folding(transforms, views),
    init = "lbw_eigen"
  )
}

#' @noRd
.lbw_lambda <- function(lam, d) {
  if (is.null(lam)) lam <- 1
  lam <- as.numeric(lam)
  if (length(lam) == 1L) lam <- rep(lam, d)
  if (length(lam) != d) {
    .gproc_stop("invalid_problem", "reference_covariance must be scalar or length d.")
  }
  if (any(!is.finite(lam) | lam <= 0)) {
    .gproc_stop("invalid_problem", "reference_covariance must be positive.")
  }
  lam
}

#' @noRd
.lbw_prepare_views <- function(data, spec, metric) {
  nms <- names(data$views)
  d <- ncol(data$views[[1L]])
  n <- length(data$global_ids)
  out <- vector("list", length(nms))
  names(out) <- nms
  landmark <- if (inherits(metric, "proc_metric")) metric$landmark else NULL
  for (nm in nms) {
    X <- as.matrix(data$views[[nm]])
    mask <- data$observed[[nm]]
    map <- data$row_map[[nm]]
    w <- as.numeric(mask)
    if (!is.null(landmark)) {
      src <- if (is.list(landmark)) landmark[[nm]] else landmark
      if (!is.null(src)) {
        src <- as.numeric(src)
        if (length(src) == 1L) {
          w <- w * src
        } else if (length(src) == length(map)) {
          w <- w * src
        } else if (length(src) == n) {
          w <- w * src[map]
        }
      }
    }
    built <- .lbw_basis(spec, X, mask)
    mu <- spec$smoothness %||% 0
    A <- crossprod(built$Phi, w * built$Phi) + mu * built$L
    out[[nm]] <- list(
      X = X,
      Phi = built$Phi,
      L = built$L,
      A = A,
      w = w,
      mask = mask,
      map = map,
      mu = mu,
      controls = built$controls,
      free_translation = built$free_translation,
      name = nm
    )
  }
  out
}

#' @noRd
.lbw_basis <- function(spec, X, mask = NULL) {
  d <- ncol(X)
  if (is.null(mask)) mask <- rep(TRUE, nrow(X))
  if (identical(spec$family, "affine") ||
      (identical(spec$family, "lbw") && identical(spec$basis, "affine"))) {
    Phi <- cbind(X, 1)
    L <- matrix(0, d + 1L, d + 1L)
    a <- c(rep(0, d), 1)
    free <- sqrt(sum((Phi %*% a - 1)^2)) < 1e-10 && sqrt(sum((L %*% a)^2)) < 1e-10
    return(list(Phi = Phi, L = L, controls = NULL, free_translation = free))
  }
  if (identical(spec$family, "tps")) {
    ctrl <- spec$control_points
    if (is.null(ctrl)) {
      ctrl <- X[mask, , drop = FALSE]
    } else {
      ctrl <- as.matrix(ctrl)
    }
    Phi <- .gproc_tps_phi(X, ctrl)
    L <- .gproc_tps_L(ctrl)
    a <- c(1, rep(0, ncol(Phi) - 1L))
    free <- sqrt(sum((Phi %*% a - 1)^2)) < 1e-8 && sqrt(sum((L %*% a)^2)) < 1e-8
    return(list(Phi = Phi, L = L, controls = ctrl, free_translation = free))
  }
  if (is.function(spec$basis)) {
    Phi <- as.matrix(spec$basis(X))
  } else {
    .gproc_stop("invalid_problem", "proc_lbw() basis must be a function or 'affine'.")
  }
  L <- spec$penalty
  if (is.function(L)) {
    L <- as.matrix(L(Phi))
  } else if (is.null(L)) {
    L <- matrix(0, ncol(Phi), ncol(Phi))
  } else {
    L <- as.matrix(L)
  }
  ones <- rep(1, nrow(Phi))
  a <- tryCatch(.gproc_psd_solve(crossprod(Phi), crossprod(Phi, ones)), error = function(e) NULL)
  free <- FALSE
  if (!is.null(a)) {
    free <- sqrt(sum((Phi %*% a - 1)^2)) < 1e-6 && sqrt(sum((L %*% a)^2)) < 1e-6
  }
  list(Phi = Phi, L = L, controls = NULL, free_translation = free)
}

#' @noRd
.gproc_tps_kernel <- function(r, d) {
  r <- pmax(as.numeric(r), 0)
  if (d == 2L) {
    out <- r * r * log(r)
    out[!is.finite(out) | r < 1e-15] <- 0
    out
  } else if (d == 3L) {
    r
  } else {
    .gproc_stop("invalid_problem", "Built-in TPS supports d = 2 or d = 3.")
  }
}

#' @noRd
.gproc_tps_phi <- function(X, controls) {
  X <- as.matrix(X)
  controls <- as.matrix(controls)
  n <- nrow(X)
  c <- nrow(controls)
  d <- ncol(X)
  U <- matrix(0, n, c)
  for (k in seq_len(c)) {
    delta <- sweep(X, 2L, controls[k, ], "-")
    r <- sqrt(rowSums(delta * delta))
    U[, k] <- .gproc_tps_kernel(r, d)
  }
  cbind(1, X, U)
}

#' @noRd
.gproc_tps_L <- function(controls) {
  controls <- as.matrix(controls)
  c <- nrow(controls)
  d <- ncol(controls)
  K <- matrix(0, c, c)
  for (i in seq_len(c)) {
    delta <- sweep(controls, 2L, controls[i, ], "-")
    r <- sqrt(rowSums(delta * delta))
    K[i, ] <- .gproc_tps_kernel(r, d)
  }
  K <- (K + t(K)) / 2
  L <- matrix(0, 1L + d + c, 1L + d + c)
  if (c) {
    idx <- (2L + d):(1L + d + c)
    L[idx, idx] <- K
  }
  L
}

#' @noRd
.lbw_form_P <- function(views, alphas, n) {
  P <- matrix(0, n, n)
  for (nm in names(views)) {
    v <- views[[nm]]
    a <- alphas[[nm]]
    Hi <- v$Phi %*% .gproc_psd_solve(v$A, t(v$Phi * v$w))
    obs <- v$map
    contrib <- diag(v$w, length(v$w)) - Hi
    P[obs, obs] <- P[obs, obs] + a * contrib
  }
  P
}

#' @noRd
.lbw_fitted <- function(spec, B, view, d) {
  if (identical(spec$family, "affine") ||
      (identical(spec$family, "lbw") && identical(spec$basis, "affine"))) {
    A <- B[seq_len(d), , drop = FALSE]
    t <- as.numeric(B[d + 1L, ])
    return(.proc_fitted(spec, R = A, s = 1, t = t, B = B))
  }
  controls <- view$controls
  apply_warp <- function(Xnew) {
    if (identical(spec$family, "tps")) {
      .gproc_tps_phi(Xnew, if (is.null(controls)) Xnew else controls) %*% B
    } else if (is.function(spec$basis)) {
      as.matrix(spec$basis(Xnew)) %*% B
    } else {
      cbind(Xnew, 1) %*% B
    }
  }
  .proc_fitted(
    spec,
    R = diag(d),
    s = 1,
    t = rep(0, d),
    B = B,
    apply_warp = apply_warp,
    extra = list(controls = controls)
  )
}

#' @noRd
.lbw_folding <- function(transforms, views) {
  dets <- list()
  for (nm in names(transforms)) {
    tr <- transforms[[nm]]
    if (identical(tr$spec$family, "affine") || is.null(tr$apply_warp)) {
      dets[[nm]] <- det(tr$R)
    } else {
      X <- views[[nm]]$X
      dets[[nm]] <- .lbw_numeric_jacobian_dets(tr, X)
    }
  }
  nonlinear <- vapply(names(transforms), function(nm) {
    is.function(transforms[[nm]]$apply_warp)
  }, logical(1))
  n_folded <- if (any(nonlinear)) {
    sum(mapply(function(z, nl) if (nl) sum(z < 0, na.rm = TRUE) else 0L, dets, nonlinear))
  } else {
    0L
  }
  list(
    jacobian_det = dets,
    n_folded = as.integer(n_folded),
    any_folded = n_folded > 0L
  )
}

#' @noRd
.lbw_numeric_jacobian_dets <- function(tr, X, eps = 1e-5) {
  d <- ncol(X)
  n <- nrow(X)
  out <- numeric(n)
  for (i in seq_len(n)) {
    J <- matrix(0, d, d)
    for (j in seq_len(d)) {
      e <- rep(0, d)
      e[[j]] <- eps
      xp <- X[i, , drop = FALSE] + e
      xm <- X[i, , drop = FALSE] - e
      yp <- apply_proc_transform(tr, xp)
      ym <- apply_proc_transform(tr, xm)
      J[, j] <- as.numeric((yp - ym) / (2 * eps))
    }
    out[[i]] <- det(J)
  }
  out
}
