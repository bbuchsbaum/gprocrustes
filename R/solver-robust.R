#' Landmark-vector IRLS (Huber / Tukey) around Gower BCD.
#'
#' Inner steps are weighted GPA. The recorded objective is the true \eqn{\rho}
#' criterion. Weights are rolled back if that criterion rises.
#'
#' @noRd
.irls_gpa <- function(data, spec, gauge, metric, alphas, control, anchor, loss) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  w0 <- .gproc_initial_landmark_weights(data, metric)
  w <- w0
  hist <- list()
  t0 <- proc.time()[["elapsed"]]
  best <- NULL
  status <- "maximum_iterations"
  inner <- control
  inner$max_iterations <- min(as.integer(control$max_iterations %||% 500L), 80L)
  inner$nstart <- 1L
  tol <- control$tolerance %||% 1e-8
  for (it in seq_len(control$max_iterations %||% 500L)) {
    met <- metric
    met$landmark <- w
    raw <- .gower_gpa(data, spec, gauge, met, alphas, inner, anchor)
    r <- .gproc_residual_norms(raw$aligned, raw$consensus)
    obj <- .gproc_rho_objective(r, alphas, loss, w0)
    q <- lapply(names(r), function(nm) {
      qi <- .gproc_irls_weight(r[[nm]], loss)
      qi[!is.finite(rowSums(raw$aligned[[nm]]))] <- 0
      qi
    })
    names(q) <- names(r)
    w_new <- lapply(names(r), function(nm) q[[nm]] * w0[[nm]])
    names(w_new) <- names(r)
    wchg <- max(mapply(function(a, b) max(abs(a - b)), w, w_new))
    if (!is.null(best) && obj > best$true_objective + 1e-10) {
      break
    }
    best <- raw
    best$true_objective <- obj
    best$landmark_weights <- w_new
    best$irls_q <- q
    hist[[it]] <- data.frame(
      iteration = it,
      objective = obj,
      relative_change = if (it == 1L) Inf else abs(hist[[it - 1L]]$objective - obj) / max(1, abs(obj)),
      weight_change = wchg,
      stationarity = raw$stationarity %||% NA_real_,
      elapsed_time = proc.time()[["elapsed"]] - t0
    )
    w <- w_new
    inner_stat <- raw$stationarity %||% Inf
    if (it > 1L && hist[[it]]$relative_change <= tol &&
        wchg <= 1e-6 && is.finite(inner_stat) && inner_stat <= sqrt(tol)) {
      status <- "converged"
      break
    }
  }
  best$history <- if (length(hist)) do.call(rbind, hist) else data.frame()
  best$objective <- best$true_objective
  best$numerical_status <- status
  best$effective_n <- lapply(w, function(wi) {
    s <- sum(wi)
    if (!(s > 0)) 0 else s^2 / sum(wi^2)
  })
  inner_stat <- best$stationarity %||% Inf
  best$optimality_status <- if (identical(status, "converged") &&
      is.finite(inner_stat) && inner_stat <= sqrt(tol)) {
    "blockwise_stationary"
  } else {
    "not_converged"
  }
  best$loss_family <- loss$family
  best
}

#' Cell-mask / anisotropic-coordinate majorization.
#'
#' Fills unobserved target coordinates with the current consensus, solves a
#' complete weighted GPA, and accepts the step only if the true observed
#' objective decreases.
#'
#' @noRd
.mm_gpa <- function(data, spec, gauge, metric, alphas, control, anchor, loss) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  chan <- .gproc_cell_channel(data, metric)
  cells <- chan$mask
  cweight <- chan$weight
  covar <- .gproc_fitting_covariance(metric)
  Q <- .gproc_coordinate_precision(covar)
  inner <- control
  inner$max_iterations <- min(as.integer(control$max_iterations %||% 500L), 40L)
  inner$nstart <- 1L
  raw <- .gower_gpa(data, spec, gauge, metric, alphas, inner, anchor)
  transforms <- raw$transforms
  Y <- .gproc_apply_all(data, transforms, n, d)
  M <- .gproc_cell_consensus(Y, alphas, cells, n, d, cweight)
  obj <- .gproc_observed_objective(Y, M, alphas, cells, Q, loss, cweight)
  hist <- list()
  t0 <- proc.time()[["elapsed"]]
  status <- "maximum_iterations"
  tol <- control$tolerance %||% 1e-8
  for (it in seq_len(control$max_iterations %||% 500L)) {
    M_surr <- M
    Tnew <- transforms
    for (nm in nms) {
      if (!is.null(anchor) && identical(nm, anchor)) next
      Targ <- M_surr
      miss <- !cells[[nm]]
      if (any(miss)) {
        fill <- Y[[nm]]
        use <- miss & is.finite(fill)
        Targ[use] <- fill[use]
      }
      map <- data$row_map[[nm]]
      X <- as.matrix(data$views[[nm]])
      row_keep <- apply(cells[[nm]][map, , drop = FALSE], 1L, any)
      if (!any(row_keep)) next
      w_loc <- .gproc_local_landmark_weights(metric$landmark, map, n, nm)
      if (!is.null(w_loc)) w_loc <- w_loc[row_keep]
      pair <- .gproc_mm_block(
        X[row_keep, , drop = FALSE],
        Targ[map[row_keep], , drop = FALSE],
        spec,
        w_loc,
        Q
      )
      Tnew[[nm]] <- pair$transform
    }
    Y_try <- .gproc_apply_all(data, Tnew, n, d)
    M_try <- .gproc_cell_consensus(Y_try, alphas, cells, n, d, cweight)
    obj_try <- .gproc_observed_objective(Y_try, M_try, alphas, cells, Q, loss, cweight)
    accepted <- obj_try <= obj + 1e-12
    if (accepted) {
      transforms <- Tnew
      Y <- Y_try
      M <- M_try
      obj <- obj_try
    }
    rel <- abs((if (length(hist)) hist[[length(hist)]]$objective else obj) - obj) / max(1, abs(obj))
    stat <- .gproc_mm_stationarity(data, transforms, M, cells, cweight, Q, spec, alphas)
    hist[[it]] <- data.frame(
      iteration = it,
      objective = obj,
      relative_change = rel,
      accepted_mm = accepted,
      stationarity = stat,
      elapsed_time = proc.time()[["elapsed"]] - t0
    )
    if (!accepted) {
      status <- "stalled"
      break
    }
    if (it > 1L && rel <= tol) {
      status <- "converged"
      break
    }
  }
  raw$transforms <- transforms
  raw$work_transforms <- transforms
  raw$aligned <- Y
  raw$consensus <- M
  raw$objective <- obj
  raw$history <- if (length(hist)) do.call(rbind, hist) else data.frame()
  raw$numerical_status <- status
  stat_final <- if (nrow(raw$history)) raw$history$stationarity[[nrow(raw$history)]] else Inf
  raw$stationarity <- stat_final
  raw$optimality_status <- if (identical(status, "converged")) {
    "first_order_stationary"
  } else {
    "not_converged"
  }
  raw$cell_masks <- cells
  raw$cell_weights <- cweight
  raw
}

#' @noRd
.gproc_local_landmark_weights <- function(landmark, map, n, view_name = NULL) {
  if (is.null(landmark)) {
    return(NULL)
  }
  src <- if (is.list(landmark)) {
    if (!is.null(view_name) && !is.null(landmark[[view_name]])) landmark[[view_name]] else NULL
  } else {
    landmark
  }
  if (is.null(src)) {
    return(NULL)
  }
  src <- as.numeric(src)
  if (length(src) == 1L) {
    return(rep(src, length(map)))
  }
  if (length(src) == length(map)) {
    return(src)
  }
  if (length(src) == n) {
    return(src[map])
  }
  NULL
}

#' @noRd
.gproc_apply_all <- function(data, transforms, n, d) {
  out <- vector("list", length(data$views))
  names(out) <- names(data$views)
  for (nm in names(data$views)) {
    Yi <- apply_proc_transform(transforms[[nm]], data$views[[nm]])
    G <- matrix(NA_real_, n, d)
    G[data$row_map[[nm]], ] <- Yi
    unobs <- data$row_map[[nm]][!data$observed[[nm]]]
    if (length(unobs)) G[unobs, ] <- NA_real_
    out[[nm]] <- G
  }
  out
}

#' @noRd
.gproc_cell_consensus <- function(Y, alphas, cells, n, d, weights = NULL) {
  num <- matrix(0, n, d)
  den <- matrix(0, n, d)
  for (nm in names(Y)) {
    Yi <- Y[[nm]]
    a <- alphas[[nm]]
    keep <- cells[[nm]] & is.finite(Yi)
    ww <- a
    if (!is.null(weights) && !is.null(weights[[nm]])) {
      W <- weights[[nm]]
      ww <- a * W
      keep <- keep & is.finite(W) & W > 0
    }
    num[keep] <- num[keep] + ww[keep] * Yi[keep]
    den[keep] <- den[keep] + ww[keep]
  }
  M <- num
  ok <- den > 0
  M[ok] <- num[ok] / den[ok]
  M[!ok] <- 0
  M
}

#' @noRd
.gproc_initial_landmark_weights <- function(data, metric) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  out <- vector("list", length(nms))
  names(out) <- nms
  for (nm in nms) {
    w <- rep(0, n)
    w[data$row_map[[nm]][data$observed[[nm]]]] <- 1
    if (!is.null(metric$landmark)) {
      src <- if (is.list(metric$landmark)) metric$landmark[[nm]] else metric$landmark
      if (!is.null(src)) {
        src <- as.numeric(src)
        if (length(src) == 1L) {
          w <- w * src
        } else if (length(src) == n) {
          w <- w * src
        } else if (length(src) == length(data$row_map[[nm]])) {
          w[data$row_map[[nm]]] <- w[data$row_map[[nm]]] * src
        }
      }
    }
    out[[nm]] <- w
  }
  out
}

#' @noRd
.gproc_residual_norms <- function(Y, M) {
  out <- vector("list", length(Y))
  names(out) <- names(Y)
  for (nm in names(Y)) {
    dlt <- Y[[nm]] - M
    r <- sqrt(rowSums(dlt * dlt))
    r[!is.finite(rowSums(Y[[nm]]))] <- 0
    out[[nm]] <- r
  }
  out
}

#' @noRd
.gproc_irls_weight <- function(r, loss) {
  r <- pmax(as.numeric(r), 0)
  if (identical(loss$family, "huber")) {
    k <- loss$k %||% 1.345
    w <- ifelse(r <= k, 1, k / pmax(r, .Machine$double.eps))
  } else if (identical(loss$family, "tukey")) {
    c <- loss$c %||% 4.685
    u <- r / c
    w <- ifelse(r < c, (1 - u^2)^2, 0)
  } else {
    w <- rep(1, length(r))
  }
  w
}

#' @noRd
.gproc_rho <- function(r, loss) {
  r <- pmax(as.numeric(r), 0)
  if (identical(loss$family, "huber")) {
    k <- loss$k %||% 1.345
    ifelse(r <= k, 0.5 * r^2, k * r - 0.5 * k^2)
  } else if (identical(loss$family, "tukey")) {
    c <- loss$c %||% 4.685
    u <- r / c
    ifelse(r < c, (c^2 / 6) * (1 - (1 - u^2)^3), c^2 / 6)
  } else {
    0.5 * r^2
  }
}

#' @noRd
.gproc_rho_objective <- function(r, alphas, loss, landmark = NULL) {
  f <- 0
  for (nm in names(r)) {
    rho <- .gproc_rho(r[[nm]], loss)
    w <- if (!is.null(landmark) && !is.null(landmark[[nm]])) landmark[[nm]] else 1
    f <- f + alphas[[nm]] * sum(w * rho)
  }
  f
}

#' @noRd
.gproc_cell_masks <- function(data, metric) {
  .gproc_cell_channel(data, metric)$mask
}

#' Logical masks and numeric cell weights. A matrix in `metric$cell` is a weight.
#'
#' @noRd
.gproc_cell_channel <- function(data, metric) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  mask <- vector("list", length(nms))
  weight <- vector("list", length(nms))
  names(mask) <- nms
  names(weight) <- nms
  raw <- if (inherits(metric, "proc_metric")) metric$cell else NULL
  for (nm in nms) {
    G <- matrix(FALSE, n, d)
    W <- matrix(0, n, d)
    loc <- data$observed[[nm]]
    G[data$row_map[[nm]][loc], ] <- TRUE
    W[G] <- 1
    src <- NULL
    if (!is.null(data$cells) && !is.null(data$cells[[nm]])) {
      src <- data$cells[[nm]]
    } else if (is.list(raw) && !is.null(raw[[nm]])) {
      src <- raw[[nm]]
    } else if (is.matrix(raw) && is.null(names(raw))) {
      src <- raw
    }
    if (!is.null(src)) {
      C <- as.matrix(src)
      if (is.logical(C)) {
        G[] <- FALSE
        if (nrow(C) == n) G <- C else G[data$row_map[[nm]], ] <- C
        W <- matrix(0, n, d)
        W[G] <- 1
      } else {
        storage.mode(C) <- "double"
        if (nrow(C) == n) {
          W <- C
        } else {
          W[] <- 0
          W[data$row_map[[nm]], ] <- C
        }
        G <- is.finite(W) & W > 0
        W[!G] <- 0
      }
    }
    mask[[nm]] <- G
    weight[[nm]] <- W
  }
  list(mask = mask, weight = weight)
}

#' @noRd
.gproc_observed_objective <- function(Y, M, alphas, cells, Q, loss, weights = NULL) {
  f <- 0
  for (nm in names(Y)) {
    dlt <- Y[[nm]] - M
    keep <- cells[[nm]]
    dlt[!keep] <- 0
    if (is.null(Q) && !is.null(weights) && !is.null(weights[[nm]])) {
      r2 <- rowSums(weights[[nm]] * dlt * dlt)
    } else {
      r2 <- .gproc_row_quadratic(dlt, Q)
      if (!is.null(weights) && !is.null(weights[[nm]])) {
        r2 <- r2 * rowMeans(weights[[nm]])
      }
    }
    r <- sqrt(pmax(r2, 0))
    f <- f + alphas[[nm]] * sum(.gproc_rho(r, if (is.null(loss)) proc_squared_l2() else loss))
  }
  f
}

#' Row-wise e^T Q e. If Q is NULL this is the Euclidean squared length.
#'
#' @noRd
.gproc_row_quadratic <- function(E, Q) {
  E <- as.matrix(E)
  if (is.null(Q)) {
    return(rowSums(E * E))
  }
  rowSums((E %*% Q) * E)
}

#' @noRd
.gproc_coordinate_precision <- function(covar) {
  if (is.null(covar) || is.null(covar$coordinate)) {
    return(NULL)
  }
  A <- covar$coordinate
  if (isTRUE(covar$as_precision)) {
    return(A)
  }
  ev <- .gproc_eigen_sym(A)
  ev$vectors %*% (ifelse(ev$values > 1e-12, 1 / ev$values, 0) * t(ev$vectors))
}

#' Pairwise block update for the MM surrogate, including anisotropic Q.
#'
#' @noRd
.gproc_mm_block <- function(X, Y, spec, w, Q) {
  pair <- procrustes(X, Y, transform = spec, weights = w)
  if (is.null(Q) || !.gproc_coordinate_anisotropic(list(coordinate = Q, as_precision = TRUE))) {
    return(pair)
  }
  .gproc_anisotropic_refine(X, Y, spec, w, Q, pair)
}

#' Descent on O(d)/SO(d) for the e^T Q e pairwise residual.
#'
#' @noRd
.gproc_anisotropic_refine <- function(X, Y, spec, w, Q, pair) {
  tr <- pair$transform
  R <- tr$R
  s <- tr$s
  t <- tr$t
  obj <- function(R, s, t) {
    E <- sweep(s * X %*% R, 2L, t, "+") - Y
    ww <- if (is.null(w)) 1 else w
    sum(ww * .gproc_row_quadratic(E, Q))
  }
  f <- obj(R, s, t)
  moms <- proc_moments(X, Y, w, center = .gproc_eliminates_translation(spec))
  A <- if (.gproc_eliminates_translation(spec)) {
    centered_crossprod(X, X, w)
  } else {
    .gproc_raw_crossprod(X, X, w)
  }
  C <- moms$C
  for (it in seq_len(25L)) {
    G <- (s^2) * A %*% R %*% Q - s * C %*% Q
    RtG <- crossprod(R, G)
    skew <- 0.5 * (RtG - t(RtG))
    R_try <- .polar_factor(R - skew, spec$group)$R
    if (.gproc_eliminates_translation(spec)) {
      t_try <- as.numeric(moms$ybar - s * (moms$xbar %*% R_try))
    } else {
      t_try <- t * 0
    }
    if (identical(spec$scaling, "isotropic")) {
      aQ <- sum(diag(crossprod(R_try, A %*% R_try %*% Q)))
      gQ <- sum(R_try * (C %*% Q))
      s_try <- if (aQ > 0) max(gQ / aQ, 0) else 0
      if (.gproc_eliminates_translation(spec)) {
        t_try <- as.numeric(moms$ybar - s_try * (moms$xbar %*% R_try))
      }
    } else {
      s_try <- s
    }
    f_try <- obj(R_try, s_try, t_try)
    if (f_try <= f + 1e-14) {
      R <- R_try
      s <- s_try
      t <- t_try
      f <- f_try
    } else {
      break
    }
  }
  tr <- .proc_fitted(spec, R, s, t)
  pair$transform <- tr
  pair$objective <- .gproc_transform_residual(X, Y, tr, if (is.null(w)) {
    rep(1, nrow(X))
  } else {
    w
  })
  pair
}

#' First-order residual of the observed-cell / anisotropic objective.
#'
#' @noRd
.gproc_mm_stationarity <- function(data, transforms, M, cells, weights, Q, spec, alphas) {
  rmax <- 0
  nms <- names(data$views)
  for (nm in nms) {
    map <- data$row_map[[nm]]
    X <- as.matrix(data$views[[nm]])
    Y <- apply_proc_transform(transforms[[nm]], X)
    E <- Y - M[map, , drop = FALSE]
    keep <- cells[[nm]][map, , drop = FALSE]
    E[!keep] <- 0
    if (!is.null(weights) && !is.null(weights[[nm]])) {
      E <- E * weights[[nm]][map, , drop = FALSE]
    }
    if (!is.null(Q)) {
      E <- E %*% Q
    }
    if (.gproc_eliminates_translation(spec)) {
      w <- rep(1, nrow(X))
      Xc <- sweep(X, 2L, weighted_centroid(X, w), "-")
      Ec <- sweep(E, 2L, weighted_centroid(E, w), "-")
      G <- crossprod(Xc, Ec)
    } else {
      G <- crossprod(X, E)
    }
    R <- transforms[[nm]]$R
    S <- crossprod(R, G)
    skew <- 0.5 * (S - t(S))
    rmax <- max(rmax, sqrt(sum(skew^2)) * (alphas[[nm]] %||% 1))
  }
  rmax
}
