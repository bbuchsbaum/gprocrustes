#' Landmark-vector IRLS (Huber / Tukey) around Gower BCD.
#'
#' Inner steps are weighted GPA. The recorded objective is the true \(\rho\)
#' criterion. Weights are rolled back if that criterion rises.
#'
#' @noRd
.irls_gpa <- function(data, spec, gauge, metric, alphas, control, anchor, loss) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  w <- .gproc_initial_landmark_weights(data, metric)
  hist <- list()
  t0 <- proc.time()[["elapsed"]]
  best <- NULL
  status <- "maximum_iterations"
  inner <- control
  inner$max_iterations <- min(as.integer(control$max_iterations %||% 500L), 80L)
  inner$nstart <- 1L
  for (it in seq_len(control$max_iterations %||% 500L)) {
    met <- metric
    met$landmark <- w
    raw <- .gower_gpa(data, spec, gauge, met, alphas, inner, anchor)
    r <- .gproc_residual_norms(raw$aligned, raw$consensus)
    obj <- .gproc_rho_objective(r, alphas, loss)
    w_new <- lapply(names(r), function(nm) {
      wi <- .gproc_irls_weight(r[[nm]], loss)
      wi[!is.finite(rowSums(raw$aligned[[nm]]))] <- 0
      wi
    })
    names(w_new) <- names(r)
    wchg <- max(mapply(function(a, b) max(abs(a - b)), w, w_new))
    if (!is.null(best) && obj > best$true_objective + 1e-10) {
      break
    }
    best <- raw
    best$true_objective <- obj
    best$landmark_weights <- w_new
    hist[[it]] <- data.frame(
      iteration = it,
      objective = obj,
      relative_change = if (it == 1L) Inf else abs(hist[[it - 1L]]$objective - obj) / max(1, abs(obj)),
      weight_change = wchg,
      elapsed_time = proc.time()[["elapsed"]] - t0
    )
    w <- w_new
    if (it > 1L && hist[[it]]$relative_change <= (control$tolerance %||% 1e-8) &&
        wchg <= 1e-6) {
      status <- "converged"
      break
    }
  }
  best$history <- if (length(hist)) do.call(rbind, hist) else data.frame()
  best$objective <- best$true_objective
  best$numerical_status <- status
  best$optimality_status <- if (identical(status, "converged")) "blockwise_stationary" else "not_converged"
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
  cells <- .gproc_cell_masks(data, metric)
  covar <- .gproc_as_covariance(metric$precision)
  inner <- control
  inner$max_iterations <- min(as.integer(control$max_iterations %||% 500L), 40L)
  inner$nstart <- 1L
  raw <- .gower_gpa(data, spec, gauge, metric, alphas, inner, anchor)
  transforms <- raw$transforms
  Y <- .gproc_apply_all(data, transforms, n, d)
  M <- .gproc_cell_consensus(Y, alphas, cells, n, d)
  obj <- .gproc_observed_objective(Y, M, alphas, cells, covar, loss)
  hist <- list()
  t0 <- proc.time()[["elapsed"]]
  status <- "maximum_iterations"
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
      pair <- procrustes(
        X[row_keep, , drop = FALSE],
        Targ[map[row_keep], , drop = FALSE],
        transform = spec,
        weights = w_loc
      )
      Tnew[[nm]] <- pair$transform
    }
    Y_try <- .gproc_apply_all(data, Tnew, n, d)
    M_try <- .gproc_cell_consensus(Y_try, alphas, cells, n, d)
    obj_try <- .gproc_observed_objective(Y_try, M_try, alphas, cells, covar, loss)
    accepted <- obj_try <= obj + 1e-12
    if (accepted) {
      transforms <- Tnew
      Y <- Y_try
      M <- M_try
      obj <- obj_try
    }
    rel <- abs((if (length(hist)) hist[[length(hist)]]$objective else obj) - obj) / max(1, abs(obj))
    hist[[it]] <- data.frame(
      iteration = it,
      objective = obj,
      relative_change = rel,
      accepted_mm = accepted,
      elapsed_time = proc.time()[["elapsed"]] - t0
    )
    if (it > 1L && rel <= (control$tolerance %||% 1e-8)) {
      status <- "converged"
      break
    }
    if (!accepted) {
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
  raw$optimality_status <- if (identical(status, "converged")) "first_order_stationary" else "not_converged"
  raw$cell_masks <- cells
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
.gproc_cell_consensus <- function(Y, alphas, cells, n, d) {
  num <- matrix(0, n, d)
  den <- matrix(0, n, d)
  for (nm in names(Y)) {
    Yi <- Y[[nm]]
    a <- alphas[[nm]]
    keep <- cells[[nm]] & is.finite(Yi)
    num[keep] <- num[keep] + a * Yi[keep]
    den[keep] <- den[keep] + a
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
.gproc_rho_objective <- function(r, alphas, loss) {
  f <- 0
  for (nm in names(r)) {
    f <- f + alphas[[nm]] * sum(.gproc_rho(r[[nm]], loss))
  }
  f
}

#' @noRd
.gproc_cell_masks <- function(data, metric) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  out <- vector("list", length(nms))
  names(out) <- nms
  raw <- metric$cell
  for (nm in nms) {
    G <- matrix(FALSE, n, d)
    loc <- data$observed[[nm]]
    G[data$row_map[[nm]][loc], ] <- TRUE
    if (!is.null(data$cells) && !is.null(data$cells[[nm]])) {
      C <- data$cells[[nm]]
      if (!is.matrix(C)) C <- as.matrix(C)
      storage.mode(C) <- "logical"
      G[] <- FALSE
      G[data$row_map[[nm]], ] <- C
    } else if (is.list(raw) && !is.null(raw[[nm]])) {
      C <- raw[[nm]]
      if (!is.matrix(C)) C <- as.matrix(C)
      storage.mode(C) <- "logical"
      G[] <- FALSE
      if (nrow(C) == n) {
        G <- C
      } else {
        G[data$row_map[[nm]], ] <- C
      }
    }
    out[[nm]] <- G
  }
  out
}

#' @noRd
.gproc_observed_objective <- function(Y, M, alphas, cells, covar, loss) {
  Wd <- if (!is.null(covar) && !is.null(covar$coordinate)) {
    A <- covar$coordinate
    if (isTRUE(covar$as_precision)) A else {
      ev <- .gproc_eigen_sym(A)
      ev$vectors %*% (ifelse(ev$values > 1e-12, 1 / ev$values, 0) * t(ev$vectors))
    }
  } else {
    NULL
  }
  f <- 0
  for (nm in names(Y)) {
    dlt <- Y[[nm]] - M
    if (!is.null(Wd)) {
      dlt <- dlt %*% Wd
    }
    keep <- cells[[nm]]
    dlt[!keep] <- 0
    r <- sqrt(rowSums(dlt * dlt))
    f <- f + alphas[[nm]] * sum(.gproc_rho(r, if (is.null(loss)) proc_squared_l2() else loss))
  }
  f
}
