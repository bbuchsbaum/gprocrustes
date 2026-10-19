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
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  masks <- weights <- list()
  raw <- metric$cell
  for (nm in names(data$views)) {
    map <- data$row_map[[nm]]
    allowed <- matrix(FALSE, n, d)
    allowed[map[data$observed[[nm]]], ] <- TRUE
    if (!is.null(data$cells[[nm]])) {
      C <- data$cells[[nm]]
      if (anyNA(C)) .gproc_stop("invalid_problem", "Cell masks cannot contain NA.")
      local <- matrix(FALSE, n, d)
      local[map, ] <- C
      allowed <- allowed & local
    }
    W <- matrix(1, n, d)
    src <- if (is.list(raw)) raw[[nm]] else raw
    if (!is.null(src)) {
      C <- as.matrix(src)
      if (ncol(C) != d || !nrow(C) %in% c(n, length(map)) ||
          any(!is.finite(C)) || any(C < 0)) {
        .gproc_stop("invalid_problem", "Cell weights must be finite nonnegative matrices matching the view or global entities.")
      }
      if (is.list(raw) && nrow(C) == length(map)) {
        W[] <- 0
        W[map, ] <- C
      } else if (nrow(C) == n) {
        W <- C * 1
      } else {
        W[] <- 0
        W[map, ] <- C
      }
    }
    W[!allowed] <- 0
    masks[[nm]] <- allowed & W > 0
    weights[[nm]] <- W
  }
  list(mask = masks, weight = weights)
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
