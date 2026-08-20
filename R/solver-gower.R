#' Generalized similarity / orthogonal GPA (Gower block-coordinate descent).
#'
#' Consensus-first BCD. Scale conventions follow math spec section 9.
#' The only allowed optimality claim is `blockwise_stationary`.
#'
#' @noRd
.gower_gpa <- function(data, spec, gauge, metric, alphas, control, anchor = NULL) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  if (any(vapply(data$views, ncol, integer(1)) != d)) {
    .gproc_stop("dimension_mismatch", "Gower GPA requires a common target dimension.")
  }
  if (!is.null(anchor) && !anchor %in% nms) {
    .gproc_stop("invalid_problem", "Unknown anchor view.")
  }

  scale_info <- .gower_scale_mode(spec, gauge)
  views <- .gower_prepare_views(data, spec, scale_info$mode)
  views <- .gower_attach_landmark_weights(views, data, metric)
  spec_step <- spec
  spec_step$translation <- isTRUE(spec$translation) && !identical(scale_info$mode, "preshape")
  if (scale_info$mode %in% c("preshape", "gower")) {
    spec_step$scaling <- "none"
    if (identical(scale_info$mode, "preshape")) {
      spec_step$family <- "orthogonal"
    }
  }

  starts <- .gower_start_names(control)
  best <- NULL
  for (init in starts) {
    cand <- .gower_run(
      views = views,
      spec_step = spec_step,
      spec = spec,
      scale_info = scale_info,
      alphas = alphas,
      n = n,
      d = d,
      nms = nms,
      control = control,
      init = init,
      anchor = anchor
    )
    if (is.null(best) || cand$objective < best$objective) {
      best <- cand
    }
  }
  best$transforms <- .gower_original_transforms(best$transforms, views, spec, d)
  Y <- .gower_aligned(views, best$work_transforms, n, d)
  best$aligned <- Y
  best$objective <- .gower_objective(Y, best$consensus, alphas, views)
  best$landmark_weights <- lapply(views, `[[`, "w")
  best
}

#' @keywords internal
.gower_scale_mode <- function(spec, gauge) {
  gscale <- if (inherits(gauge, "proc_gauge")) gauge$scale else "none"
  isotropic <- identical(spec$scaling, "isotropic")
  if (identical(gscale, "preshape")) {
    return(list(mode = "preshape", overridden = FALSE,
                note = "Each centered configuration is unit-Frobenius; rotations only."))
  }
  if (!isotropic) {
    return(list(mode = "none", overridden = FALSE, note = NULL))
  }
  if (identical(gscale, "fixed_consensus")) {
    return(list(mode = "fixed_consensus", overridden = FALSE,
                note = "Consensus Frobenius norm is held at 1."))
  }
  overridden <- identical(gscale, "none")
  list(
    mode = "gower",
    overridden = overridden,
    note = if (overridden) {
      "Isotropic scales with gauge='none' would collapse; using the Gower energy constraint."
    } else {
      "Ten Berge / Gower collective scaling with a global energy constraint."
    }
  )
}

#' @keywords internal
.gower_start_names <- function(control) {
  init <- control$init %||% "sequential"
  init <- .gproc_match_arg(init, c("sequential", "medoid", "spectral"), "init")
  nstart <- as.integer(control$nstart %||% 1L)
  if (nstart < 1L) nstart <- 1L
  starts <- init
  if (nstart >= 2L) {
    starts <- unique(c(init, "medoid", "sequential"))
  }
  if (nstart >= 3L) {
    starts <- unique(c(starts, "spectral"))
  }
  starts
}

#' @keywords internal
.gower_prepare_views <- function(data, spec, scale_mode) {
  nms <- names(data$views)
  out <- vector("list", length(nms))
  names(out) <- nms
  d <- ncol(data$views[[1L]])
  for (nm in nms) {
    X <- as.matrix(data$views[[nm]])
    mask <- data$observed[[nm]]
    map <- data$row_map[[nm]]
    pre_s <- 1
    pre_t <- rep(0, d)
    X_work <- X
    if (identical(scale_mode, "preshape")) {
      xbar <- weighted_centroid(X[mask, , drop = FALSE], w = NULL)
      X_work <- sweep(X, 2L, xbar, "-")
      pre_t <- -xbar
    }
    if (identical(scale_mode, "preshape")) {
      nrm <- sqrt(sum(X_work[mask, , drop = FALSE]^2))
      if (nrm > 0) {
        X_work <- X_work / nrm
        pre_s <- 1 / nrm
        pre_t <- pre_t / nrm
      }
    }
    energy <- sum(X_work[mask, , drop = FALSE]^2)
    out[[nm]] <- list(
      X = X_work,
      X_orig = X,
      mask = mask,
      map = map,
      name = nm,
      pre_s = pre_s,
      pre_t = pre_t,
      energy = energy,
      w = as.numeric(mask)
    )
  }
  out
}

#' @noRd
.gower_attach_landmark_weights <- function(views, data, metric) {
  nms <- names(views)
  n <- length(data$global_ids)
  raw <- if (inherits(metric, "proc_metric")) metric$landmark else NULL
  for (nm in nms) {
    v <- views[[nm]]
    w <- as.numeric(v$mask)
    if (is.null(raw)) {
      views[[nm]]$w <- w
      next
    }
    src <- if (is.list(raw)) raw[[nm]] else raw
    if (is.null(src)) {
      views[[nm]]$w <- w
      next
    }
    src <- as.numeric(src)
    if (length(src) == 1L) {
      w <- w * src
    } else if (length(src) == length(v$map)) {
      w <- w * src
    } else if (length(src) == n) {
      w <- w * src[v$map]
    } else {
      .gproc_stop("invalid_problem", sprintf("Landmark weights for '%s' have the wrong length.", nm))
    }
    views[[nm]]$w <- w
  }
  views
}

#' @keywords internal
.gower_original_transforms <- function(work, views, spec, d) {
  out <- work
  for (nm in names(work)) {
    pre <- .proc_fitted(spec, diag(d), views[[nm]]$pre_s, views[[nm]]$pre_t)
    tr <- compose_proc_transform(pre, work[[nm]])
    tr$spec <- spec
    out[[nm]] <- tr
  }
  out
}

#' @keywords internal
.gower_run <- function(views, spec_step, spec, scale_info, alphas, n, d, nms,
                       control, init, anchor) {
  t0 <- proc.time()[["elapsed"]]
  state <- .gower_initialize(views, spec_step, alphas, n, d, nms, init, anchor)
  Y <- .gower_aligned(views, state$transforms, n, d)
  M <- .gower_consensus(Y, alphas, n, d, views)
  state <- .gower_apply_scale(state, views, Y, M, alphas, scale_info$mode, n, d)
  Y <- .gower_aligned(views, state$transforms, n, d)
  M <- .gower_consensus(Y, alphas, n, d, views)
  state$transforms <- .gower_refit_translations(views, M, state$transforms, spec_step)
  Y <- .gower_aligned(views, state$transforms, n, d)
  M <- .gower_consensus(Y, alphas, n, d, views)
  obj <- .gower_objective(Y, M, alphas, views)
  hist <- list()
  M_prev <- NULL
  maxit <- control$max_iterations %||% 500L
  tol <- control$tolerance %||% 1e-8
  accelerate <- isTRUE(control$accelerate %||% TRUE)
  status <- "maximum_iterations"

  for (it in seq_len(maxit)) {
    old_obj <- obj
    old_M <- M
    old_state <- state

    state$transforms <- .gower_fit_transforms(views, M, spec_step, state$transforms, anchor)
    Y <- .gower_aligned(views, state$transforms, n, d)
    M <- .gower_consensus(Y, alphas, n, d, views)
    state <- .gower_apply_scale(state, views, Y, M, alphas, scale_info$mode, n, d)
    Y <- .gower_aligned(views, state$transforms, n, d)
    M <- .gower_consensus(Y, alphas, n, d, views)
    state$transforms <- .gower_refit_translations(views, M, state$transforms, spec_step)
    Y <- .gower_aligned(views, state$transforms, n, d)
    M <- .gower_consensus(Y, alphas, n, d, views)
    obj <- .gower_objective(Y, M, alphas, views)
    accepted_acc <- FALSE

    if (accelerate && !is.null(M_prev) && obj <= old_obj + 1e-12) {
      extra <- M + (M - old_M)
      state_try <- state
      state_try$transforms <- .gower_fit_transforms(views, extra, spec_step, state$transforms, anchor)
      Y_try <- .gower_aligned(views, state_try$transforms, n, d)
      M_try <- .gower_consensus(Y_try, alphas, n, d, views)
      state_try <- .gower_apply_scale(state_try, views, Y_try, M_try, alphas, scale_info$mode, n, d)
      Y_try <- .gower_aligned(views, state_try$transforms, n, d)
      M_try <- .gower_consensus(Y_try, alphas, n, d, views)
      state_try$transforms <- .gower_refit_translations(views, M_try, state_try$transforms, spec_step)
      Y_try <- .gower_aligned(views, state_try$transforms, n, d)
      M_try <- .gower_consensus(Y_try, alphas, n, d, views)
      obj_try <- .gower_objective(Y_try, M_try, alphas, views)
      if (obj_try <= obj + 1e-12) {
        state <- state_try
        Y <- Y_try
        M <- M_try
        obj <- obj_try
        accepted_acc <- TRUE
      }
    }

    if (obj > old_obj + 1e-10) {
      state <- old_state
      M <- old_M
      Y <- .gower_aligned(views, state$transforms, n, d)
      obj <- old_obj
      accepted_acc <- FALSE
    }

    rel <- abs(old_obj - obj) / max(1, abs(old_obj))
    cchg <- .gower_consensus_change(old_M, M)
    stat <- .gower_stationarity(views, M, state$transforms, spec_step)
    orth <- .gower_orthogonality(state$transforms)
    hist[[it]] <- data.frame(
      iteration = it,
      objective = obj,
      relative_change = rel,
      consensus_change = cchg,
      stationarity = stat,
      orthogonality_residual = orth,
      elapsed_time = proc.time()[["elapsed"]] - t0,
      accepted_acceleration = accepted_acc
    )
    M_prev <- old_M
    if (rel <= tol && cchg <= sqrt(tol) && stat <= sqrt(tol)) {
      status <- "converged"
      break
    }
  }

  hist_df <- if (length(hist)) do.call(rbind, hist) else data.frame()
  list(
    transforms = state$transforms,
    work_transforms = state$transforms,
    consensus = M,
    aligned = Y,
    objective = obj,
    history = hist_df,
    numerical_status = status,
    optimality_status = {
      stat_final <- if (nrow(hist_df)) hist_df$stationarity[[nrow(hist_df)]] else Inf
      if (identical(status, "converged") && is.finite(stat_final) && stat_final <= sqrt(tol)) {
        "blockwise_stationary"
      } else {
        "not_converged"
      }
    },
    init = init,
    scale_mode = scale_info$mode,
    scale_note = scale_info$note,
    scale_overridden = isTRUE(scale_info$overridden),
    stationarity = if (nrow(hist_df)) hist_df$stationarity[[nrow(hist_df)]] else NA_real_
  )
}

#' @keywords internal
.gower_initialize <- function(views, spec_step, alphas, n, d, nms, init, anchor) {
  if (!is.null(anchor)) {
    transforms <- .gower_align_to_view(views, spec_step, anchor, n)
    return(list(transforms = transforms))
  }
  if (identical(init, "medoid")) {
    med <- .gower_medoid(views, spec_step)
    return(list(transforms = .gower_align_to_view(views, spec_step, med, n)))
  }
  if (identical(init, "spectral")) {
    spec_ok <- all(vapply(views, function(v) all(v$mask) && length(v$map) == n, logical(1)))
    if (spec_ok) {
      return(list(transforms = .gower_spectral_init(views, spec_step, d)))
    }
  }
  list(transforms = .gower_sequential_init(views, spec_step, alphas, n, d))
}

#' @keywords internal
.gower_sequential_init <- function(views, spec_step, alphas, n, d) {
  nms <- names(views)
  transforms <- lapply(views, function(v) proc_identity(spec_step, d))
  names(transforms) <- nms
  Y <- .gower_aligned(views, transforms, n, d)
  M <- .gower_consensus(Y, alphas, n, d, views)
  for (nm in nms[-1L]) {
    transforms[[nm]] <- .gower_fit_one(views[[nm]], M, spec_step, transforms[[nm]])
    Y <- .gower_aligned(views, transforms, n, d)
    M <- .gower_consensus(Y, alphas, n, d, views)
  }
  transforms
}

#' @keywords internal
.gower_medoid <- function(views, spec_step) {
  nms <- names(views)
  k <- length(nms)
  D <- matrix(0, k, k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i == j) next
      D[i, j] <- .gower_pair_distance(views[[i]], views[[j]], spec_step)
    }
  }
  nms[[which.min(rowSums(D))]]
}

#' @keywords internal
.gower_pair_distance <- function(vi, vj, spec_step) {
  shared <- intersect(vi$map[vi$mask], vj$map[vj$mask])
  if (!length(shared)) {
    return(Inf)
  }
  ii <- match(shared, vi$map)
  jj <- match(shared, vj$map)
  pair <- procrustes(vi$X[ii, , drop = FALSE], vj$X[jj, , drop = FALSE], transform = spec_step)
  pair$objective
}

#' @keywords internal
.gower_align_to_view <- function(views, spec_step, ref_name, n) {
  ref <- views[[ref_name]]
  d <- ncol(ref$X)
  M <- matrix(0, n, d)
  M[ref$map[ref$mask], ] <- ref$X[ref$mask, , drop = FALSE]
  transforms <- lapply(views, function(v) proc_identity(spec_step, d))
  names(transforms) <- names(views)
  for (nm in names(views)) {
    if (identical(nm, ref_name)) next
    transforms[[nm]] <- .gower_fit_one(views[[nm]], M, spec_step, transforms[[nm]])
  }
  transforms
}

#' @keywords internal
.gower_spectral_init <- function(views, spec_step, d) {
  k <- length(views)
  C <- matrix(0, k * d, k * d)
  Xs <- lapply(views, function(v) v$X)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      blk <- as.matrix(centered_crossprod(Xs[[i]], Xs[[j]]))
      ii <- ((i - 1L) * d + 1L):(i * d)
      jj <- ((j - 1L) * d + 1L):(j * d)
      C[ii, jj] <- blk
    }
  }
  ev <- .gproc_eigen_sym(C)
  U <- ev$vectors[, seq_len(d), drop = FALSE]
  transforms <- vector("list", k)
  names(transforms) <- names(views)
  for (i in seq_len(k)) {
    Ui <- U[((i - 1L) * d + 1L):(i * d), , drop = FALSE]
    R <- .polar_factor(Ui, spec_step$group)$R
    transforms[[i]] <- .proc_fitted(spec_step, R, 1, rep(0, d))
  }
  transforms
}

#' @keywords internal
.gower_fit_transforms <- function(views, M, spec_step, transforms, anchor) {
  for (nm in names(views)) {
    if (!is.null(anchor) && identical(nm, anchor)) next
    transforms[[nm]] <- .gower_fit_one(views[[nm]], M, spec_step, transforms[[nm]])
  }
  transforms
}

#' @keywords internal
#' After a collective scale step, restore the exact translation block update.
#'
#' @noRd
.gower_refit_translations <- function(views, M, transforms, spec_step) {
  if (!.gproc_eliminates_translation(spec_step)) {
    return(transforms)
  }
  for (nm in names(views)) {
    view <- views[[nm]]
    idx <- view$map[view$mask]
    if (!length(idx)) next
    Y <- M[idx, , drop = FALSE]
    keep <- rowSums(is.finite(Y)) == ncol(Y)
    if (!any(keep)) next
    X <- view$X[view$mask, , drop = FALSE][keep, , drop = FALSE]
    Y <- Y[keep, , drop = FALSE]
    w <- if (!is.null(view$w)) view$w[view$mask][keep] else NULL
    if (!is.null(w) && sum(w) <= 0) next
    tr <- transforms[[nm]]
    xbar <- weighted_centroid(X, w)
    ybar <- weighted_centroid(Y, w)
    tr$t <- as.numeric(ybar - tr$s * (xbar %*% tr$R))
    transforms[[nm]] <- tr
  }
  transforms
}

#' @keywords internal
.gower_fit_one <- function(view, M, spec_step, current) {
  idx <- view$map[view$mask]
  if (!length(idx)) {
    return(current)
  }
  X <- view$X[view$mask, , drop = FALSE]
  Y <- M[idx, , drop = FALSE]
  keep <- rowSums(is.finite(Y)) == ncol(Y)
  if (!any(keep)) {
    return(current)
  }
  w <- if (!is.null(view$w)) view$w[view$mask][keep] else NULL
  if (!is.null(w) && sum(w) <= 0) {
    return(current)
  }
  pair <- procrustes(X[keep, , drop = FALSE], Y[keep, , drop = FALSE],
                     transform = spec_step, weights = w)
  pair$transform
}

#' @keywords internal
.gower_aligned <- function(views, transforms, n, d) {
  out <- vector("list", length(views))
  names(out) <- names(views)
  for (nm in names(views)) {
    Yi <- apply_proc_transform(transforms[[nm]], views[[nm]]$X)
    G <- matrix(NA_real_, n, d)
    obs <- views[[nm]]$map[views[[nm]]$mask]
    G[obs, ] <- Yi[views[[nm]]$mask, , drop = FALSE]
    out[[nm]] <- G
  }
  out
}

#' @keywords internal
.gower_consensus <- function(Y, alphas, n, d, views = NULL) {
  num <- matrix(0, n, d)
  den <- rep(0, n)
  for (nm in names(Y)) {
    Yi <- Y[[nm]]
    a <- alphas[[nm]]
    obs <- is.finite(Yi[, 1L])
    if (!any(obs)) next
    w <- rep(1, n)
    if (!is.null(views) && !is.null(views[[nm]]$w)) {
      w[views[[nm]]$map] <- views[[nm]]$w
    }
    ww <- a * w[obs]
    num[obs, ] <- num[obs, , drop = FALSE] + ww * Yi[obs, , drop = FALSE]
    den[obs] <- den[obs] + ww
  }
  M <- num
  ok <- den > 0
  M[ok, ] <- num[ok, , drop = FALSE] / den[ok]
  M[!ok, ] <- 0
  M
}

#' @keywords internal
.gower_objective <- function(Y, M, alphas, views = NULL) {
  f <- 0
  for (nm in names(Y)) {
    Yi <- Y[[nm]]
    obs <- is.finite(Yi[, 1L])
    if (!any(obs)) next
    w <- rep(1, sum(obs))
    if (!is.null(views) && !is.null(views[[nm]]$w)) {
      wg <- rep(0, length(obs))
      wg[views[[nm]]$map] <- views[[nm]]$w
      w <- wg[obs]
    }
    dlt <- Yi[obs, , drop = FALSE] - M[obs, , drop = FALSE]
    f <- f + alphas[[nm]] * sum(w * rowSums(dlt * dlt))
  }
  f
}

#' @keywords internal
.gower_apply_scale <- function(state, views, Y, M, alphas, scale_mode, n, d) {
  if (identical(scale_mode, "gower")) {
    return(.gower_tenberge(state, views, alphas, n, d))
  }
  if (identical(scale_mode, "fixed_consensus")) {
    nrm <- sqrt(sum(M^2))
    if (nrm > 0) {
      fac <- 1 / nrm
      for (nm in names(state$transforms)) {
        state$transforms[[nm]]$s <- state$transforms[[nm]]$s * fac
        state$transforms[[nm]]$t <- state$transforms[[nm]]$t * fac
      }
    }
  }
  state
}

#' Ten Berge / Gower collective scaling: Gs = lambda Ds on unscaled rotations.
#'
#' @noRd
.gower_tenberge <- function(state, views, alphas, n, d) {
  nms <- names(views)
  k <- length(nms)
  Z <- vector("list", k)
  names(Z) <- nms
  for (nm in nms) {
    tr <- state$transforms[[nm]]
    Zi <- views[[nm]]$X %*% tr$R
    G <- matrix(NA_real_, n, d)
    G[views[[nm]]$map[views[[nm]]$mask], ] <- Zi[views[[nm]]$mask, , drop = FALSE]
    Z[[nm]] <- G
  }
  Gmat <- matrix(0, k, k)
  D <- rep(0, k)
  a <- as.numeric(alphas[nms])
  names(a) <- nms
  for (i in seq_len(k)) {
    Zi <- Z[[i]]
    obs_i <- is.finite(Zi[, 1L])
    D[i] <- a[[i]] * sum(Zi[obs_i, , drop = FALSE]^2)
    for (j in seq_len(k)) {
      Zj <- Z[[j]]
      obs <- obs_i & is.finite(Zj[, 1L])
      if (any(obs)) {
        Gmat[i, j] <- a[[i]] * a[[j]] * sum(Zi[obs, , drop = FALSE] * Zj[obs, , drop = FALSE])
      }
    }
  }
  if (k == 1L || any(!is.finite(D)) || all(D <= 0)) {
    return(state)
  }
  keep <- D > 0
  if (sum(keep) < 2L) {
    return(state)
  }
  Ds <- D[keep]
  Gs <- Gmat[keep, keep, drop = FALSE]
  Dinvs <- 1 / sqrt(Ds)
  A <- t(t(Gs * Dinvs) * Dinvs)
  ev <- .gproc_eigen_sym(A)
  v <- ev$vectors[, 1L]
  s_keep <- abs(v) * Dinvs
  c0 <- sum(vapply(views[nms[keep]], `[[`, numeric(1), "energy") * a[keep])
  energy <- sum(s_keep^2 * Ds)
  if (energy > 0 && c0 > 0) {
    s_keep <- s_keep * sqrt(c0 / energy)
  }
  s <- rep(1, k)
  names(s) <- nms
  s[keep] <- s_keep
  for (nm in nms) {
    old_s <- state$transforms[[nm]]$s
    fac <- if (is.finite(old_s) && old_s > 0) unname(s[[nm]]) / old_s else unname(s[[nm]])
    state$transforms[[nm]]$s <- unname(s[[nm]])
    state$transforms[[nm]]$t <- state$transforms[[nm]]$t * fac
  }
  state
}

#' @keywords internal
.gower_consensus_change <- function(M_old, M_new) {
  ok <- is.finite(M_old[, 1L]) & is.finite(M_new[, 1L])
  if (!any(ok)) {
    return(0)
  }
  pair <- tryCatch(
    procrustes(M_old[ok, , drop = FALSE], M_new[ok, , drop = FALSE],
               transform = proc_orthogonal("O")),
    error = function(e) NULL
  )
  if (is.null(pair)) {
    return(sqrt(sum((M_new[ok, , drop = FALSE] - M_old[ok, , drop = FALSE])^2)))
  }
  aligned <- apply_proc_transform(pair$transform, M_old[ok, , drop = FALSE])
  sqrt(sum((aligned - M_new[ok, , drop = FALSE])^2)) / max(1, sqrt(sum(M_new[ok, ]^2)))
}

#' @keywords internal
.gower_stationarity <- function(views, M, transforms, spec_step = NULL) {
  rmax <- 0
  translate <- .gproc_eliminates_translation(spec_step)
  scale_on <- identical(spec_step$scaling, "isotropic")
  for (nm in names(views)) {
    view <- views[[nm]]
    idx <- view$map[view$mask]
    if (!length(idx)) next
    Y <- M[idx, , drop = FALSE]
    keep <- rowSums(is.finite(Y)) == ncol(Y)
    if (!any(keep)) next
    X <- view$X[view$mask, , drop = FALSE][keep, , drop = FALSE]
    Y <- Y[keep, , drop = FALSE]
    w <- if (!is.null(view$w)) view$w[view$mask][keep] else NULL
    if (!is.null(w) && sum(w) <= 0) next
    moms <- proc_moments(X, Y, w, center = translate)
    C <- moms$C
    R <- transforms[[nm]]$R
    s <- transforms[[nm]]$s
    S <- crossprod(R, C)
    skew <- 0.5 * (S - t(S))
    rmax <- max(rmax, sqrt(sum(skew^2)))
    if (translate) {
      Yhat <- apply_proc_transform(transforms[[nm]], X)
      r_t <- sqrt(sum((weighted_centroid(Yhat, w) - weighted_centroid(Y, w))^2))
      rmax <- max(rmax, r_t)
    }
    if (scale_on) {
      gamma <- sum(R * C)
      rmax <- max(rmax, abs(gamma - s * moms$a))
    }
  }
  rmax
}

#' @keywords internal
.gower_orthogonality <- function(transforms) {
  max(vapply(transforms, function(tr) {
    constraint_residual(tr)$orthogonality
  }, numeric(1)))
}
