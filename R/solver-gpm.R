#' Matrix-free generalized power method for orthogonal GOPP.
#'
#' Never forms \eqn{C=B^\top B}. A step is \eqn{Z=\sum_j\alpha_j X_j R_j},
#' \eqn{B_i=\alpha_i X_i^\top Z}, then \eqn{R_i\leftarrow\mathrm{polar}(B_i)}.
#' Row masks use the \eqn{D^\dagger} operator of math section 18.
#'
#' @noRd
.gpm_gpa <- function(data, spec, gauge, metric, alphas, control, anchor = NULL) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  d <- ncol(data$views[[1L]])
  if (any(vapply(data$views, ncol, integer(1)) != d)) {
    .gproc_stop("dimension_mismatch", "GPM requires a common target dimension.")
  }
  if (!is.null(anchor) && !anchor %in% nms) {
    .gproc_stop("invalid_problem", "Unknown anchor view.")
  }
  views <- .gpm_prepare_views(data)
  complete <- all(vapply(views, function(v) all(v$mask) && length(v$map) == n, logical(1)))
  backend <- .gpm_backend(control, length(nms) * d)
  requested_init <- control$init %||% "spectral"
  starts <- .gpm_start_names(control)
  candidates <- list()
  for (init in starts) {
    candidates[[init]] <- .gpm_run(
      views, spec, alphas, control, anchor, n, d, nms, complete, backend, init
    )
  }

  best_name <- names(candidates)[[which.min(vapply(candidates, `[[`, numeric(1), "objective"))]]
  best <- candidates[[best_name]]
  restart_reason <- NULL
  if (.gpm_zero_consensus(best, alphas)) {
    restart_reason <- paste(
      sprintf("%s initialization converged to a zero-consensus positive-objective point;", best_name),
      "tried a deterministic alternative"
    )
    remaining <- setdiff(c("spectral", "medoid", "sequential"), names(candidates))
    for (init in remaining) {
      candidates[[init]] <- .gpm_run(
        views, spec, alphas, control, anchor, n, d, nms, complete, backend, init
      )
      best_name <- names(candidates)[[which.min(vapply(candidates, `[[`, numeric(1), "objective"))]]
      best <- candidates[[best_name]]
      if (!.gpm_zero_consensus(best, alphas)) {
        break
      }
    }
  }

  objectives <- vapply(candidates, `[[`, numeric(1), "objective")
  best_name <- names(candidates)[[which.min(objectives)]]
  best <- candidates[[best_name]]
  best$requested_init <- requested_init
  best$starts <- names(candidates)
  best$start_objectives <- objectives
  best$degenerate_starts <- names(candidates)[vapply(
    candidates, .gpm_zero_consensus, logical(1), alphas = alphas
  )]
  best$restart_reason <- restart_reason
  best
}

#' @noRd
.gpm_start_names <- function(control) {
  init <- control$init %||% "spectral"
  init <- .gproc_match_arg(init, c("sequential", "medoid", "spectral"), "init")
  nstart <- as.integer(control$nstart %||% 1L)
  if (!is.finite(nstart) || nstart < 1L) nstart <- 1L
  nstart <- min(nstart, 3L)
  utils::head(unique(c(init, "medoid", "spectral", "sequential")), nstart)
}

#' One GPM run from one deterministic initializer.
#'
#' @noRd
.gpm_run <- function(views, spec, alphas, control, anchor, n, d, nms,
                     complete, backend, init) {
  Rs <- .gpm_initialize(
    views, spec, alphas, n, d, nms, init, complete, backend, anchor
  )
  t0 <- proc.time()[["elapsed"]]
  maxit <- control$max_iterations %||% 500L
  tol <- control$tolerance %||% 1e-8
  status <- "maximum_iterations"
  hist <- list()
  B <- .gpm_apply_C(views, alphas, Rs, n, complete)
  obj <- .gpm_objective(views, Rs, alphas, n, d)
  for (it in seq_len(maxit)) {
    old_Rs <- Rs
    old_obj <- obj
    Rs <- .gpm_polar_blocks(B, spec$group, anchor)
    B <- .gpm_apply_C(views, alphas, Rs, n, complete)
    obj <- .gpm_objective(views, Rs, alphas, n, d)
    stat <- .gpm_stationarity(Rs, B)
    orth <- .gpm_orthogonality(Rs)
    gchg <- .gpm_gauge_change(old_Rs, Rs)
    rel <- abs(old_obj - obj) / max(1, abs(old_obj))
    hist[[it]] <- data.frame(
      iteration = it,
      objective = obj,
      relative_change = rel,
      consensus_change = gchg,
      stationarity = stat,
      orthogonality_residual = orth,
      elapsed_time = proc.time()[["elapsed"]] - t0
    )
    if (rel <= tol && gchg <= sqrt(tol) && stat <= sqrt(tol)) {
      status <- "converged"
      break
    }
  }
  hist_df <- if (length(hist)) do.call(rbind, hist) else data.frame()
  Y <- .gpm_aligned(views, Rs, n, d)
  M <- .gower_consensus(Y, alphas, n, d)
  transforms <- .gpm_original_transforms(Rs, views, spec, d)
  cert <- .gpm_maybe_certificate(views, alphas, Rs, B, spec, control, complete, backend, n, d)
  opt <- if (!identical(status, "converged")) {
    "not_converged"
  } else if (identical(cert$status, "certified_global")) {
    "certified_global"
  } else {
    "first_order_stationary"
  }
  list(
    transforms = transforms,
    work_transforms = transforms,
    consensus = M,
    aligned = Y,
    objective = .gower_objective(Y, M, alphas),
    history = hist_df,
    numerical_status = status,
    optimality_status = opt,
    init = init,
    stationarity = if (nrow(hist_df)) hist_df$stationarity[[nrow(hist_df)]] else NA_real_,
    certificate = cert,
    backend = backend,
    complete = complete
  )
}

#' Detect a nontrivial fit whose consensus has collapsed relative to its data energy.
#'
#' The ratio is scale-free. All-zero inputs are legitimate and do not trigger a
#' restart; a collapsed nonzero problem has essentially all energy in residuals.
#'
#' @noRd
.gpm_zero_consensus <- function(candidate, alphas) {
  total_energy <- sum(vapply(names(candidate$aligned), function(nm) {
    alphas[[nm]] * sum(candidate$aligned[[nm]]^2)
  }, numeric(1)))
  if (!is.finite(total_energy) || total_energy <= .Machine$double.xmin) {
    return(FALSE)
  }
  consensus_energy <- sum(alphas) * sum(candidate$consensus^2)
  tol <- 256 * .Machine$double.eps
  is.finite(consensus_energy) &&
    consensus_energy <= tol * total_energy &&
    candidate$objective > tol * total_energy
}

#' @noRd
.gpm_backend <- function(control, kd) {
  req <- control$backend %||% "auto"
  thresh <- control$dense_block_threshold %||% 96L
  if (identical(req, "matrix_free")) {
    return("matrix_free")
  }
  if (identical(req, "dense") || (identical(req, "auto") && kd <= thresh)) {
    return("dense")
  }
  "matrix_free"
}

#' @noRd
.gpm_prepare_views <- function(data) {
  nms <- names(data$views)
  n <- length(data$global_ids)
  out <- vector("list", length(nms))
  names(out) <- nms
  for (nm in nms) {
    X <- as.matrix(data$views[[nm]])
    mask <- data$observed[[nm]]
    map <- data$row_map[[nm]]
    xbar <- weighted_centroid(X[mask, , drop = FALSE], w = NULL)
    Xc <- sweep(X, 2L, xbar, "-")
    out[[nm]] <- list(
      X = Xc,
      mask = mask,
      map = map,
      name = nm,
      pre_t = -xbar,
      n_global = n
    )
  }
  out
}

#' @noRd
.gpm_initialize <- function(views, spec, alphas, n, d, nms, init, complete, backend, anchor) {
  if (!is.null(anchor)) {
    Rs <- lapply(views, function(v) diag(d))
    names(Rs) <- nms
    ref <- views[[anchor]]
    M <- matrix(0, n, d)
    M[ref$map[ref$mask], ] <- ref$X[ref$mask, , drop = FALSE]
    for (nm in nms) {
      if (identical(nm, anchor)) next
      B <- .gpm_fit_cross(views[[nm]], M)
      Rs[[nm]] <- .polar_factor(B, spec$group)$R
    }
    return(Rs)
  }
  if (identical(init, "medoid")) {
    spec_step <- proc_orthogonal(spec$group)
    med <- .gower_medoid(views, spec_step)
    tr <- .gower_align_to_view(views, spec_step, med, n)
    return(lapply(tr, `[[`, "R"))
  }
  if (identical(init, "sequential")) {
    spec_step <- proc_orthogonal(spec$group)
    tr <- .gower_sequential_init(views, spec_step, alphas, n, d)
    return(lapply(tr, `[[`, "R"))
  }
  .gpm_spectral_init(views, spec, alphas, n, d, complete, backend)
}

#' Leading d-eigenspace of the implicit GOPP operator, then block polar factors.
#'
#' @noRd
.gpm_spectral_init <- function(views, spec, alphas, n, d, complete, backend) {
  nms <- names(views)
  k <- length(nms)
  if (complete && identical(backend, "dense")) {
    C <- .gpm_form_C(views, alphas, d, n)
    ev <- .gproc_eigen_sym(C)
    U <- ev$vectors[, seq_len(d), drop = FALSE]
  } else {
    U <- matrix(0, k * d, d)
    U[seq_len(d), ] <- diag(d)
    for (it in seq_len(40L)) {
      blocks <- .gpm_unstack(U, nms, d)
      B <- .gpm_apply_C(views, alphas, blocks, n, complete)
      W <- .gpm_stack(B)
      qrW <- qr(W)
      U_new <- qr.Q(qrW)[, seq_len(d), drop = FALSE]
      gap <- sqrt(sum((U_new - U)^2))
      U <- U_new
      if (gap < 1e-10) break
    }
  }
  blocks <- .gpm_unstack(U, nms, d)
  out <- vector("list", k)
  names(out) <- nms
  for (nm in nms) {
    out[[nm]] <- .polar_factor(blocks[[nm]], spec$group)$R
  }
  out
}

#' @noRd
.gpm_form_C <- function(views, alphas, d, n) {
  nms <- names(views)
  k <- length(nms)
  Xs <- lapply(views, function(v) {
    G <- matrix(0, n, d)
    G[v$map, ] <- v$X
    G
  })
  C <- matrix(0, k * d, k * d)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      blk <- as.matrix(crossprod(Xs[[i]], Xs[[j]]))
      ii <- ((i - 1L) * d + 1L):(i * d)
      jj <- ((j - 1L) * d + 1L):(j * d)
      C[ii, jj] <- alphas[[i]] * alphas[[j]] * blk
    }
  }
  (C + t(C)) / 2
}

#' @noRd
.gpm_apply_C <- function(views, alphas, Rs, n, complete) {
  nms <- names(views)
  d <- ncol(Rs[[1L]])
  p <- d
  Z <- matrix(0, n, p)
  den <- rep(0, n)
  for (nm in nms) {
    v <- views[[nm]]
    Yi <- v$X %*% Rs[[nm]]
    a <- alphas[[nm]]
    obs <- v$map[v$mask]
    Z[obs, ] <- Z[obs, , drop = FALSE] + a * Yi[v$mask, , drop = FALSE]
    den[obs] <- den[obs] + a
  }
  if (!complete) {
    ok <- den > 0
    Z[ok, ] <- Z[ok, , drop = FALSE] / den[ok]
    Z[!ok, ] <- 0
  }
  B <- vector("list", length(nms))
  names(B) <- nms
  for (nm in nms) {
    v <- views[[nm]]
    loc <- Z[v$map, , drop = FALSE]
    if (!all(v$mask)) {
      loc[!v$mask, ] <- 0
    }
    Bi <- alphas[[nm]] * crossprod(v$X, loc)
    B[[nm]] <- Bi
  }
  B
}

#' @noRd
.gpm_polar_blocks <- function(B, group, anchor) {
  out <- B
  for (nm in names(B)) {
    if (!is.null(anchor) && identical(nm, anchor)) {
      out[[nm]] <- diag(nrow(B[[nm]]))
    } else {
      out[[nm]] <- .polar_factor(B[[nm]], group)$R
    }
  }
  out
}

#' @noRd
.gpm_fit_cross <- function(view, M) {
  idx <- view$map[view$mask]
  X <- view$X[view$mask, , drop = FALSE]
  Y <- M[idx, , drop = FALSE]
  keep <- rowSums(is.finite(Y)) == ncol(Y)
  if (!any(keep)) {
    return(diag(ncol(view$X)))
  }
  as.matrix(crossprod(X[keep, , drop = FALSE], Y[keep, , drop = FALSE]))
}

#' @noRd
.gpm_aligned <- function(views, Rs, n, d) {
  transforms <- lapply(Rs, function(R) {
    .proc_fitted(proc_orthogonal("O"), R, 1, rep(0, ncol(R)))
  })
  names(transforms) <- names(views)
  .gower_aligned(views, transforms, n, d)
}

#' @noRd
.gpm_objective <- function(views, Rs, alphas, n, d) {
  Y <- .gpm_aligned(views, Rs, n, d)
  M <- .gower_consensus(Y, alphas, n, d)
  .gower_objective(Y, M, alphas)
}

#' @noRd
.gpm_stationarity <- function(Rs, B) {
  rmax <- 0
  for (nm in names(Rs)) {
    S <- crossprod(Rs[[nm]], B[[nm]])
    skew <- 0.5 * (S - t(S))
    rmax <- max(rmax, sqrt(sum(skew^2)))
  }
  rmax
}

#' @noRd
.gpm_orthogonality <- function(Rs) {
  max(vapply(Rs, function(R) {
    sqrt(sum((crossprod(R) - diag(nrow(R)))^2))
  }, numeric(1)))
}

#' Ling \eqn{d_F(S,T)=\min_Q\|S-TQ\|_F}.
#'
#' @noRd
.gpm_gauge_change <- function(Rs_old, Rs_new) {
  Sold <- .gpm_stack(Rs_old)
  Snew <- .gpm_stack(Rs_new)
  C <- crossprod(Sold, Snew)
  Q <- .polar_factor(C, "O")$R
  sqrt(sum((Snew - Sold %*% Q)^2)) / max(1, sqrt(sum(Snew^2)))
}

#' @noRd
.gpm_stack <- function(blocks) {
  do.call(rbind, blocks[names(blocks)])
}

#' @noRd
.gpm_unstack <- function(S, nms, d) {
  out <- vector("list", length(nms))
  names(out) <- nms
  for (i in seq_along(nms)) {
    out[[i]] <- S[((i - 1L) * d + 1L):(i * d), , drop = FALSE]
  }
  out
}

#' @noRd
.gpm_original_transforms <- function(Rs, views, spec, d) {
  out <- vector("list", length(Rs))
  names(out) <- names(Rs)
  for (nm in names(Rs)) {
    work <- .proc_fitted(spec, Rs[[nm]], 1, rep(0, d))
    pre <- .proc_fitted(spec, diag(d), 1, views[[nm]]$pre_t)
    tr <- compose_proc_transform(pre, work)
    tr$spec <- spec
    out[[nm]] <- tr
  }
  out
}

#' @noRd
.gpm_maybe_certificate <- function(views, alphas, Rs, B, spec, control, complete,
                                   backend, n, d) {
  mode <- control$certify %||% "auto"
  if (identical(mode, "never")) {
    return(list(
      status = "unavailable",
      reason = "certify='never'"
    ))
  }
  if (!identical(spec$group, "O")) {
    return(list(
      status = "unavailable",
      reason = "Ling's dual certificate is for the O(d) GOPP trace problem; it is not reused for SO(d)."
    ))
  }
  .gpm_certificate(views, alphas, Rs, B, complete, backend, n, d)
}

#' Dual slack \eqn{\Lambda-C} at a GPM candidate (math section 17).
#'
#' @noRd
.gpm_certificate <- function(views, alphas, Rs, B, complete, backend, n, d) {
  nms <- names(Rs)
  k <- length(nms)
  Lambda <- lapply(B, .gpm_sym_sqrt_bbt)
  names(Lambda) <- nms
  CS <- .gpm_stack(B)
  LS <- .gpm_stack(lapply(nms, function(nm) Lambda[[nm]] %*% Rs[[nm]]))
  r_dual <- sqrt(sum((CS - LS)^2)) / (1 + sqrt(sum(CS^2)))
  eigs <- .gpm_dual_eigs(views, alphas, Lambda, Rs, complete, backend, n, d)
  ev <- sort(as.numeric(eigs$values))
  lam_min <- ev[[1L]]
  lam_next <- if (length(ev) >= d + 1L) ev[[d + 1L]] else NA_real_
  scale <- max(1, max(abs(ev)))
  tol_e <- 1e-7 * scale
  nullity <- as.integer(sum(abs(ev) <= tol_e))
  bound_ok <- isTRUE(eigs$validated_lower_bound)
  psd <- bound_ok && lam_min >= -tol_e
  unique_mod <- bound_ok && is.finite(lam_next) && lam_next > tol_e
  certified <- isTRUE(psd) && r_dual <= 1e-6
  list(
    status = if (certified) "certified_global" else "not_certified",
    r_dual = unname(r_dual),
    lambda_min = unname(lam_min),
    lambda_d_plus_1 = unname(lam_next),
    uniqueness_modulo_O = unique_mod,
    nullity = nullity,
    expected_nullity = as.integer(d),
    dual_spectrum = eigs$method,
    reason = if (certified) {
      "CS = Lambda S and a numerical lower bound on lambda_min(Lambda - C) is nonnegative."
    } else if (!bound_ok) {
      paste(
        "matrix-free Ritz values are an upper bound on lambda_min(Lambda - C),",
        "not a validated lower bound; certificate unavailable"
      )
    } else if (!psd) {
      sprintf("estimated smallest dual eigenvalue = %.6g; stationary candidate only", lam_min)
    } else {
      sprintf("dual stationarity residual r_dual = %.6g is above tolerance", r_dual)
    }
  )
}

#' Lambda_i is the symmetric square root of B_i B_i^T.
#'
#' @noRd
.gpm_sym_sqrt_bbt <- function(B) {
  sv <- .gproc_svd(B, nu = nrow(B), nv = 0L)
  sv$u %*% (sv$d * t(sv$u))
}

#' Smallest eigenvalues of \eqn{\Lambda-C}, forming the block Gram only when dense.
#'
#' @noRd
.gpm_dual_eigs <- function(views, alphas, Lambda, Rs, complete, backend, n, d) {
  nms <- names(views)
  k <- length(nms)
  kd <- k * d
  if (complete && identical(backend, "dense") && kd <= 256L) {
    C <- .gpm_form_C(views, alphas, d, n)
    Lbig <- matrix(0, kd, kd)
    for (i in seq_len(k)) {
      ii <- ((i - 1L) * d + 1L):(i * d)
      Lbig[ii, ii] <- Lambda[[i]]
    }
    ev <- .gproc_eigen_sym(Lbig - C)$values
    return(list(values = ev, validated_lower_bound = TRUE, method = "dense_full"))
  }
  ev <- .gpm_dual_eigs_matrix_free(views, alphas, Lambda, complete, n, d)
  list(values = ev, validated_lower_bound = FALSE, method = "ritz")
}

#' @noRd
.gpm_dual_eigs_matrix_free <- function(views, alphas, Lambda, complete, n, d) {
  nms <- names(views)
  k <- length(nms)
  kd <- k * d
  apply_L <- function(V) {
    blocks <- .gpm_unstack(V, nms, d)
    Cv <- .gpm_apply_C(views, alphas, blocks, n, complete)
    out <- vector("list", k)
    names(out) <- nms
    for (nm in nms) {
      out[[nm]] <- Lambda[[nm]] %*% blocks[[nm]] - Cv[[nm]]
    }
    .gpm_stack(out)
  }
  # Power iteration on a shifted operator for the smallest algebraic eigenvalues.
  nrm <- 0
  v <- stats::rnorm(kd)
  v <- v / sqrt(sum(v^2))
  for (it in seq_len(8L)) {
    w <- as.numeric(apply_L(matrix(v, kd, 1L)))
    nrm <- sqrt(sum(w^2))
    if (nrm > 0) v <- w / nrm
  }
  shift <- nrm + 1
  U <- matrix(stats::rnorm(kd * min(d + 3L, kd)), kd, min(d + 3L, kd))
  U <- qr.Q(qr(U))
  for (it in seq_len(40L)) {
    W <- shift * U - apply_L(U)
    U <- qr.Q(qr(W))
  }
  G <- crossprod(U, apply_L(U))
  .gproc_eigen_sym(G)$values
}
