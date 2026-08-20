#' @export
print.gpa_fit <- function(x, ...) {
  ds <- x$data_summary
  niter <- if (is.data.frame(x$history) && nrow(x$history)) max(x$history$iteration) else 0L
  cat(sprintf(
    paste0(
      "Generalized Procrustes fit\n",
      "  Configurations: %d\n",
      "  Consensus entities: %d\n",
      "  Target dimension: %d\n",
      "  Transformation: %s, scaling=%s, translation=%s\n",
      "  Correspondence: %s\n",
      "  Overlap graph: %s\n",
      "  Solver: %s\n",
      "  Iterations: %s\n",
      "  Numerical status: %s\n",
      "  Optimality: %s\n",
      "  Certificate: %s\n",
      "  Aligned data: %s\n"
    ),
    ds$n_views,
    ds$n_entities,
    ds$target_dimension,
    x$problem$transform$group,
    x$problem$transform$scaling,
    x$problem$transform$translation,
    ds$correspondence,
    if (ds$overlap$n_components == 1L) "connected" else sprintf("%d components", ds$overlap$n_components),
    x$solver,
    niter,
    x$numerical_status,
    x$optimality_status,
    x$certificate$status %||% "unavailable",
    x$keep_aligned
  ))
  invisible(x)
}

#' @export
summary.gpa_fit <- function(object, ...) {
  print(object)
  cat(sprintf("  Objective: %.6g\n", object$objective))
  cat(sprintf("  Effective rank: %s\n", object$rank$effective_rank))
  invisible(object)
}

#' Consensus matrix.
#'
#' @param object A `gpa_fit`.
#' @param gauge `"native"` for the optimizer's representative,
#'   `"canonical"` for a principal-axis presentation, or a reference
#'   `gpa_fit` / matrix to align to. Canonicalization is display, not estimation.
#' @param ... Unused.
#' @return The consensus matrix in the requested gauge.
#' @export
consensus <- function(object, gauge = "native", ...) {
  UseMethod("consensus")
}

#' @export
consensus.gpa_fit <- function(object, gauge = "native", ...) {
  M <- object$consensus
  if (inherits(gauge, "gpa_fit") || is.matrix(gauge)) {
    ref <- if (inherits(gauge, "gpa_fit")) gauge$consensus else as.matrix(gauge)
    grp <- if (identical(object$problem$transform$group, "SO")) "SO" else "O"
    pair <- procrustes(M, ref, transform = proc_orthogonal(grp))
    return(M %*% pair$transform$R)
  }
  gauge <- as.character(gauge)[1L]
  if (identical(gauge, "native")) {
    return(M)
  }
  if (identical(gauge, "canonical")) {
    return(.gproc_canonical_consensus(M, object$problem$transform$group))
  }
  .gproc_stop("invalid_problem", "gauge must be 'native', 'canonical', or a reference.")
}

#' Fitted transforms.
#' @param object A `gpa_fit`.
#' @param ... Unused.
#' @export
transformations <- function(object, ...) {
  UseMethod("transformations")
}

#' @export
transformations.gpa_fit <- function(object, ...) {
  object$transformations
}

#' Lazy aligned views.
#' @param object A `gpa_fit`.
#' @param ... Unused.
#' @export
aligned <- function(object, ...) {
  UseMethod("aligned")
}

#' @export
aligned.gpa_fit <- function(object, ...) {
  store <- object$aligned_store
  if (identical(object$keep_aligned, "materialize")) {
    lapply(store, function(v) {
      if (inherits(v, "proc_aligned_view")) as.matrix(v) else v
    })
  } else {
    store
  }
}

#' @export
fitted.gpa_fit <- function(object, ...) {
  .gproc_fitted_global(object)
}

#' @export
residuals.gpa_fit <- function(object, level = c("landmark", "configuration"), ...) {
  level <- match.arg(level)
  M <- consensus(object)
  Y <- .gproc_fitted_global(object)
  res <- lapply(Y, function(Yi) Yi - M)
  if (identical(level, "configuration")) {
    return(vapply(res, function(E) sqrt(sum(E^2, na.rm = TRUE)), numeric(1)))
  }
  res
}

#' Scatter lazy aligned views onto the consensus index.
#' @noRd
.gproc_fitted_global <- function(object) {
  n <- nrow(object$consensus)
  d <- ncol(object$consensus)
  lapply(object$aligned_store, function(v) {
    if (inherits(v, "proc_aligned_view")) {
      return(.gproc_scatter_aligned(v, n, d))
    }
    if (is.matrix(v) && nrow(v) == n) {
      return(v)
    }
    v
  })
}

#' @keywords internal
.gproc_scatter_aligned <- function(view, n, d) {
  Yi <- as.matrix(view)
  G <- matrix(NA_real_, n, d)
  map <- view$row_map
  if (is.null(map)) map <- seq_len(nrow(Yi))
  G[map, ] <- Yi
  if (!is.null(view$mask) && !all(view$mask)) {
    G[map[!view$mask], ] <- NA_real_
  }
  G
}

#' @export
coef.gpa_fit <- function(object, ...) {
  lapply(transformations(object), function(tr) {
    list(R = tr$R, s = tr$s, t = tr$t)
  })
}

#' Align new data to a fitted consensus.
#'
#' @param object A `gpa_fit`.
#' @param newdata A matrix or `proc_data` with one view.
#' @param ... Unused.
#' @export
predict.gpa_fit <- function(object, newdata, ...) {
  spec <- object$problem$transform
  if (inherits(newdata, "proc_data")) {
    newdata <- newdata$views[[1L]]
  }
  newdata <- as.matrix(newdata)
  if (spec$family %in% c("affine", "lbw", "tps")) {
    M <- consensus(object)
    if (nrow(newdata) != nrow(M)) {
      .gproc_stop("invalid_problem", "predict() for LBW needs newdata with n_entities rows.")
    }
    built <- .lbw_basis(spec, newdata)
    w <- rep(1, nrow(newdata))
    A <- crossprod(built$Phi, w * built$Phi) + (spec$smoothness %||% 0) * built$L
    B <- .gproc_psd_solve(A, crossprod(built$Phi, w * M))
    return(built$Phi %*% B)
  }
  pair <- procrustes(newdata, consensus(object), transform = spec)
  apply_proc_transform(pair$transform, newdata)
}

#' Rank, overlap, and stationarity diagnostics.
#' @param object A fit.
#' @param ... Unused.
#' @export
diagnose <- function(object, ...) {
  UseMethod("diagnose")
}

#' @export
diagnose.gpa_fit <- function(object, ...) {
  list(
    rank = object$rank,
    overlap = object$missing_support,
    constraint = lapply(transformations(object), constraint_residual),
    stationarity = object$rank$stationarity %||% NA_real_,
    numerical_status = object$numerical_status,
    optimality_status = object$optimality_status,
    certificate = object$certificate,
    scale_mode = object$gauge$scale_mode %||% object$gauge$spec$scale,
    folding = object$folding,
    reference_covariance = object$reference_covariance,
    tied_axis_blocks = object$gauge$tied_axis_blocks %||% tied_axis_blocks(object$consensus)
  )
}

#' Datum-space error \eqn{\sum_i\|X_i-T_i^{-1}(P_iM)\|^2}.
#'
#' Distinct from the reference-space objective. Unavailable when the
#' transform has no exact inverse (math section 34).
#'
#' @param object A `gpa_fit`.
#' @export
datum_space_error <- function(object) {
  if (!inherits(object, "gpa_fit")) {
    .gproc_stop("invalid_problem", "datum_space_error() requires a gpa_fit.")
  }
  spec <- object$problem$transform
  if (spec$family %in% c("tps", "lbw") && !identical(spec$family, "affine")) {
    .gproc_stop(
      "inverse_unavailable",
      "Datum-space error needs an exact inverse; TPS / general LBW do not have one."
    )
  }
  M <- consensus(object)
  trs <- transformations(object)
  av <- object$aligned_store
  raw <- object$raw_data
  f <- 0
  for (nm in names(trs)) {
    if (!is.null(raw) && inherits(raw, "proc_data")) {
      X <- raw$views[[nm]]
      map <- raw$row_map[[nm]]
    } else if (inherits(av[[nm]], "proc_aligned_view")) {
      X <- av[[nm]]$x
      map <- av[[nm]]$row_map
    } else {
      X <- av[[nm]]
      map <- seq_len(nrow(X))
    }
    inv <- inverse_proc_transform(trs[[nm]])
    pred <- apply_proc_transform(inv, M[map, , drop = FALSE])
    f <- f + sum((as.matrix(X) - pred)^2)
  }
  f
}

#' Held-out CV over a smoothness grid for LBW / TPS.
#'
#' Training residual always falls as \eqn{\mu\downarrow}. Use this, not \eqn{L_r},
#' to choose bending energy.
#'
#' @param data `proc_data` or a list of matrices.
#' @param transform An LBW / TPS / affine spec.
#' @param smoothness Candidate \eqn{\mu} values.
#' @param folds Landmark folds.
#' @param ... Passed to `gpa()`.
#' @export
tune_smoothness <- function(data,
                            transform,
                            smoothness = c(0.1, 1, 10),
                            folds = 3L,
                            ...) {
  spec <- as_proc_transform(transform)
  scores <- vapply(smoothness, function(mu) {
    spec$smoothness <- mu
    fit <- gpa(data, transform = spec, ...)
    cross_validate(fit, folds = folds)$total
  }, numeric(1))
  structure(
    list(
      smoothness = as.numeric(smoothness),
      scores = scores,
      best = smoothness[[which.min(scores)]],
      note = "Held-out landmark discrepancy after gauge alignment; not the training residual."
    ),
    class = "proc_cv"
  )
}

#' @export
diagnose.proc_pair_fit <- function(object, ...) {
  list(
    rank = object$rank,
    singular_values = object$singular_values,
    transform_unique = object$transform_unique,
    unidentified_subspace_dimension = object$unidentified_subspace_dimension,
    constraint = constraint_residual(object$transform)
  )
}

#' Optimality certificate (pairwise: unavailable by construction).
#' @param object A fit.
#' @param ... Unused.
#' @export
certify <- function(object, ...) {
  UseMethod("certify")
}

#' @export
certify.gpa_fit <- function(object, ...) {
  object$certificate
}

#' @export
certify.proc_pair_fit <- function(object, ...) {
  list(
    status = "unavailable",
    reason = "Pairwise polar factor is an exact closed form; Ling's dual certificate is for GOPP."
  )
}

#' Energy decomposition of aligned configurations.
#'
#' @param object A `gpa_fit`, or a list of aligned matrices when `consensus` is supplied.
#' @param consensus Optional consensus matrix.
#' @param weights Configuration weights.
#' @param ... Unused.
#' @export
decompose <- function(object, consensus = NULL, weights = NULL, ...) {
  UseMethod("decompose")
}

#' @export
decompose.gpa_fit <- function(object, ...) {
  object$decomposition
}

#' @export
decompose.default <- function(object, consensus = NULL, weights = NULL, ...) {
  decompose_energy(object, consensus, weights)
}

#' @keywords internal
decompose_energy <- function(aligned_mats, M, weights) {
  nms <- names(aligned_mats)
  if (is.null(weights)) weights <- rep(1, length(aligned_mats))
  weights <- as.numeric(weights)
  names(weights) <- nms
  by_view <- mapply(function(Y, a) {
    a * sum((Y - M)^2, na.rm = TRUE)
  }, aligned_mats, as.list(weights))
  names(by_view) <- nms
  resid <- sum(by_view)
  total <- sum(mapply(function(Y, a) a * sum(Y^2, na.rm = TRUE), aligned_mats, as.list(weights)))
  row_w <- Reduce(`+`, Map(function(Y, a) {
    a * as.numeric(is.finite(Y[, 1L]))
  }, aligned_mats, as.list(weights)))
  cons <- sum(row_w * rowSums(M^2))
  by_entity <- Reduce(`+`, Map(function(Y, a) {
    dlt <- Y - M
    dlt[!is.finite(dlt)] <- 0
    a * rowSums(dlt^2)
  }, aligned_mats, as.list(weights)))
  list(
    total_energy = unname(total),
    consensus_energy = unname(cons),
    residual_energy = unname(resid),
    by_configuration = by_view,
    by_entity = by_entity,
    kind = "energy"
  )
}

#' Canonical presentation of a fit (principal axes of the consensus).
#'
#' This is a presentation, not part of estimation.
#'
#' @param object A `gpa_fit`.
#' @param ... Unused.
#' @export
canonicalize <- function(object, ...) {
  UseMethod("canonicalize")
}

#' @export
canonicalize.gpa_fit <- function(object, ...) {
  M <- object$consensus
  Mc <- .gproc_canonical_consensus(M, object$problem$transform$group)
  pair <- procrustes(M, Mc, transform = proc_orthogonal(
    if (identical(object$problem$transform$group, "SO")) "SO" else "O"
  ))
  out <- align_gauge(object, pair$transform$R)
  out$gauge$tied_axis_blocks <- tied_axis_blocks(Mc)
  out$gauge$canonicalization <- "principal_axes"
  out
}

#' Right-multiply a fit by a common orthogonal `Q`.
#'
#' @param x A `gpa_fit` or matrix.
#' @param reference An orthogonal matrix, or a `gpa_fit` whose consensus is the target.
#' @export
align_gauge <- function(x, reference) {
  UseMethod("align_gauge")
}

#' @export
align_gauge.gpa_fit <- function(x, reference) {
  grp <- if (identical(x$problem$transform$group, "SO")) "SO" else "O"
  if (inherits(reference, "gpa_fit")) {
    pair <- procrustes(
      consensus(x),
      consensus(reference),
      transform = proc_orthogonal(grp)
    )
    Q <- pair$transform$R
  } else {
    Q <- as.matrix(reference)
  }
  .gproc_check_gauge_Q(Q, x)
  x$consensus <- x$consensus %*% Q
  x$transformations <- lapply(x$transformations, function(tr) {
    tr$R <- tr$R %*% Q
    tr$t <- as.numeric(tr$t %*% Q)
    if (!is.null(tr$B)) {
      tr$B <- tr$B %*% Q
      if (is.function(tr$apply_warp) || tr$spec$family %in% c("tps", "lbw", "affine")) {
        d <- ncol(tr$B)
        tr <- .lbw_fitted(tr$spec, tr$B, list(controls = tr$extra$controls), d)
      }
    }
    tr
  })
  x$aligned_store <- lapply(names(x$aligned_store), function(nm) {
    v <- x$aligned_store[[nm]]
    if (inherits(v, "proc_aligned_view")) {
      v$transform <- x$transformations[[nm]]
    } else if (is.matrix(v)) {
      v <- v %*% Q
    }
    v
  })
  names(x$aligned_store) <- names(x$transformations)
  x
}

#' @noRd
.gproc_check_gauge_Q <- function(Q, fit) {
  d <- ncol(fit$consensus)
  if (!is.matrix(Q) || nrow(Q) != d || ncol(Q) != d) {
    .gproc_stop("invalid_problem", "Gauge matrix Q must be d by d.")
  }
  if (sqrt(sum((crossprod(Q) - diag(d))^2)) > 1e-8) {
    .gproc_stop("invalid_problem", "Gauge matrix Q must be orthogonal.")
  }
  grp <- fit$problem$transform$group
  if (identical(grp, "SO") && det(Q) < 0) {
    .gproc_stop("invalid_problem", "An SO(d) fit cannot be right-multiplied by a reflection.")
  }
  if (fit$problem$transform$family %in% c("affine", "lbw", "tps")) {
    lam <- fit$reference_covariance
    if (!is.null(lam) && length(unique(round(as.numeric(lam), 10))) > 1L) {
      L <- diag(as.numeric(lam), length(lam), length(lam))
      if (sqrt(sum((crossprod(Q, L %*% Q) - L)^2)) > 1e-8) {
        .gproc_stop(
          "invalid_problem",
          "Gauge Q must stabilize the declared reference covariance Lambda."
        )
      }
    }
  }
  invisible(Q)
}

#' @export
align_gauge.default <- function(x, reference) {
  as.matrix(x) %*% as.matrix(reference)
}
