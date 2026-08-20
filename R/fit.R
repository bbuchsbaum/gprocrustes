#' @export
print.gpa_fit <- function(x, ...) {
  ds <- x$data_summary
  miss <- NA_real_
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
      "  Numerical status: %s\n",
      "  Optimality: %s\n",
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
    x$numerical_status,
    x$optimality_status,
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
#' @param object A `gpa_fit`.
#' @param ... Unused.
#' @export
consensus <- function(object, ...) {
  UseMethod("consensus")
}

#' @export
consensus.gpa_fit <- function(object, ...) {
  object$consensus
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
  object$aligned_store
}

#' @export
fitted.gpa_fit <- function(object, ...) {
  lapply(aligned(object), as.matrix)
}

#' @export
residuals.gpa_fit <- function(object, level = c("landmark", "configuration"), ...) {
  level <- match.arg(level)
  M <- consensus(object)
  Y <- fitted(object)
  res <- lapply(Y, function(Yi) Yi - M)
  if (identical(level, "configuration")) {
    return(vapply(res, function(E) sqrt(sum(E^2)), numeric(1)))
  }
  res
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
    numerical_status = object$numerical_status,
    optimality_status = object$optimality_status
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
decompose <- function(object, ...) {
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
  A <- sum(weights)
  total <- sum(mapply(function(Y, a) a * sum(Y^2), aligned_mats, as.list(weights)))
  cons <- A * sum(M^2)
  by_view <- mapply(function(Y, a) a * sum((Y - M)^2), aligned_mats, as.list(weights))
  names(by_view) <- nms
  resid <- sum(by_view)
  by_entity <- Reduce(`+`, Map(function(Y, a) {
    a * rowSums((Y - M)^2)
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
  M <- consensus(object)
  sv <- svd(M, nu = ncol(M), nv = ncol(M))
  Q <- sv$v
  if (det(Q) < 0 && identical(object$problem$transform$group, "SO")) {
    Q[, ncol(Q)] <- -Q[, ncol(Q)]
  }
  align_gauge(object, Q)
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
  if (inherits(reference, "gpa_fit")) {
    pair <- procrustes(consensus(x), consensus(reference), transform = proc_orthogonal("O"))
    Q <- pair$transform$R
  } else {
    Q <- as.matrix(reference)
  }
  x$consensus <- x$consensus %*% Q
  x$transformations <- lapply(x$transformations, function(tr) {
    tr$R <- tr$R %*% Q
    tr$t <- as.numeric(tr$t %*% Q)
    tr
  })
  if (identical(x$keep_aligned, "lazy")) {
    x$aligned_store <- lapply(names(x$aligned_store), function(nm) {
      v <- x$aligned_store[[nm]]
      v$transform <- x$transformations[[nm]]
      v
    })
    names(x$aligned_store) <- names(x$transformations)
  }
  x
}

#' @export
align_gauge.default <- function(x, reference) {
  as.matrix(x) %*% as.matrix(reference)
}
