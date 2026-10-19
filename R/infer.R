#' Declared statistical model for inference after fitting.
#'
#' \eqn{\Sigma_M} is the model covariance. It is not used by `gpa()` unless the
#' user also put it in `proc_metric(precision=)` as a superimposition metric.
#'
#' @param covariance A `proc_covariance`, typically `kind = "model"`.
#' @param resampling `"entities"` (paired row bootstrap) or `"configurations"`.
#' @param use_fitting_metric If `TRUE`, analytic / tangent \eqn{Q} is \eqn{\Sigma_S}
#'   from the fit. The default is to use `covariance` only.
#' @export
proc_shape_model <- function(covariance = proc_covariance(kind = "model"),
                             resampling = c("entities", "configurations"),
                             use_fitting_metric = FALSE) {
  resampling <- match.arg(resampling)
  if (!inherits(covariance, "proc_covariance")) {
    covariance <- proc_covariance(kind = "model")
  }
  structure(
    list(
      covariance = covariance,
      resampling = resampling,
      use_fitting_metric = isTRUE(use_fitting_metric)
    ),
    class = "proc_shape_model"
  )
}

#' Statistical inference for a fitted GPA problem.
#'
#' Fitting and inference are separate. Bootstrap replicates are gauge-aligned
#' to the original consensus before any coordinatewise summary (math section 39).
#' Entity bootstrap uses frequency weights for resampled entity contributions,
#' then evaluates each fitted transform on the fixed original landmarks; it
#' estimates conditional alignment uncertainty at those landmarks.
#' Analytic standard errors are a declared tangent-space approximation.
#'
#' @details Bootstrap currently supports complete squared-L2 orthogonal fits
#'   and unscaled similarities, with optional configuration and landmark weights.
#'   Bootstrap requests with robust losses, cell masks, precision metrics,
#'   scale constraints, or anchors are rejected explicitly. Configuration draws
#'   are refitted even when every sampled configuration is the same.
#'   Permutation inference tests smaller-than-null fitted loss by permuting
#'   coordinate values relative to fixed entity IDs; it assumes exchangeable rows.
#'
#' @param object A `gpa_fit`.
#' @param model A `proc_shape_model`.
#' @param method `"bootstrap"`, `"permutation"`, or `"analytic"`.
#' @param n_resamples Number of bootstrap or permutation draws.
#' @param ... Unused.
#' @return A `proc_inference` object.
#' @export
infer <- function(object,
                  model = proc_shape_model(),
                  method = c("bootstrap", "permutation", "analytic"),
                  n_resamples = 100L,
                  ...) {
  UseMethod("infer")
}

#' @export
infer.gpa_fit <- function(object,
                          model = proc_shape_model(),
                          method = c("bootstrap", "permutation", "analytic"),
                          n_resamples = 100L,
                          ...) {
  if (!inherits(model, "proc_shape_model")) {
    model <- proc_shape_model()
  }
  method <- match.arg(method)
  if (length(n_resamples) != 1L) {
    .gproc_stop("invalid_problem", "n_resamples must be a finite integer.")
  }
  n_input <- suppressWarnings(as.numeric(n_resamples)[1L])
  if (!is.finite(n_input) || n_input > .Machine$integer.max ||
      n_input < -.Machine$integer.max || n_input != as.integer(n_input)) {
    .gproc_stop("invalid_problem", "n_resamples must be a finite integer.")
  }
  n_resamples <- as.integer(n_input)
  min_draws <- if (identical(method, "bootstrap")) 2L else if (identical(method, "permutation")) 1L else 0L
  if (n_resamples < min_draws) {
    .gproc_stop(
      "invalid_problem",
      sprintf("n_resamples must be at least %d for %s inference.", min_draws, method)
    )
  }
  switch(
    method,
    bootstrap = .gproc_infer_bootstrap(object, model, n_resamples),
    permutation = .gproc_infer_permutation(object, model, n_resamples),
    analytic = .gproc_infer_analytic(object, model)
  )
}

#' Landmark-group cross-validation of a GPA fit (math section 35, ordinary GPA).
#'
#' Held-out landmarks are predicted from a refit that never saw them. The
#' cross-validation consensus is gauge-aligned to the full consensus before
#' the held-out discrepancy is scored.
#'
#' @param object A `gpa_fit`.
#' @param folds Number of landmark folds, or a list of integer index vectors.
#' @param level `"landmark"` (Bai-Bartoli held-out entities) or
#'   `"configuration"` (leave one view out).
#' @export
cross_validate <- function(object, folds = 5L, level = c("landmark", "configuration")) {
  if (!inherits(object, "gpa_fit")) {
    .gproc_stop("invalid_problem", "cross_validate() requires a gpa_fit.")
  }
  level <- match.arg(level)
  data <- .gproc_data_from_fit(object)
  if (identical(level, "configuration")) {
    return(.gproc_cv_configurations(object, data))
  }
  n <- length(data$global_ids)
  if (!isTRUE(object$problem$complete_covering)) {
    .gproc_stop("invalid_problem", "Landmark CV currently requires a complete covering.")
  }
  if (is.numeric(folds) && length(folds) == 1L) {
    k <- max(2L, as.integer(folds))
    perm <- sample.int(n)
    folds <- split(perm, rep(seq_len(k), length.out = n))
  }
  M0 <- consensus(object)
  spec <- object$problem$transform
  scores <- numeric(length(folds))
  names(scores) <- paste0("fold", seq_along(folds))
  for (g in seq_along(folds)) {
    hold <- sort(unique(as.integer(folds[[g]])))
    keep <- setdiff(seq_len(n), hold)
    if (!length(keep) || !length(hold)) {
      scores[[g]] <- NA_real_
      next
    }
    dat_g <- .gproc_subset_entities(data, keep)
    fit_g <- .gproc_refit(object, dat_g)
    pair <- procrustes(
      consensus(fit_g),
      M0[keep, , drop = FALSE],
      transform = proc_orthogonal(if (identical(spec$group, "SO")) "SO" else "O")
    )
    Q <- pair$transform$R
    pred <- 0
    for (nm in names(data$views)) {
      Xh <- data$views[[nm]][hold, , drop = FALSE]
      Yh <- apply_proc_transform(fit_g$transformations[[nm]], Xh) %*% Q
      pred <- pred + sum((Yh - M0[hold, , drop = FALSE])^2)
    }
    scores[[g]] <- pred
  }
  structure(
    list(
      scores = scores,
      total = sum(scores, na.rm = TRUE),
      folds = folds,
      level = "landmark",
      kind = "held_out_landmarks",
      note = paste(
        "Held-out discrepancy after gauge-aligning each training consensus",
        "to the full-data consensus. Not a claim about warp smoothness until LBW/TPS exists."
      )
    ),
    class = "proc_cv"
  )
}

#' @export
print.proc_inference <- function(x, ...) {
  cat(sprintf(
    "GPA inference (%s)\n  Resamples: %s\n  Gauge: %s\n  Assumptions: %s\n",
    x$method,
    x$n_resamples %||% "none",
    x$gauge,
    paste(x$assumptions, collapse = "; ")
  ))
  invisible(x)
}

#' @export
print.proc_cv <- function(x, ...) {
  cat(sprintf("Landmark CV total discrepancy: %.6g\n", x$total))
  invisible(x)
}

#' @noRd
.gproc_inference_covariance <- function(object, model) {
  if (is.null(model)) {
    return(NULL)
  }
  if (isTRUE(model$use_fitting_metric)) {
    return(.gproc_fitting_covariance(object$problem$metric))
  }
  covar <- model$covariance
  if (is.null(covar) || !inherits(covar, "proc_covariance")) {
    return(NULL)
  }
  covar
}

#' @noRd
.gproc_cv_configurations <- function(object, data) {
  nms <- names(data$views)
  if (length(nms) < 3L) {
    .gproc_stop("invalid_problem", "Configuration CV needs at least three views.")
  }
  scores <- numeric(length(nms))
  names(scores) <- nms
  for (nm in nms) {
    rest <- setdiff(nms, nm)
    dat_g <- proc_data(
      data$views[rest],
      observed = data$observed[rest],
      cells = if (!is.null(data$cells)) data$cells[rest] else NULL
    )
    fit_g <- .gproc_refit(object, dat_g)
    pred <- predict(fit_g, data$views[[nm]])
    scores[[nm]] <- sum((pred - consensus(fit_g))^2)
  }
  structure(
    list(
      scores = scores,
      total = sum(scores),
      folds = nms,
      level = "configuration",
      kind = "leave_one_configuration_out",
      note = paste(
        "Each view is aligned to the consensus of the others.",
        "Not a claim about warp smoothness until LBW/TPS exists."
      )
    ),
    class = "proc_cv"
  )
}

#' @noRd
.gproc_data_from_fit <- function(object) {
  if (!is.null(object$raw_data) && inherits(object$raw_data, "proc_data")) {
    return(object$raw_data)
  }
  av <- object$aligned_store
  views <- lapply(av, function(v) {
    if (inherits(v, "proc_aligned_view")) v$x else as.matrix(v)
  })
  observed <- lapply(av, function(v) {
    if (inherits(v, "proc_aligned_view")) v$mask else NULL
  })
  cells <- NULL
  if (!is.null(object$cell_masks)) {
    cells <- lapply(names(views), function(nm) {
      C <- object$cell_masks[[nm]]
      map <- if (inherits(av[[nm]], "proc_aligned_view")) av[[nm]]$row_map else seq_len(nrow(views[[nm]]))
      if (is.null(C)) NULL else C[map, , drop = FALSE]
    })
    names(cells) <- names(views)
  }
  ids <- NULL
  if (!isTRUE(object$problem$complete_covering)) {
    ids <- lapply(av, function(v) {
      if (inherits(v, "proc_aligned_view") && !is.null(v$row_map)) as.character(v$row_map) else NULL
    })
  }
  proc_data(views, ids = ids, observed = observed, cells = cells)
}

#' @noRd
.gproc_subset_entities <- function(data, keep) {
  views <- lapply(data$views, function(X) X[keep, , drop = FALSE])
  observed <- lapply(data$observed, function(m) m[keep])
  cells <- if (!is.null(data$cells)) {
    lapply(data$cells, function(C) if (is.null(C)) NULL else C[keep, , drop = FALSE])
  } else {
    NULL
  }
  proc_data(views, observed = observed, cells = cells)
}

#' @noRd
.gproc_refit <- function(object, data, metric = NULL) {
  spec <- object$problem$transform
  gauge <- object$gauge$spec
  if (!inherits(gauge, "proc_gauge")) gauge <- proc_gauge()
  metric <- metric %||% object$problem$metric
  if (!inherits(metric, "proc_metric")) metric <- proc_metric()
  loss <- object$problem$loss
  if (!inherits(loss, "proc_loss")) loss <- proc_squared_l2()
  gpa(
    data,
    transform = spec,
    gauge = gauge,
    metric = metric,
    loss = loss,
    solver = "auto",
    control = gpa_control(
      tolerance = 1e-10,
      max_iterations = 500L,
      certify = "never",
      accelerate = FALSE,
      init = "sequential"
    ),
    anchor = object$problem$anchor
  )
}

#' @noRd
.gproc_align_consensus <- function(M_b, M0, spec) {
  grp <- if (identical(spec$group, "SO")) "SO" else "O"
  translate <- .gproc_eliminates_translation(spec)
  pair <- procrustes(
    M_b,
    M0,
    transform = proc_similarity(group = grp, translation = translate, scaling = "none")
  )
  list(M = apply_proc_transform(pair$transform, M_b), R = pair$transform$R,
       t = pair$transform$t)
}

#' @noRd
.gproc_bootstrap_metric <- function(object, data, counts) {
  metric <- object$problem$metric
  if (!inherits(metric, "proc_metric")) metric <- proc_metric()
  n <- length(data$global_ids)
  metric$landmark <- lapply(names(data$views), function(nm) {
    map <- data$row_map[[nm]]
    base <- .gproc_local_landmark_weights(metric$landmark, map, n, nm)
    if (is.null(base)) base <- rep(1, length(map))
    as.numeric(base) * counts[map]
  })
  names(metric$landmark) <- names(data$views)
  metric
}

#' @noRd
.gproc_full_original_consensus <- function(object, data, fit_b) {
  metric <- object$problem$metric
  if (!inherits(metric, "proc_metric")) metric <- proc_metric()
  alphas <- .gproc_config_weights(metric$configuration, names(data$views))
  n <- length(data$global_ids)
  d <- ncol(consensus(fit_b))
  num <- matrix(0, n, d)
  den <- numeric(n)
  for (nm in names(data$views)) {
    map <- data$row_map[[nm]]
    Y <- apply_proc_transform(transformations(fit_b)[[nm]], data$views[[nm]])
    keep <- data$observed[[nm]]
    ww <- .gproc_local_landmark_weights(metric$landmark, map, n, nm)
    if (is.null(ww)) ww <- rep(1, length(map))
    ww <- alphas[[nm]] * as.numeric(ww)
    use <- keep & is.finite(ww) & ww > 0
    if (!any(use)) next
    gi <- map[use]
    num[gi, ] <- num[gi, , drop = FALSE] + Y[use, , drop = FALSE] * ww[use]
    den[gi] <- den[gi] + ww[use]
  }
  out <- num
  use <- den > 0
  out[use, ] <- num[use, , drop = FALSE] / den[use]
  out[!use, ] <- 0
  out
}

#' @noRd
.gproc_bootstrap_supported <- function(object, data, configurations = FALSE) {
  compiled <- object$problem
  if (!compiled$transform$family %in% c("orthogonal", "similarity")) {
    .gproc_stop(
      "invalid_problem",
      "Bootstrap currently supports orthogonal or unscaled similarity transforms only."
    )
  }
  if (!identical(compiled$loss$family, "squared_l2")) {
    .gproc_stop("invalid_problem", "Bootstrap currently supports squared-L2 fits only.")
  }
  if (isTRUE(compiled$cell_masked) || isTRUE(compiled$cell_metric)) {
    .gproc_stop("invalid_problem", "Bootstrap with cell masks or coordinate precision is not implemented.")
  }
  if (!is.null(compiled$metric$precision)) {
    .gproc_stop("invalid_problem", "Bootstrap with a precision metric is not implemented.")
  }
  if (!identical(compiled$transform$scaling, "none") ||
      !identical(compiled$gauge$scale, "none")) {
    .gproc_stop("invalid_problem", "Bootstrap with scale constraints is not implemented.")
  }
  if (!isTRUE(compiled$complete_covering)) {
    .gproc_stop("invalid_problem", "Entity bootstrap currently requires a complete covering.")
  }
  if (!is.null(compiled$anchor)) {
    .gproc_stop("invalid_problem", "Bootstrap with an anchor is not implemented.")
  }
  invisible(data)
}

#' @noRd
.gproc_resampled_configurations <- function(data, take, names_b) {
  proc_data(
    stats::setNames(data$views[take], names_b),
    ids = stats::setNames(data$ids[take], names_b),
    observed = stats::setNames(data$observed[take], names_b),
    cells = stats::setNames(data$cells[take], names_b),
    metadata = if (is.null(data$metadata)) NULL else stats::setNames(data$metadata[take], names_b)
  )
}

#' @noRd
.gproc_resampled_configuration_metric <- function(object, take, names_b) {
  metric <- object$problem$metric
  if (!inherits(metric, "proc_metric")) metric <- proc_metric()
  alpha <- object$weights$configuration %||% .gproc_config_weights(metric$configuration, names(object$raw_data$views))
  metric$configuration <- unname(alpha[take])
  if (is.list(metric$landmark)) {
    metric$landmark <- stats::setNames(metric$landmark[take], names_b)
  }
  metric
}

#' @noRd
.gproc_infer_bootstrap <- function(object, model, B) {
  data <- .gproc_data_from_fit(object)
  n <- length(data$global_ids)
  spec <- object$problem$transform
  M0 <- consensus(object)
  reps <- array(NA_real_, c(nrow(M0), ncol(M0), B))
  if (identical(model$resampling, "configurations")) {
    .gproc_bootstrap_supported(object, data, configurations = TRUE)
    nms <- names(data$views)
    for (b in seq_len(B)) {
      take <- sample.int(length(nms), length(nms), replace = TRUE)
      names_b <- paste0("boot", seq_along(take))
      dat_b <- .gproc_resampled_configurations(data, take, names_b)
      metric_b <- .gproc_resampled_configuration_metric(object, take, names_b)
      fit_b <- .gproc_refit(object, dat_b, metric_b)
      reps[, , b] <- .gproc_align_consensus(consensus(fit_b), M0, spec)$M
    }
    assumption <- "configurations resampled with replacement and refit; then gauge-aligned"
  } else {
    .gproc_bootstrap_supported(object, data)
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      counts <- tabulate(idx, n)
      metric_b <- .gproc_bootstrap_metric(object, data, counts)
      fit_b <- .gproc_refit(object, data, metric_b)
      M_b <- .gproc_full_original_consensus(object, data, fit_b)
      reps[, , b] <- if (is.null(object$problem$anchor)) {
        .gproc_align_consensus(M_b, M0, spec)$M
      } else {
        M_b
      }
    }
    assumption <- paste(
      "entities resampled through multiplicity weights on the original data; then gauge-aligned",
      "conditional alignment uncertainty at fixed original landmarks"
    )
  }
  mu <- apply(reps, c(1L, 2L), mean)
  se <- apply(reps, c(1L, 2L), stats::sd)
  gauge_note <- if (is.null(object$problem$anchor)) {
    "each replicate aligned to the original consensus"
  } else {
    "anchor fixes the replicate frame; no common gauge alignment"
  }
  structure(
    list(
      method = "bootstrap",
      model = model,
      estimate = M0,
      aligned_mean = mu,
      se = se,
      replicates = reps,
      n_resamples = B,
      gauge = gauge_note,
      assumptions = c(
        assumption,
        "coordinatewise summaries are after gauge alignment only",
        "repeated consensus eigenvalues make individual axes unstable"
      )
    ),
    class = "proc_inference"
  )
}

#' @noRd
.gproc_infer_permutation <- function(object, model, B) {
  data <- .gproc_data_from_fit(object)
  n <- length(data$global_ids)
  if (!isTRUE(object$problem$complete_covering)) {
    .gproc_stop("invalid_problem", "Permutation currently requires a complete covering.")
  }
  nms <- names(data$views)
  T0 <- object$objective
  Tperm <- numeric(B)
  for (b in seq_len(B)) {
    views_p <- data$views
    for (nm in nms[-1L]) {
      X <- views_p[[nm]][sample.int(n), , drop = FALSE]
      rownames(X) <- rownames(views_p[[nm]])
      views_p[[nm]] <- X
    }
    dat_p <- proc_data(
      views_p,
      ids = data$ids,
      observed = data$observed,
      cells = data$cells,
      metadata = data$metadata
    )
    fit_p <- .gproc_refit(object, dat_p)
    Tperm[[b]] <- fit_p$objective
  }
  p <- (1 + sum(Tperm <= T0)) / (1 + B)
  structure(
    list(
      method = "permutation",
      model = model,
      estimate = consensus(object),
      statistic = T0,
      null_statistics = Tperm,
      p_value = p,
      n_resamples = B,
      gauge = "not applicable; compares fitted objective",
      assumptions = c(
        "independent value permutations destroy correspondence while IDs remain fixed",
        "the test statistic is the fitted objective after a refit",
        "this is not a test of a Gaussian shape model"
      )
    ),
    class = "proc_inference"
  )
}

#' @noRd
.gproc_infer_analytic <- function(object, model) {
  M <- consensus(object)
  spec <- object$problem$transform
  tan <- tangent_coordinates(object, model = model)
  n <- nrow(M)
  d <- ncol(M)
  nuis <- attr(tan, "nuisance_dimension") %||% 0L
  rss <- sum(vapply(tan, function(Z) sum(Z^2), numeric(1)))
  k <- length(tan)
  shape_dim <- max(1L, n * d - as.integer(nuis))
  df <- max(1, (k - 1L) * shape_dim)
  sigma2 <- rss / df
  se <- matrix(sqrt(pmax(sigma2, 0)), n, d)
  structure(
    list(
      method = "analytic",
      model = model,
      estimate = M,
      se = se,
      sigma2 = sigma2,
      df = df,
      tangent = tan,
      n_resamples = NULL,
      gauge = "horizontal projection at the fitted consensus",
      assumptions = c(
        "local tangent chart; shape space is not affine",
        "residuals after P_H are treated as isotropic Gaussian in that chart",
        "this is a declared approximation, not a certified sampling distribution",
        if (isTRUE(model$use_fitting_metric))
          "Q is the superimposition metric (user declaration)"
        else
          "Q is the model covariance, not automatically Sigma_S"
      )
    ),
    class = "proc_inference"
  )
}
