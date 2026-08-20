#' Declared statistical model for inference after fitting.
#'
#' \(\Sigma_M\) is the model covariance. It is not used by `gpa()` unless the
#' user also put it in `proc_metric(precision=)` as a superimposition metric.
#'
#' @param covariance A `proc_covariance`, typically `kind = "model"`.
#' @param resampling `"entities"` (paired row bootstrap) or `"configurations"`.
#' @param use_fitting_metric If `TRUE`, analytic / tangent \(Q\) is \(\Sigma_S\)
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
#' to the original consensus before any coordinatewise summary (math §39).
#' Analytic standard errors are a declared tangent-space approximation.
#'
#' @param object A `gpa_fit`.
#' @param model A `proc_shape_model`.
#' @param method `"bootstrap"`, `"permutation"`, or `"analytic"`.
#' @param n_resamples Number of bootstrap or permutation draws.
#' @param ... Unused.
#' @return A `proc_inference` object.
#' @export
infer <- function(object, ...) {
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
  n_resamples <- as.integer(n_resamples)[1L]
  switch(
    method,
    bootstrap = .gproc_infer_bootstrap(object, model, n_resamples),
    permutation = .gproc_infer_permutation(object, model, n_resamples),
    analytic = .gproc_infer_analytic(object, model)
  )
}

#' Landmark-group cross-validation of a GPA fit (math §35, ordinary GPA).
#'
#' Held-out landmarks are predicted from a refit that never saw them. The
#' cross-validation consensus is gauge-aligned to the full consensus before
#' the held-out discrepancy is scored.
#'
#' @param object A `gpa_fit`.
#' @param folds Number of landmark folds, or a list of integer index vectors.
#' @param level `"landmark"` (Bai–Bartoli held-out entities) or
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
  av <- aligned(object)
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
.gproc_refit <- function(object, data) {
  spec <- object$problem$transform
  gauge <- object$gauge$spec
  if (!inherits(gauge, "proc_gauge")) gauge <- proc_gauge()
  metric <- object$problem$metric
  if (!inherits(metric, "proc_metric")) metric <- proc_metric()
  loss <- object$problem$loss
  if (!inherits(loss, "proc_loss")) loss <- proc_squared_l2()
  solver <- object$solver
  if (!solver %in% c("gpm", "gower_bcd", "irls", "mm", "pairwise_polar", "auto")) {
    solver <- "auto"
  }
  gpa(
    data,
    transform = spec,
    gauge = gauge,
    metric = metric,
    loss = loss,
    solver = solver,
    control = gpa_control(
      tolerance = 1e-6,
      max_iterations = 80L,
      certify = "never",
      accelerate = FALSE,
      init = "sequential"
    )
  )
}

#' @noRd
.gproc_align_consensus <- function(M_b, M0, spec) {
  grp <- if (identical(spec$group, "SO")) "SO" else "O"
  pair <- procrustes(M_b, M0, transform = proc_orthogonal(grp))
  list(M = M_b %*% pair$transform$R, R = pair$transform$R)
}

#' @noRd
.gproc_infer_bootstrap <- function(object, model, B) {
  data <- .gproc_data_from_fit(object)
  n <- length(data$global_ids)
  spec <- object$problem$transform
  M0 <- consensus(object)
  if (!isTRUE(object$problem$complete_covering)) {
    .gproc_stop("invalid_problem", "Entity bootstrap currently requires a complete covering.")
  }
  reps <- array(NA_real_, c(nrow(M0), ncol(M0), B))
  if (identical(model$resampling, "configurations")) {
    nms <- names(data$views)
    for (b in seq_len(B)) {
      take <- sample(nms, length(nms), replace = TRUE)
      views_b <- data$views[take]
      names(views_b) <- paste0("boot", seq_along(views_b))
      if (length(unique(take)) < 2L) {
        reps[, , b] <- M0
        next
      }
      fit_b <- .gproc_refit(object, proc_data(views_b))
      reps[, , b] <- .gproc_align_consensus(consensus(fit_b), M0, spec)$M
    }
    assumption <- "configurations resampled with replacement; then gauge-aligned"
  } else {
    for (b in seq_len(B)) {
      idx <- sample.int(n, n, replace = TRUE)
      dat_b <- .gproc_subset_entities(data, idx)
      fit_b <- .gproc_refit(object, dat_b)
      reps[, , b] <- .gproc_align_consensus(consensus(fit_b), M0, spec)$M
    }
    assumption <- "entities resampled together across views; then gauge-aligned"
  }
  mu <- apply(reps, c(1L, 2L), mean)
  se <- apply(reps, c(1L, 2L), stats::sd)
  structure(
    list(
      method = "bootstrap",
      model = model,
      estimate = M0,
      aligned_mean = mu,
      se = se,
      replicates = reps,
      n_resamples = B,
      gauge = "each replicate polar-aligned to the original consensus",
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
  T0 <- object$decomposition$residual_energy %||% object$objective
  Tperm <- numeric(B)
  for (b in seq_len(B)) {
    views_p <- data$views
    for (nm in nms[-1L]) {
      views_p[[nm]] <- views_p[[nm]][sample.int(n), , drop = FALSE]
    }
    fit_p <- .gproc_refit(object, proc_data(views_p))
    Tperm[[b]] <- fit_p$decomposition$residual_energy %||% fit_p$objective
  }
  p <- (1 + sum(Tperm >= T0)) / (1 + B)
  structure(
    list(
      method = "permutation",
      model = model,
      estimate = consensus(object),
      statistic = T0,
      null_statistics = Tperm,
      p_value = p,
      n_resamples = B,
      gauge = "not applicable; compares residual energy",
      assumptions = c(
        "independent row permutations destroy correspondence in every view but the first",
        "the test statistic is residual energy after a refit",
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
          "Q is the model covariance, not automatically Σ_S"
      )
    ),
    class = "proc_inference"
  )
}
