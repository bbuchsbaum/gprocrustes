#' Fit a Procrustes / GPA problem.
#'
#' Two complete views with squared \(L_2\) use the exact pairwise kernel.
#' Three or more views are not implemented until the Gower solver lands.
#'
#' @param data A `proc_data` object, a list of matrices, or a single matrix
#'   (used as the source when `...` is unused).
#' @param transform Transform specification or shortcut.
#' @param gauge Gauge specification.
#' @param metric Weight / metric channels.
#' @param loss Residual loss. Only squared \(L_2\) is used by the pairwise kernel.
#' @param solver `"auto"` selects from the capability table.
#' @param control A `gpa_control` object.
#' @param anchor Optional view name treated as a fixed reference.
#' @param allow_disconnected If `TRUE`, disconnected overlap graphs are not fatal.
#' @return A `gpa_fit`.
#' @export
gpa <- function(data,
                transform = "similarity",
                gauge = proc_gauge(scale = "none", orientation = "free"),
                metric = proc_metric(),
                loss = proc_squared_l2(),
                solver = "auto",
                control = gpa_control(),
                anchor = NULL,
                allow_disconnected = FALSE) {
  spec <- as_proc_transform(transform)
  if (!inherits(gauge, "proc_gauge")) gauge <- proc_gauge()
  if (!inherits(metric, "proc_metric")) metric <- proc_metric()
  if (!inherits(loss, "proc_loss")) loss <- proc_squared_l2()
  if (!inherits(control, "gpa_control")) control <- gpa_control()

  data <- .gproc_as_data(data)
  ov <- overlap_graph(data)
  if (ov$n_components > 1L && !isTRUE(allow_disconnected)) {
    .gproc_stop(
      "disconnected_overlap_graph",
      sprintf("Overlap graph has %d connected components.", ov$n_components)
    )
  }
  alphas <- .gproc_config_weights(metric$configuration, names(data$views))
  if (sum(alphas) <= 0) {
    .gproc_stop("zero_total_configuration_weight", "All configuration weights are zero.")
  }

  compiled <- compile_proc_problem(data, spec, gauge, metric, loss)
  plan <- explain_solver(compiled)
  if (length(data$views) == 2L && identical(plan$solver, "pairwise_polar")) {
    fit_pair <- .gpa_pairwise(data, spec, metric, anchor)
    return(.gpa_from_pairwise(fit_pair, data, spec, gauge, metric, loss, compiled, plan, ov, alphas, anchor, control))
  }
  if (length(data$views) == 2L && identical(spec$family, "signed_permutation")) {
    fit_pair <- .gpa_pairwise(data, spec, metric, anchor)
    return(.gpa_from_pairwise(fit_pair, data, spec, gauge, metric, loss, compiled, plan, ov, alphas, anchor, control))
  }
  .gproc_stop(
    "invalid_problem",
    paste(
      "This compiled problem is not implemented yet.",
      "Milestone 1 ships the pairwise kernel only.",
      sprintf("Selected solver would be: %s.", plan$solver)
    )
  )
}

#' @keywords internal
.gproc_as_data <- function(data) {
  if (inherits(data, "proc_data")) {
    return(data)
  }
  proc_data(data)
}

#' @keywords internal
.gproc_config_weights <- function(alpha, nms) {
  if (is.null(alpha)) {
    return(stats::setNames(rep(1, length(nms)), nms))
  }
  alpha <- as.numeric(alpha)
  if (length(alpha) == 1L) alpha <- rep(alpha, length(nms))
  if (length(alpha) != length(nms)) {
    .gproc_stop("invalid_problem", "configuration weights must match the number of views.")
  }
  stats::setNames(alpha, nms)
}

#' @keywords internal
.gpa_pairwise <- function(data, spec, metric, anchor) {
  nms <- names(data$views)
  src <- nms[[1L]]
  tgt <- nms[[2L]]
  if (!is.null(anchor)) {
    if (!anchor %in% nms) {
      .gproc_stop("invalid_problem", "Unknown anchor view.")
    }
    tgt <- anchor
    src <- setdiff(nms, anchor)
    src <- src[[1L]]
  }
  X <- data$views[[src]]
  Y <- data$views[[tgt]]
  w <- .gproc_landmark_weights(metric$landmark, data, src, tgt)
  mask <- data$observed[[src]] & data$observed[[tgt]]
  if (!is.null(w)) {
    w <- w[mask]
  } else {
    w <- rep(1, sum(mask))
  }
  pair <- procrustes(X[mask, , drop = FALSE], Y[mask, , drop = FALSE],
                     transform = spec, weights = w)
  pair$source_name <- src
  pair$target_name <- tgt
  pair
}

#' @keywords internal
.gproc_landmark_weights <- function(landmark, data, src, tgt) {
  if (is.null(landmark)) {
    return(NULL)
  }
  if (is.list(landmark)) {
    w <- landmark[[src]]
    if (is.null(w)) w <- landmark[[tgt]]
    return(w)
  }
  as.numeric(landmark)
}

#' @keywords internal
.gpa_from_pairwise <- function(pair, data, spec, gauge, metric, loss, compiled, plan, ov, alphas, anchor, control) {
  nms <- names(data$views)
  transforms <- lapply(nms, function(nm) proc_identity(spec, ncol(data$views[[nm]])))
  names(transforms) <- nms
  transforms[[pair$source_name]] <- pair$transform
  aligned_lazy <- lapply(nms, function(nm) {
    proc_aligned_view(
      x = data$views[[nm]],
      transform = transforms[[nm]],
      row_map = data$row_map[[nm]],
      mask = data$observed[[nm]]
    )
  })
  names(aligned_lazy) <- nms
  mats <- lapply(aligned_lazy, as.matrix)
  A <- sum(alphas)
  M <- Reduce(`+`, Map(function(a, Y) a * Y, as.list(alphas), mats)) / A
  energy <- decompose_energy(mats, M, alphas)
  structure(
    list(
      problem = compiled,
      data_summary = list(
        n_views = length(nms),
        view_names = nms,
        n_entities = length(data$global_ids),
        target_dimension = ncol(M),
        correspondence = data$correspondence,
        overlap = ov
      ),
      consensus = M,
      transformations = transforms,
      gauge = list(spec = gauge, freedoms = .gproc_gauge_freedoms(spec, ncol(M), pair)),
      objective = pair$objective,
      decomposition = energy,
      history = data.frame(
        iteration = 0L,
        objective = pair$objective,
        relative_change = 0,
        accepted_acceleration = FALSE
      ),
      rank = list(
        effective_rank = pair$rank,
        singular_values = pair$singular_values,
        transform_unique = pair$transform_unique,
        objective_unique = pair$objective_unique,
        unidentified_subspace_dimension = pair$unidentified_subspace_dimension
      ),
      missing_support = ov,
      weights = list(configuration = alphas, landmark = metric$landmark),
      solver = plan$solver,
      backend = "base+Matrix",
      timings = list(),
      numerical_status = "converged",
      optimality_status = pair$optimality_status,
      certificate = list(status = "unavailable"),
      warnings = list(),
      aligned_store = if (identical(control$keep_aligned, "lazy")) aligned_lazy else mats,
      keep_aligned = control$keep_aligned,
      call = match.call(),
      version = .gproc_pkg_version()
    ),
    class = "gpa_fit"
  )
}

#' @keywords internal
.gproc_gauge_freedoms <- function(spec, d, pair) {
  list(
    common_rotation = identical(spec$family, "orthogonal") || identical(spec$scaling, "none"),
    common_translation = isTRUE(spec$translation),
    common_scale = identical(spec$scaling, "isotropic"),
    unidentified_subspace_dimension = pair$unidentified_subspace_dimension
  )
}

#' Compile a problem into planner input.
#'
#' @param data `proc_data`.
#' @param transform Transform spec.
#' @param gauge Gauge spec.
#' @param metric Metric.
#' @param loss Loss.
#' @export
compile_proc_problem <- function(data,
                                 transform = proc_similarity(),
                                 gauge = proc_gauge(scale = "none"),
                                 metric = proc_metric(),
                                 loss = proc_squared_l2()) {
  spec <- as_proc_transform(transform)
  data <- .gproc_as_data(data)
  structure(
    list(
      n_views = length(data$views),
      dimensions = vapply(data$views, ncol, integer(1)),
      transform = spec,
      gauge = gauge,
      metric = metric,
      loss = if (inherits(loss, "proc_loss")) loss else proc_squared_l2(),
      complete = all(vapply(data$observed, all, logical(1))),
      row_masked = any(vapply(data$observed, function(m) !all(m), logical(1))),
      cell_metric = !is.null(metric$cell) || !is.null(metric$precision)
    ),
    class = "proc_compiled_problem"
  )
}

#' Explain the selected solver for a compiled problem or fit.
#'
#' @param x A compiled problem or `gpa_fit`.
#' @export
explain_solver <- function(x) {
  if (inherits(x, "gpa_fit")) {
    return(explain_solver(x$problem))
  }
  spec <- x$transform
  pairwise_ok <- x$n_views == 2L &&
    identical(x$loss$family, "squared_l2") &&
    !isTRUE(x$cell_metric) &&
    spec$family %in% c("orthogonal", "similarity", "signed_permutation") &&
    length(unique(x$dimensions)) == 1L
  if (pairwise_ok) {
    reasons <- c(
      sprintf("transformation is %s", spec$group),
      "loss is squared L2",
      "two configurations",
      "no coordinate-specific precision"
    )
    plan <- c(
      "weighted sufficient-statistic moments",
      "thin SVD / polar factor of the cross-covariance",
      if (identical(spec$group, "SO")) "determinant correction for SO(d)",
      if (identical(spec$scaling, "isotropic")) "positive isotropic scale",
      if (isTRUE(spec$translation)) "analytic translation elimination"
    )
    return(structure(
      list(
        solver = if (identical(spec$family, "signed_permutation")) "signed_permutation" else "pairwise_polar",
        reasons = reasons,
        plan = plan,
        allowed_claim = "exact_closed_form"
      ),
      class = "proc_solver_plan"
    ))
  }
  structure(
    list(
      solver = "unimplemented",
      reasons = "no Milestone-1 kernel matches this compilation",
      plan = "see docs/spec/01-solver-guarantee-matrix.md",
      allowed_claim = "not_applicable"
    ),
    class = "proc_solver_plan"
  )
}

#' @export
print.proc_solver_plan <- function(x, ...) {
  cat(sprintf("Selected solver: %s\n", x$solver))
  cat("Reasons:\n")
  for (r in x$reasons) cat(sprintf("  %s\n", r))
  cat("Execution plan:\n")
  for (r in x$plan) cat(sprintf("  %s\n", r))
  cat(sprintf("Allowed optimality claim: %s\n", x$allowed_claim))
  invisible(x)
}
