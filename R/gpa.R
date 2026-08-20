#' Fit a Procrustes / GPA problem.
#'
#' Unanchored orthogonal problems use the matrix-free generalized power
#' method. Similarity problems use consensus-first Gower BCD. Two views
#' with an explicit `anchor` (or a signed-permutation group) use the exact
#' pairwise kernel.
#'
#' @param data A `proc_data` object, a list of matrices, or a single matrix
#'   (used as the source when `...` is unused).
#' @param transform Transform specification or shortcut.
#' @param gauge Gauge specification.
#' @param metric Weight / metric channels.
#' @param loss Residual loss. Huber / Tukey use landmark-vector IRLS.
#' @param solver `"auto"` selects from the capability table. `"gpm"`,
#'   `"gower_bcd"`, `"irls"`, `"mm"`, `"lbw_eigen"`, and `"pairwise_polar"`
#'   force engines.
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
  metric_use <- .gproc_apply_goodall_metric(metric, data)
  use_pairwise <- .gpa_use_pairwise(data, spec, solver, plan, anchor)
  if (use_pairwise) {
    pair_plan <- .gpa_pairwise_plan(spec)
    fit_pair <- .gpa_pairwise(data, spec, metric_use, anchor)
    return(.gpa_from_pairwise(fit_pair, data, spec, gauge, metric_use, loss, compiled,
                              pair_plan, ov, alphas, anchor, control))
  }
  if (identical(solver, "irls") || (identical(solver, "auto") && identical(plan$solver, "irls"))) {
    raw <- .irls_gpa(data, spec, gauge, metric_use, alphas, control, anchor, loss)
    plan_use <- if (identical(plan$solver, "irls")) plan else {
      structure(list(solver = "irls", reasons = "solver='irls'",
                     plan = "landmark-vector IRLS around Gower BCD",
                     allowed_claim = "blockwise_stationary"),
                class = "proc_solver_plan")
    }
    fit <- .gpa_from_gower(raw, data, spec, gauge, metric_use, loss, compiled, plan_use,
                           ov, alphas, anchor, control)
    if (!is.null(raw$landmark_weights)) {
      fit$weights$robust <- raw$landmark_weights
      fit$weights$landmark <- raw$landmark_weights
    }
    return(fit)
  }
  if (identical(solver, "mm") || (identical(solver, "auto") && identical(plan$solver, "mm"))) {
    raw <- .mm_gpa(data, spec, gauge, metric_use, alphas, control, anchor, loss)
    plan_use <- if (identical(plan$solver, "mm")) plan else {
      structure(list(solver = "mm", reasons = "solver='mm'",
                     plan = "majorization with filled unobserved coordinates",
                     allowed_claim = "first_order_stationary"),
                class = "proc_solver_plan")
    }
    fit <- .gpa_from_gower(raw, data, spec, gauge, metric_use, loss, compiled, plan_use,
                           ov, alphas, anchor, control)
    fit$cell_masks <- raw$cell_masks
    return(fit)
  }
  if (identical(solver, "lbw_eigen") || (identical(solver, "auto") && identical(plan$solver, "lbw_eigen"))) {
    raw <- .lbw_gpa(data, spec, gauge, metric_use, alphas, control, anchor)
    plan_use <- if (identical(plan$solver, "lbw_eigen")) plan else {
      structure(
        list(
          solver = "lbw_eigen",
          reasons = "solver='lbw_eigen'",
          plan = "variable projection plus constrained eigenproblem",
          allowed_claim = "exact_closed_form for the stated LBW formulation"
        ),
        class = "proc_solver_plan"
      )
    }
    return(.gpa_from_lbw(raw, data, spec, gauge, metric_use, loss, compiled, plan_use,
                         ov, alphas, anchor, control))
  }
  if (identical(solver, "gpm") || (identical(solver, "auto") && identical(plan$solver, "gpm"))) {
    if (!identical(plan$solver, "gpm")) {
      .gproc_stop("invalid_problem", "gpm does not apply to this compiled problem.")
    }
    raw <- .gpm_gpa(data, spec, gauge, metric_use, alphas, control, anchor)
    fit <- .gpa_from_gpm(raw, data, spec, gauge, metric_use, loss, compiled, plan,
                         ov, alphas, anchor, control)
    if (identical(gauge$orientation, "principal")) {
      fit <- canonicalize(fit)
    }
    return(fit)
  }
  if (identical(solver, "gower_bcd") || (identical(solver, "auto") && identical(plan$solver, "gower_bcd"))) {
    raw <- .gower_gpa(data, spec, gauge, metric_use, alphas, control, anchor)
    gower_plan <- if (identical(plan$solver, "gower_bcd")) plan else explain_solver(compiled)
    if (!identical(gower_plan$solver, "gower_bcd") && identical(solver, "gower_bcd")) {
      # Orthogonal problems prefer GPM; Gower BCD remains a legitimate forced engine.
      gower_plan <- structure(
        list(
          solver = "gower_bcd",
          reasons = c("solver='gower_bcd' requested", sprintf("transformation is %s", spec$group)),
          plan = "block-coordinate descent against the consensus",
          allowed_claim = "blockwise_stationary"
        ),
        class = "proc_solver_plan"
      )
    }
    fit <- .gpa_from_gower(raw, data, spec, gauge, metric_use, loss, compiled, gower_plan,
                           ov, alphas, anchor, control)
    if (identical(gauge$orientation, "principal")) {
      fit <- canonicalize(fit)
    }
    return(fit)
  }
  .gproc_stop(
    "invalid_problem",
    paste(
      "This compiled problem is not implemented yet.",
      sprintf("Selected solver would be: %s.", plan$solver)
    )
  )
}

#' @keywords internal
.gpa_use_pairwise <- function(data, spec, solver, plan, anchor) {
  k <- length(data$views)
  if (identical(solver, "pairwise_polar") || identical(solver, "signed_permutation")) {
    if (k != 2L) {
      .gproc_stop("invalid_problem", "The pairwise kernel requires exactly two views.")
    }
    return(TRUE)
  }
  if (k == 2L && identical(spec$family, "signed_permutation")) {
    return(TRUE)
  }
  if (k == 2L && !is.null(anchor) && plan$solver %in% c("gower_bcd", "pairwise_polar", "gpm")) {
    return(TRUE)
  }
  FALSE
}

#' @keywords internal
.gpa_pairwise_plan <- function(spec) {
  structure(
    list(
      solver = if (identical(spec$family, "signed_permutation")) "signed_permutation" else "pairwise_polar",
      reasons = c(
        sprintf("transformation is %s", spec$group),
        "loss is squared L2",
        "two configurations",
        "explicit anchor or signed-permutation kernel"
      ),
      plan = c(
        "weighted sufficient-statistic moments",
        "thin SVD / polar factor of the cross-covariance",
        if (identical(spec$group, "SO")) "determinant correction for SO(d)",
        if (identical(spec$scaling, "isotropic")) "positive isotropic scale",
        if (isTRUE(spec$translation)) "analytic translation elimination"
      ),
      allowed_claim = "exact_closed_form"
    ),
    class = "proc_solver_plan"
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
      solver_plan = plan,
      call = match.call(),
      version = .gproc_pkg_version()
    ),
    class = "gpa_fit"
  )
}

#' @keywords internal
.gpa_from_gower <- function(raw, data, spec, gauge, metric, loss, compiled, plan,
                            ov, alphas, anchor, control) {
  nms <- names(data$views)
  transforms <- raw$transforms
  aligned_lazy <- lapply(nms, function(nm) {
    proc_aligned_view(
      x = data$views[[nm]],
      transform = transforms[[nm]],
      row_map = data$row_map[[nm]],
      mask = data$observed[[nm]]
    )
  })
  names(aligned_lazy) <- nms
  global <- raw$aligned
  if (is.null(global) || !length(global)) {
    global <- lapply(aligned_lazy, function(v) {
      .gproc_scatter_aligned(v, length(data$global_ids), ncol(raw$consensus))
    })
  }
  energy <- decompose_energy(global, raw$consensus, alphas)
  sv <- .gproc_svd(raw$consensus, nu = 0L, nv = min(5L, ncol(raw$consensus)))
  warnings <- list()
  if (isTRUE(raw$scale_overridden) && !is.null(raw$scale_note)) {
    warnings <- c(warnings, list(list(code = "scale_constraint_applied", message = raw$scale_note)))
  }
  landmark_w <- raw$landmark_weights %||% metric$landmark
  fit <- structure(
    list(
      problem = compiled,
      data_summary = list(
        n_views = length(nms),
        view_names = nms,
        n_entities = length(data$global_ids),
        target_dimension = ncol(raw$consensus),
        correspondence = data$correspondence,
        overlap = ov
      ),
      consensus = raw$consensus,
      transformations = transforms,
      gauge = list(
        spec = gauge,
        scale_mode = raw$scale_mode,
        scale_note = raw$scale_note,
        freedoms = .gproc_gauge_freedoms(spec, ncol(raw$consensus), NULL)
      ),
      objective = raw$objective,
      decomposition = energy,
      history = raw$history,
      rank = list(
        effective_rank = sum(sv$d > 1e-10 * max(sv$d, 0)),
        singular_values = sv$d,
        transform_unique = NA,
        objective_unique = NA,
        unidentified_subspace_dimension = NA_integer_,
        stationarity = raw$stationarity
      ),
      missing_support = ov,
      weights = list(configuration = alphas, landmark = metric$landmark),
      solver = plan$solver,
      backend = "base+Matrix",
      timings = list(
        elapsed = if (nrow(raw$history)) raw$history$elapsed_time[[nrow(raw$history)]] else NA_real_,
        iterations = nrow(raw$history)
      ),
      numerical_status = raw$numerical_status,
      optimality_status = raw$optimality_status,
      certificate = list(
        status = "unavailable",
        reason = "Gower BCD is monotone / blockwise stationary, not a Ling dual certificate."
      ),
      warnings = warnings,
      aligned_store = if (identical(control$keep_aligned, "lazy")) aligned_lazy else lapply(aligned_lazy, as.matrix),
      keep_aligned = control$keep_aligned,
      solver_plan = plan,
      init = raw$init,
      call = NULL,
      version = .gproc_pkg_version()
    ),
    class = "gpa_fit"
  )
  fit$weights$landmark <- landmark_w
  if (identical(plan$solver, "lbw_eigen")) {
    fit$certificate$reason <- paste(
      "Global eigen-solution for the constrained reference-space LBW formulation;",
      "not a claim that deformable GPA is globally solved."
    )
    fit$backend <- "eigencore"
  }
  if (identical(plan$solver, "irls")) {
    fit$certificate$reason <- paste(
      "IRLS is monotone in the landmark-vector rho criterion;",
      "Ling's dual certificate is for squared L2 GOPP."
    )
  } else if (identical(plan$solver, "mm")) {
    fit$certificate$reason <- paste(
      "Cell-mask / anisotropic MM is first-order stationary;",
      "it is not an SVD or Ling certificate."
    )
  }
  fit
}

#' @keywords internal
.gpa_from_lbw <- function(raw, data, spec, gauge, metric, loss, compiled, plan,
                          ov, alphas, anchor, control) {
  fit <- .gpa_from_gower(raw, data, spec, gauge, metric, loss, compiled, plan,
                         ov, alphas, anchor, control)
  fit$certificate <- list(
    status = if (identical(raw$optimality_status, "exact_closed_form")) {
      "eigen_global_for_formulation"
    } else {
      "unavailable"
    },
    reason = if (identical(raw$optimality_status, "exact_closed_form")) {
      paste(
        "Global eigen-solution for the constrained reference-space LBW formulation;",
        "not a claim that deformable GPA is globally solved."
      )
    } else {
      "Free-translation residual is large; first-order claim only."
    },
    eigen_gap = raw$eigen_gap,
    free_translation = raw$free_translation
  )
  fit$reference_covariance <- raw$reference_covariance
  fit$folding <- raw$folding
  fit$data_objective <- raw$data_objective
  fit$penalty <- raw$penalty
  fit$backend <- "eigencore"
  fit$optimality_status <- raw$optimality_status
  fit
}

#' @keywords internal
.gpa_from_gpm <- function(raw, data, spec, gauge, metric, loss, compiled, plan,
                          ov, alphas, anchor, control) {
  fit <- .gpa_from_gower(raw, data, spec, gauge, metric, loss, compiled, plan,
                         ov, alphas, anchor, control)
  fit$certificate <- raw$certificate
  fit$backend <- raw$backend %||% "matrix_free"
  fit$optimality_status <- raw$optimality_status
  fit
}

#' @keywords internal
.gproc_gauge_freedoms <- function(spec, d, pair = NULL) {
  list(
    common_rotation = identical(spec$family, "orthogonal") || identical(spec$scaling, "none") ||
      identical(spec$family, "similarity"),
    common_translation = isTRUE(spec$translation),
    common_scale = identical(spec$scaling, "isotropic"),
    unidentified_subspace_dimension = if (is.null(pair)) NA_integer_ else pair$unidentified_subspace_dimension
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
  n_ent <- length(data$global_ids)
  n_rows <- vapply(data$views, nrow, integer(1))
  complete <- all(vapply(data$observed, all, logical(1)))
  structure(
    list(
      n_views = length(data$views),
      dimensions = vapply(data$views, ncol, integer(1)),
      n_rows = n_rows,
      n_entities = n_ent,
      transform = spec,
      gauge = gauge,
      metric = metric,
      loss = if (inherits(loss, "proc_loss")) loss else proc_squared_l2(),
      complete = complete,
      complete_covering = complete && all(n_rows == n_ent),
      row_masked = any(vapply(data$observed, function(m) !all(m), logical(1))),
      cell_masked = !is.null(data$cells) && any(vapply(data$cells, function(C) {
        !is.null(C) && !all(C)
      }, logical(1))),
      cell_metric = !is.null(metric$cell) ||
        .gproc_coordinate_anisotropic(.gproc_fitting_covariance(metric)),
      landmark_weighted = .gproc_landmark_weighted(metric, n_ent)
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
    if (!is.null(x$solver_plan)) {
      return(x$solver_plan)
    }
    return(explain_solver(x$problem))
  }
  spec <- x$transform
  if (spec$family %in% c("affine", "lbw", "tps") &&
      !identical(x$loss$family, "squared_l2")) {
    return(structure(
      list(
        solver = "unimplemented",
        reasons = "LBW eigen-solver is for squared reference-space L2 only",
        plan = "see docs/spec/01-solver-guarantee-matrix.md",
        allowed_claim = "not_applicable"
      ),
      class = "proc_solver_plan"
    ))
  }
  if (isTRUE(x$cell_masked) || isTRUE(x$cell_metric)) {
    return(structure(
      list(
        solver = "mm",
        reasons = c(
          if (isTRUE(x$cell_masked)) "cell masks differ by coordinate" else "coordinate-specific precision",
          "missing cells are not missing rows",
          "no ordinary SVD solution",
          if (identical(x$loss$family, "huber") || identical(x$loss$family, "tukey"))
            sprintf("landmark-vector %s is evaluated on observed cells", x$loss$family)
        ),
        plan = c(
          "fill unobserved target coordinates from the current residual",
          "solve a complete surrogate",
          "accept only a decrease of the true observed-data objective"
        ),
        allowed_claim = "first_order_stationary"
      ),
      class = "proc_solver_plan"
    ))
  }
  if (spec$family %in% c("affine", "lbw", "tps") &&
      identical(x$loss$family, "squared_l2") &&
      !isTRUE(x$cell_masked) &&
      !isTRUE(x$cell_metric)) {
    return(structure(
      list(
        solver = "lbw_eigen",
        reasons = c(
          sprintf("transformation is a linear-basis warp (%s)", spec$family),
          "loss is squared L2 in the reference space",
          "M^T M = Lambda is an explicit reference-covariance constraint",
          "free-translation will be checked at solve time"
        ),
        plan = c(
          "build Phi and L for each view",
          "eliminate B_i by a penalized smoother (no explicit inverse)",
          "form P = sum_i (I - H_i) or the partial-shape operator",
          "take the bottom d eigenvectors after the translation null",
          "scale axes by the declared Lambda^{1/2}"
        ),
        allowed_claim = "exact_closed_form only for this constrained reference-space formulation"
      ),
      class = "proc_solver_plan"
    ))
  }
  if (identical(x$loss$family, "huber") || identical(x$loss$family, "tukey")) {
    return(structure(
      list(
        solver = "irls",
        reasons = c(
          sprintf("loss is landmark-vector %s", x$loss$family),
          "IRLS weights act on residual vector lengths",
          if (isTRUE(x$loss$convex)) "Huber is convex" else "Tukey is nonconvex"
        ),
        plan = c(
          "weighted Gower BCD with current IRLS weights",
          "update w_ij = q(||Y_ij - M_j||)",
          "accept only a decrease of the true rho objective"
        ),
        allowed_claim = "blockwise_stationary"
      ),
      class = "proc_solver_plan"
    ))
  }
  signed_pair <- x$n_views == 2L &&
    identical(spec$family, "signed_permutation") &&
    identical(x$loss$family, "squared_l2") &&
    !isTRUE(x$cell_metric) &&
    length(unique(x$dimensions)) == 1L
  if (signed_pair) {
    return(.gpa_pairwise_plan(spec))
  }
  gpm_ok <- x$n_views >= 2L &&
    identical(x$loss$family, "squared_l2") &&
    !isTRUE(x$cell_metric) &&
    !isTRUE(x$landmark_weighted) &&
    identical(spec$family, "orthogonal") &&
    length(unique(x$dimensions)) == 1L &&
    isTRUE(x$complete_covering)
  if (gpm_ok) {
    return(structure(
      list(
        solver = "gpm",
        reasons = c(
          sprintf("transformation is %s", spec$group),
          "loss is squared L2",
          "configurations are complete after implicit centering",
          "no coordinate-specific precision",
          sprintf("%d configurations; GOPP / generalized power method", x$n_views)
        ),
        plan = c(
          "spectral initialization via the implicit operator C",
          "matrix-free GPM: Z = sum_j X_j R_j, B_i = X_i^T Z",
          "blockwise d x d polar projections",
          if (identical(spec$group, "O")) "dual certificate after a fixed point" else
            "no O(d) dual certificate for SO(d)"
        ),
        allowed_claim = "first_order_stationary; certified global only if dual test succeeds"
      ),
      class = "proc_solver_plan"
    ))
  }
  gower_ok <- x$n_views >= 2L &&
    identical(x$loss$family, "squared_l2") &&
    !isTRUE(x$cell_metric) &&
    spec$family %in% c("orthogonal", "similarity") &&
    length(unique(x$dimensions)) == 1L
  if (gower_ok) {
    return(structure(
      list(
        solver = "gower_bcd",
        reasons = c(
          sprintf("transformation is %s", spec$group),
          "loss is squared L2",
          sprintf("%d configurations; consensus-first Gower BCD", x$n_views),
          "no coordinate-specific precision",
          if (isTRUE(x$row_masked)) "row masks are handled by observed-set means" else "rows are complete"
        ),
        plan = c(
          "initialize by sequential, medoid, or spectral polar factors",
          "exact pairwise Procrustes of each view against the current consensus",
          "consensus is the weighted mean of aligned observed rows",
          "Gower / Ten Berge, preshape, or fixed-consensus scale convention",
          "safeguarded consensus extrapolation with rollback if F rises"
        ),
        allowed_claim = "blockwise_stationary"
      ),
      class = "proc_solver_plan"
    ))
  }
  structure(
    list(
      solver = "unimplemented",
      reasons = "no implemented kernel matches this compilation",
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
