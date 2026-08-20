#' Tidy data for a GPA diagnostic plot.
#'
#' @param object A `gpa_fit`.
#' @param type One of `"overlay"`, `"residuals"`, `"convergence"`,
#'   `"decomposition"`, `"influence"`, `"support"`, `"weights"`,
#'   `"uncertainty"`.
#' @param ... Unused.
#' @return A data frame. For large overlays (`n > 5000`), a random sample
#'   of entities is returned and `sampled` / `n_total` attributes are set.
#' @export
plot_data <- function(object, type = "overlay", ...) {
  UseMethod("plot_data")
}

#' @export
plot_data.gpa_fit <- function(object,
                              type = c("overlay", "residuals", "convergence", "decomposition",
                                       "influence", "support", "weights", "uncertainty",
                                       "deformation", "cv"),
                              ...) {
  type <- match.arg(type)
  switch(
    type,
    overlay = .gproc_plot_data_overlay(object),
    residuals = .gproc_plot_data_residuals(object),
    convergence = .gproc_plot_data_convergence(object),
    decomposition = .gproc_plot_data_decomposition(object),
    influence = .gproc_plot_data_influence(object),
    support = .gproc_plot_data_support(object),
    weights = .gproc_plot_data_weights(object),
    uncertainty = .gproc_plot_data_uncertainty(object, ...),
    deformation = .gproc_plot_data_deformation(object),
    cv = .gproc_plot_data_cv(object)
  )
}

#' ggplot2 diagnostics for a `gpa_fit` or `proc_inference`.
#'
#' @param object A `gpa_fit` or `proc_inference`.
#' @param type Plot type; see `plot_data()`.
#' @param ... Unused.
#' @return A ggplot object.
#' @export
autoplot.gpa_fit <- function(object,
                             type = c("overlay", "residuals", "convergence", "decomposition",
                                      "influence", "support", "weights", "uncertainty",
                                      "deformation", "cv"),
                             ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .gproc_stop(
      "package_missing",
      "ggplot2 is required for autoplot(); install it or use plot_data()."
    )
  }
  type <- match.arg(type)
  dat <- plot_data(object, type = type)
  switch(
    type,
    overlay = .gproc_autoplot_overlay(dat),
    residuals = .gproc_autoplot_residuals(dat),
    convergence = .gproc_autoplot_convergence(dat),
    decomposition = .gproc_autoplot_decomposition(dat),
    influence = .gproc_autoplot_influence(dat),
    support = .gproc_autoplot_support(dat),
    weights = .gproc_autoplot_weights(dat),
    uncertainty = .gproc_autoplot_uncertainty(dat),
    deformation = .gproc_autoplot_deformation(dat),
    cv = .gproc_autoplot_cv(dat)
  )
}

#' @export
plot.gpa_fit <- function(x,
                         type = c("overlay", "residuals", "convergence", "decomposition",
                                  "influence", "support", "weights", "uncertainty"),
                         ...) {
  print(autoplot.gpa_fit(x, type = type, ...))
  invisible(x)
}

#' @keywords internal
.gproc_plot_data_overlay <- function(object) {
  Y <- .gproc_fitted_global(object)
  M <- consensus(object)
  d <- ncol(M)
  rot <- diag(d)
  note <- NULL
  coords <- M
  if (d > 2L) {
    sv <- .gproc_svd(M, nu = min(2L, nrow(M)), nv = min(2L, d))
    rot <- sv$v[, 1:2, drop = FALSE]
    coords <- M %*% rot
    note <- "First two principal axes of the consensus (presentation only)."
  } else if (d == 1L) {
    coords <- cbind(M, 0)
  }
  sample_n <- 5000L
  idx <- seq_len(nrow(M))
  sampled <- FALSE
  if (length(idx) > sample_n) {
    idx <- sort(sample.int(nrow(M), sample_n))
    sampled <- TRUE
  }
  rows <- lapply(names(Y), function(nm) {
    Yi <- Y[[nm]]
    if (d > 2L) {
      Yi <- Yi %*% rot
    } else if (d == 1L) {
      Yi <- cbind(Yi, 0)
    }
    data.frame(
      entity = idx,
      view = nm,
      x = Yi[idx, 1L],
      y = Yi[idx, 2L],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  cons <- data.frame(
    entity = idx,
    view = "consensus",
    x = coords[idx, 1L],
    y = coords[idx, 2L],
    stringsAsFactors = FALSE
  )
  out <- rbind(out, cons)
  attr(out, "sampled") <- sampled
  attr(out, "n_total") <- nrow(M)
  attr(out, "note") <- note
  out
}

#' @keywords internal
.gproc_plot_data_residuals <- function(object) {
  res <- residuals(object, level = "landmark")
  rows <- lapply(names(res), function(nm) {
    E <- res[[nm]]
    data.frame(
      entity = seq_len(nrow(E)),
      view = nm,
      residual = sqrt(rowSums(E^2)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @keywords internal
.gproc_plot_data_convergence <- function(object) {
  h <- object$history
  if (!is.data.frame(h) || !nrow(h)) {
    return(data.frame(
      iteration = integer(),
      objective = numeric(),
      relative_change = numeric(),
      stationarity = numeric(),
      accepted_acceleration = logical()
    ))
  }
  h
}

#' @keywords internal
.gproc_plot_data_decomposition <- function(object) {
  dec <- decompose(object)
  data.frame(
    configuration = names(dec$by_configuration),
    residual_energy = as.numeric(dec$by_configuration),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.gproc_autoplot_overlay <- function(dat) {
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], colour = .data[["view"]])) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(
      title = "Aligned configurations and consensus",
      x = "Dimension 1",
      y = "Dimension 2",
      colour = "View"
    )
  note <- attr(dat, "note")
  if (!is.null(note)) {
    p <- p + ggplot2::labs(subtitle = note)
  }
  if (isTRUE(attr(dat, "sampled"))) {
    p <- p + ggplot2::labs(
      caption = sprintf("Sampled %d of %d entities.", nrow(dat), attr(dat, "n_total"))
    )
  }
  p
}

#' @keywords internal
.gproc_autoplot_residuals <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["entity"]], y = .data[["residual"]], colour = .data[["view"]])) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(
      title = "Landmark residual norms",
      x = "Entity",
      y = "||Y_i - M||",
      colour = "View"
    )
}

#' @keywords internal
.gproc_autoplot_decomposition <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["configuration"]], y = .data[["residual_energy"]])) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "Residual energy by configuration",
      x = "Configuration",
      y = "Residual energy"
    )
}

#' @keywords internal
.gproc_view_weights <- function(object) {
  w <- object$weights$robust
  if (is.null(w)) {
    w <- object$weights$landmark
  }
  nms <- object$data_summary$view_names
  n <- object$data_summary$n_entities
  out <- vector("list", length(nms))
  names(out) <- nms
  if (is.null(w)) {
    for (nm in nms) out[[nm]] <- rep(1, n)
    return(out)
  }
  if (!is.list(w)) {
    w <- as.numeric(w)
    for (nm in nms) {
      out[[nm]] <- if (length(w) == 1L) rep(w, n) else if (length(w) == n) w else rep(1, n)
    }
    return(out)
  }
  for (nm in nms) {
    src <- w[[nm]]
    if (is.null(src)) {
      out[[nm]] <- rep(1, n)
    } else {
      src <- as.numeric(src)
      out[[nm]] <- if (length(src) == 1L) {
        rep(src, n)
      } else if (length(src) == n) {
        src
      } else {
        c(src, rep(NA_real_, max(0L, n - length(src))))[seq_len(n)]
      }
    }
  }
  out
}

#' @keywords internal
.gproc_plot_data_influence <- function(object) {
  res <- .gproc_plot_data_residuals(object)
  ww <- .gproc_view_weights(object)
  res$weight <- mapply(function(view, entity) {
    w <- ww[[as.character(view)]]
    if (is.null(w) || entity > length(w)) NA_real_ else w[[entity]]
  }, res$view, res$entity)
  res$influence <- res$residual * res$weight
  res
}

#' @keywords internal
.gproc_plot_data_support <- function(object) {
  Y <- .gproc_fitted_global(object)
  n <- nrow(object$consensus)
  d <- ncol(object$consensus)
  n_obs <- integer(n)
  for (Yi in Y) {
    n_obs <- n_obs + as.integer(is.finite(rowSums(Yi)))
  }
  n_cells <- if (!is.null(object$cell_masks)) {
    tot <- integer(n)
    for (C in object$cell_masks) {
      if (!is.null(C)) tot <- tot + as.integer(rowSums(C))
    }
    tot
  } else {
    n_obs * d
  }
  data.frame(
    entity = seq_len(n),
    n_observed_views = n_obs,
    n_observed_cells = n_cells,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.gproc_plot_data_weights <- function(object) {
  ww <- .gproc_view_weights(object)
  n <- object$data_summary$n_entities
  rows <- lapply(names(ww), function(nm) {
    data.frame(
      entity = seq_len(n),
      view = nm,
      weight = ww[[nm]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @keywords internal
.gproc_autoplot_convergence <- function(dat) {
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[["iteration"]], y = .data[["objective"]])) +
    ggplot2::geom_line()
  if ("accepted_acceleration" %in% names(dat)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data[["accepted_acceleration"]]))
  } else {
    p <- p + ggplot2::geom_point()
  }
  p + ggplot2::labs(
    title = "GPA objective",
    x = "Iteration",
    y = "F",
    shape = "Accepted acceleration"
  )
}

#' @keywords internal
.gproc_autoplot_influence <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["entity"]], y = .data[["influence"]], colour = .data[["view"]])) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(
      title = "Landmark influence (residual x weight)",
      x = "Entity",
      y = "r_ij w_ij",
      colour = "View"
    )
}

#' @keywords internal
.gproc_autoplot_support <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["entity"]], y = .data[["n_observed_views"]])) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "Observation support",
      x = "Entity",
      y = "Observed views"
    )
}

#' @keywords internal
.gproc_autoplot_weights <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["entity"]], y = .data[["weight"]], colour = .data[["view"]])) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(
      title = "Landmark weights",
      x = "Entity",
      y = "w",
      colour = "View"
    )
}

#' @keywords internal
.gproc_plot_data_uncertainty <- function(object, inference = NULL, ...) {
  if (is.null(inference)) {
    inference <- object$inference
  }
  if (inherits(object, "proc_inference")) {
    inference <- object
  }
  if (is.null(inference) || !inherits(inference, "proc_inference")) {
    .gproc_stop(
      "invalid_problem",
      "Uncertainty plots need infer() output via inference= or object$inference."
    )
  }
  if (is.null(inference$se)) {
    .gproc_stop("invalid_problem", "This inference object has no coordinatewise se.")
  }
  M <- inference$estimate
  se <- inference$se
  d <- ncol(M)
  data.frame(
    entity = seq_len(nrow(M)),
    x = M[, 1L],
    y = if (d >= 2L) M[, 2L] else 0,
    se_x = se[, 1L],
    se_y = if (d >= 2L) se[, 2L] else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' @export
plot_data.proc_inference <- function(object, type = "uncertainty", ...) {
  .gproc_plot_data_uncertainty(object, inference = object)
}

#' @rdname autoplot.gpa_fit
#' @export
autoplot.proc_inference <- function(object, type = "uncertainty", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .gproc_stop(
      "package_missing",
      "ggplot2 is required for autoplot(); install it or use plot_data()."
    )
  }
  .gproc_autoplot_uncertainty(plot_data(object, type = type))
}

#' @keywords internal
.gproc_autoplot_uncertainty <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data[["y"]] - .data[["se_y"]], ymax = .data[["y"]] + .data[["se_y"]]),
      width = 0
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data[["x"]] - .data[["se_x"]], xmax = .data[["x"]] + .data[["se_x"]]),
      height = 0
    ) +
    ggplot2::labs(
      title = "Gauge-aligned consensus uncertainty",
      x = "Dimension 1",
      y = "Dimension 2"
    )
}

#' @keywords internal
.gproc_plot_data_deformation <- function(object) {
  fold <- object$folding
  if (is.null(fold)) {
    return(data.frame(view = character(), jacobian_det = numeric(), folded = logical()))
  }
  rows <- lapply(names(fold$jacobian_det), function(nm) {
    z <- fold$jacobian_det[[nm]]
    data.frame(
      view = nm,
      entity = seq_along(z),
      jacobian_det = as.numeric(z),
      folded = as.numeric(z) < 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' @keywords internal
.gproc_plot_data_cv <- function(object) {
  cv <- object$cv
  if (is.null(cv) && !is.null(object$smoothness_cv)) {
    cv <- object$smoothness_cv
  }
  if (is.null(cv)) {
    .gproc_stop("invalid_problem", "No cross-validation result is stored on this fit.")
  }
  data.frame(
    smoothness = cv$smoothness %||% seq_along(cv$scores),
    score = as.numeric(cv$scores),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.gproc_autoplot_deformation <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["entity"]], y = .data[["jacobian_det"]], colour = .data[["view"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(
      title = "Warp Jacobian determinant",
      x = "Entity",
      y = "det J",
      colour = "View"
    )
}

#' @keywords internal
.gproc_autoplot_cv <- function(dat) {
  ggplot2::ggplot(dat, ggplot2::aes(x = .data[["smoothness"]], y = .data[["score"]])) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "Held-out landmark CV",
      x = "smoothness",
      y = "CV score"
    )
}
