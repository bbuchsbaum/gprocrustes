#' Transformation-level parameters of a GPA fit.
#'
#' One row per configuration. Columns are quantities already stored on the
#' fitted transforms: scale, translation, and `constraint_residual()`.
#'
#' @param x A `gpa_fit`.
#' @param ... Unused.
#' @return A data frame.
#' @export
#' @examples
#' views <- list(
#'   A = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
#'   B = matrix(c(0.1, 0, 1.1, 0.1, 1, 1.1, 0, 0.9), 4, 2, byrow = TRUE)
#' )
#' tidy.gpa_fit(gpa(views, transform = proc_orthogonal("O")))
tidy.gpa_fit <- function(x, ...) {
  trs <- transformations(x)
  nms <- names(trs)
  if (is.null(nms)) {
    nms <- as.character(seq_along(trs))
  }
  rows <- lapply(seq_along(trs), function(i) {
    tr <- trs[[i]]
    cons <- constraint_residual(tr)
    tvec <- as.numeric(tr$t)
    out <- data.frame(
      configuration = nms[[i]],
      family = tr$spec$family,
      group = tr$spec$group %||% NA_character_,
      s = as.numeric(tr$s)[1L],
      orthogonality = cons$orthogonality,
      determinant = cons$determinant,
      stringsAsFactors = FALSE
    )
    if (length(tvec)) {
      tdf <- as.data.frame(as.list(tvec), optional = TRUE)
      names(tdf) <- paste0("t", seq_along(tvec))
      out <- cbind(out, tdf)
    }
    out
  })
  do.call(rbind, rows)
}

#' Fit-level diagnostics of a GPA fit.
#'
#' @param x A `gpa_fit`.
#' @param ... Unused.
#' @return A one-row data frame.
#' @export
#' @examples
#' views <- list(
#'   A = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
#'   B = matrix(c(0.1, 0, 1.1, 0.1, 1, 1.1, 0, 0.9), 4, 2, byrow = TRUE)
#' )
#' glance.gpa_fit(gpa(views, transform = proc_orthogonal("O")))
glance.gpa_fit <- function(x, ...) {
  ds <- x$data_summary
  niter <- if (is.data.frame(x$history) && nrow(x$history)) {
    max(x$history$iteration)
  } else {
    0L
  }
  dec <- x$decomposition
  data.frame(
    n_configurations = ds$n_views,
    n_entities = ds$n_entities,
    dimension = ds$target_dimension,
    objective = unname(x$objective),
    residual_energy = unname(dec$residual_energy %||% NA_real_),
    consensus_energy = unname(dec$consensus_energy %||% NA_real_),
    solver = x$solver,
    numerical_status = x$numerical_status,
    optimality_status = x$optimality_status,
    certificate = x$certificate$status %||% "unavailable",
    effective_rank = x$rank$effective_rank %||% NA_real_,
    niter = as.integer(niter),
    stringsAsFactors = FALSE
  )
}

#' Entity- or configuration-level residuals and weights.
#'
#' @param x A `gpa_fit`.
#' @param level `"landmark"` for one row per entity and configuration, or
#'   `"configuration"` for one row per view.
#' @param ... Unused.
#' @return A data frame.
#' @export
#' @examples
#' views <- list(
#'   A = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
#'   B = matrix(c(0.1, 0, 1.1, 0.1, 1, 1.1, 0, 0.9), 4, 2, byrow = TRUE)
#' )
#' augment.gpa_fit(gpa(views, transform = proc_orthogonal("O")))
augment.gpa_fit <- function(x, level = c("landmark", "configuration"), ...) {
  level <- match.arg(level)
  if (identical(level, "configuration")) {
    return(.gproc_augment_configuration(x))
  }
  .gproc_augment_landmark(x)
}

#' @noRd
.gproc_augment_landmark <- function(x) {
  res <- residuals(x, level = "landmark")
  ww <- .gproc_view_weights(x)
  rows <- lapply(names(res), function(nm) {
    E <- as.matrix(res[[nm]])
    w <- ww[[nm]]
    if (is.null(w)) {
      w <- rep(1, nrow(E))
    }
    out <- data.frame(
      entity = seq_len(nrow(E)),
      configuration = nm,
      residual = sqrt(rowSums(E^2, na.rm = TRUE)),
      weight = as.numeric(w)[seq_len(nrow(E))],
      stringsAsFactors = FALSE
    )
    if (ncol(E)) {
      rdf <- as.data.frame(E, optional = TRUE)
      names(rdf) <- paste0("r", seq_len(ncol(E)))
      out <- cbind(out, rdf)
    }
    out
  })
  do.call(rbind, rows)
}

#' @noRd
.gproc_augment_configuration <- function(x) {
  r <- residuals(x, level = "configuration")
  nms <- names(r)
  if (is.null(nms)) {
    nms <- as.character(seq_along(r))
  }
  alpha <- x$weights$configuration
  if (is.null(alpha)) {
    w <- rep(1, length(r))
  } else {
    w <- as.numeric(alpha)
    if (!is.null(names(alpha))) {
      w <- as.numeric(alpha[nms])
    }
    if (length(w) == 1L) {
      w <- rep(w, length(r))
    }
  }
  data.frame(
    configuration = nms,
    residual = as.numeric(r),
    weight = w[seq_along(r)],
    stringsAsFactors = FALSE
  )
}
