#' Gauge specification.
#'
#' @param scale Scale convention: `"gower"`, `"preshape"`, `"fixed_consensus"`, or `"none"`.
#' @param orientation `"free"` during estimation or `"principal"` for a presentation.
#' @export
proc_gauge <- function(scale = c("gower", "preshape", "fixed_consensus", "none"),
                       orientation = c("free", "principal")) {
  structure(
    list(scale = match.arg(scale), orientation = match.arg(orientation)),
    class = "proc_gauge"
  )
}

#' Superimposition metric / weight channels.
#'
#' @param configuration Configuration weights \eqn{\alpha_i}.
#' @param landmark Landmark / entity weights.
#' @param cell Coordinate-specific weights.
#' @param precision Precision operator \eqn{Q}.
#' @export
proc_metric <- function(configuration = NULL,
                        landmark = NULL,
                        cell = NULL,
                        precision = NULL) {
  structure(
    list(
      configuration = configuration,
      landmark = landmark,
      cell = cell,
      precision = precision
    ),
    class = "proc_metric"
  )
}

#' @rdname proc_metric
#' @export
proc_weights <- function(configuration = NULL,
                         landmark = NULL,
                         cell = NULL,
                         precision = NULL) {
  proc_metric(configuration, landmark, cell, precision)
}

#' Squared Frobenius loss.
#' @export
proc_squared_l2 <- function() {
  structure(list(family = "squared_l2"), class = c("proc_squared_l2", "proc_loss"))
}

#' Landmark-vector Huber loss.
#'
#' @param level Residual level; only `"landmark"` is supported.
#' @param k Huber threshold.
#' @export
proc_huber <- function(level = "landmark", k = 1.345) {
  if (!identical(level, "landmark")) {
    .gproc_stop("invalid_problem", "Huber loss must act on landmark residual vectors.")
  }
  structure(
    list(family = "huber", level = level, k = k, convex = TRUE),
    class = c("proc_huber", "proc_loss")
  )
}

#' Landmark-vector Tukey bisquare loss (nonconvex).
#'
#' @param level Residual level; only `"landmark"` is supported.
#' @param c Tuning constant.
#' @export
proc_tukey <- function(level = "landmark", c = 4.685) {
  if (!identical(level, "landmark")) {
    .gproc_stop("invalid_problem", "Tukey loss must act on landmark residual vectors.")
  }
  structure(
    list(family = "tukey", level = level, c = c, convex = FALSE),
    class = c("proc_tukey", "proc_loss")
  )
}

#' Solver control.
#'
#' @param tolerance Convergence tolerance.
#' @param max_iterations Maximum iterations for iterative solvers.
#' @param keep_aligned `"lazy"` or `"materialize"`.
#' @param certify `"auto"`, `"never"`, or `"always"`.
#' @param init Initialization for iterative GPA: `"sequential"`, `"medoid"`,
#'   or `"spectral"`. `NULL` selects the engine default: spectral for GPM and
#'   sequential for Gower BCD.
#' @param accelerate If `TRUE`, try a safeguarded extrapolated consensus step
#'   and roll it back if the true objective rises.
#' @param nstart Number of distinct deterministic starts, up to three. The
#'   requested or engine-default initialization is tried first, followed by
#'   the remaining medoid, spectral, and sequential alternatives. The lowest
#'   exact objective is kept.
#' @param backend `"auto"` forms the small block Gram when \eqn{Kd} is below
#'   `dense_block_threshold`; `"matrix_free"` never does; `"dense"` always does.
#' @param dense_block_threshold Maximum \eqn{Kd} for a dense block-Gram eigenstep.
#' @export
gpa_control <- function(tolerance = 1e-8,
                        max_iterations = 500L,
                        keep_aligned = c("lazy", "materialize"),
                        certify = c("auto", "never", "always"),
                        init = NULL,
                        accelerate = TRUE,
                        nstart = 1L,
                        backend = c("auto", "matrix_free", "dense"),
                        dense_block_threshold = 96L) {
  structure(
    list(
      tolerance = tolerance,
      max_iterations = as.integer(max_iterations),
      keep_aligned = match.arg(keep_aligned),
      certify = match.arg(certify),
      init = if (is.null(init)) {
        NULL
      } else {
        match.arg(init, c("sequential", "medoid", "spectral"))
      },
      accelerate = isTRUE(accelerate),
      nstart = as.integer(nstart)[1L],
      backend = match.arg(backend),
      dense_block_threshold = as.integer(dense_block_threshold)[1L]
    ),
    class = "gpa_control"
  )
}
