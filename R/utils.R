#' @keywords internal
.gproc_stop <- function(code, message, ...) {
  cond <- structure(
    list(message = message, call = NULL, error_code = code, extra = list(...)),
    class = c(code, "gprocrustes_error", "error", "condition")
  )
  stop(cond)
}

#' @keywords internal
.gproc_is_matrixlike <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

#' @keywords internal
.gproc_as_numeric_matrix <- function(x, name = "matrix") {
  if (inherits(x, "Matrix")) {
    return(x)
  }
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  x
}

#' @keywords internal
.gproc_check_finite <- function(x, name = "matrix") {
  vals <- if (inherits(x, "Matrix")) x@x else as.vector(x)
  if (any(!is.finite(vals))) {
    .gproc_stop("nonfinite_values", sprintf("Nonfinite values in %s.", name))
  }
  invisible(TRUE)
}

#' @keywords internal
.gproc_nrow <- function(x) nrow(x)

#' @keywords internal
.gproc_ncol <- function(x) ncol(x)

#' @keywords internal
.gproc_ones <- function(n) {
  matrix(1, n, 1L)
}

#' @keywords internal
.gproc_match_arg <- function(arg, choices, name) {
  if (length(arg) != 1L || !is.character(arg)) {
    .gproc_stop("invalid_problem", sprintf("`%s` must be one of: %s.", name, paste(choices, collapse = ", ")))
  }
  if (!arg %in% choices) {
    .gproc_stop("invalid_problem", sprintf("`%s` must be one of: %s.", name, paste(choices, collapse = ", ")))
  }
  arg
}

#' @keywords internal
.gproc_pkg_version <- function() {
  as.character(utils::packageVersion("gprocrustes"))
}
