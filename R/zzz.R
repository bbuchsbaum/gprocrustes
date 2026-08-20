#' @keywords internal
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    registerS3method("autoplot", "gpa_fit", autoplot.gpa_fit, envir = asNamespace("ggplot2"))
  }
}
