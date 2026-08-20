#' @keywords internal
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    registerS3method("autoplot", "gpa_fit", autoplot.gpa_fit, envir = asNamespace("ggplot2"))
    registerS3method("autoplot", "proc_inference", autoplot.proc_inference,
                     envir = asNamespace("ggplot2"))
  }
  if (requireNamespace("broom", quietly = TRUE)) {
    registerS3method("tidy", "gpa_fit", tidy.gpa_fit, envir = asNamespace("broom"))
    registerS3method("glance", "gpa_fit", glance.gpa_fit, envir = asNamespace("broom"))
    registerS3method("augment", "gpa_fit", augment.gpa_fit, envir = asNamespace("broom"))
  }
}
