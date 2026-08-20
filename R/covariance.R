#' Superimposition or model covariance.
#'
#' \(\Sigma_S\) is the metric used to fit. \(\Sigma_M\) is reserved for
#' inference and is not used by the solvers. A Kronecker structure
#' \(\Sigma_N\otimes\Sigma_d\) is the small set the package implements.
#' If \(\Sigma_d\propto I\), entity whitening or landmark weights recover
#' the exact SVD path. A genuine coordinate metric does not.
#'
#' @param entity Landmark / entity factor \(\Sigma_N\) or its precision.
#' @param coordinate Coordinate factor \(\Sigma_d\) or its precision.
#' @param kind `"superimposition"` (\(\Sigma_S\)) or `"model"` (\(\Sigma_M\)).
#' @param as_precision If `TRUE`, `entity` / `coordinate` are precisions.
#' @export
proc_covariance <- function(entity = NULL,
                            coordinate = NULL,
                            kind = c("superimposition", "model"),
                            as_precision = FALSE) {
  kind <- match.arg(kind)
  if (!is.null(entity)) entity <- as.matrix(entity)
  if (!is.null(coordinate)) coordinate <- as.matrix(coordinate)
  structure(
    list(
      entity = entity,
      coordinate = coordinate,
      kind = kind,
      as_precision = isTRUE(as_precision)
    ),
    class = "proc_covariance"
  )
}

#' @keywords internal
.gproc_as_covariance <- function(precision) {
  if (is.null(precision)) {
    return(NULL)
  }
  if (inherits(precision, "proc_covariance")) {
    return(precision)
  }
  if (is.matrix(precision) || is.numeric(precision)) {
    A <- as.matrix(precision)
    if (nrow(A) == ncol(A)) {
      return(proc_covariance(coordinate = A, as_precision = TRUE))
    }
  }
  .gproc_stop("invalid_problem", "precision must be a proc_covariance or a square matrix.")
}

#' @keywords internal
.gproc_coordinate_anisotropic <- function(covar) {
  if (is.null(covar) || is.null(covar$coordinate)) {
    return(FALSE)
  }
  A <- covar$coordinate
  off <- A - diag(diag(A), nrow(A))
  spread <- diff(range(diag(A)))
  any(abs(off) > 1e-12 * max(1, max(abs(A)))) || spread > 1e-12 * max(1, max(abs(diag(A))))
}

#' @noRd
.gproc_fitting_covariance <- function(metric) {
  if (is.null(metric) || is.null(metric$precision)) {
    return(NULL)
  }
  covar <- .gproc_as_covariance(metric$precision)
  if (identical(covar$kind, "model")) {
    return(NULL)
  }
  covar
}

#' @keywords internal
.gproc_entity_landmark_weights <- function(covar, n) {
  A <- covar$entity
  if (is.null(A)) {
    return(NULL)
  }
  A <- as.matrix(A)
  if (nrow(A) != n || ncol(A) != n) {
    .gproc_stop("invalid_problem", "Entity covariance must be n_entities x n_entities.")
  }
  A <- (A + t(A)) / 2
  off <- A - diag(diag(A), n)
  if (any(abs(off) > 1e-12 * max(1, max(abs(A))))) {
    .gproc_stop(
      "invalid_problem",
      "Off-diagonal entity covariance is not implemented; supply a diagonal entity factor."
    )
  }
  d <- diag(A)
  if (any(!is.finite(d) | d <= 0)) {
    .gproc_stop("invalid_problem", "Entity diagonal must be positive and finite.")
  }
  if (isTRUE(covar$as_precision)) d else 1 / d
}

#' @keywords internal
.gproc_combine_landmark_weights <- function(landmark, extra, nms) {
  if (is.null(landmark)) {
    return(extra)
  }
  if (is.list(landmark)) {
    out <- landmark
    for (nm in nms) {
      src <- landmark[[nm]]
      if (is.null(src)) {
        out[[nm]] <- extra
      } else {
        src <- as.numeric(src)
        if (length(src) == 1L) {
          out[[nm]] <- extra * src
        } else if (length(src) == length(extra)) {
          out[[nm]] <- extra * src
        } else {
          out[[nm]] <- src
        }
      }
    }
    return(out)
  }
  landmark <- as.numeric(landmark)
  if (length(landmark) == 1L) {
    return(extra * landmark)
  }
  if (length(landmark) == length(extra)) {
    return(extra * landmark)
  }
  landmark
}

#' @noRd
.gproc_apply_goodall_metric <- function(metric, data) {
  if (!inherits(metric, "proc_metric")) {
    return(metric)
  }
  covar <- .gproc_fitting_covariance(metric)
  if (is.null(covar)) {
    if (!is.null(metric$precision) && inherits(metric$precision, "proc_covariance") &&
        identical(metric$precision$kind, "model")) {
      metric$precision <- NULL
    }
    return(metric)
  }
  if (!is.null(covar$entity)) {
    w <- .gproc_entity_landmark_weights(covar, length(data$global_ids))
    metric$landmark <- .gproc_combine_landmark_weights(metric$landmark, w, names(data$views))
  }
  if (.gproc_coordinate_anisotropic(covar)) {
    return(metric)
  }
  metric$precision <- NULL
  metric
}

#' @keywords internal
.gproc_landmark_weighted <- function(metric, n) {
  if (!is.null(metric$landmark)) {
    src <- metric$landmark
    if (is.list(src)) {
      return(any(vapply(src, function(v) {
        !is.null(v) && any(abs(as.numeric(v) - 1) > 1e-12)
      }, logical(1))))
    }
    v <- as.numeric(src)
    if (!length(v)) {
      return(FALSE)
    }
    return(any(abs(v - 1) > 1e-12))
  }
  covar <- .gproc_fitting_covariance(metric)
  if (is.null(covar) || is.null(covar$entity) || .gproc_coordinate_anisotropic(covar)) {
    return(FALSE)
  }
  w <- .gproc_entity_landmark_weights(covar, n)
  any(abs(w - 1) > 1e-12)
}
