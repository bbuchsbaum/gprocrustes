#' Lazy aligned view: original data plus a fitted transform.
#'
#' Materializes only when subset or converted to a matrix. A general rotation
#' makes sparse columns dense; that is reported, not hidden.
#'
#' @param x Original matrix-like object.
#' @param transform Fitted transform.
#' @param row_map Optional local-to-global row map.
#' @param mask Optional observed-row mask.
#' @param densifies_sparse Logical; `TRUE` if applying `R` fills in zeros.
#' @export
proc_aligned_view <- function(x,
                              transform,
                              row_map = NULL,
                              mask = NULL,
                              densifies_sparse = NULL) {
  if (is.null(densifies_sparse)) {
    densifies_sparse <- inherits(x, "Matrix") && !isTRUE(transform$s == 1 &&
      isTRUE(all.equal(transform$R, diag(nrow(transform$R)), tolerance = 1e-14)) &&
      all(transform$t == 0))
  }
  structure(
    list(
      x = x,
      transform = transform,
      row_map = row_map,
      mask = mask,
      densifies_sparse = isTRUE(densifies_sparse)
    ),
    class = "proc_aligned_view"
  )
}

#' @export
dim.proc_aligned_view <- function(x) dim(x$x)

#' @export
`[.proc_aligned_view` <- function(x, i, j, drop = TRUE) {
  if (missing(i)) {
    i <- seq_len(nrow(x$x))
  }
  raw <- x$x[i, , drop = FALSE]
  out <- apply_proc_transform(x$transform, raw)
  if (!missing(j)) {
    out <- out[, j, drop = drop]
  } else if (isTRUE(drop) && nrow(out) == 1L) {
    out <- drop(out)
  }
  out
}

#' @export
as.matrix.proc_aligned_view <- function(x, ...) {
  apply_proc_transform(x$transform, x$x)
}

#' Apply a function to a materialized aligned view.
#'
#' @param x A `proc_aligned_view`.
#' @param FUN A function of a matrix.
#' @param ... Passed to `FUN`.
#' @export
block_apply <- function(x, FUN, ...) {
  FUN(as.matrix(x), ...)
}

#' @export
print.proc_aligned_view <- function(x, ...) {
  cat(sprintf(
    "Lazy aligned view %s, densifies_sparse=%s\n",
    paste(dim(x), collapse = "x"),
    x$densifies_sparse
  ))
  invisible(x)
}
