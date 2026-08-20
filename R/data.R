#' Construct a multi-view Procrustes data object.
#'
#' @param views Named list of entity-by-dimension matrix-like objects.
#' @param ids Optional named list of entity identifiers, one vector per view.
#' @param observed Optional named list of logical row masks (`TRUE` = observed).
#' @param metadata Optional named list of per-view metadata.
#' @return A `proc_data` object.
#' @export
proc_data <- function(views, ids = NULL, observed = NULL, metadata = NULL) {
  if (is.matrix(views) || inherits(views, "Matrix")) {
    views <- list(view1 = views)
  }
  if (!is.list(views) || length(views) < 1L) {
    .gproc_stop("invalid_problem", "`views` must be a nonempty list of matrices.")
  }
  if (is.null(names(views)) || any(!nzchar(names(views)))) {
    names(views) <- paste0("view", seq_along(views))
  }
  if (anyDuplicated(names(views))) {
    .gproc_stop("invalid_problem", "View names must be unique.")
  }
  views <- lapply(views, .gproc_as_numeric_matrix)
  for (nm in names(views)) {
    .gproc_check_finite(views[[nm]], nm)
  }
  ids <- .gproc_align_named(ids, names(views), "ids")
  observed <- .gproc_align_named(observed, names(views), "observed")
  maps <- .gproc_correspondence(views, ids)
  masks <- .gproc_row_masks(views, observed)
  structure(
    list(
      views = views,
      ids = maps$ids,
      row_map = maps$row_map,
      global_ids = maps$global_ids,
      observed = masks,
      correspondence = maps$assumption,
      metadata = metadata
    ),
    class = "proc_data"
  )
}

#' @export
print.proc_data <- function(x, ...) {
  dims <- vapply(x$views, function(v) paste(dim(v), collapse = "x"), character(1))
  cat(sprintf(
    "Procrustes data: %d views, %d consensus entities, correspondence=%s\n",
    length(x$views), length(x$global_ids), x$correspondence
  ))
  for (nm in names(x$views)) {
    cat(sprintf("  %s: %s\n", nm, dims[[nm]]))
  }
  invisible(x)
}

#' @keywords internal
.gproc_align_named <- function(x, nms, label) {
  if (is.null(x)) {
    return(stats::setNames(vector("list", length(nms)), nms))
  }
  if (!is.list(x)) {
    .gproc_stop("invalid_problem", sprintf("`%s` must be a named list.", label))
  }
  if (is.null(names(x))) {
    if (length(x) != length(nms)) {
      .gproc_stop("invalid_problem", sprintf("`%s` must be named or have one entry per view.", label))
    }
    names(x) <- nms
  }
  out <- stats::setNames(vector("list", length(nms)), nms)
  extra <- setdiff(names(x), nms)
  if (length(extra)) {
    .gproc_stop("invalid_problem", sprintf("`%s` has unknown views: %s.", label, paste(extra, collapse = ", ")))
  }
  for (nm in intersect(names(x), nms)) {
    out[[nm]] <- x[[nm]]
  }
  out
}

#' @keywords internal
.gproc_correspondence <- function(views, ids) {
  nms <- names(views)
  resolved <- vector("list", length(nms))
  names(resolved) <- nms
  have <- stats::setNames(rep(FALSE, length(nms)), nms)
  for (nm in nms) {
    id <- ids[[nm]]
    if (is.null(id)) {
      rn <- rownames(views[[nm]])
      if (!is.null(rn)) id <- rn
    }
    if (!is.null(id)) {
      id <- as.character(id)
      if (length(id) != .gproc_nrow(views[[nm]])) {
        .gproc_stop("invalid_problem", sprintf("ids for '%s' must have length nrow(view).", nm))
      }
      if (anyDuplicated(id)) {
        .gproc_stop("duplicate_entity_ids", sprintf("Duplicate entity identifiers in '%s'.", nm))
      }
      resolved[[nm]] <- id
      have[[nm]] <- TRUE
    }
  }
  if (length(have) && all(have)) {
    global <- unique(unlist(resolved, use.names = FALSE))
    row_map <- lapply(resolved, function(id) match(id, global))
    return(list(ids = resolved, row_map = row_map, global_ids = global, assumption = "explicit_ids"))
  }
  if (any(have)) {
    .gproc_stop(
      "invalid_problem",
      "Partial row names or ids require an explicit map for every view."
    )
  }
  nr <- vapply(views, .gproc_nrow, integer(1))
  if (length(unique(nr)) != 1L) {
    .gproc_stop(
      "invalid_problem",
      "Unnamed views of unequal size cannot use positional correspondence."
    )
  }
  global <- as.character(seq_len(nr[[1L]]))
  resolved <- lapply(views, function(v) global)
  row_map <- lapply(resolved, function(id) match(id, global))
  list(ids = resolved, row_map = row_map, global_ids = global, assumption = "positional")
}

#' @keywords internal
.gproc_row_masks <- function(views, observed) {
  out <- vector("list", length(views))
  names(out) <- names(views)
  for (nm in names(views)) {
    m <- observed[[nm]]
    n <- .gproc_nrow(views[[nm]])
    if (is.null(m)) {
      out[[nm]] <- rep(TRUE, n)
    } else {
      m <- as.logical(m)
      if (length(m) != n) {
        .gproc_stop("invalid_problem", sprintf("Row mask for '%s' has the wrong length.", nm))
      }
      if (anyNA(m)) {
        .gproc_stop("invalid_problem", "Row masks cannot contain NA.")
      }
      out[[nm]] <- m
    }
  }
  out
}

#' Overlap graph of configurations.
#'
#' @param data A `proc_data` object.
#' @param min_overlap Minimum shared observed entities for an edge.
#' @return A list with adjacency metadata and connected components.
#' @export
overlap_graph <- function(data, min_overlap = 1L) {
  nms <- names(data$views)
  k <- length(nms)
  shared <- matrix(0L, k, k, dimnames = list(nms, nms))
  for (nm_i in nms) {
    for (nm_j in nms) {
      if (identical(nm_i, nm_j)) next
      ii <- data$ids[[nm_i]][data$observed[[nm_i]]]
      jj <- data$ids[[nm_j]][data$observed[[nm_j]]]
      shared[nm_i, nm_j] <- length(intersect(ii, jj))
    }
  }
  adj <- shared >= min_overlap
  diag(adj) <- TRUE
  comp <- .gproc_components(adj)
  list(
    names = nms,
    overlap_count = shared,
    adjacency = adj,
    components = comp$membership,
    n_components = comp$n
  )
}

#' @keywords internal
.gproc_components <- function(adj) {
  k <- nrow(adj)
  seen <- rep(FALSE, k)
  memb <- integer(k)
  n <- 0L
  for (i in seq_len(k)) {
    if (seen[i]) next
    n <- n + 1L
    stack <- i
    while (length(stack)) {
      v <- stack[[1L]]
      stack <- stack[-1L]
      if (seen[v]) next
      seen[v] <- TRUE
      memb[v] <- n
      nbr <- which(adj[v, ])
      stack <- c(stack, nbr[!seen[nbr]])
    }
  }
  list(membership = memb, n = n)
}
