#' Construct a multi-view Procrustes data object.
#'
#' @param views Named list of entity-by-dimension matrix-like objects.
#' @param ids Optional named list of entity identifiers, one vector per view.
#' @param observed Optional named list of logical row masks (`TRUE` = observed).
#' @param cells Optional named list of logical entity-by-dimension masks.
#'   A missing cell is not a missing row.
#' @param metadata Optional named list of per-view metadata.
#' @return A `proc_data` object.
#' @export
proc_data <- function(views, ids = NULL, observed = NULL, cells = NULL, metadata = NULL) {
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
  cells <- .gproc_align_named(cells, names(views), "cells")
  maps <- .gproc_correspondence(views, ids)
  masks <- .gproc_row_masks(views, observed)
  cell_masks <- .gproc_cell_local(views, cells)
  for (nm in names(views)) {
    if (!is.null(cell_masks[[nm]])) {
      masks[[nm]] <- masks[[nm]] & apply(cell_masks[[nm]], 1L, any)
    }
  }
  structure(
    list(
      views = views,
      ids = maps$ids,
      row_map = maps$row_map,
      global_ids = maps$global_ids,
      observed = masks,
      cells = cell_masks,
      correspondence = maps$assumption,
      metadata = metadata
    ),
    class = "proc_data"
  )
}

#' Domain adapter: a 3-way entity-by-dimension-by-view array.
#'
#' @param x An array with dimensions \eqn{(n, d, K)} or a list of matrices.
#' @param view_names Optional view names.
#' @export
proc_from_array <- function(x, view_names = NULL) {
  if (is.list(x) && !is.array(x)) {
    if (!is.null(view_names)) names(x) <- view_names
    return(proc_data(x))
  }
  if (length(dim(x)) != 3L) {
    .gproc_stop("invalid_problem", "proc_from_array() expects an n x d x K array.")
  }
  k <- dim(x)[[3L]]
  if (is.null(view_names)) {
    view_names <- paste0("view", seq_len(k))
  }
  views <- lapply(seq_len(k), function(i) {
    m <- x[, , i, drop = FALSE]
    dim(m) <- dim(x)[1:2]
    m
  })
  names(views) <- view_names
  proc_data(views)
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

#' @keywords internal
.gproc_cell_local <- function(views, cells) {
  out <- vector("list", length(views))
  names(out) <- names(views)
  for (nm in names(views)) {
    C <- cells[[nm]]
    if (is.null(C)) {
      next
    }
    C <- as.matrix(C)
    storage.mode(C) <- "logical"
    if (nrow(C) != .gproc_nrow(views[[nm]]) || ncol(C) != .gproc_ncol(views[[nm]])) {
      .gproc_stop("invalid_problem", sprintf("Cell mask for '%s' has the wrong dimension.", nm))
    }
    out[[nm]] <- C
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
  ranks <- matrix(NA_integer_, k, k, dimnames = list(nms, nms))
  for (nm_i in nms) {
    for (nm_j in nms) {
      ii <- data$ids[[nm_i]][data$observed[[nm_i]]]
      jj <- data$ids[[nm_j]][data$observed[[nm_j]]]
      common <- intersect(ii, jj)
      shared[nm_i, nm_j] <- length(common)
      if (length(common) && !is.null(data$views[[nm_i]])) {
        ia <- match(common, data$ids[[nm_i]])
        ja <- match(common, data$ids[[nm_j]])
        C <- centered_crossprod(
          data$views[[nm_i]][ia, , drop = FALSE],
          data$views[[nm_j]][ja, , drop = FALSE]
        )
        sig <- .gproc_svd(C, nu = 0L, nv = 0L)$d
        s1 <- if (length(sig)) max(sig[[1L]], 0) else 0
        ranks[nm_i, nm_j] <- if (s1 <= 0) 0L else as.integer(sum(sig > 1e-10 * s1))
      }
    }
  }
  adj <- shared >= min_overlap
  diag(adj) <- TRUE
  comp <- .gproc_components(adj)
  support <- integer(length(data$global_ids))
  names(support) <- data$global_ids
  for (nm in nms) {
    ids <- data$ids[[nm]][data$observed[[nm]]]
    support[ids] <- support[ids] + 1L
  }
  edge_rank <- ranks
  diag(edge_rank) <- NA_integer_
  min_rank <- suppressWarnings(min(edge_rank, na.rm = TRUE))
  if (!is.finite(min_rank)) min_rank <- NA_integer_
  d <- ncol(data$views[[1L]])
  list(
    names = nms,
    overlap_count = shared,
    overlap_rank = ranks,
    min_overlap_rank = as.integer(min_rank),
    adjacency = adj,
    components = comp$membership,
    n_components = comp$n,
    n_unobserved = as.integer(sum(support == 0L)),
    n_singleton = as.integer(sum(support == 1L)),
    relative_transforms_estimable = comp$n == 1L && is.finite(min_rank) && min_rank >= 1L,
    unique_relative_rotation = comp$n == 1L && is.finite(min_rank) && min_rank >= d
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
