conformance_dir <- function() {
  d <- system.file("conformance", package = "gprocrustes")
  if (!nzchar(d) || !dir.exists(d)) {
    d <- file.path(testthat::test_path("..", ".."), "inst", "conformance")
  }
  d
}

load_fixture <- function(id) {
  path <- file.path(conformance_dir(), "fixtures", paste0(id, ".json"))
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt <- gsub("\\bInfinity\\b", "null", txt)
  txt <- gsub("\\bNaN\\b", "null", txt)
  jsonlite::fromJSON(txt, simplifyVector = TRUE)
}

list_fixtures <- function(kind = NULL) {
  files <- list.files(file.path(conformance_dir(), "fixtures"), pattern = "\\.json$",
                      full.names = TRUE)
  ids <- sub("\\.json$", "", basename(files))
  if (is.null(kind)) {
    return(ids)
  }
  keep <- vapply(ids, function(id) {
    fx <- tryCatch(load_fixture(id), error = function(e) NULL)
    !is.null(fx) && identical(fx$kind, kind)
  }, logical(1))
  ids[keep]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

read_matrix <- function(obj, base = conformance_dir()) {
  if (is.null(obj)) {
    return(NULL)
  }
  if (is.list(obj) && !is.null(obj$path)) {
    raw <- readBin(
      file.path(base, obj$path),
      what = "double",
      n = as.integer(obj$nrow) * as.integer(obj$ncol),
      size = 8L,
      endian = "little"
    )
    return(matrix(raw, as.integer(obj$nrow), as.integer(obj$ncol),
                  byrow = identical(obj$order, "row-major")))
  }
  if (is.matrix(obj)) {
    return(obj)
  }
  if (is.numeric(obj)) {
    return(as.matrix(obj))
  }
  if (is.list(obj) && length(obj)) {
    return(do.call(rbind, lapply(obj, as.numeric)))
  }
  as.matrix(obj)
}

expect_close <- function(actual, expected, atol = 1e-12, rtol = 0, info = NULL) {
  if (is.null(expected)) {
    return(invisible())
  }
  testthat::expect_equal(actual, expected, tolerance = max(atol, rtol), info = info)
}

rotation_distance <- function(R, S) {
  sqrt(sum((R - S)^2))
}

fixture_view_list <- function(fx) {
  views <- fx$problem$views
  if (is.null(views)) {
    return(NULL)
  }
  parse_one <- function(name, matrix, ids = NULL, row_mask = NULL, cell_mask = NULL) {
    empty <- function(x) is.null(x) || (is.list(x) && !length(x))
    list(
      name = as.character(name),
      matrix = read_matrix(matrix),
      ids = if (empty(ids)) NULL else as.character(unlist(ids)),
      row_mask = if (empty(row_mask)) NULL else as.logical(unlist(row_mask)),
      cell_mask = if (empty(cell_mask)) NULL else {
        C <- read_matrix(cell_mask)
        storage.mode(C) <- "logical"
        C
      }
    )
  }
  if (is.data.frame(views)) {
    lapply(seq_len(nrow(views)), function(i) {
      parse_one(
        views$name[[i]],
        views$matrix[[i]],
        if ("ids" %in% names(views)) views$ids[[i]],
        if ("row_mask" %in% names(views)) views$row_mask[[i]],
        if ("cell_mask" %in% names(views)) views$cell_mask[[i]]
      )
    })
  } else {
    lapply(views, function(v) {
      parse_one(v$name, v$matrix, v$ids, v$row_mask, v$cell_mask)
    })
  }
}

fixture_proc_data <- function(fx) {
  vs <- fixture_view_list(fx)
  if (is.null(vs) || !length(vs)) {
    return(NULL)
  }
  mats <- lapply(vs, `[[`, "matrix")
  names(mats) <- vapply(vs, `[[`, character(1), "name")
  pick <- function(field) {
    out <- lapply(vs, `[[`, field)
    names(out) <- names(mats)
    if (all(vapply(out, is.null, logical(1)))) NULL else out
  }
  proc_data(mats, ids = pick("ids"), observed = pick("row_mask"), cells = pick("cell_mask"))
}

run_pairwise_fixture <- function(fx) {
  spec <- as_proc_transform(fx$problem$transform)
  X <- read_matrix(fx$problem$source)
  Y <- read_matrix(fx$problem$target)
  w <- fx$problem$row_weights
  if (is.null(X) || is.null(Y)) {
    views <- fx$problem$views
    if (is.null(views)) {
      return(NULL)
    }
    # views may be a data.frame-like parsed list
    return(NULL)
  }
  procrustes(X, Y, transform = spec, weights = w)
}
