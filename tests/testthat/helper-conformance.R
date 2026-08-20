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
