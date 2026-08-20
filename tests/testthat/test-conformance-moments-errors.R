test_that("centered cross-products match moment fixtures", {
  for (id in list_fixtures("moments")) {
    fx <- load_fixture(id)
    X <- read_matrix(fx$problem$source)
    Y <- read_matrix(fx$problem$target)
    w <- fx$problem$row_weights
    C <- centered_crossprod(X, Y, w)
    Cexp <- read_matrix(fx$expected$crossproduct)
    expect_equal(C, Cexp, tolerance = fx$expected$objective_atol %||% 1e-12, info = id)
    expect_equal(centered_trace(X, w), fx$expected$a, tolerance = 1e-12, info = paste(id, "a"))
    expect_equal(centered_trace(Y, w), fx$expected$b, tolerance = 1e-12, info = paste(id, "b"))
    Xs <- Matrix::Matrix(X, sparse = TRUE)
    Ys <- Matrix::Matrix(Y, sparse = TRUE)
    expect_equal(centered_crossprod(Xs, Ys, w), Cexp, tolerance = 1e-12, info = paste(id, "sparse"))
  }
})

test_that("error fixtures raise the documented codes", {
  expect_error(
    procrustes(matrix(c(1, NA, 0, 1), 2, 2), diag(2)),
    class = "nonfinite_values"
  )
  expect_error(
    procrustes(matrix(c(1, Inf, 0, 1), 2, 2), diag(2)),
    class = "nonfinite_values"
  )
  expect_error(
    proc_data(list(A = diag(2)), ids = list(A = c("a", "a"))),
    class = "duplicate_entity_ids"
  )
  expect_error(
    procrustes(matrix(rnorm(10), 5, 2), matrix(rnorm(15), 5, 3)),
    class = "dimension_mismatch"
  )
  dat <- proc_data(
    list(A = matrix(rnorm(8), 4, 2), B = matrix(rnorm(8), 4, 2)),
    ids = list(A = c("a", "b", "c", "d"), B = c("e", "f", "g", "h"))
  )
  expect_error(gpa(dat, transform = "O"), class = "disconnected_overlap_graph")
  expect_error(
    gpa(list(matrix(rnorm(8), 4, 2), matrix(rnorm(8), 4, 2)),
        transform = "O",
        metric = proc_metric(configuration = c(0, 0))),
    class = "zero_total_configuration_weight"
  )
})

test_that("overlap graph reports connectivity", {
  dat <- proc_data(
    list(
      A = matrix(rnorm(6), 3, 2),
      B = matrix(rnorm(6), 3, 2),
      C = matrix(rnorm(6), 3, 2)
    ),
    ids = list(A = c("a", "b", "c"), B = c("c", "d", "e"), C = c("e", "f", "a"))
  )
  ov <- overlap_graph(dat)
  expect_equal(ov$n_components, 1L)
})
