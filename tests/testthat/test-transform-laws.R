test_that("identity, compose, and inverse round-trip", {
  spec <- proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic")
  X <- matrix(rnorm(20), 10, 2)
  Y <- X %*% matrix(c(0, -1, 1, 0), 2, 2) * 1.3 + matrix(c(0.4, -0.2), 10, 2, byrow = TRUE)
  pair <- procrustes(X, Y, transform = spec)
  tr <- pair$transform
  id <- proc_identity(spec, 2)
  expect_equal(apply_proc_transform(id, X), X)
  expect_equal(
    apply_proc_transform(compose_proc_transform(tr, inverse_proc_transform(tr)), X),
    X,
    tolerance = 1e-10
  )
  expect_equal(
    apply_proc_transform(tr, X),
    apply_proc_transform(compose_proc_transform(id, tr), X)
  )
})

test_that("O preserves norms and SO has det +1", {
  X <- matrix(rnorm(24), 8, 3)
  Q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  Y <- X %*% Q
  so <- procrustes(X, Y, transform = proc_orthogonal("SO"))
  expect_equal(det(so$transform$R), 1, tolerance = 1e-10)
  expect_equal(sqrt(sum((X %*% so$transform$R)^2)), sqrt(sum(X^2)), tolerance = 1e-10)
  Href <- diag(3)
  Href[3, 3] <- -1
  o <- procrustes(X, X %*% Href, transform = proc_orthogonal("O"))
  expect_lt(o$objective, 1e-10)
  so_ref <- procrustes(X, X %*% Href, transform = proc_orthogonal("SO"))
  expect_gt(so_ref$objective, o$objective - 1e-12)
  expect_equal(det(so_ref$transform$R), 1, tolerance = 1e-10)
})

test_that("gauge right-action leaves the pairwise residual unchanged", {
  X <- matrix(rnorm(16), 8, 2)
  Y <- X %*% matrix(c(cos(0.4), sin(0.4), -sin(0.4), cos(0.4)), 2, 2)
  Q <- matrix(c(cos(1.1), sin(1.1), -sin(1.1), cos(1.1)), 2, 2)
  a <- procrustes(X, Y, transform = "O")
  b <- procrustes(X %*% Q, Y %*% Q, transform = "O")
  expect_equal(a$objective, b$objective, tolerance = 1e-12)
})

test_that("energy identity holds for aligned views", {
  A1 <- matrix(rnorm(10), 5, 2)
  A2 <- A1 + 0.1 * matrix(rnorm(10), 5, 2)
  A3 <- A1 + 0.2 * matrix(rnorm(10), 5, 2)
  M <- (A1 + A2 + A3) / 3
  dec <- decompose_energy(list(A1 = A1, A2 = A2, A3 = A3), M, c(1, 1, 1))
  expect_equal(dec$total_energy, dec$consensus_energy + dec$residual_energy, tolerance = 1e-12)
})

test_that("lazy aligned views match materialized transforms", {
  X <- matrix(rnorm(20), 10, 2)
  Y <- X %*% matrix(c(0, -1, 1, 0), 2, 2)
  fit <- gpa(list(source = X, target = Y), transform = proc_orthogonal("SO"))
  lazy <- aligned(fit)[["source"]]
  expect_s3_class(lazy, "proc_aligned_view")
  expect_equal(as.matrix(lazy), apply_proc_transform(transformations(fit)[["source"]], X))
  expect_equal(lazy[1:3, ], apply_proc_transform(transformations(fit)[["source"]], X[1:3, , drop = FALSE]))
})

test_that("d = 1 distinguishes O and SO", {
  X <- matrix(c(1, 2, 3, 4), 4, 1)
  o <- procrustes(X, -X, transform = "O")
  so <- procrustes(X, -X, transform = "SO")
  expect_equal(o$transform$R[1, 1], -1)
  expect_equal(so$transform$R[1, 1], 1)
  expect_equal(o$objective, 0, tolerance = 1e-12)
  expect_gt(so$objective, 1)
})
