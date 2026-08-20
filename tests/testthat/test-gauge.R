test_that("quotient_distance is zero after a pure gauge rotation", {
  set.seed(1)
  X <- scale(matrix(rnorm(24), 12, 2), scale = FALSE)
  Q <- matrix(c(0, -1, 1, 0), 2, 2)
  expect_lt(quotient_distance(X, X %*% Q), 1e-10)
  expect_gt(quotient_distance(X, X %*% diag(c(1.4, 0.7))), 0.1)
})

test_that("principal angles and subspace distance see a shared plane", {
  U <- diag(3)[, 1:2]
  V <- cbind(c(0, 1, 0), c(1, 0, 0))
  ang <- principal_angles(U, V)
  expect_equal(length(ang), 2L)
  expect_lt(max(ang), 1e-10)
  expect_lt(subspace_distance(U, V), 1e-10)
})

test_that("horizontal plus vertical recovers the perturbation", {
  set.seed(2)
  M <- scale(matrix(rnorm(20), 10, 2), scale = FALSE)
  Z <- matrix(rnorm(20), 10, 2)
  H <- horizontal_component(Z, M)
  V <- vertical_component(Z, M)
  expect_equal(H + V, Z, tolerance = 1e-8)
  B <- gprocrustes:::.gproc_nuisance_basis(M, proc_orthogonal("O"))
  expect_lt(max(abs(crossprod(B, as.vector(H)))), 1e-6)
})

test_that("consensus gauge is presentation, not a new estimate", {
  views <- list(
    A = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
    B = matrix(c(0.1, 0, 1.1, 0.1, 1, 1.1, 0, 0.9), 4, 2, byrow = TRUE),
    C = matrix(c(-0.1, 0.1, 0.9, 0, 1.1, 0.9, 0.1, 1), 4, 2, byrow = TRUE)
  )
  fit <- gpa(views, transform = proc_orthogonal("O"))
  native <- consensus(fit, gauge = "native")
  can <- consensus(fit, gauge = "canonical")
  expect_equal(native, fit$consensus)
  expect_lt(quotient_distance(native, can), 1e-8)
  expect_equal(fit$objective, canonicalize(fit)$objective, tolerance = 1e-12)
  blocks <- tied_axis_blocks(native)
  expect_true(is.list(blocks))
})
