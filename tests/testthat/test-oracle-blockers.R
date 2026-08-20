rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

test_that("matrix-free GPM never emits certified_global from Ritz values", {
  set.seed(2)
  X <- scale(matrix(rnorm(28), 14, 2), scale = FALSE)
  views <- list(A = X, B = X %*% rot2(0.8), C = X %*% rot2(-1.1), D = X %*% rot2(0.25))
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(backend = "matrix_free", init = "spectral", certify = "auto")
  )
  cert <- certify(fit)
  expect_false(identical(cert$status, "certified_global"))
  expect_equal(cert$status, "not_certified")
  expect_match(cert$reason, "Ritz|lower bound")
})

test_that("a missed negative eigendirection makes a Ritz value look PSD", {
  A <- diag(c(2, 1, -0.05))
  U <- cbind(c(1, 0, 0), c(0, 1, 0))
  ritz <- eigen(crossprod(U, A %*% U), symmetric = TRUE)$values
  ev <- eigen(A, symmetric = TRUE)$values
  expect_gt(min(ritz), 0)
  expect_lt(min(ev), 0)
})

test_that("pairwise orthogonal stores t-star so the transform reproduces the objective", {
  X <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  Y <- X + matrix(c(0.4, -0.25), 4, 2, byrow = TRUE)
  pair <- procrustes(X, Y, transform = proc_orthogonal("O"))
  expect_lt(pair$objective, 1e-12)
  got <- apply_proc_transform(pair$transform, X)
  expect_equal(got, Y, tolerance = 1e-12)
  expect_equal(pair$transform$t, c(0.4, -0.25), tolerance = 1e-12)
  raw <- procrustes(X, Y, transform = proc_orthogonal("O"), center = FALSE)
  expect_gt(raw$objective, 0.1)
  expect_equal(
    sum((apply_proc_transform(raw$transform, X) - Y)^2),
    raw$objective,
    tolerance = 1e-12
  )
})

test_that("SO uniqueness is not the O(d) rank==d rule", {
  X <- matrix(c(1, 0, 0, 1, 0, 0), 3, 2, byrow = TRUE)
  Y <- X
  so <- procrustes(X, Y, transform = proc_orthogonal("SO"))
  expect_true(so$transform_unique)
  expect_equal(so$unidentified_subspace_dimension, 0L)
  C <- matrix(c(1, 0, 0, 0), 2, 2)
  pol_so <- polar_factor(C, group = "SO")
  pol_o <- polar_factor(C, group = "O")
  expect_true(pol_so$transform_unique)
  expect_false(pol_o$transform_unique)
  expect_equal(pol_o$unidentified, 1L)
})

test_that("signed-permutation ties are reported", {
  X <- matrix(1, 2, 2)
  Y <- matrix(1, 2, 2)
  pair <- match_components(X, Y)
  expect_false(pair$transform_unique)
})

test_that("pairwise gpa uses row_map incidence, not local masks", {
  A <- matrix(c(0, 0, 1, 0, 1, 1), 3, 2, byrow = TRUE)
  B <- A[c(3, 1, 2), ] %*% rot2(0.4)
  dat <- proc_data(
    list(A = A, B = B),
    ids = list(A = c("a", "b", "c"), B = c("c", "a", "b"))
  )
  fit <- gpa(dat, transform = proc_orthogonal("O"), solver = "pairwise_polar", anchor = "A")
  expect_equal(fit$solver, "pairwise_polar")
  Ya <- as.matrix(aligned(fit)[["A"]])
  Yb <- as.matrix(aligned(fit)[["B"]])
  expect_lt(max(abs(Ya - Yb[match(c("a", "b", "c"), c("c", "a", "b")), ])), 1e-8)
  M <- consensus(fit)
  expect_equal(nrow(M), 3L)
})

test_that("compose does not type a reflection as SO", {
  so <- proc_identity(proc_orthogonal("SO"), 2)
  Href <- diag(c(1, -1))
  X <- matrix(c(1, 0, 0, 1, 0.2, 0.3), 3, 2, byrow = TRUE)
  o <- procrustes(X, X %*% Href, transform = proc_orthogonal("O"))$transform
  tr <- compose_proc_transform(so, o)
  expect_equal(tr$spec$group, "O")
  expect_lt(det(tr$R), 0)
})

test_that("forced inapplicable solvers are rejected", {
  views <- list(
    A = scale(matrix(rnorm(20), 10, 2), scale = FALSE),
    B = scale(matrix(rnorm(20), 10, 2), scale = FALSE),
    C = scale(matrix(rnorm(20), 10, 2), scale = FALSE)
  )
  expect_error(
    gpa(views, transform = proc_orthogonal("O"), loss = proc_huber(), solver = "pairwise_polar"),
    class = "invalid_problem"
  )
  expect_error(
    gpa(views, transform = proc_orthogonal("O"), loss = proc_huber(), solver = "lbw_eigen"),
    class = "invalid_problem"
  )
  expect_error(
    gpa(views, transform = proc_orthogonal("O"), loss = proc_huber(), solver = "gower_bcd"),
    class = "invalid_problem"
  )
  expect_error(
    gpa(views, gauge = list(scale = "none"), transform = proc_orthogonal("O")),
    class = "invalid_problem"
  )
})

test_that("IRLS multiplies robust q(r) by the declared landmark weights", {
  set.seed(4)
  n <- 16
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  views <- list(A = X, B = X %*% rot2(0.4), C = X %*% rot2(-0.3))
  views$C[1, ] <- views$C[1, ] + c(12, 10)
  w <- c(2, rep(1, n - 1))
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    metric = proc_metric(landmark = w),
    loss = proc_huber(k = 0.8),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-6, max_iterations = 40)
  )
  expect_equal(fit$solver, "irls")
  expect_true(!is.null(fit$weights$effective_n))
  expect_gt(fit$weights$landmark$A[[2L]] / fit$weights$landmark$A[[1L]], 0.2)
})

test_that("anisotropic MM evaluates e' Q e, not e' Q^2 e", {
  E <- matrix(c(1, 0, 0, 2), 2, 2, byrow = TRUE)
  Q <- diag(c(1, 4))
  got <- gprocrustes:::.gproc_row_quadratic(E, Q)
  expect_equal(got, c(1, 16))
  wrong <- rowSums((E %*% Q) * (E %*% Q))
  expect_false(isTRUE(all.equal(got, wrong)))
})

test_that("a rejected MM proposal is not first-order stationarity", {
  X <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  cells <- list(
    A = matrix(TRUE, 4, 2),
    B = matrix(TRUE, 4, 2)
  )
  dat <- proc_data(list(A = X, B = X %*% rot2(0.2)), cells = cells)
  fit <- gpa(
    dat,
    transform = proc_orthogonal("O"),
    control = gpa_control(accelerate = FALSE, max_iterations = 1L, tolerance = 1e-20)
  )
  expect_false(identical(fit$optimality_status, "first_order_stationary") &&
                 identical(fit$numerical_status, "stalled"))
})

test_that("LBW reduced operator includes the left weight factor", {
  Phi <- matrix(c(1, 0, 1, 1, 0, 1), 3, 2, byrow = TRUE)
  w <- c(1, 0.2, 0)
  A <- crossprod(Phi, w * Phi) + 0.1 * diag(2)
  Hi_wrong <- Phi %*% gprocrustes:::.gproc_psd_solve(A, t(Phi * w))
  Hi_right <- (Phi * w) %*% gprocrustes:::.gproc_psd_solve(A, t(Phi * w))
  expect_gt(max(abs(Hi_wrong - Hi_right)), 1e-8)
  expect_lt(max(abs(Hi_right - t(Hi_right))), 1e-12)
})

test_that("align_gauge refuses a reflection on an SO fit", {
  set.seed(1)
  X <- scale(matrix(rnorm(24), 12, 2), scale = FALSE)
  views <- list(A = X, B = X %*% rot2(0.5), C = X %*% rot2(-0.4))
  fit <- gpa(views, transform = proc_orthogonal("SO"))
  Q <- diag(c(1, -1))
  expect_error(align_gauge(fit, Q), class = "invalid_problem")
})

test_that("materialized aligned views stay consistent after canonicalize", {
  set.seed(3)
  X <- scale(matrix(rnorm(20), 10, 2), scale = FALSE)
  views <- list(A = X, B = X %*% rot2(0.7), C = X %*% rot2(-0.2))
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(keep_aligned = "materialize")
  )
  can <- canonicalize(fit)
  Ya <- aligned(can)[["A"]]
  expect_true(is.matrix(Ya))
  expect_equal(Ya, apply_proc_transform(transformations(can)[["A"]], views$A), tolerance = 1e-10)
})
