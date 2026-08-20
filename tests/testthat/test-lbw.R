rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

affine_views <- function(n = 12, seed = 5) {
  set.seed(seed)
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  list(
    A = X,
    B = X %*% matrix(c(1.2, 0.1, 0.2, 0.9), 2, 2) + matrix(c(0.4, -0.3), n, 2, byrow = TRUE),
    C = X %*% matrix(c(0.8, -0.2, 0.1, 1.1), 2, 2) + matrix(c(-0.2, 0.5), n, 2, byrow = TRUE)
  )
}

test_that("explain_solver selects the LBW eigen-engine for affine", {
  compiled <- compile_proc_problem(
    proc_data(affine_views()),
    transform = proc_affine(reference_covariance = 1)
  )
  plan <- explain_solver(compiled)
  expect_equal(plan$solver, "lbw_eigen")
  expect_match(plan$allowed_claim, "exact_closed_form")
})

test_that("noiseless affine GPA recovers a common shape under declared Lambda", {
  views <- affine_views()
  fit <- gpa(views, transform = proc_affine(reference_covariance = 1))
  expect_equal(fit$solver, "lbw_eigen")
  expect_equal(fit$optimality_status, "exact_closed_form")
  expect_equal(fit$backend, "eigencore")
  expect_true(isTRUE(fit$certificate$free_translation))
  expect_false(identical(fit$certificate$status, "certified_global"))
  M <- consensus(fit)
  expect_equal(unname(diag(crossprod(M))), c(1, 1), tolerance = 1e-6)
  expect_lt(fit$data_objective, 1e-6)
  expect_false(isTRUE(fit$folding$any_folded))
})

test_that("reference covariance is explicit and changes the consensus scale", {
  views <- affine_views()
  a <- gpa(views, transform = proc_affine(reference_covariance = 1))
  b <- gpa(views, transform = proc_affine(reference_covariance = c(4, 1)))
  expect_equal(a$reference_covariance, c(1, 1))
  expect_equal(b$reference_covariance, c(4, 1))
  expect_gt(abs(sum(consensus(a)^2) - sum(consensus(b)^2)), 0.5)
})

test_that("affine inverse exists and TPS inverse is refused", {
  views <- affine_views(n = 10)
  fit <- gpa(views, transform = proc_affine())
  inv <- inverse_proc_transform(transformations(fit)[["A"]])
  X <- views$A
  Y <- apply_proc_transform(transformations(fit)[["A"]], X)
  back <- apply_proc_transform(inv, Y)
  expect_equal(as.matrix(back), as.matrix(X), tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(is.finite(datum_space_error(fit)))
  tps <- gpa(views, transform = proc_tps(smoothness = 1))
  expect_error(
    inverse_proc_transform(transformations(tps)[["A"]]),
    class = "inverse_unavailable"
  )
  expect_error(datum_space_error(tps), class = "inverse_unavailable")
})

test_that("TPS with large smoothness approaches affine", {
  views <- affine_views(n = 10)
  aff <- gpa(views, transform = proc_affine())
  tps <- gpa(views, transform = proc_tps(smoothness = 1e6))
  expect_equal(tps$solver, "lbw_eigen")
  gap <- procrustes(consensus(tps), consensus(aff), transform = proc_orthogonal("O"))$objective
  expect_lt(gap, 0.2)
})

test_that("partial rows still compile to the LBW eigen-engine", {
  views <- affine_views(n = 10)
  mask <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  dat <- proc_data(
    views,
    observed = list(A = rep(TRUE, 10), B = mask, C = rep(TRUE, 10))
  )
  fit <- gpa(dat, transform = proc_affine())
  expect_equal(fit$solver, "lbw_eigen")
  expect_true(is.finite(fit$objective))
})

test_that("tune_smoothness uses held-out landmarks, not training residual", {
  set.seed(2)
  views <- affine_views(n = 9)
  views$C <- views$C + matrix(rnorm(18, sd = 0.15), 9, 2)
  cv <- tune_smoothness(
    views,
    transform = proc_tps(),
    smoothness = c(0.2, 5),
    folds = 2L
  )
  expect_equal(length(cv$scores), 2L)
  expect_true(all(is.finite(cv$scores)))
  expect_true(cv$best %in% c(0.2, 5))
})

test_that("deformation plot data reports Jacobian signs", {
  fit <- gpa(affine_views(n = 8), transform = proc_affine())
  dat <- plot_data(fit, type = "deformation")
  expect_true("jacobian_det" %in% names(dat))
  expect_true(all(is.finite(dat$jacobian_det)))
})
