rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

infer_views <- function(n = 12, seed = 7, noise = 0.05) {
  set.seed(seed)
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  list(
    A = X + matrix(rnorm(n * 2, sd = noise), n, 2),
    B = X %*% rot2(0.4) + matrix(rnorm(n * 2, sd = noise), n, 2),
    C = X %*% rot2(-0.5) + matrix(rnorm(n * 2, sd = noise), n, 2)
  )
}

test_that("tangent residuals are orthogonal to the nuisance basis", {
  views <- infer_views(noise = 0.1)
  fit <- gpa(views, transform = proc_orthogonal("O"))
  tan <- tangent_coordinates(fit)
  M <- consensus(fit)
  B <- gprocrustes:::.gproc_nuisance_basis(M, fit$problem$transform)
  for (Z in tan) {
    expect_lt(max(abs(crossprod(B, as.vector(Z)))), 1e-6)
  }
  Z0 <- tangent_project(views$A, fit, as_residual = FALSE)
  expect_equal(dim(Z0), dim(M))
})

test_that("bootstrap summaries use gauge-aligned replicates only", {
  views <- infer_views()
  fit <- gpa(views, transform = proc_orthogonal("O"))
  inf <- infer(fit, method = "bootstrap", n_resamples = 8)
  expect_equal(inf$method, "bootstrap")
  expect_match(inf$gauge, "aligned")
  expect_equal(dim(inf$replicates)[[3]], 8L)
  expect_equal(dim(inf$se), dim(consensus(fit)))
  expect_true(all(is.finite(inf$se)))
  raw_gap <- sqrt(sum((inf$replicates[, , 1] - consensus(fit))^2))
  expect_true(is.finite(raw_gap))
})

test_that("permutation p-value is in (0, 1] and declares the null", {
  views <- infer_views(noise = 0.02)
  fit <- gpa(views, transform = proc_orthogonal("O"))
  inf <- infer(fit, method = "permutation", n_resamples = 8)
  expect_equal(inf$method, "permutation")
  expect_true(inf$p_value > 0 && inf$p_value <= 1)
  expect_true(any(grepl("correspondence", inf$assumptions)))
})

test_that("analytic inference uses model covariance, not Σ_S, by default", {
  views <- infer_views()
  fit <- gpa(views, transform = proc_orthogonal("O"))
  inf <- infer(fit, method = "analytic")
  expect_equal(inf$method, "analytic")
  expect_true(any(grepl("not automatically", inf$assumptions)))
  expect_false(identical(inf$optimality_status, "certified_global"))
  expect_true(inf$sigma2 > 0)
  expect_equal(dim(inf$se), dim(consensus(fit)))
})

test_that("landmark CV scores held-out discrepancy after gauge alignment", {
  views <- infer_views(n = 10)
  fit <- gpa(views, transform = proc_orthogonal("O"))
  set.seed(3)
  cv <- cross_validate(fit, folds = 2L)
  expect_s3_class(cv, "proc_cv")
  expect_true(is.finite(cv$total))
  expect_equal(length(cv$scores), 2L)
  cvk <- cross_validate(fit, level = "configuration")
  expect_equal(cvk$level, "configuration")
  expect_equal(length(cvk$scores), 3L)
  expect_true(is.finite(cvk$total))
})

test_that("proc_from_array is a complete-covering proc_data", {
  arr <- array(rnorm(5 * 2 * 3), c(5, 2, 3))
  dat <- proc_from_array(arr, view_names = c("A", "B", "C"))
  expect_s3_class(dat, "proc_data")
  expect_equal(names(dat$views), c("A", "B", "C"))
  expect_equal(nrow(dat$views$A), 5L)
})

test_that("uncertainty plot data comes from infer()", {
  views <- infer_views()
  fit <- gpa(views, transform = proc_orthogonal("O"))
  inf <- infer(fit, method = "analytic")
  dat <- plot_data(inf)
  expect_true(all(c("x", "y", "se_x", "se_y") %in% names(dat)))
  fit$inference <- inf
  dat2 <- plot_data(fit, type = "uncertainty")
  expect_equal(nrow(dat2), nrow(consensus(fit)))
})
