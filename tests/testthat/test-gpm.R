rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

gpm_views <- function(n = 14, seed = 2) {
  set.seed(seed)
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  list(A = X, B = X %*% rot2(0.8), C = X %*% rot2(-1.1), D = X %*% rot2(0.25))
}

test_that("noiseless O(d) GPM is certified global", {
  views <- gpm_views()
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(init = "spectral", certify = "auto", tolerance = 1e-12)
  )
  expect_equal(fit$solver, "gpm")
  expect_lt(fit$objective, 1e-10)
  cert <- certify(fit)
  expect_equal(cert$status, "certified_global")
  expect_equal(fit$optimality_status, "certified_global")
  expect_true(is.finite(cert$r_dual))
  expect_lt(cert$r_dual, 1e-6)
  expect_true(cert$lambda_min > -1e-6)
  expect_true(cert$uniqueness_modulo_O)
  expect_equal(cert$expected_nullity, 2L)
})

test_that("SO(d) GPM does not reuse the O(d) certificate", {
  views <- gpm_views()
  fit <- gpa(views, transform = proc_orthogonal("SO"))
  expect_equal(fit$solver, "gpm")
  cert <- certify(fit)
  expect_equal(cert$status, "unavailable")
  expect_match(cert$reason, "SO")
  expect_true(fit$optimality_status %in% c("first_order_stationary", "not_converged"))
  expect_false(identical(fit$optimality_status, "certified_global"))
})

test_that("certify='never' leaves a GPM fit uncertified", {
  views <- gpm_views()
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(certify = "never")
  )
  expect_equal(certify(fit)$status, "unavailable")
  expect_equal(fit$optimality_status, "first_order_stationary")
})

test_that("matrix-free CS matches the formed block Gram", {
  views <- gpm_views(n = 10, seed = 9)
  dat <- proc_data(views)
  prep <- gprocrustes:::.gpm_prepare_views(dat)
  nms <- names(prep)
  d <- 2L
  n <- 10L
  alphas <- stats::setNames(rep(1, 4), nms)
  Rs <- lapply(prep, function(v) rot2(0.2))
  names(Rs) <- nms
  B <- gprocrustes:::.gpm_apply_C(prep, alphas, Rs, n, complete = TRUE)
  C <- gprocrustes:::.gpm_form_C(prep, alphas, d, n)
  S <- gprocrustes:::.gpm_stack(Rs)
  CS <- C %*% S
  got <- gprocrustes:::.gpm_stack(B)
  expect_equal(got, CS, tolerance = 1e-10)
})

test_that("GPM history reports stationarity and orthogonality", {
  fit <- gpa(gpm_views(), transform = proc_orthogonal("O"))
  h <- fit$history
  expect_true(all(c("stationarity", "orthogonality_residual", "consensus_change") %in% names(h)))
  expect_lt(utils::tail(h$stationarity, 1), 1e-6)
  expect_lt(utils::tail(h$orthogonality_residual, 1), 1e-10)
})

test_that("row-masked orthogonal problems fall back to Gower BCD", {
  set.seed(4)
  n <- 12
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  mask <- rep(TRUE, n)
  mask[c(2, 11)] <- FALSE
  dat <- proc_data(
    list(A = X, B = X %*% rot2(0.4), C = X %*% rot2(-0.7)),
    observed = list(A = rep(TRUE, n), B = mask, C = rep(TRUE, n))
  )
  fit <- gpa(dat, transform = proc_orthogonal("O"))
  expect_equal(fit$solver, "gower_bcd")
  expect_lt(fit$objective, 1e-8)
})

test_that("forced gower_bcd remains available on orthogonal problems", {
  fit <- gpa(
    gpm_views(),
    transform = proc_orthogonal("O"),
    solver = "gower_bcd",
    control = gpa_control(accelerate = FALSE)
  )
  expect_equal(fit$solver, "gower_bcd")
  expect_equal(fit$optimality_status, "blockwise_stationary")
  expect_equal(certify(fit)$status, "unavailable")
})

test_that("overlap graph reports ranks and singleton entities", {
  dat <- proc_data(
    list(
      A = matrix(rnorm(10), 5, 2),
      B = matrix(rnorm(8), 4, 2),
      C = matrix(rnorm(8), 4, 2)
    ),
    ids = list(
      A = c("a", "b", "c", "d", "e"),
      B = c("a", "b", "c", "d"),
      C = c("c", "d", "e", "f")
    )
  )
  ov <- overlap_graph(dat)
  expect_equal(ov$n_components, 1L)
  expect_true(is.matrix(ov$overlap_rank))
  expect_true(ov$n_singleton >= 1L)
  expect_true(ov$relative_transforms_estimable)
})

test_that("certificate fixtures document the required fields", {
  skip_if_not_installed("jsonlite")
  fx <- load_fixture("cert-gpm-fixed-point-shape")
  expect_true(all(c(
    "report r_dual = ||CS-Lambda S|| / (1+||CS||)",
    "report a lower bound on lambda_min(Lambda-C)",
    "uniqueness if lambda_{d+1}(Lambda-C)>0"
  ) %in% fx$expected$laws))
  so <- load_fixture("cert-od-unavailable-default")
  expect_true("do not reuse for SO(d), robust losses, or cell masks" %in% so$expected$laws)
})

test_that("matrix_free backend does not change a noiseless recovery", {
  views <- gpm_views(n = 8, seed = 6)
  dense <- gpa(views, transform = proc_orthogonal("O"),
               control = gpa_control(backend = "dense", init = "spectral"))
  mf <- gpa(views, transform = proc_orthogonal("O"),
            control = gpa_control(backend = "matrix_free", init = "spectral"))
  expect_equal(mf$backend, "matrix_free")
  expect_equal(dense$objective, mf$objective, tolerance = 1e-8)
})
