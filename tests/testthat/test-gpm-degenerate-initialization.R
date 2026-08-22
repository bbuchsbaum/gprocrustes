.ez1_views <- function(scale = 1, epsilon = 0) {
  set.seed(20260822)
  X <- t(matrix(rnorm(24), 3)) * scale
  if (epsilon == 0) {
    Y <- -X
  } else {
    set.seed(20260823)
    Y <- -X + epsilon * matrix(rnorm(length(X)), nrow(X), ncol(X)) * scale
  }
  list(a = X, b = Y)
}

.ez1_warning_codes <- function(fit) {
  vapply(fit$warnings, `[[`, character(1), "code")
}

.ez1_rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

test_that("auto dispatch sends compatible two-view GPA to the exact oracle", {
  views <- .ez1_views()
  spec <- proc_orthogonal("O")
  pair <- procrustes(views$a, views$b, transform = spec)
  fit <- gpa(views, transform = spec)
  plan <- explain_solver(compile_proc_problem(proc_data(views), spec))

  expect_equal(plan$solver, "pairwise_polar")
  expect_equal(fit$solver, "pairwise_polar")
  expect_equal(fit$optimality_status, "exact_closed_form")
  expect_equal(fit$numerical_status, "converged")
  expect_lte(abs(fit$objective - pair$objective / 2), 1e-12 * sum(views$a^2))
  expect_lte(fit$objective, 1e-20 * sum(views$a^2))
  expect_gt(sum(consensus(fit)^2), 0.5 * sum(views$a^2))
  expect_lte(max(abs(fitted(fit)$a - fitted(fit)$b)), 1e-12 * max(abs(views$a)))
})

test_that("GPM defaults to spectral while Gower keeps its sequential default", {
  views2 <- .ez1_views()
  views4 <- list(a = views2$a, b = views2$b, c = views2$a, d = views2$b)
  gpm <- gpa(views4, transform = proc_orthogonal("O"))

  expect_null(gpa_control()$init)
  expect_equal(gpm$solver, "gpm")
  expect_equal(gpm$requested_init, "spectral")
  expect_equal(gpm$init, "spectral")
  expect_equal(gpm$starts, "spectral")
  expect_lte(gpm$objective, 1e-20 * sum(views2$a^2))
  expect_gt(sum(consensus(gpm)^2), 0.5 * sum(views2$a^2))
  expect_true(all(diff(gpm$history$objective) <= 1e-10 * sum(views2$a^2)))
  expect_lte(utils::tail(gpm$history$stationarity, 1), 1e-8 * sum(views2$a^2))

  set.seed(11)
  Z <- scale(matrix(rnorm(24), 12, 2), scale = FALSE)
  similarity_views <- list(
    a = Z,
    b = Z %*% .ez1_rot2(0.6),
    c = Z %*% .ez1_rot2(-0.4)
  )
  gower <- gpa(similarity_views, transform = proc_similarity())
  expect_equal(gower$solver, "gower_bcd")
  expect_equal(gower$init, "sequential")
})

test_that("GPM escapes a zero-consensus sequential start across scales", {
  for (data_scale in c(1e-6, 1, 1e6)) {
    views <- .ez1_views(scale = data_scale)
    energy <- sum(vapply(views, function(x) sum(x^2), numeric(1)))
    centered_energy <- sum(vapply(views, function(x) {
      sum(scale(x, scale = FALSE)^2)
    }, numeric(1)))
    fit <- gpa(
      views,
      transform = proc_orthogonal("O"),
      solver = "gpm",
      control = gpa_control(init = "sequential", nstart = 1L)
    )

    expect_equal(fit$requested_init, "sequential")
    expect_equal(fit$init, "spectral")
    expect_equal(fit$starts, c("sequential", "spectral"))
    expect_true("sequential" %in% fit$degenerate_starts)
    expect_true("zero_consensus_restart" %in% .ez1_warning_codes(fit))
    expect_gt(fit$start_objectives[["sequential"]], 0.9 * centered_energy)
    expect_equal(fit$objective, min(fit$start_objectives))
    expect_lte(fit$objective, 1e-12 * energy)
    expect_equal(fit$numerical_status, "converged")
    expect_true(
      fit$optimality_status %in% c("first_order_stationary", "certified_global")
    )
    expect_lte(utils::tail(fit$history$stationarity, 1), 1e-8 * energy)
  }
})

test_that("GPM nstart runs distinct starts and keeps the best exact objective", {
  views <- .ez1_views()
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    solver = "gpm",
    control = gpa_control(init = "sequential", nstart = 3L)
  )

  expect_equal(fit$requested_init, "sequential")
  expect_equal(fit$starts, c("sequential", "medoid", "spectral"))
  expect_named(fit$start_objectives, fit$starts)
  expect_gt(fit$start_objectives[["sequential"]], 1)
  expect_lte(fit$start_objectives[["medoid"]], 1e-20 * sum(views$a^2))
  expect_lte(fit$start_objectives[["spectral"]], 1e-20 * sum(views$a^2))
  expect_equal(fit$objective, min(fit$start_objectives))
  expect_equal(fit$objective, fit$start_objectives[[fit$init]])
})

test_that("near-antipodal sequential GPM does not trigger a false restart", {
  views <- .ez1_views(epsilon = 1e-12)
  oracle <- gpa(views, transform = proc_orthogonal("O"), solver = "pairwise_polar")
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    solver = "gpm",
    control = gpa_control(init = "sequential", nstart = 1L, certify = "never")
  )
  energy <- sum(vapply(views, function(x) sum(x^2), numeric(1)))

  expect_equal(fit$requested_init, "sequential")
  expect_equal(fit$init, "sequential")
  expect_equal(fit$starts, "sequential")
  expect_length(fit$degenerate_starts, 0L)
  expect_false("zero_consensus_restart" %in% .ez1_warning_codes(fit))
  expect_lte(abs(fit$objective - oracle$objective), 1e-10 * energy)
  expect_gt(sum(consensus(fit)^2), 0.5 * sum(views$a^2))
  expect_true(all(diff(fit$history$objective) <= 1e-10 * energy))
})

test_that("zero-energy data are not mistaken for a collapsed nonzero problem", {
  views <- list(
    a = matrix(0, 6, 3),
    b = matrix(0, 6, 3),
    c = matrix(0, 6, 3)
  )
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    solver = "gpm",
    control = gpa_control(init = "sequential", nstart = 1L)
  )

  expect_equal(fit$objective, 0)
  expect_equal(sum(consensus(fit)^2), 0)
  expect_equal(fit$init, "sequential")
  expect_equal(fit$starts, "sequential")
  expect_length(fit$degenerate_starts, 0L)
  expect_false("zero_consensus_restart" %in% .ez1_warning_codes(fit))
})
