rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

robust_views <- function(n = 18, seed = 4, outlier = TRUE) {
  set.seed(seed)
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  views <- list(
    A = X,
    B = X %*% rot2(0.45),
    C = X %*% rot2(-0.55)
  )
  if (isTRUE(outlier)) {
    views$C[1, ] <- views$C[1, ] + c(18, 16)
  }
  views
}

test_that("Huber IRLS down-weights a landmark-vector outlier", {
  views <- robust_views()
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    loss = proc_huber(k = 0.75),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-6, max_iterations = 40)
  )
  expect_equal(fit$solver, "irls")
  expect_equal(fit$optimality_status, "blockwise_stationary")
  expect_equal(certify(fit)$status, "unavailable")
  expect_true(all(diff(fit$history$objective) <= 1e-8))
  w <- fit$weights$robust$C
  expect_length(w, nrow(views$C))
  expect_lt(w[[1L]], min(w[-1L]) * 0.6)
})

test_that("Tukey IRLS can zero a landmark-vector outlier", {
  views <- robust_views()
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    loss = proc_tukey(c = 1.2),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-6, max_iterations = 40)
  )
  expect_equal(fit$solver, "irls")
  expect_equal(fit$optimality_status, "blockwise_stationary")
  expect_equal(fit$weights$robust$C[[1L]], 0)
  expect_true(all(fit$weights$robust$C[-1L] > 0))
})

test_that("explain_solver selects irls for landmark Huber and never GPM", {
  compiled <- compile_proc_problem(
    proc_data(robust_views(outlier = FALSE)),
    transform = proc_orthogonal("O"),
    loss = proc_huber()
  )
  plan <- explain_solver(compiled)
  expect_equal(plan$solver, "irls")
  expect_equal(plan$allowed_claim, "blockwise_stationary")
  expect_error(
    gpa(
      robust_views(outlier = FALSE),
      transform = proc_orthogonal("O"),
      loss = proc_huber(),
      solver = "gpm"
    ),
    class = "invalid_problem"
  )
})

test_that("cell masks are not missing rows and use MM", {
  X <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  B <- X %*% rot2(0.3)
  C <- X %*% rot2(-0.25)
  C_bad <- C
  C_bad[2, 1] <- 80
  cells <- list(
    A = matrix(TRUE, 4, 2),
    B = matrix(TRUE, 4, 2),
    C = matrix(TRUE, 4, 2)
  )
  cells$C[2, 1] <- FALSE
  dat <- proc_data(list(A = X, B = B, C = C_bad), cells = cells)
  expect_true(dat$observed$C[[2L]])
  expect_false(dat$cells$C[2, 1])
  expect_true(dat$cells$C[2, 2])

  fit <- gpa(
    dat,
    transform = proc_orthogonal("O"),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-8, max_iterations = 40)
  )
  expect_equal(fit$solver, "mm")
  expect_equal(fit$optimality_status, "first_order_stationary")
  expect_equal(certify(fit)$status, "unavailable")
  expect_true(all(diff(fit$history$objective) <= 1e-8))

  clean <- gpa(
    list(A = X, B = B, C = C),
    transform = proc_orthogonal("O")
  )
  dirty <- gpa(
    list(A = X, B = B, C = C_bad),
    transform = proc_orthogonal("O")
  )
  to_clean <- function(fit) {
    procrustes(consensus(fit), consensus(clean), transform = proc_orthogonal("O"))$objective
  }
  expect_lt(to_clean(fit), to_clean(dirty))
  expect_lt(to_clean(fit), 1e-4)
})

test_that("anisotropic coordinate precision is first-order stationary, not SVD", {
  views <- robust_views(outlier = FALSE)
  prec <- proc_covariance(
    coordinate = diag(c(1, 12)),
    kind = "superimposition",
    as_precision = TRUE
  )
  compiled <- compile_proc_problem(
    proc_data(views),
    transform = proc_orthogonal("O"),
    metric = proc_metric(precision = prec)
  )
  expect_true(compiled$cell_metric)
  expect_equal(explain_solver(compiled)$solver, "mm")
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    metric = proc_metric(precision = prec),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-6, max_iterations = 30)
  )
  expect_equal(fit$solver, "mm")
  expect_equal(fit$optimality_status, "first_order_stationary")
  expect_false(identical(fit$optimality_status, "exact_closed_form"))
  expect_equal(certify(fit)$status, "unavailable")
})

test_that("model covariance is unused by the solver", {
  views <- robust_views(outlier = FALSE, n = 12)
  base <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(tolerance = 1e-10)
  )
  model <- proc_covariance(
    coordinate = diag(c(1, 20)),
    kind = "model",
    as_precision = TRUE
  )
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    metric = proc_metric(precision = model),
    control = gpa_control(tolerance = 1e-10)
  )
  expect_equal(fit$solver, base$solver)
  expect_equal(fit$objective, base$objective, tolerance = 1e-8)
})

test_that("diagonal entity precision is the Goodall landmark-weight path", {
  views <- robust_views(outlier = FALSE, n = 12)
  w <- c(rep(1, 11), 0.05)
  ent <- proc_covariance(entity = diag(w), as_precision = TRUE, kind = "superimposition")
  compiled <- compile_proc_problem(
    proc_data(views),
    transform = proc_orthogonal("O"),
    metric = proc_metric(precision = ent)
  )
  expect_true(compiled$landmark_weighted)
  expect_false(compiled$cell_metric)
  expect_equal(explain_solver(compiled)$solver, "gower_bcd")
  ctrl <- gpa_control(accelerate = FALSE, tolerance = 1e-8)
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    metric = proc_metric(precision = ent),
    control = ctrl
  )
  expect_equal(fit$solver, "gower_bcd")
  direct <- gpa(
    views,
    transform = proc_orthogonal("O"),
    metric = proc_metric(landmark = w),
    control = ctrl
  )
  expect_equal(fit$objective, direct$objective, tolerance = 1e-8)
  expect_lt(
    procrustes(
      consensus(fit),
      consensus(direct),
      transform = proc_orthogonal("O")
    )$objective,
    1e-8
  )
})

test_that("influence, support, and weight plot data exist", {
  views <- robust_views()
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    loss = proc_huber(k = 0.8),
    control = gpa_control(accelerate = FALSE, max_iterations = 20)
  )
  inf <- plot_data(fit, type = "influence")
  expect_true(all(c("entity", "view", "residual", "weight", "influence") %in% names(inf)))
  sup <- plot_data(fit, type = "support")
  expect_equal(nrow(sup), nrow(views$A))
  expect_true(all(sup$n_observed_views == 3L))
  wt <- plot_data(fit, type = "weights")
  expect_true("weight" %in% names(wt))
})
