square_views <- function() {
  list(
    A = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
    B = matrix(c(0.1, 0, 1.1, 0.1, 1, 1.1, 0, 0.9), 4, 2, byrow = TRUE),
    C = matrix(c(-0.1, 0.1, 0.9, 0, 1.1, 0.9, 0.1, 1), 4, 2, byrow = TRUE)
  )
}

test_that("tidy/glance/augment report existing fit quantities", {
  fit <- gpa(square_views(), transform = proc_orthogonal("O"))
  td <- tidy.gpa_fit(fit)
  expect_equal(nrow(td), 3L)
  expect_true(all(c("configuration", "family", "s", "orthogonality", "t1", "t2") %in% names(td)))
  expect_equal(td$configuration, names(transformations(fit)))
  expect_equal(td$s, unname(vapply(transformations(fit), function(tr) tr$s, numeric(1))))

  gl <- glance.gpa_fit(fit)
  expect_equal(nrow(gl), 1L)
  expect_equal(gl$n_configurations, 3L)
  expect_equal(gl$n_entities, 4L)
  expect_equal(gl$objective, fit$objective, tolerance = 1e-12)
  expect_equal(gl$solver, fit$solver)
  expect_equal(gl$optimality_status, fit$optimality_status)

  aug <- augment.gpa_fit(fit)
  expect_equal(nrow(aug), 12L)
  expect_true(all(c("entity", "configuration", "residual", "weight", "r1", "r2") %in% names(aug)))
  res <- residuals(fit, level = "landmark")
  expect_equal(aug$residual[aug$configuration == "A"], sqrt(rowSums(res$A^2)))

  cfg <- augment.gpa_fit(fit, level = "configuration")
  expect_equal(cfg$residual, unname(residuals(fit, level = "configuration")))
})

test_that("Huber weights appear in augment()", {
  set.seed(4)
  X <- scale(matrix(rnorm(36), 18, 2), scale = FALSE)
  views <- list(A = X, B = X, C = X)
  views$C[1, ] <- views$C[1, ] + c(18, 16)
  fit <- gpa(
    views,
    transform = proc_orthogonal("O"),
    loss = proc_huber(k = 0.75),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-6, max_iterations = 40)
  )
  aug <- augment.gpa_fit(fit)
  wC <- aug$weight[aug$configuration == "C"]
  expect_equal(wC, fit$weights$robust$C)
  expect_lt(wC[[1L]], min(wC[-1L]))
})

test_that("broom generics dispatch when broom is installed", {
  skip_if_not_installed("broom")
  fit <- gpa(square_views(), transform = proc_orthogonal("O"))
  expect_equal(broom::tidy(fit), tidy.gpa_fit(fit))
  expect_equal(broom::glance(fit), glance.gpa_fit(fit))
  expect_equal(broom::augment(fit), augment.gpa_fit(fit))
})
