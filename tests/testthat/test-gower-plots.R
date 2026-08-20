test_that("plot_data returns tidy frames for the four Gower diagnostics", {
  views <- list(
    A = matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE),
    B = matrix(c(0.1, -0.1, 1.1, 0, 1, 1.1, -0.1, 1), 4, 2, byrow = TRUE),
    C = matrix(c(-0.1, 0.1, 0.9, 0.1, 1.1, 0.9, 0, 1.1), 4, 2, byrow = TRUE)
  )
  fit <- gpa(views, transform = proc_orthogonal("O"))
  ov <- plot_data(fit, type = "overlay")
  expect_true(all(c("entity", "view", "x", "y") %in% names(ov)))
  expect_true("consensus" %in% ov$view)
  rs <- plot_data(fit, type = "residuals")
  expect_true(all(c("entity", "view", "residual") %in% names(rs)))
  cv <- plot_data(fit, type = "convergence")
  expect_true("objective" %in% names(cv))
  dec <- plot_data(fit, type = "decomposition")
  expect_equal(nrow(dec), 3L)
})

test_that("autoplot returns ggplot objects when ggplot2 is installed", {
  skip_if_not_installed("ggplot2")
  views <- list(
    A = matrix(rnorm(16), 8, 2),
    B = matrix(rnorm(16), 8, 2),
    C = matrix(rnorm(16), 8, 2)
  )
  fit <- gpa(views, transform = proc_orthogonal("O"))
  for (tp in c("overlay", "residuals", "convergence", "decomposition",
               "influence", "support", "weights")) {
    p <- autoplot.gpa_fit(fit, type = tp)
    expect_s3_class(p, "ggplot")
  }
})
