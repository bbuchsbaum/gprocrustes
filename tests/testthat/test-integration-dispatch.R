test_that("automatic two-view routing preserves the declared consensus problem", {
  set.seed(391)
  X <- matrix(rnorm(24), 12, 2)
  views <- list(A = X, B = 2 * X + matrix(rnorm(24, sd = .1), 12, 2))
  cases <- list(
    list(transform = proc_orthogonal("O"), gauge = proc_gauge(scale = "preshape"),
         metric = proc_metric()),
    list(transform = proc_similarity("O"), gauge = proc_gauge(scale = "gower"),
         metric = proc_metric()),
    list(transform = proc_similarity("O"), gauge = proc_gauge(scale = "fixed_consensus"),
         metric = proc_metric()),
    list(transform = proc_orthogonal("O"), gauge = proc_gauge(),
         metric = proc_metric(landmark = list(A = seq_len(12), B = rev(seq_len(12)))))
  )
  for (args in cases) {
    auto <- do.call(gpa, c(list(data = views), args))
    reference <- do.call(gpa, c(list(data = views, solver = "gower_bcd"), args))
    expect_equal(auto$solver, "gower_bcd")
    expect_equal(auto$objective, reference$objective, tolerance = 1e-10)
    expect_equal(consensus(auto), consensus(reference), tolerance = 1e-10)
    if (args$gauge$scale == "preshape") {
      expect_equal(vapply(fitted(auto), function(Y) sum(Y^2), numeric(1)),
                   c(A = 1, B = 1), tolerance = 1e-10)
    }
    if (args$gauge$scale == "fixed_consensus") {
      expect_equal(sum(consensus(auto)^2), 1, tolerance = 1e-10)
    }
  }
})
