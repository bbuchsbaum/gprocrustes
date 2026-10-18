test_that("orthogonal preshape solves the normalized problem on every requested path", {
  set.seed(71)
  X <- scale(matrix(rnorm(24), 12, 2), scale = FALSE)
  views <- list(A = X, B = 2 * X, C = 3 * X)
  for (solver in c("auto", "gower_bcd")) {
    fit <- gpa(views, transform = proc_orthogonal("O"),
               gauge = proc_gauge(scale = "preshape"), solver = solver)
    expect_equal(vapply(fitted(fit), function(Y) sqrt(sum(Y^2)), numeric(1)),
                 c(A = 1, B = 1, C = 1), tolerance = 1e-10)
    expect_lt(fit$objective, 1e-20)
    expect_equal(unname(consensus(fit) %*% t(consensus(fit))),
                 unname(tcrossprod(X) / sum(X^2)), tolerance = 1e-10)
  }
  expect_error(gpa(views, transform = "O", gauge = proc_gauge(scale = "preshape"),
                   solver = "gpm"), class = "invalid_problem")
  expect_error(gpa(views[1:2], transform = "O", gauge = proc_gauge(scale = "preshape"),
                   solver = "pairwise_polar"), class = "invalid_problem")
})

test_that("preshape is invariant to view translation, positive scale and row order", {
  set.seed(10)
  X <- matrix(rnorm(30), 15, 2)
  v <- list(A = X, B = sweep(4 * X, 2, c(12, -8), "+"), C = 2 * X)
  f <- gpa(v, transform = "SO", gauge = proc_gauge(scale = "preshape"))
  expect_lt(f$objective, 1e-20)
  expect_equal(colSums(consensus(f)), c(0, 0), tolerance = 1e-12)
  zero <- lapply(v, function(Y) Y * 0)
  expect_error(gpa(zero, transform = "SO", gauge = proc_gauge(scale = "preshape")),
               class = "invalid_problem")
})

test_that("centering cannot override a prohibited similarity translation", {
  X<-matrix(c(0,0,1,0,0,1),3,2,byrow=TRUE)
  expect_error(procrustes(X,X+1,proc_similarity("SO",translation=FALSE),center=TRUE),
               class="invalid_problem")
})
