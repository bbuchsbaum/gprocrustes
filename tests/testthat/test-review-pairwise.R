signed_permutation_oracle <- function(X, Y, w) {
  d <- ncol(X)
  permutations <- function(values) {
    if (length(values) <= 1L) return(list(values))
    unlist(lapply(seq_along(values), function(i) {
      lapply(permutations(values[-i]), function(rest) c(values[i], rest))
    }), recursive = FALSE)
  }
  sign_vectors <- function(d) {
    as.matrix(expand.grid(rep(list(c(-1, 1)), d)))
  }
  best <- list(objective = Inf)
  for (p in permutations(seq_len(d))) {
    for (sgn in split(sign_vectors(d), row(sign_vectors(d)))) {
      R <- matrix(0, d, d)
      R[cbind(p, seq_len(d))] <- as.numeric(sgn)
      objective <- sum(w * rowSums((X %*% R - Y)^2))
      if (objective < best$objective) {
        best <- list(R = R, objective = objective)
      }
    }
  }
  best
}

raw_similarity_oracle <- function(X, Y, w) {
  C <- crossprod(X, w * Y)
  decomp <- svd(C)
  R <- decomp$u %*% t(decomp$v)
  a <- sum(w * rowSums(X^2))
  s <- max(sum(R * C) / a, 0)
  list(R = R, s = s, objective = sum(w * rowSums((s * X %*% R - Y)^2)))
}

test_that("signed permutations use the declared source-to-target direction", {
  X3 <- matrix(c(
    1.3, -0.2, 0.8,
    -0.7, 1.1, 0.5,
    0.4, 0.9, -1.2,
    1.7, -0.5, 0.3,
    -1.1, 0.6, 1.4
  ), 5, 3, byrow = TRUE)
  X4 <- cbind(X3, c(-0.4, 1.2, 0.7, -1.5, 0.2))
  cases <- list(
    list(X = X3, p = c(2L, 3L, 1L), signs = c(-1, 1, -1), w = c(1, 3, 2, 4, 1)),
    list(X = X4, p = c(2L, 3L, 4L, 1L), signs = c(1, -1, 1, -1), w = c(2, 1, 4, 1, 3))
  )
  for (case in cases) {
    d <- ncol(case$X)
    R_true <- matrix(0, d, d)
    R_true[cbind(case$p, seq_len(d))] <- case$signs
    Y <- case$X %*% R_true + 1e-5 * matrix(seq_len(length(case$X)), ncol = d)
    oracle <- signed_permutation_oracle(case$X, Y, case$w)
    fit <- procrustes(case$X, Y, proc_signed_permutation(), weights = case$w)
    expect_equal(fit$transform$R, oracle$R, tolerance = 1e-12)
    expect_equal(fit$objective, oracle$objective, tolerance = 1e-12)
  }
})

test_that("translation-prohibited similarity uses raw weighted moments", {
  X <- matrix(c(-1, 0, 2, 1, 0, 3, 4, -2, -3, 1), 5, 2, byrow = TRUE)
  R <- matrix(c(cos(0.6), sin(0.6), -sin(0.6), cos(0.6)), 2, 2)
  Y <- 1.4 * X %*% R + matrix(c(7, -5), 5, 2, byrow = TRUE)
  w <- c(1, 4, 2, 3, 1)
  oracle <- raw_similarity_oracle(X, Y, w)
  fit <- procrustes(X, Y, proc_similarity(group = "O", translation = FALSE), weights = w)
  expect_equal(fit$transform$t, c(0, 0))
  expect_equal(fit$transform$R, oracle$R, tolerance = 1e-12)
  expect_equal(fit$transform$s, oracle$s, tolerance = 1e-12)
  expect_equal(fit$objective, oracle$objective, tolerance = 1e-10)
})

test_that("translation-enabled similarity is invariant to independent shifts", {
  X <- matrix(c(-1, 0, 2, 1, 0, 3, 4, -2, -3, 1), 5, 2, byrow = TRUE)
  R <- matrix(c(cos(0.6), sin(0.6), -sin(0.6), cos(0.6)), 2, 2)
  Y <- 1.4 * X %*% R + matrix(c(0.3, -0.8), 5, 2, byrow = TRUE)
  w <- c(1, 4, 2, 3, 1)
  spec <- proc_similarity(group = "O", translation = TRUE)
  base <- procrustes(X, Y, spec, weights = w)
  shifted <- procrustes(
    sweep(X, 2, c(4, -7), "+"),
    sweep(Y, 2, c(-3, 5), "+"),
    spec,
    weights = w
  )
  expect_equal(shifted$transform$R, base$transform$R, tolerance = 1e-12)
  expect_equal(shifted$transform$s, base$transform$s, tolerance = 1e-12)
  expect_equal(shifted$objective, base$objective, tolerance = 1e-12)
})

test_that("orthogonal fits retain their documented centered default", {
  X <- matrix(c(-1, 0, 2, 1, 0, 3, 4, -2, -3, 1), 5, 2, byrow = TRUE)
  Y <- sweep(X, 2, c(2.5, -4), "+")
  fit <- procrustes(X, Y, proc_orthogonal("O"))
  expect_lt(fit$objective, 1e-12)
  expect_equal(fit$transform$t, c(2.5, -4), tolerance = 1e-12)
})

test_that("pairwise polar rejects LBW-only transformation families", {
  X <- matrix(seq_len(12), 6, 2)
  for (spec in list(proc_affine(), proc_lbw("affine"), proc_tps())) {
    expect_error(procrustes(X, X, spec), class = "invalid_problem")
  }
})
