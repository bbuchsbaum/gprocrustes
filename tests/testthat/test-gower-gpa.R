rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

make_views <- function(n = 12, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * 2), n, 2)
  X <- scale(X, scale = FALSE)
  list(
    A = X,
    B = X %*% rot2(0.7) * 1.4 + matrix(c(0.3, -0.2), n, 2, byrow = TRUE),
    C = X %*% rot2(-0.4) * 0.8 + matrix(c(-0.1, 0.5), n, 2, byrow = TRUE)
  )
}

test_that("K>=3 orthogonal GPA recovers a common shape up to gauge", {
  set.seed(11)
  n <- 16
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  views <- list(
    A = X,
    B = X %*% rot2(0.9),
    C = X %*% rot2(-1.2)
  )
  fit <- gpa(
    views,
    transform = proc_orthogonal("SO"),
    gauge = proc_gauge(scale = "none"),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-10)
  )
  expect_equal(fit$solver, "gpm")
  expect_true(fit$optimality_status %in% c("first_order_stationary", "certified_global"))
  expect_match(explain_solver(fit)$allowed_claim, "first_order_stationary")
  expect_lt(fit$objective, 1e-8)
  aligned_a <- as.matrix(aligned(fit)[["A"]])
  pair <- procrustes(aligned_a, X, transform = proc_orthogonal("O"))
  expect_lt(pair$objective, 1e-8)
})

test_that("Gower history is monotone and has the documented columns", {
  views <- make_views()
  fit <- gpa(
    views,
    transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
    gauge = proc_gauge(scale = "gower"),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-10)
  )
  h <- fit$history
  expect_true(is.data.frame(h))
  expect_true(nrow(h) >= 1)
  expect_true(all(c(
    "iteration", "objective", "relative_change", "consensus_change",
    "stationarity", "orthogonality_residual", "elapsed_time",
    "accepted_acceleration"
  ) %in% names(h)))
  expect_true(all(diff(h$objective) <= 1e-8))
  expect_equal(fit$numerical_status, "converged")
  expect_lt(utils::tail(h$stationarity, 1), 1e-4)
})

test_that("energy identity holds after Gower GPA", {
  views <- make_views()
  fit <- gpa(views, transform = proc_orthogonal("O"), gauge = proc_gauge(scale = "none"))
  dec <- decompose(fit)
  expect_equal(dec$kind, "energy")
  expect_equal(dec$total_energy, dec$consensus_energy + dec$residual_energy, tolerance = 1e-10)
  expect_equal(dec$residual_energy, fit$objective, tolerance = 1e-10)
})

test_that("pairwise cluster identity matches Gower equation (2)", {
  views <- make_views(n = 8, seed = 4)
  fit <- gpa(views, transform = proc_orthogonal("O"))
  Y <- fitted(fit)
  nms <- names(Y)
  pair <- 0
  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      if (i >= j) next
      pair <- pair + sum((Y[[i]] - Y[[j]])^2)
    }
  }
  k <- length(nms)
  expect_equal(pair, k * fit$objective, tolerance = 1e-10)
})

test_that("row masks keep unobserved consensus rows from that view", {
  set.seed(3)
  n <- 10
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  views <- list(A = X, B = X %*% rot2(0.5), C = X %*% rot2(-0.3))
  mask_c <- rep(TRUE, n)
  mask_c[c(1, 2, 9)] <- FALSE
  dat <- proc_data(
    views,
    observed = list(A = rep(TRUE, n), B = rep(TRUE, n), C = mask_c)
  )
  fit <- gpa(dat, transform = proc_orthogonal("SO"), gauge = proc_gauge(scale = "none"))
  expect_equal(fit$solver, "gower_bcd")
  expect_lt(fit$objective, 1e-6)
  Yc <- fitted(fit)[["C"]]
  expect_true(all(is.na(Yc[c(1, 2, 9), 1])))
  expect_true(all(is.finite(Yc[mask_c, 1])))
})

test_that("partial correspondence via ids still aligns shared landmarks", {
  set.seed(5)
  X <- scale(matrix(rnorm(12), 6, 2), scale = FALSE)
  dat <- proc_data(
    list(
      A = X[1:5, ],
      B = X[2:6, ] %*% rot2(0.6),
      C = X[c(1, 3, 5, 6), ] %*% rot2(-0.25)
    ),
    ids = list(
      A = c("a", "b", "c", "d", "e"),
      B = c("b", "c", "d", "e", "f"),
      C = c("a", "c", "e", "f")
    )
  )
  fit <- gpa(dat, transform = proc_orthogonal("O"), gauge = proc_gauge(scale = "none"))
  expect_equal(fit$data_summary$overlap$n_components, 1L)
  expect_lt(fit$objective, 1e-6)
  expect_equal(nrow(consensus(fit)), 6L)
})

test_that("isotropic + gauge none records the Gower constraint and does not collapse", {
  views <- make_views()
  fit <- gpa(
    views,
    transform = proc_similarity(group = "O", translation = TRUE, scaling = "isotropic"),
    gauge = proc_gauge(scale = "none")
  )
  expect_equal(fit$gauge$scale_mode, "gower")
  expect_true(length(fit$warnings) >= 1)
  expect_gt(sqrt(sum(consensus(fit)^2)), 0.1)
  expect_true(all(vapply(transformations(fit), function(tr) tr$s, numeric(1)) > 0))
})

test_that("preshape GPA rotates only after unit-Frobenius centering", {
  views <- make_views()
  fit <- gpa(
    views,
    transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
    gauge = proc_gauge(scale = "preshape")
  )
  expect_equal(fit$gauge$scale_mode, "preshape")
  Y <- fitted(fit)
  for (nm in names(Y)) {
    expect_equal(sqrt(sum(Y[[nm]]^2)), 1, tolerance = 1e-8)
  }
})

test_that("fixed-consensus scale holds ||M||_F = 1", {
  views <- make_views()
  fit <- gpa(
    views,
    transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
    gauge = proc_gauge(scale = "fixed_consensus")
  )
  expect_equal(fit$gauge$scale_mode, "fixed_consensus")
  expect_equal(sqrt(sum(consensus(fit)^2)), 1, tolerance = 1e-8)
})

test_that("aligned configurations share a centroid when translation is free", {
  views <- make_views()
  fit <- gpa(
    views,
    transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
    gauge = proc_gauge(scale = "gower")
  )
  Y <- fitted(fit)
  cents <- lapply(Y, function(Yi) colMeans(Yi))
  expect_equal(cents[[1]], cents[[2]], tolerance = 1e-8)
  expect_equal(cents[[1]], cents[[3]], tolerance = 1e-8)
  expect_equal(colMeans(consensus(fit)), cents[[1]], tolerance = 1e-8)
})

test_that("K=2 with an explicit anchor stays on the pairwise kernel", {
  views <- make_views()[1:2]
  fit <- gpa(
    views,
    transform = proc_orthogonal("SO"),
    anchor = "A"
  )
  expect_equal(fit$solver, "pairwise_polar")
  expect_equal(fit$optimality_status, "exact_closed_form")
  expect_equal(transformations(fit)[["A"]]$R, diag(2), tolerance = 1e-12)
})

test_that("K=2 without anchor is consensus-first Gower BCD", {
  views <- make_views()[1:2]
  fit <- gpa(views, transform = proc_orthogonal("SO"))
  expect_equal(fit$solver, "gpm")
  expect_true(fit$optimality_status %in% c("first_order_stationary", "certified_global"))
})

test_that("inits sequential, medoid, and spectral all converge on complete data", {
  set.seed(8)
  X <- scale(matrix(rnorm(20), 10, 2), scale = FALSE)
  views <- list(A = X, B = X %*% rot2(0.8), C = X %*% rot2(-0.5))
  objs <- vapply(c("sequential", "medoid", "spectral"), function(init) {
    fit <- gpa(
      views,
      transform = proc_orthogonal("O"),
      control = gpa_control(init = init, accelerate = FALSE, tolerance = 1e-10)
    )
    expect_equal(fit$numerical_status, "converged")
    fit$objective
  }, numeric(1))
  expect_true(all(is.finite(objs)))
  expect_true(all(objs < 1e-8))
})

test_that("safeguarded acceleration never increases the recorded objective", {
  views <- make_views()
  fit <- gpa(
    views,
    transform = proc_similarity(),
    gauge = proc_gauge(scale = "gower"),
    control = gpa_control(accelerate = TRUE, tolerance = 1e-10)
  )
  expect_true(all(diff(fit$history$objective) <= 1e-8))
})

test_that("configuration weights change the consensus", {
  views <- make_views()
  eq <- gpa(views, transform = proc_orthogonal("O"))
  w <- gpa(
    views,
    transform = proc_orthogonal("O"),
    metric = proc_metric(configuration = c(A = 10, B = 0.1, C = 0.1))
  )
  expect_gt(sqrt(sum((consensus(eq) - consensus(w))^2)), 1e-6)
})

test_that("canonicalize is a presentation and leaves the objective unchanged", {
  views <- make_views()
  fit <- gpa(views, transform = proc_orthogonal("O"))
  can <- canonicalize(fit)
  expect_equal(can$objective, fit$objective, tolerance = 1e-12)
  M <- consensus(can)
  expect_lt(abs(M[1, 2] * M[2, 1]), 1 + sqrt(sum(M^2)))
})

test_that("explain_solver selects gpm for multi-view orthogonal L2", {
  compiled <- compile_proc_problem(
    proc_data(make_views()),
    transform = proc_orthogonal("SO")
  )
  plan <- explain_solver(compiled)
  expect_equal(plan$solver, "gpm")
  expect_match(plan$allowed_claim, "first_order_stationary")
})

test_that("explain_solver selects gower_bcd for multi-view similarity", {
  compiled <- compile_proc_problem(
    proc_data(make_views()),
    transform = proc_similarity()
  )
  plan <- explain_solver(compiled)
  expect_equal(plan$solver, "gower_bcd")
  expect_equal(plan$allowed_claim, "blockwise_stationary")
})

test_that("nstart records a finite best objective and spectral falls back on masks", {
  views <- make_views()
  multi <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(nstart = 3L, accelerate = FALSE)
  )
  expect_true(is.finite(multi$objective))
  n <- 8
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  mask <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
  dat <- proc_data(
    list(A = X, B = X %*% rot2(0.3), C = X %*% rot2(-0.2)),
    observed = list(A = rep(TRUE, n), B = mask, C = rep(TRUE, n))
  )
  spec_fit <- gpa(
    dat,
    transform = proc_orthogonal("O"),
    control = gpa_control(init = "spectral", accelerate = FALSE)
  )
  expect_equal(spec_fit$init, "spectral")
  expect_equal(spec_fit$numerical_status, "converged")
})

test_that("certificate remains unavailable for Gower BCD", {
  fit <- gpa(make_views(), transform = proc_similarity(), gauge = proc_gauge(scale = "gower"))
  cert <- certify(fit)
  expect_equal(cert$status, "unavailable")
  expect_false(identical(fit$optimality_status, "certified_global"))
})

test_that("gower-identity fixture satisfies the pairwise-cluster law after fitting", {
  skip_if_not_installed("jsonlite")
  fx <- load_fixture("gower-identity-pairwise-clusters")
  mats <- lapply(fx$problem$views$matrix, function(m) {
    if (is.list(m)) do.call(rbind, lapply(m, as.numeric)) else as.matrix(m)
  })
  if (is.null(names(mats)) || !length(names(mats))) {
    nms <- fx$problem$views$name
    if (is.null(nms)) nms <- paste0("X", seq_along(mats))
    names(mats) <- nms
  }
  # jsonlite may parse views as a data.frame
  if (is.data.frame(fx$problem$views)) {
    nms <- fx$problem$views$name
    mats <- lapply(seq_along(nms), function(i) {
      m <- fx$problem$views$matrix[[i]]
      if (is.list(m)) do.call(rbind, m) else matrix(m, ncol = 2, byrow = TRUE)
    })
    names(mats) <- nms
  }
  fit <- gpa(mats, transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
             gauge = proc_gauge(scale = "gower"))
  Y <- fitted(fit)
  pair <- 0
  nms <- names(Y)
  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      if (i >= j) next
      pair <- pair + sum((Y[[i]] - Y[[j]])^2)
    }
  }
  expect_equal(pair, length(nms) * fit$objective, tolerance = 1e-8)
  cents <- lapply(Y, colMeans)
  expect_equal(cents[[1]], cents[[2]], tolerance = 1e-8)
})

test_that("law-gower-energy fixture matches the algebraic identity", {
  skip_if_not_installed("jsonlite")
  fx <- load_fixture("law-gower-energy")
  views <- fx$problem$views
  if (is.data.frame(views)) {
    mats <- lapply(seq_len(nrow(views)), function(i) {
      m <- views$matrix[[i]]
      if (is.list(m)) do.call(rbind, lapply(m, as.numeric)) else as.matrix(m)
    })
    names(mats) <- views$name
  } else {
    mats <- lapply(views, function(v) {
      m <- v$matrix
      if (is.list(m)) do.call(rbind, lapply(m, as.numeric)) else as.matrix(m)
    })
    names(mats) <- vapply(views, `[[`, character(1), "name")
  }
  w <- fx$problem$configuration_weights
  M <- Reduce(`+`, mats) / length(mats)
  dec <- decompose_energy(mats, M, w)
  expect_equal(dec$total_energy, fx$expected$a, tolerance = 1e-12)
  expect_equal(dec$consensus_energy, fx$expected$b, tolerance = 1e-12)
  expect_equal(dec$residual_energy, fx$expected$objective, tolerance = 1e-12)
  expect_equal(dec$total_energy, dec$consensus_energy + dec$residual_energy, tolerance = 1e-12)
})

test_that("gower-1975-published-history transcribes Tables 2 and 5", {
  skip_if_not_installed("jsonlite")
  fx <- load_fixture("gower-1975-published-history")
  expect_equal(fx$citation, "Gower, J. C. (1975). Psychometrika 40:33-51, Tables 2 and 5.")
  nms <- fx$problem$views$name
  expect_equal(nms, c("judge_1", "judge_2", "judge_3"))
  mats <- lapply(seq_along(nms), function(i) {
    m <- fx$problem$views$matrix[[i]]
    if (is.list(m)) do.call(rbind, m) else as.matrix(m)
  })
  names(mats) <- nms
  expect_equal(dim(mats$judge_1), c(9L, 7L))
  expect_equal(mats$judge_1[1, ], c(47, 44, 49, 38, 35, 40, 40))
  expect_equal(mats$judge_2[8, 1], 5)
  expect_equal(mats$judge_3[8, ], c(5, 95, 95, 3, 20, 2, 24))
  sr <- fx$expected$published_sr
  expect_equal(sr$initial, 0.661438)
  expect_equal(utils::tail(sr$with_scaling$after_scaling, 1), 0.598842)
  expect_equal(utils::tail(sr$without_scaling$after_rotation, 1), 0.657137)
  fit <- gpa(
    mats,
    transform = proc_similarity(group = "SO", translation = TRUE, scaling = "isotropic"),
    gauge = proc_gauge(scale = "gower")
  )
  expect_true(all(diff(fit$history$objective) <= 1e-8))
  expect_equal(certify(fit)$status, "unavailable")
  expect_gt(sqrt(sum(consensus(fit)^2)), 0)
})

test_that("gower-1975-algorithm-laws remain documented on a live fit", {
  skip_if_not_installed("jsonlite")
  fx <- load_fixture("gower-1975-algorithm-laws")
  expect_true(all(c(
    "the residual is monotone nonincreasing",
    "a global energy constraint prevents the collapse s_i=0, M=0",
    "monotone decrease does not prove a global optimum"
  ) %in% fx$expected$laws))
  fit <- gpa(make_views(), transform = "similarity", gauge = proc_gauge(scale = "gower"))
  expect_true(all(diff(fit$history$objective) <= 1e-8))
  expect_gt(sqrt(sum(consensus(fit)^2)), 0)
  expect_equal(certify(fit)$status, "unavailable")
})
