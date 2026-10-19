review_infer_data <- function(named = FALSE, noise = 0) {
  set.seed(71)
  X <- scale(matrix(rnorm(24), 12, 2), scale = FALSE)
  turn <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  v <- list(A = X, B = X %*% turn(.5), C = X %*% turn(-.7))
  v <- lapply(v, function(Y) Y + matrix(rnorm(24, sd = noise), 12, 2))
  if (named) v <- lapply(v, function(Y) { rownames(Y) <- paste0("p", 1:12); Y })
  v
}

test_that("entity bootstrap evaluates every replicate on the original landmark index", {
  a <- gpa(review_infer_data(), transform = "O")
  b <- gpa(review_infer_data(TRUE), transform = "O")
  set.seed(5)
  ba <- infer(a, method = "bootstrap", n_resamples = 8)
  set.seed(5)
  bb <- infer(b, method = "bootstrap", n_resamples = 8)
  expect_lt(max(ba$se), 1e-8)
  expect_equal(bb$replicates, ba$replicates, tolerance = 1e-8)
  for (i in 1:8) {
    expect_equal(unname(ba$replicates[, , i]), unname(consensus(a)), tolerance = 1e-8)
  }
  expect_match(paste(ba$assumptions, collapse = " "), "original|fixed")
})

test_that("entity bootstrap equals an independently constructed multiplicity-weight refit", {
  v <- review_infer_data(TRUE, noise = .1)
  fit <- gpa(v, transform = "O")
  set.seed(45)
  draws <- replicate(3, sample.int(12, 12, replace = TRUE))
  set.seed(45)
  boot <- infer(fit, method = "bootstrap", n_resamples = 3)
  for (b in 1:3) {
    # Duplicate rows have exactly the same LS contribution as integer frequency weights.
    w <- tabulate(draws[, b], 12)
    ref <- gpa(v, transform = "O", metric = proc_metric(landmark = w),
               solver = "gower_bcd", control = gpa_control(tolerance = 1e-10,
                 max_iterations = 500, certify = "never", accelerate = FALSE))
    Y <- lapply(names(v), function(nm) apply_proc_transform(transformations(ref)[[nm]], v[[nm]]))
    M <- Reduce(`+`, Y) / length(Y)
    # Independent polar alignment, including the common translation nuisance.
    Mc <- sweep(M, 2, colMeans(M), "-")
    target <- consensus(fit)
    Tc <- sweep(target, 2, colMeans(target), "-")
    sv <- svd(crossprod(Mc, Tc))
    expected <- sweep(Mc %*% sv$u %*% t(sv$v), 2, colMeans(target), "+")
    expect_equal(unname(boot$replicates[, , b]), unname(expected), tolerance = 2e-6)
  }
})

test_that("permutations break named correspondence and use the lower residual tail", {
  a <- gpa(review_infer_data(), transform = "O")
  b <- gpa(review_infer_data(TRUE), transform = "O")
  set.seed(5)
  pa <- infer(a, method = "permutation", n_resamples = 12)
  set.seed(5)
  pb <- infer(b, method = "permutation", n_resamples = 12)
  expect_gt(min(pa$null_statistics), 1)
  expect_equal(pb$null_statistics, pa$null_statistics, tolerance = 1e-7)
  expect_equal(pa$p_value, 1 / 13)
  expect_equal(pb$p_value, 1 / 13)
  expect_equal(pa$p_value, (1 + sum(pa$null_statistics <= pa$statistic)) / 13)
})

test_that("resampling rejects invalid draw counts and preserves weighted model channels", {
  v <- review_infer_data(TRUE, noise = .05)
  f <- gpa(v, transform = "O", metric = proc_metric(landmark = seq_len(12)))
  expect_error(infer(f, method = "bootstrap", n_resamples = 1), class = "invalid_problem")
  expect_error(infer(f, method = "permutation", n_resamples = 0), class = "invalid_problem")
  expect_error(infer(f, method = "permutation", n_resamples = 1e20), class = "invalid_problem")
  anchored <- gpa(v, transform="O", anchor="A")
  expect_error(infer(anchored,method="bootstrap",n_resamples=3),class="invalid_problem")
  set.seed(7)
  expect_s3_class(infer(f, method = "bootstrap", n_resamples = 3), "proc_inference")
  set.seed(7)
  p <- infer(f, method = "permutation", n_resamples = 3)
  expect_equal(p$statistic, f$objective)
})

test_that("entity bootstrap respects local landmark weights after row reordering", {
  v <- review_infer_data(TRUE, noise = .08)
  ids <- lapply(v, rownames)
  global_w <- c(2, 5, 1, 4, 3, 6, 2, 7, 1, 5, 4, 3)
  local_w <- lapply(ids, function(id) global_w[match(id, ids[[1L]])])
  metric <- proc_metric(landmark = stats::setNames(local_w, names(v)))
  named_fit <- gpa(v, transform = "O", metric = metric)

  order_b <- c(4L, 1L, 10L, 2L, 12L, 7L, 5L, 11L, 3L, 9L, 6L, 8L)
  order_c <- c(8L, 3L, 11L, 5L, 1L, 9L, 4L, 12L, 6L, 2L, 10L, 7L)
  take <- list(A = seq_len(12L), B = order_b, C = order_c)
  shuffled <- Map(function(X, ii) {
    out <- X[ii, , drop = FALSE]
    rownames(out) <- NULL
    out
  }, v, take)
  ids_shuffled <- Map(function(id, ii) id[ii], ids, take)
  weights_shuffled <- Map(function(id) global_w[match(id, ids[[1L]])], ids_shuffled)
  explicit_fit <- gpa(
    proc_data(shuffled, ids = ids_shuffled),
    transform = "O",
    metric = proc_metric(landmark = stats::setNames(weights_shuffled, names(shuffled)))
  )

  set.seed(19)
  named_boot <- infer(named_fit, method = "bootstrap", n_resamples = 3)
  set.seed(19)
  explicit_boot <- infer(explicit_fit, method = "bootstrap", n_resamples = 3)
  expect_equal(consensus(explicit_fit), consensus(named_fit), tolerance = 1e-8)
  expect_equal(explicit_boot$replicates, named_boot$replicates, tolerance = 2e-6)
})

test_that("configuration bootstrap refits duplicate-only configuration draws", {
  v <- review_infer_data(TRUE, noise = .05)
  fit <- gpa(v, transform = "O")
  data <- gprocrustes:::.gproc_data_from_fit(fit)
  names_b <- c("boot1", "boot2", "boot3")
  duplicate <- gprocrustes:::.gproc_resampled_configurations(data, rep(1L, 3L), names_b)
  metric <- gprocrustes:::.gproc_resampled_configuration_metric(fit, rep(1L, 3L), names_b)
  refit <- gprocrustes:::.gproc_refit(fit, duplicate, metric)
  expect_true(all(is.finite(consensus(refit))))
  expect_equal(names(transformations(refit)), names_b)
})

test_that("bootstrap rejects unsupported transformation families explicitly", {
  X <- matrix(c(-1, 0, 2, 1, 0, 3, 4, -2, -3, 1), 5, 2, byrow = TRUE)
  specs <- list(proc_signed_permutation(), proc_affine(), proc_lbw("affine"), proc_tps())
  for (spec in specs) {
    fit <- gpa(list(A = X, B = X), transform = spec)
    expect_error(infer(fit, method = "bootstrap", n_resamples = 2), class = "invalid_problem")
  }
})
