classic_tps_kernel <- function(r, d) {
  r <- pmax(as.numeric(r), 0)
  if (d == 2L) {
    out <- r^2 * log(r)
    out[!is.finite(out) | r < 1e-15] <- 0
    return(out)
  }
  if (d == 3L) return(-r)
  stop("oracle supports only 2D and 3D")
}

classic_tps_oracle <- function(X, Y, Xnew) {
  d <- ncol(X)
  controls <- X
  P <- cbind(1, X)
  K <- vapply(seq_len(nrow(controls)), function(k) {
    delta <- sweep(controls, 2L, controls[k, ], "-")
    classic_tps_kernel(sqrt(rowSums(delta * delta)), d)
  }, numeric(nrow(controls)))
  K <- (K + t(K)) / 2
  saddle <- rbind(
    cbind(K, P),
    cbind(t(P), matrix(0, ncol(P), ncol(P)))
  )
  coef <- solve(saddle, rbind(Y, matrix(0, ncol(P), ncol(Y))))
  Unew <- vapply(seq_len(nrow(controls)), function(k) {
    delta <- sweep(Xnew, 2L, controls[k, ], "-")
    classic_tps_kernel(sqrt(rowSums(delta * delta)), d)
  }, numeric(nrow(Xnew)))
  Unew %*% coef[seq_len(nrow(controls)), , drop = FALSE] +
    cbind(1, Xnew) %*% coef[(nrow(controls) + 1L):nrow(coef), , drop = FALSE]
}

test_that("constrained TPS basis agrees with the augmented saddle oracle in 2D and 3D", {
  for (d in 2:3) {
    X <- if (d == 2L) {
      rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(0.25, 0.6), c(0.8, 0.3))
    } else {
      rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1),
            c(1, 1, 0), c(1, 0, 1), c(0, 1, 1), c(0.3, 0.4, 0.7))
    }
    Y <- cbind(X[, 1L]^2 + 0.2 * X[, d], X[, 1L] - 0.4 * X[, 2L])
    Xnew <- X + 0.07
    built <- gprocrustes:::.lbw_basis(proc_tps(control_points = X, smoothness = 0), X)
    B <- gprocrustes:::.gproc_psd_solve(crossprod(built$Phi), crossprod(built$Phi, Y))
    fitted <- built$Phi %*% B
    tr <- gprocrustes:::.lbw_fitted(
      proc_tps(control_points = X, smoothness = 0), B,
      list(controls = X, tps_qc = built$tps_qc), d
    )
    off_training <- apply_proc_transform(tr, Xnew)
    oracle <- classic_tps_oracle(X, Y, Xnew)
    expect_equal(fitted, Y, tolerance = 1e-8)
    expect_equal(off_training, oracle, tolerance = 1e-7)
    expect_gte(sum(B * (built$L %*% B)), -1e-10)
  }
})

test_that("TPS bending penalty is PSD and vanishes only for affine represented functions", {
  for (d in 2:3) {
    X <- if (d == 2L) {
      rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(0.2, 0.7), c(0.8, 0.3))
    } else {
      rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1),
            c(1, 1, 0), c(1, 0, 1), c(0, 1, 1), c(0.3, 0.4, 0.7))
    }
    built <- gprocrustes:::.lbw_basis(proc_tps(control_points = X), X)
    vals <- gprocrustes:::.gproc_eigen_sym(built$L)$values
    expect_gte(min(vals), -1e-10)
    expect_equal(sum(abs(vals) < 1e-10), ncol(X) + 1L)

    affine <- cbind(1, X) %*% seq_len(ncol(X) + 1L)
    nonlinear <- gprocrustes:::.gproc_tps_phi(X, X)[, ncol(cbind(1, X)) + 1L]
    Ba <- gprocrustes:::.gproc_psd_solve(crossprod(built$Phi), crossprod(built$Phi, affine))
    Bn <- gprocrustes:::.gproc_psd_solve(crossprod(built$Phi), crossprod(built$Phi, nonlinear))
    expect_lt(sum(Ba * (built$L %*% Ba)), 1e-10)
    expect_gt(sum(Bn * (built$L %*% Bn)), 1e-8)
  }
})

test_that("LBW rejects invalid penalties, invalid smoothness, and anchored eigenproblems", {
  X <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(0.3, 0.7))
  bad_penalties <- list(
    diag(c(0, -1, 0)),
    matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0), 3, 3),
    diag(c(0, NA_real_, 0))
  )
  for (L in bad_penalties) {
    bad <- proc_lbw(function(z) cbind(1, z), penalty = L, smoothness = 1)
    expect_error(gpa(list(A = X, B = X + 0.1), transform = bad), class = "invalid_problem")
  }
  expect_error(gpa(list(A = X, B = X + 0.1), transform = proc_tps(smoothness = -1)),
               class = "invalid_problem")
  for (spec in list(proc_affine(), proc_lbw("affine"), proc_tps())) {
    for (solver in c("auto", "lbw_eigen")) {
      expect_error(gpa(list(A = X, B = X + 0.1), transform = spec,
                       anchor = "A", solver = solver), class = "invalid_problem")
    }
  }
})

test_that("high TPS smoothness approaches affine for non-affinely related data", {
  X <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(0.15, 0.65), c(0.8, 0.35))
  bend <- cbind(0.25 * X[, 1L]^2, -0.3 * X[, 2L]^2)
  views <- list(A = X, B = X + bend)
  affine <- gpa(views, transform = proc_affine())
  tps <- gpa(views, transform = proc_tps(smoothness = 1e8))
  gap <- procrustes(consensus(tps), consensus(affine), transform = proc_orthogonal("O"))$objective
  expect_lt(gap, 1e-4)
})

test_that("TPS prediction rebuilds a constrained basis for new data and gauge alignment retains fitted controls", {
  X <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(0.2, 0.7), c(0.8, 0.3))
  views <- list(A = X, B = X + cbind(0.15 * X[, 1L]^2, -0.1 * X[, 2L]^2))
  fit <- gpa(views, transform = proc_tps(smoothness = 0))
  newdata <- X + 0.04
  built <- gprocrustes:::.lbw_basis(fit$problem$transform, newdata)
  B <- gprocrustes:::.gproc_psd_solve(
    crossprod(built$Phi),
    crossprod(built$Phi, consensus(fit))
  )
  expect_equal(predict(fit, newdata), built$Phi %*% B, tolerance = 1e-10)
  tr_new <- gprocrustes:::.lbw_fitted(
    fit$problem$transform, B,
    list(controls = newdata, tps_qc = built$tps_qc), 2L
  )
  expect_equal(
    apply_proc_transform(tr_new, newdata + 0.03),
    classic_tps_oracle(newdata, consensus(fit), newdata + 0.03),
    tolerance = 1e-7
  )

  tr <- transformations(fit)[["A"]]
  Q <- matrix(c(0, 1, -1, 0), 2, 2)
  moved <- align_gauge(fit, Q)
  expect_equal(
    apply_proc_transform(transformations(moved)[["A"]], newdata),
    apply_proc_transform(tr, newdata) %*% Q,
    tolerance = 1e-10
  )
})
