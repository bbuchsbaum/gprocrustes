rot2 <- function(theta) {
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
}

campaign_views <- function(n = 12, seed = 21) {
  set.seed(seed)
  X <- scale(matrix(rnorm(n * 2), n, 2), scale = FALSE)
  list(
    A = X,
    B = X %*% rot2(0.55),
    C = X %*% rot2(-0.4)
  )
}

test_that("every fixture kind required by the conformance spec is present", {
  skip_if_not_installed("jsonlite")
  kinds <- c(
    "pairwise", "moments", "law", "error",
    "gower_history", "identifiability", "certificate"
  )
  counts <- vapply(kinds, function(k) length(list_fixtures(k)), integer(1))
  expect_true(all(counts > 0), info = paste(names(counts), counts, collapse = ", "))
  expect_true("gower-1975-published-history" %in% list_fixtures("gower_history"))
  expect_true("moments-dense-sparse-agree" %in% list_fixtures("moments"))
  expect_true("law-gauge-objective" %in% list_fixtures("law"))
  expect_true("cellmask-not-svd" %in% list_fixtures("law"))
  expect_true("id-disconnected-overlap" %in% list_fixtures("identifiability"))
  expect_true("cert-od-unavailable-default" %in% list_fixtures("certificate"))
})

test_that("identifiability fixtures report the documented overlap graph", {
  skip_if_not_installed("jsonlite")
  for (id in list_fixtures("identifiability")) {
    fx <- load_fixture(id)
    exp <- fx$expected
    if (identical(exp$numerical_status, "invalid_problem") ||
        identical(exp$numerical_status, "error")) {
      dat <- fixture_proc_data(fx)
      if (identical(id, "id-disconnected-overlap")) {
        expect_equal(overlap_graph(dat)$n_components, as.integer(exp$connected_components %||% 2L))
        expect_error(gpa(dat, transform = "O"), class = exp$error_code %||% "disconnected_overlap_graph")
      } else {
        expect_error(
          if (is.null(dat)) stop("no views") else gpa(dat, transform = "O"),
          class = exp$error_code %||% "error",
          info = id
        )
      }
      next
    }
    dat <- fixture_proc_data(fx)
    ov <- overlap_graph(dat)
    if (!is.null(exp$connected_components)) {
      expect_equal(ov$n_components, as.integer(exp$connected_components), info = id)
    }
    if (identical(id, "id-entity-observed-nowhere")) {
      expect_false("z" %in% dat$global_ids)
      expect_equal(nrow(gpa(dat, transform = proc_orthogonal("O"))$consensus), 2L)
    }
    if (identical(id, "id-connected-partial")) {
      expect_equal(ov$n_components, 1L)
      fit <- gpa(dat, transform = proc_orthogonal("O"), gauge = proc_gauge(scale = "none"))
      expect_equal(nrow(consensus(fit)), 6L)
    }
  }
})

test_that("error fixtures raise the documented codes from JSON", {
  skip_if_not_installed("jsonlite")
  for (id in list_fixtures("error")) {
    fx <- load_fixture(id)
    code <- fx$expected$error_code %||% "error"
    if (identical(id, "err-nan")) {
      expect_error(procrustes(matrix(c(1, NA, 0, 1), 2, 2), diag(2)), class = code)
    } else if (identical(id, "err-inf")) {
      expect_error(procrustes(matrix(c(1, Inf, 0, 1), 2, 2), diag(2)), class = code)
    } else if (identical(id, "err-dimension-mismatch")) {
      expect_error(
        procrustes(matrix(rnorm(10), 5, 2), matrix(rnorm(15), 5, 3)),
        class = code
      )
    } else if (identical(id, "err-duplicate-ids")) {
      expect_error(fixture_proc_data(fx), class = code)
    } else if (identical(id, "err-zero-config-weight-all")) {
      expect_error(
        gpa(
          list(matrix(rnorm(8), 4, 2), matrix(rnorm(8), 4, 2)),
          transform = "O",
          metric = proc_metric(configuration = c(0, 0))
        ),
        class = code
      )
    }
  }
})

test_that("gauge and missingness law fixtures execute their identities", {
  skip_if_not_installed("jsonlite")
  fx <- load_fixture("law-gauge-objective")
  vs <- fixture_view_list(fx)
  mats <- stats::setNames(lapply(vs, `[[`, "matrix"), vapply(vs, `[[`, character(1), "name"))
  a <- procrustes(mats$X, mats$Y, transform = proc_orthogonal("O"))
  b <- procrustes(mats$XQ, mats$YQ, transform = proc_orthogonal("O"))
  expect_equal(a$objective, b$objective, tolerance = fx$expected$objective_atol %||% 1e-12)

  so <- load_fixture("law-so-det-plus-one")
  expect_true("det(R)=+1 for every proper-rotation fit" %in% so$expected$laws)
  X <- matrix(c(1, 0, 0, 1, 0.3, -0.2), 3, 2, byrow = TRUE)
  Href <- diag(c(1, -1))
  pair <- procrustes(X, X %*% Href, transform = proc_orthogonal("SO"))
  expect_equal(det(pair$transform$R), 1, tolerance = 1e-12)

  ord <- load_fixture("law-entity-order")
  X <- read_matrix(ord$problem$source)
  Y <- read_matrix(ord$problem$target)
  perm <- c(3L, 1L, 4L, 2L, 8L, 5L, 7L, 6L)
  perm <- perm[perm <= nrow(X)]
  p1 <- procrustes(X, Y, transform = proc_orthogonal("O"))
  p2 <- procrustes(X[perm, , drop = FALSE], Y[perm, , drop = FALSE], transform = proc_orthogonal("O"))
  expect_equal(p1$objective, p2$objective, tolerance = 1e-12)

  rowfx <- load_fixture("rowmask-consensus-mean")
  dat <- fixture_proc_data(rowfx)
  fit <- gpa(dat, transform = proc_orthogonal("O"), gauge = proc_gauge(scale = "none"))
  M <- consensus(fit)
  Y <- fitted(fit)
  expect_equal(M[1, ], colMeans(rbind(Y$A[1, ], Y$B[1, ])), tolerance = 1e-8)
  expect_equal(M[2, ], Y$A[2, ], tolerance = 1e-8)
  expect_equal(M[3, ], Y$B[3, ], tolerance = 1e-8)
  expect_true(all(is.na(Y$A[3, ])))
  expect_true(all(is.na(Y$B[2, ])))

  cellfx <- load_fixture("cellmask-not-svd")
  expect_true(all(c(
    "missing rows != missing arbitrary cells",
    "ordinary polar factor is not an exact block update"
  ) %in% cellfx$expected$laws))
  X <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  cells <- list(
    A = matrix(c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE), 4, 2, byrow = TRUE),
    B = matrix(TRUE, 4, 2)
  )
  compiled <- compile_proc_problem(
    proc_data(list(A = X, B = X %*% rot2(0.2)), cells = cells),
    transform = proc_orthogonal("O")
  )
  plan <- explain_solver(compiled)
  expect_equal(plan$solver, "mm")
  expect_equal(plan$allowed_claim, "first_order_stationary")
})

test_that("certificate fixtures keep Ling off SO, robust, and cell-mask paths", {
  skip_if_not_installed("jsonlite")
  shape <- load_fixture("cert-gpm-fixed-point-shape")
  expect_true(all(c(
    "report r_dual = ||CS-Lambda S|| / (1+||CS||)",
    "report a lower bound on lambda_min(Lambda-C)"
  ) %in% shape$expected$laws))
  so <- load_fixture("cert-od-unavailable-default")
  expect_true("do not reuse for SO(d), robust losses, or cell masks" %in% so$expected$laws)

  views <- campaign_views()
  gpm <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(init = "spectral", certify = "auto", tolerance = 1e-12)
  )
  cert <- certify(gpm)
  expect_true(cert$status %in% c("certified_global", "not_certified"))
  if (identical(cert$status, "certified_global")) {
    expect_true(is.finite(cert$r_dual))
    expect_true(is.finite(cert$lambda_min))
    expect_gte(cert$lambda_min, -1e-6)
  }
  expect_equal(certify(gpa(views, transform = proc_orthogonal("SO")))$status, "unavailable")
  expect_equal(
    certify(gpa(views, transform = proc_orthogonal("O"), loss = proc_huber()))$status,
    "unavailable"
  )
  cells <- lapply(views, function(v) matrix(TRUE, nrow(v), ncol(v)))
  cells$C[2, 1] <- FALSE
  mm <- gpa(
    proc_data(views, cells = cells),
    transform = proc_orthogonal("O"),
    control = gpa_control(accelerate = FALSE, max_iterations = 25)
  )
  expect_equal(mm$solver, "mm")
  expect_equal(certify(mm)$status, "unavailable")
  expect_false(identical(mm$optimality_status, "certified_global"))
})

test_that("lazy and materialized aligned stores agree after gauge presentation", {
  views <- campaign_views()
  lazy <- gpa(views, transform = proc_orthogonal("O"), control = gpa_control(keep_aligned = "lazy"))
  mat <- gpa(views, transform = proc_orthogonal("O"), control = gpa_control(keep_aligned = "materialize"))
  expect_equal(lazy$keep_aligned, "lazy")
  expect_equal(mat$keep_aligned, "materialize")
  expect_true(inherits(aligned(lazy)[["A"]], "proc_aligned_view"))
  expect_true(is.matrix(aligned(mat)[["A"]]))
  expect_equal(as.matrix(aligned(lazy)[["A"]]), aligned(mat)[["A"]], tolerance = 1e-10)
  expect_equal(lazy$objective, mat$objective, tolerance = 1e-12)
  expect_equal(fitted(lazy)[["B"]], fitted(mat)[["B"]], tolerance = 1e-10)
  expect_equal(residuals(lazy)[["C"]], residuals(mat)[["C"]], tolerance = 1e-10)
  expect_true(!is.null(lazy$raw_data))
  expect_true(!is.null(mat$raw_data))
  can_lazy <- canonicalize(lazy)
  can_mat <- canonicalize(mat)
  expect_equal(can_lazy$objective, lazy$objective, tolerance = 1e-12)
  expect_equal(
    as.matrix(aligned(can_lazy)[["A"]]),
    apply_proc_transform(transformations(can_lazy)[["A"]], views$A),
    tolerance = 1e-10
  )
  expect_equal(
    aligned(can_mat)[["A"]],
    apply_proc_transform(transformations(can_mat)[["A"]], views$A),
    tolerance = 1e-10
  )
})

test_that("dense and sparse moments agree for raw and centered pairwise kernels", {
  set.seed(8)
  X <- matrix(rnorm(24), 8, 3)
  Y <- X %*% matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 1), 3, 3) + 0.2
  w <- c(1, 0.5, 2, 1, 0, 1, 1.5, 1)
  Xs <- Matrix::Matrix(X, sparse = TRUE)
  Ys <- Matrix::Matrix(Y, sparse = TRUE)
  for (center in c(TRUE, FALSE)) {
    dense <- proc_moments(X, Y, w, center = center)
    sparse <- proc_moments(Xs, Ys, w, center = center)
    expect_equal(sparse$C, dense$C, tolerance = 1e-12)
    expect_equal(sparse$a, dense$a, tolerance = 1e-12)
    expect_equal(sparse$b, dense$b, tolerance = 1e-12)
    expect_equal(
      centered_crossprod(Xs, Ys, w),
      centered_crossprod(X, Y, w),
      tolerance = 1e-12
    )
  }
  pair_d <- procrustes(X, Y, transform = proc_orthogonal("O"), weights = w)
  pair_s <- procrustes(Xs, Ys, transform = proc_orthogonal("O"), weights = w)
  expect_equal(pair_s$objective, pair_d$objective, tolerance = 1e-10)
  expect_equal(pair_s$transform$R, pair_d$transform$R, tolerance = 1e-10)
})

test_that("iterative optimality claims match the recorded stationarity residual", {
  views <- campaign_views()
  check_status <- function(fit, allowed, cert = NULL) {
    expect_true(fit$optimality_status %in% allowed, info = fit$solver)
    expect_true(fit$numerical_status %in% c("converged", "stalled", "not_converged"))
    if (identical(fit$numerical_status, "converged") &&
        fit$optimality_status %in% c(
          "blockwise_stationary", "first_order_stationary", "certified_global"
        )) {
      stat <- if (!is.null(fit$history) && nrow(fit$history)) {
        utils::tail(fit$history$stationarity, 1)
      } else {
        fit$stationarity
      }
      expect_true(is.finite(stat), info = fit$solver)
      expect_lt(stat, 1e-4)
    }
    if (identical(fit$optimality_status, "certified_global")) {
      expect_equal(fit$numerical_status, "converged")
      expect_equal(certify(fit)$status, "certified_global")
    }
    if (!is.null(cert)) {
      expect_equal(certify(fit)$status, cert)
    }
    expect_false(
      identical(fit$numerical_status, "stalled") &&
        identical(fit$optimality_status, "first_order_stationary")
    )
    expect_false(
      identical(fit$numerical_status, "stalled") &&
        identical(fit$optimality_status, "certified_global")
    )
  }

  gpm <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(init = "spectral", certify = "auto", tolerance = 1e-12)
  )
  check_status(gpm, c("first_order_stationary", "certified_global"))

  so <- gpa(views, transform = proc_orthogonal("SO"))
  check_status(so, c("first_order_stationary", "not_converged"), cert = "unavailable")

  gower <- gpa(
    views,
    transform = proc_similarity(),
    gauge = proc_gauge(scale = "gower"),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-10)
  )
  check_status(gower, "blockwise_stationary", cert = "unavailable")

  irls <- gpa(
    views,
    transform = proc_orthogonal("O"),
    loss = proc_huber(k = 1.2),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-6, max_iterations = 40)
  )
  check_status(irls, "blockwise_stationary", cert = "unavailable")

  cells <- lapply(views, function(v) matrix(TRUE, nrow(v), ncol(v)))
  cells$C[2, 1] <- FALSE
  mm <- gpa(
    proc_data(views, cells = cells),
    transform = proc_orthogonal("O"),
    control = gpa_control(accelerate = FALSE, tolerance = 1e-8, max_iterations = 40)
  )
  expect_equal(mm$solver, "mm")
  check_status(mm, c("first_order_stationary", "not_converged"), cert = "unavailable")

  mf <- gpa(
    views,
    transform = proc_orthogonal("O"),
    control = gpa_control(backend = "matrix_free", init = "spectral", certify = "auto")
  )
  expect_false(identical(certify(mf)$status, "certified_global"))
})

test_that("frozen 1.0 exports and S3 methods are present and ICP/OT are not", {
  ns <- getNamespaceExports("gprocrustes")
  required <- c(
    "procrustes", "gpa", "proc_data", "proc_from_array",
    "compile_proc_problem", "explain_solver",
    "aligned", "block_apply",
    "proc_orthogonal", "proc_similarity", "proc_signed_permutation",
    "proc_affine", "proc_lbw", "proc_tps",
    "proc_metric", "proc_weights", "proc_covariance",
    "proc_squared_l2", "proc_huber", "proc_tukey",
    "proc_gauge", "gpa_control",
    "consensus", "transformations",
    "align_gauge", "canonicalize", "diagnose", "certify", "decompose",
    "quotient_distance", "principal_angles", "subspace_distance",
    "horizontal_component", "vertical_component", "tied_axis_blocks",
    "tangent_coordinates", "tangent_project",
    "apply_proc_transform", "inverse_proc_transform", "compose_proc_transform",
    "datum_space_error",
    "infer", "proc_shape_model", "cross_validate", "tune_smoothness",
    "tidy.gpa_fit", "glance.gpa_fit", "augment.gpa_fit",
    "plot_data", "autoplot.gpa_fit"
  )
  missing <- setdiff(required, ns)
  expect_equal(missing, character(0))
  forbidden <- c(
    "icp", "iterative_closest_point", "optimal_transport", "ot_align",
    "hungarian", "slide_semilandmarks", "proc_pca", "proc_cca", "proc_nmf"
  )
  expect_equal(intersect(tolower(ns), forbidden), character(0))

  s3 <- utils::methods(class = "gpa_fit")
  s3_names <- sub("\\.gpa_fit$", "", as.character(s3))
  expect_true(all(c(
    "predict", "fitted", "residuals", "coef", "aligned", "consensus",
    "diagnose", "certify", "decompose", "canonicalize", "align_gauge"
  ) %in% s3_names))
  expect_true(is.function(getExportedValue("gprocrustes", "tidy.gpa_fit")))
  expect_true(is.function(getExportedValue("gprocrustes", "glance.gpa_fit")))
  expect_true(is.function(getExportedValue("gprocrustes", "augment.gpa_fit")))
})
