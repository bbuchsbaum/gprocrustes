test_that("pairwise fixtures match the closed-form kernel", {
  ids <- list_fixtures("pairwise")
  skip_if(!length(ids), "no pairwise fixtures")
  for (id in ids) {
    fx <- load_fixture(id)
    exp <- fx$expected
    if (identical(exp$numerical_status, "error") ||
        identical(exp$numerical_status, "invalid_problem")) {
      expect_error(run_pairwise_fixture(fx), class = exp$error_code %||% "error", info = id)
      next
    }
    if (is.null(fx$problem$source) || is.null(fx$problem$target)) {
      next
    }
    pair <- run_pairwise_fixture(fx)
    expect_equal(pair$optimality_status, exp[["optimality_status"]] %||% "exact_closed_form")
    obj_exp <- exp[["objective"]]
    if (!is.null(obj_exp) && is.numeric(obj_exp)) {
      atol <- exp[["objective_atol"]] %||% 1e-10
      scale <- max(1, abs(pair$moments$a), abs(pair$moments$b), abs(obj_exp))
      expect_true(
        abs(pair$objective - obj_exp) <= max(atol, 1e-10 * scale, 1e-8),
        label = sprintf("%s objective actual=%.6g expected=%.6g", id, pair$objective, obj_exp)
      )
    }
    if (!is.null(exp[["scale"]])) {
      expect_equal(pair$transform$s, exp[["scale"]], tolerance = 1e-8)
    }
    if (!is.null(exp[["translation"]]) && isTRUE(pair$transform$spec$translation)) {
      expect_equal(as.numeric(pair$transform$t), as.numeric(exp[["translation"]]),
                   tolerance = 1e-8)
    }
    if (!is.null(exp[["rotation"]])) {
      Rexp <- read_matrix(exp[["rotation"]])
      dR <- rotation_distance(pair$transform$R, Rexp)
      compare <- exp[["rotation_compare"]] %||% "polar"
      if (identical(compare, "exact")) {
        expect_true(dR < 1e-8, label = paste(id, "R exact", dR))
      } else {
        expect_true(dR < 1e-7, label = paste(id, "R polar", dR))
      }
    }
    if (!is.null(exp[["effective_rank"]])) {
      expect_equal(pair$rank, as.integer(exp[["effective_rank"]]))
    }
    if (!is.null(exp[["transform_unique"]])) {
      expect_equal(pair$transform_unique, isTRUE(exp[["transform_unique"]]))
    }
    if (!is.null(exp[["unidentified_subspace_dimension"]])) {
      expect_equal(
        pair$unidentified_subspace_dimension,
        as.integer(exp[["unidentified_subspace_dimension"]])
      )
    }
    if (is.numeric(exp[["gamma"]])) {
      expect_equal(pair$gamma, exp[["gamma"]], tolerance = 1e-8)
    }
    if (is.numeric(exp[["a"]])) {
      expect_equal(pair$moments$a, exp[["a"]], tolerance = 1e-8)
    }
    if (is.numeric(exp[["b"]])) {
      expect_equal(pair$moments$b, exp[["b"]], tolerance = 1e-8)
    }
  }
})
