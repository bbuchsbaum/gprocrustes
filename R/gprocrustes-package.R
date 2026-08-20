#' gprocrustes: a solver-transparent Procrustes engine
#'
#' Pairwise and generalized Procrustes analysis with explicit transformation
#' groups, gauges, missingness, and optimality claims. Unanchored multi-view
#' problems use Gower block-coordinate descent; two-view anchored problems use
#' the exact pairwise kernel.
#'
#' The mathematics in `docs/spec/00-mathematics.md` is the primary specification.
#'
#' @keywords internal
#' @importFrom stats coef fitted predict residuals setNames
#' @importFrom utils packageVersion
"_PACKAGE"

utils::globalVariables(".data")
