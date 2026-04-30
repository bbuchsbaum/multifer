#' Run the oneblock sequential deflation engine
#'
#' Implements the Vitale P3 test ladder on (oneblock, variance) via
#' the shared \code{\link{ladder_driver}}. The observed test statistic at step
#' \code{a} is the relative tail-ratio
#'
#'   T_a = lambda_1(E_a) / ||E_a||_F^2
#'
#' where \code{E_a} is the deflated residual matrix at step \code{a}.
#' This is mathematically identical to the Part 1 section 1 formula
#' \code{lambda_a(X) / sum(lambda_q(X) for q >= a)}: after (a-1)
#' deflations the first singular value of the residual IS the a-th
#' singular value of the original matrix.
#'
#' The null generator permutes columns of \code{E_a} and evaluates the
#' collapsed Vitale P3 null, \eqn{T_a^*(M_b) = \lambda_a(M_b) / \sum_{q \ge a} \lambda_q(M_b)},
#' on the randomized residual \code{M_b}. No projected matrix is formed
#' and no second SVD is needed; only the first \code{a} singular values
#' of \code{M_b} are required.
#'
#' The compiled ladder plan now uses the Phase 1.5 partial-SVD helpers:
#' \code{top_singular_values()} for observed and null statistics, and
#' \code{top_svd()} for rank-1 deflation. These helpers use
#' \pkg{RSpectra} when available and fall back to \code{base::svd()},
#' preserving the exact ladder semantics while avoiding full decompositions
#' on large residual matrices.
#'
#' @param recipe A compiled \code{multifer_infer_recipe}. Must have
#'   geometry \code{"oneblock"} and relation \code{"variance"}.
#' @param X Numeric matrix, n x p. Centered internally.
#' @param B Per-rung cap on Monte Carlo draws. Default \code{1000}.
#' @param B_total Optional integer global Monte Carlo budget shared
#'   across ladder rungs. Defaults to `B * max_steps`.
#' @param batch_size Positive integer, Besag-Clifford batch size within
#'   each rung. Default `32L`.
#' @param alpha Significance threshold. Default \code{0.05}.
#' @param max_steps Maximum number of ladder rungs. Default
#'   \code{min(nrow(X), ncol(X)) - 1}, capped at 50 for Phase 1.
#' @param seed Integer or NULL.
#' @param auto_subspace Logical. When `TRUE` (default), near-tied
#'   roots are automatically bundled into a subspace unit via
#'   [form_units()] `group_near_ties = TRUE`.
#' @param tie_threshold Positive numeric. Relative-gap threshold used
#'   when `auto_subspace = TRUE`. Default `0.01`.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{units}{Result of \code{\link{form_units}} on the full observed
#'       roots vector, with the \code{selected} flag derived from the
#'       ladder results.}
#'     \item{component_tests}{data.frame with one row per tested rung:
#'       columns \code{step}, \code{observed_stat}, \code{p_value},
#'       \code{mc_se}, \code{r}, \code{B}, \code{selected}.}
#'     \item{roots_observed}{Numeric vector. Full vector of squared
#'       singular values \code{svd(X_centered)$d^2}.}
#'     \item{ladder_result}{Raw output of \code{\link{ladder_driver}},
#'       retained for diagnostics.}
#'   }
#'
#' @export
run_oneblock_ladder <- function(recipe,
                                X,
                                B             = 1000L,
                                B_total       = NULL,
                                batch_size    = 32L,
                                alpha         = 0.05,
                                max_steps     = NULL,
                                seed          = NULL,
                                auto_subspace = TRUE,
                                tie_threshold = 0.01) {

  plan <- compile_oneblock_ladder_plan(recipe, X, max_steps = max_steps)

  ## --- run the ladder ---------------------------------------------------------

  ladder_result <- ladder_driver(
    observed_stat_fn = plan$observed_stat_fn,
    null_stat_fn     = plan$null_stat_fn,
    deflate_fn       = plan$deflate_fn,
    initial_data     = plan$initial_data,
    max_steps        = plan$max_steps,
    B                = B,
    B_total          = B_total,
    batch_size       = batch_size,
    alpha            = alpha,
    seed             = seed
  )

  .ladder_plan_result(
    plan,
    ladder_result,
    auto_subspace = auto_subspace,
    tie_threshold = tie_threshold
  )
}
