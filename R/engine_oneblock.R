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
#' Phase 1.5 TODO: replace full \code{base::svd()} calls in
#' \code{observed_stat_fn} and \code{null_stat_fn} with a partial SVD
#' (e.g. \code{RSpectra::svds(, k = 1)}) once the interface stabilises.
#' For Phase 1 (refit-first) the full decomposition is acceptable.
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

  ## --- validate recipe --------------------------------------------------------

  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "oneblock") {
    stop(
      paste0("`recipe` geometry must be \"oneblock\"; got \"",
             recipe$shape$geometry$kind, "\"."),
      call. = FALSE
    )
  }
  if (recipe$shape$relation$kind != "variance") {
    stop(
      paste0("`recipe` relation must be \"variance\"; got \"",
             recipe$shape$relation$kind, "\"."),
      call. = FALSE
    )
  }

  ## --- validate X -------------------------------------------------------------

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("`X` must not contain NA, NaN, or Inf.", call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 2L) {
    stop("`X` must have at least 2 rows and 2 columns.", call. = FALSE)
  }

  ## --- center X once ----------------------------------------------------------

  center <- colMeans(X)
  Xc     <- sweep(X, 2L, center, "-")

  ## --- full observed roots (for form_units) -----------------------------------

  roots_observed <- cached_svd(Xc)$d^2
  zero_tol <- max(1, sum(Xc^2)) * .Machine$double.eps

  ## --- max_steps default ------------------------------------------------------

  if (is.null(max_steps)) {
    max_steps <- min(min(dim(Xc)) - 1L, 50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  ## --- callbacks --------------------------------------------------------------

  # observed_stat_fn: top lambda of current residual / ||residual||_F^2.
  # This IS the Part 1 section 1 simplification applied to the deflated
  # residual E_a. The "first" singular value of E_a is lambda_a(X).
  observed_stat_fn <- function(step, data) {
    # Only the top-1 singular value is needed for the numerator; the
    # Frobenius norm gives the full tail-sum cheaply.
    s_top <- top_singular_values(data, 1L)[1L]
    total <- sum(data * data)
    if (total <= zero_tol) return(0)
    (s_top * s_top) / total
  }

  # null_stat_fn: column-permute the residual, then apply the collapsed
  # Vitale P3 null for rung `step`:
  #   lambda_step(M_b) / sum_{q >= step} lambda_q(M_b)
  # where M_b is the randomized residual. This is algebraically identical
  # to the paper's projected-null statistic without forming the projected
  # matrix or taking a second SVD.
  #
  # The denominator is `sum_{q >= step} s_q^2 = ||M_b||_F^2 - sum_{q < step} s_q^2`.
  # Column permutation preserves column sums of squares, so ||M_b||_F^2
  # equals ||data||_F^2 -- computed once per draw (O(np), negligible vs the
  # SVD). Only the top `step` singular values are needed, which opens the
  # door to a partial SVD (RSpectra::svds) on large problems.
  null_stat_fn <- function(step, data) {
    perm <- base::apply(data, 2L, base::sample)
    s    <- top_singular_values(perm, step)
    if (length(s) < step) return(0)
    s2   <- s * s
    total_F2 <- sum(data * data)
    tail_total <- total_F2 - sum(s2[seq_len(step - 1L)])
    if (tail_total <= zero_tol) return(0)
    s2[step] / tail_total
  }

  # deflate_fn: remove the rank-1 approximation from the top component.
  # Uses the top-1 partial SVD: data - d[1] * u[,1] %*% t(v[,1]).
  deflate_fn <- function(step, data) {
    sv <- top_svd(data, 1L)
    resid <- data - sv$d[1L] * (sv$u[, 1L, drop = FALSE] %*% t(sv$v[, 1L, drop = FALSE]))
    if (sum(resid^2) <= zero_tol) {
      return(matrix(0, nrow = nrow(data), ncol = ncol(data)))
    }
    resid
  }

  ## --- run the ladder ---------------------------------------------------------

  ladder_result <- ladder_driver(
    observed_stat_fn = observed_stat_fn,
    null_stat_fn     = null_stat_fn,
    deflate_fn       = deflate_fn,
    initial_data     = Xc,
    max_steps        = max_steps,
    B                = B,
    B_total          = B_total,
    batch_size       = batch_size,
    alpha            = alpha,
    seed             = seed
  )

  ## --- build selected vector --------------------------------------------------

  rejected_through <- ladder_result$rejected_through
  n_roots          <- length(roots_observed)

  # Roots 1..rejected_through are selected; the rest are not.
  selected <- logical(n_roots)
  if (rejected_through >= 1L) {
    selected[seq_len(rejected_through)] <- TRUE
  }

  ## --- form_units on the full observed roots ----------------------------------

  units <- form_units(roots_observed,
                      selected        = selected,
                      group_near_ties = isTRUE(auto_subspace),
                      tie_threshold   = tie_threshold)

  ## --- component_tests data.frame ---------------------------------------------

  sr <- ladder_result$step_results
  n_tested <- length(sr)

  component_tests <- data.frame(
    step          = vapply(sr, function(x) x$step,          integer(1L)),
    observed_stat = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value       = vapply(sr, function(x) x$p_value,       double(1L)),
    mc_se         = vapply(sr, function(x) x$mc_se,         double(1L)),
    r             = vapply(sr, function(x) x$r,             integer(1L)),
    B             = vapply(sr, function(x) x$B,             integer(1L)),
    selected      = vapply(sr, function(x) x$selected,      logical(1L)),
    stringsAsFactors = FALSE
  )

  ## --- return -----------------------------------------------------------------

  list(
    units           = units,
    component_tests = component_tests,
    roots_observed  = roots_observed,
    ladder_result   = ladder_result,
    labels          = list(
      statistic = "Vitale P3 tail-ratio on latent variance roots",
      null      = "column permutation of residual matrix",
      estimand  = "latent variance roots"
    )
  )
}
