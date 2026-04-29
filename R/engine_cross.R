# CCA / nuisance design helpers
#
# The exact cross-covariance core-space helpers
# (.cross_covariance_core*) live in R/core_space.R. This file keeps
# only the CCA-specific permutation and nuisance-design helpers below,
# plus the run_cross_ladder engine.

.restricted_row_permutation <- function(groups) {
  groups <- as.factor(groups)
  idx <- seq_along(groups)
  unsplit(lapply(split(idx, groups), function(g) {
    if (length(g) <= 1L) {
      return(g)
    }
    sample(g, length(g))
  }), groups)
}

.theil_keep_indices <- function(groups, n_keep) {
  groups <- as.factor(groups)
  idx_by_group <- split(seq_along(groups), groups)
  keep <- integer(0)
  pos <- rep(1L, length(idx_by_group))

  while (length(keep) < n_keep) {
    progressed <- FALSE
    for (g in seq_along(idx_by_group)) {
      if (pos[g] <= length(idx_by_group[[g]])) {
        keep <- c(keep, idx_by_group[[g]][pos[g]])
        pos[g] <- pos[g] + 1L
        progressed <- TRUE
        if (length(keep) == n_keep) break
      }
    }
    if (!progressed) {
      stop("Unable to construct a Theil selection with the requested reduced size.",
           call. = FALSE)
    }
  }

  sort(keep)
}

.nuisance_residual_basis <- function(Z, n, groups = NULL, tol = 1e-10) {
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("For nuisance-adjusted correlation designs, `design$Z` must be a numeric matrix.",
         call. = FALSE)
  }
  if (nrow(Z) != n) {
    stop("For nuisance-adjusted correlation designs, `design$Z` must have one row per observation.",
         call. = FALSE)
  }
  qz <- qr(Z, tol = tol)
  r <- qz$rank
  if (r >= n) {
    return(list(Q = matrix(0, nrow = n, ncol = 0L), groups = NULL))
  }

  if (is.null(groups)) {
    q_full <- qr.Q(qz, complete = TRUE)
    return(list(Q = q_full[, seq.int(r + 1L, n), drop = FALSE], groups = NULL))
  }

  groups <- as.factor(groups)
  if (length(groups) != n) {
    stop("For nuisance-adjusted structured designs, `groups` must have one value per observation.",
         call. = FALSE)
  }
  if (anyNA(groups)) {
    stop("For nuisance-adjusted structured designs, `groups` must not contain NA.",
         call. = FALSE)
  }

  n_keep <- n - r
  keep_idx <- .theil_keep_indices(groups, n_keep)
  qz_full <- qr.Q(qz, complete = FALSE)
  Rz <- diag(n) - qz_full %*% t(qz_full)
  K <- Rz[keep_idx, keep_idx, drop = FALSE]
  eig <- eigen((K + t(K)) / 2, symmetric = TRUE)
  vals <- pmax(eig$values, 0)
  if (any(vals <= tol)) {
    stop("Structured nuisance basis is rank-deficient for the chosen selection matrix.",
         call. = FALSE)
  }
  K_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(eig$vectors)
  Q <- Rz[, keep_idx, drop = FALSE] %*% K_inv_sqrt
  list(Q = Q, groups = droplevels(groups[keep_idx]))
}

#' Run the cross sequential deflation engine
#'
#' Implements the test ladder for (cross, covariance) and
#' (cross, correlation) recipes via the shared
#' \code{\link{ladder_driver}}. Both relations share the same scaffold
#' but use different per-step statistics:
#'
#' \itemize{
#'   \item covariance: squared leading singular value of the centered
#'     cross-product \code{X^t Y} on the deflated residual.
#'   \item correlation: squared leading singular value of the
#'     orthonormal cross-product \code{Q_x^t Q_y} on the deflated
#'     residual, where \code{Q_x} and \code{Q_y} are the thin QR
#'     factors of the centered residual blocks.
#' }
#'
#' For covariance-mode cross problems, the top singular value of the
#' deflated cross-block matrix at step \code{a} is the \code{a}-th
#' latent root of the original problem, so root-strength testing on the
#' deflated problem gives the intended stepwise ladder.
#'
#' For correlation-mode cross problems, valid multi-root testing requires
#' stepwise residualization in the whitened CCA space. `multifer`
#' implements that ladder for the plain paired-row design, and for
#' nuisance-adjusted designs it first projects both blocks into an
#' exchangeability-preserving residual basis before running the same
#' stepwise ladder. More complex designs still fall back to the
#' conservative first-root cap.
#'
#' Phase 1 is refit-first; the Part 2 section 9 core-space update is
#' a Phase 1.5 optimization.
#'
#' @param recipe A compiled \code{multifer_infer_recipe} with geometry
#'   \code{"cross"} and relation \code{"covariance"} or
#'   \code{"correlation"}.
#' @param X Numeric matrix, n x p.
#' @param Y Numeric matrix, n x q. Must have the same number of rows
#'   as \code{X}.
#' @param B Per-rung cap on Monte Carlo draws.
#' @param B_total Optional integer global Monte Carlo budget shared
#'   across ladder rungs. Defaults to `B * max_steps`.
#' @param batch_size Positive integer, Besag-Clifford batch size within
#'   each rung. Default `32L`.
#' @param alpha Significance threshold.
#' @param max_steps Maximum ladder rungs. Default
#'   \code{min(nrow, ncol(X), ncol(Y)) - 1}, capped at 50.
#' @param seed Integer or NULL.
#' @param cross_rank_cap Reserved opt-in for an approximate rank-capped
#'   core path on very large problems. Default `NULL` (exact only, no
#'   approximation). Setting a positive integer enables a truncated
#'   block-SVD pre-projection whose top-1 statistic is biased and
#'   whose use is NOT paper-faithful; it exists as a screening knob
#'   only.
#' @param auto_subspace Logical. When `TRUE` (default), near-tied
#'   roots are automatically bundled into a subspace unit via
#'   [form_units()] `group_near_ties = TRUE`.
#' @param tie_threshold Positive numeric. Relative-gap threshold used
#'   when `auto_subspace = TRUE`. Default `0.01`.
#'
#' @return A list with \code{units}, \code{component_tests},
#'   \code{roots_observed}, and \code{ladder_result}, mirroring
#'   \code{\link{run_oneblock_ladder}}.
#'
#' @export
run_cross_ladder <- function(recipe,
                             X,
                             Y,
                             B              = 1000L,
                             B_total        = NULL,
                             batch_size     = 32L,
                             alpha          = 0.05,
                             max_steps      = NULL,
                             seed           = NULL,
                             cross_rank_cap = NULL,
                             auto_subspace  = TRUE,
                             tie_threshold  = 0.01) {

  plan <- compile_cross_ladder_plan(
    recipe = recipe,
    X = X,
    Y = Y,
    max_steps = max_steps,
    cross_rank_cap = cross_rank_cap
  )

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
