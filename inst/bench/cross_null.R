# Benchmark suite: cross_null
#
# Purpose: null calibration for cross-block shapes.
# Strong within-block structure but ZERO cross-block association.
# Catches the failure mode where naive permutation CCA inflates
# error beyond the first root.
#
# This file is the canonical locked-suite artifact for Phase 0.
# The generator function and targets list are also compiled into
# the package as R/bench_generators.R so they are available as
# ordinary package functions.
#
# Generator: bench_cross_null(n, p_x, p_y, within_rank_x, within_rank_y, seed)
# Targets:   bench_cross_null_targets

.random_orthonormal <- function(nrow, ncol) {
  M <- matrix(stats::rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
  qr.Q(qr(M))
}

#' @rdname bench_generators
#' @export
bench_cross_null <- function(n = 200, p_x = 40, p_y = 30,
                              within_rank_x = 5, within_rank_y = 5,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  .make_low_rank_block <- function(nr, nc, rk) {
    rk <- min(rk, nr, nc)
    U  <- .random_orthonormal(nr, rk)
    V  <- .random_orthonormal(nc, rk)
    sv <- seq(from = 3, by = -0.4, length.out = rk)
    signal <- U %*% diag(sv, nrow = rk, ncol = rk) %*% t(V)
    noise  <- matrix(stats::rnorm(nr * nc), nrow = nr, ncol = nc)
    signal + noise
  }

  X <- .make_low_rank_block(n, p_x, within_rank_x)
  Y <- .make_low_rank_block(n, p_y, within_rank_y)

  list(
    X    = X,
    Y    = Y,
    meta = list(n             = n,
                p_x           = p_x,
                p_y           = p_y,
                within_rank_x = within_rank_x,
                within_rank_y = within_rank_y,
                seed          = seed,
                true_cross_rank = 0L)
  )
}

#' @rdname bench_generators
#' @export
bench_cross_null_targets <- list(
  over_selection_rate      = 0.05,
  alpha                    = c(0.01, 0.05, 0.10),
  measures                 = c("over_selection_rate_per_root",
                                "p_value_calibration_ks"),
  failure_modes_detected   = c("naive_permutation_inflation_beyond_first_root"),
  locked                   = TRUE
)
