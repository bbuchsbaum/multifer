# Benchmark suite: oneblock_null
#
# Purpose: null calibration for one-block shapes.
# Pure noise, no real signal. Measures over-selection rate and
# p-value calibration at alpha in {0.01, 0.05, 0.10}.
#
# This file is the canonical locked-suite artifact for Phase 0.
# The generator function and targets list are also compiled into
# the package as R/bench_generators.R so they are available as
# ordinary package functions.
#
# Generator: bench_oneblock_null(n, p, noise, seed)
# Targets:   bench_oneblock_null_targets

#' @rdname bench_generators
#' @export
bench_oneblock_null <- function(n = 200, p = 50,
                                 noise = c("gaussian", "heavy_tailed",
                                           "heteroscedastic"),
                                 seed = NULL) {
  noise <- match.arg(noise)
  if (!is.null(seed)) set.seed(seed)

  X <- switch(noise,
    gaussian = matrix(stats::rnorm(n * p), nrow = n, ncol = p),
    heavy_tailed = matrix(stats::rt(n * p, df = 5), nrow = n, ncol = p),
    heteroscedastic = {
      col_sds <- sqrt(stats::rexp(p, rate = 1))
      raw <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
      sweep(raw, 2, col_sds, `*`)
    }
  )

  list(
    X    = X,
    meta = list(n = n, p = p, noise = noise, seed = seed, true_rank = 0L)
  )
}

#' @rdname bench_generators
#' @export
bench_oneblock_null_targets <- list(
  over_selection_rate = 0.05,
  alpha               = c(0.01, 0.05, 0.10),
  measures            = c("over_selection_rate", "p_value_calibration_ks"),
  locked              = TRUE
)
