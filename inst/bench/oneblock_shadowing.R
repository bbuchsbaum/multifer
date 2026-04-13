# Benchmark suite: oneblock_shadowing
#
# Purpose: power under the Vitale shadowing case.
# Strong early roots hide weaker late roots. This is where deflation
# matters most. The weakest components (lowest entries in root_profile)
# are the "shadowed" targets.
#
# This file is the canonical locked-suite artifact for Phase 0.
# The generator function and targets list are also compiled into
# the package as R/bench_generators.R so they are available as
# ordinary package functions.
#
# Generator: bench_oneblock_shadowing(n, p, root_profile, noise, seed)
# Targets:   bench_oneblock_shadowing_targets

.random_orthonormal <- function(nrow, ncol) {
  M <- matrix(stats::rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
  qr.Q(qr(M))
}

#' @rdname bench_generators
#' @export
bench_oneblock_shadowing <- function(n = 200, p = 50,
                                      root_profile = c(8, 6, 4, 2, 1.2, 1.05),
                                      noise = c("gaussian", "heteroscedastic"),
                                      seed = NULL) {
  noise <- match.arg(noise)
  if (!is.null(seed)) set.seed(seed)

  k <- length(root_profile)
  if (k < 1L) stop("`root_profile` must have length >= 1.", call. = FALSE)
  if (k > min(n, p)) {
    stop(sprintf(
      "`root_profile` length (%d) exceeds min(n, p) = %d.",
      k, min(n, p)
    ), call. = FALSE)
  }

  U <- .random_orthonormal(n, k)
  V <- .random_orthonormal(p, k)

  signal <- U %*% diag(root_profile, nrow = k, ncol = k) %*% t(V)

  noise_mat <- switch(noise,
    gaussian = matrix(stats::rnorm(n * p), nrow = n, ncol = p),
    heteroscedastic = {
      col_sds <- sqrt(stats::rexp(p, rate = 1))
      raw <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
      sweep(raw, 2, col_sds, `*`)
    }
  )

  X <- signal + noise_mat

  list(
    X    = X,
    meta = list(n = n, p = p, noise = noise, seed = seed,
                true_rank    = k,
                root_profile = root_profile)
  )
}

#' @rdname bench_generators
#' @export
bench_oneblock_shadowing_targets <- list(
  recover_weakest = c(0.8, 0.95),
  alpha           = 0.05,
  measures        = c("power_per_component", "rank_recovery_rate",
                      "shadow_detection_rate"),
  locked          = TRUE
)
