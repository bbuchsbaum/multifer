# Benchmark suite: speed_agreement
#
# Purpose: compare the (future) core-update path against the refit path.
# Measures (a) wall-clock speedup and (b) decision agreement within
# Monte Carlo error. The comparison logic lives in the engine, not here.
#
# Size specifications:
#   medium: n = 500,  p_x = 100, p_y = 80,  true_rank = 4
#   large:  n = 2000, p_x = 500, p_y = 400, true_rank = 6
#
# This file is the canonical locked-suite artifact for Phase 0.
# The generator function and targets list are also compiled into
# the package as R/bench_generators.R so they are available as
# ordinary package functions.
#
# Generator: bench_speed_agreement(size, seed)
# Targets:   bench_speed_agreement_targets

.random_orthonormal <- function(nrow, ncol) {
  M <- matrix(stats::rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
  qr.Q(qr(M))
}

#' @rdname bench_generators
#' @export
bench_speed_agreement <- function(size = c("medium", "large"), seed = NULL) {
  size <- match.arg(size)
  if (!is.null(seed)) set.seed(seed)

  dims <- switch(size,
    medium = list(n = 500L,  p_x = 100L, p_y = 80L,  true_rank = 4L),
    large  = list(n = 2000L, p_x = 500L, p_y = 400L, true_rank = 6L)
  )
  n  <- dims$n
  px <- dims$p_x
  py <- dims$p_y
  k  <- dims$true_rank

  Ux <- .random_orthonormal(n, k)
  Uy <- .random_orthonormal(n, k)
  Vx <- .random_orthonormal(px, k)
  Vy <- .random_orthonormal(py, k)
  sv <- seq(from = 5, by = -0.5, length.out = k)

  X <- Ux %*% diag(sv, nrow = k, ncol = k) %*% t(Vx) +
       matrix(stats::rnorm(n * px), nrow = n, ncol = px)
  Y <- Uy %*% diag(sv, nrow = k, ncol = k) %*% t(Vy) +
       matrix(stats::rnorm(n * py), nrow = n, ncol = py)

  list(
    X    = X,
    Y    = Y,
    meta = list(size      = size,
                n         = n,
                p_x       = px,
                p_y       = py,
                true_rank = k,
                seed      = seed)
  )
}

#' @rdname bench_generators
#' @export
bench_speed_agreement_targets <- list(
  decision_agreement  = "within_mc_error",
  speedup_medium_min  = 5,
  speedup_large_min   = 10,
  measures            = c("wall_time_refit", "wall_time_core_update",
                           "decisions_agree", "mc_error_bound"),
  locked              = TRUE
)
