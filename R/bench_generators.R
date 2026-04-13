#' Benchmark generators and locked target specifications
#'
#' This file contains four benchmark suite generators and their
#' corresponding locked target lists. The generators produce synthetic
#' datasets; the target lists define what a future inference engine must
#' achieve. These are canonical Phase 0 artifacts: the target numbers
#' are frozen and must not be changed without a protocol revision.
#'
#' Suite overview:
#' - `bench_oneblock_null`: null calibration for one-block shapes
#' - `bench_oneblock_shadowing`: power under the Vitale shadowing case
#' - `bench_cross_null`: null calibration for cross-block shapes
#' - `bench_speed_agreement`: core-update vs. refit speedup + agreement
#'
#' @name bench_generators
NULL

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.random_orthonormal <- function(nrow, ncol, rng_state = NULL) {
  # Returns a random orthonormal nrow x ncol matrix via QR decomposition.
  M <- matrix(stats::rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
  qr.Q(qr(M))
}

# ---------------------------------------------------------------------------
# Suite 1: oneblock_null
# ---------------------------------------------------------------------------

#' Generator: one-block null calibration
#'
#' Generates a pure-noise dataset for testing over-selection rate and
#' p-value calibration in one-block inference. There is no real signal;
#' the true rank is zero.
#'
#' @param n Integer. Number of observations (rows). Default 200.
#' @param p Integer. Number of variables (columns). Default 50.
#' @param noise Character scalar. Noise model: `"gaussian"` (standard
#'   normal), `"heavy_tailed"` (t with df = 5), or `"heteroscedastic"`
#'   (Gaussian with per-column variances drawn from an exponential).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{X}{Numeric matrix of dimensions n x p.}
#'     \item{meta}{List with fields `n`, `p`, `noise`, `seed`,
#'       `true_rank` (always 0L).}
#'   }
#'
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

#' Locked targets: one-block null calibration
#'
#' Defines what a future inference engine must achieve on data from
#' [bench_oneblock_null()]. The `locked` flag is `TRUE`; these numbers
#' are frozen.
#'
#' @format A list.
#' @export
bench_oneblock_null_targets <- list(
  over_selection_rate = 0.05,
  alpha               = c(0.01, 0.05, 0.10),
  measures            = c("over_selection_rate", "p_value_calibration_ks"),
  locked              = TRUE
)

# ---------------------------------------------------------------------------
# Suite 2: oneblock_shadowing
# ---------------------------------------------------------------------------

#' Generator: one-block shadowing (Vitale shadowing case)
#'
#' Generates a dataset where strong early roots hide weaker late roots.
#' This is the primary stress test for deflation: the weakest components
#' (lowest entries in `root_profile`) are the "shadowed" targets.
#'
#' Data construction: `X = U diag(root_profile) V^T + noise_matrix`
#' where `U` is a random orthonormal `n x k` matrix, `V` is a random
#' orthonormal `p x k` matrix, and `k = length(root_profile)`.
#'
#' @param n Integer. Number of observations (rows). Default 200.
#' @param p Integer. Number of variables (columns). Default 50.
#' @param root_profile Numeric vector. Signal-to-noise ratio of each
#'   true component. Default `c(8, 6, 4, 2, 1.2, 1.05)`.
#' @param noise Character scalar. `"gaussian"` or `"heteroscedastic"`.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{X}{Numeric matrix of dimensions n x p.}
#'     \item{meta}{List with fields `n`, `p`, `noise`, `seed`,
#'       `true_rank` (`length(root_profile)`), `root_profile`.}
#'   }
#'
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

#' Locked targets: one-block shadowing
#'
#' Defines what a future inference engine must achieve on data from
#' [bench_oneblock_shadowing()]. The `locked` flag is `TRUE`; these
#' numbers are frozen.
#'
#' @format A list.
#' @export
bench_oneblock_shadowing_targets <- list(
  recover_weakest = c(0.8, 0.95),
  alpha           = 0.05,
  measures        = c("power_per_component", "rank_recovery_rate",
                      "shadow_detection_rate"),
  locked          = TRUE
)

# ---------------------------------------------------------------------------
# Suite 3: cross_null
# ---------------------------------------------------------------------------

#' Generator: cross-block null calibration
#'
#' Generates two blocks X and Y that each have their own low-rank
#' within-block structure but are drawn INDEPENDENTLY, so there is no
#' true cross-block association. The true cross rank is zero.
#'
#' This catches the failure mode where naive permutation CCA inflates
#' error beyond the first root.
#'
#' @param n Integer. Number of paired observations (rows). Default 200.
#' @param p_x Integer. Number of variables in X. Default 40.
#' @param p_y Integer. Number of variables in Y. Default 30.
#' @param within_rank_x Integer. Rank of within-block structure in X.
#'   Default 5.
#' @param within_rank_y Integer. Rank of within-block structure in Y.
#'   Default 5.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{X}{Numeric matrix of dimensions n x p_x.}
#'     \item{Y}{Numeric matrix of dimensions n x p_y.}
#'     \item{meta}{List with fields `n`, `p_x`, `p_y`,
#'       `within_rank_x`, `within_rank_y`, `seed`,
#'       `true_cross_rank` (always 0L).}
#'   }
#'
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

#' Locked targets: cross-block null calibration
#'
#' Defines what a future inference engine must achieve on data from
#' [bench_cross_null()]. The `locked` flag is `TRUE`; these numbers
#' are frozen.
#'
#' @format A list.
#' @export
bench_cross_null_targets <- list(
  over_selection_rate      = 0.05,
  alpha                    = c(0.01, 0.05, 0.10),
  measures                 = c("over_selection_rate_per_root",
                                "p_value_calibration_ks"),
  failure_modes_detected   = c("naive_permutation_inflation_beyond_first_root"),
  locked                   = TRUE
)

# ---------------------------------------------------------------------------
# Suite 4: speed_agreement
# ---------------------------------------------------------------------------

#' Generator: speed and decision-agreement benchmark
#'
#' Generates a cross-block dataset at "medium" or "large" scale for
#' benchmarking wall-clock speedup of the core-update path relative to
#' the full refit path. The statistical subtlety is minimal; scale is
#' the point.
#'
#' Size specifications:
#' - `"medium"`: n = 500, p_x = 100, p_y = 80, true_rank = 4
#' - `"large"`:  n = 2000, p_x = 500, p_y = 400, true_rank = 6
#'
#' @param size Character scalar. `"medium"` or `"large"`.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{X}{Numeric matrix of dimensions n x p_x.}
#'     \item{Y}{Numeric matrix of dimensions n x p_y.}
#'     \item{meta}{List with fields `size`, `n`, `p_x`, `p_y`,
#'       `true_rank`, `seed`.}
#'   }
#'
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

#' Locked targets: speed and decision-agreement benchmark
#'
#' Defines what a future inference engine must achieve on data from
#' [bench_speed_agreement()]. The `locked` flag is `TRUE`; these numbers
#' are frozen.
#'
#' @format A list.
#' @export
bench_speed_agreement_targets <- list(
  decision_agreement  = "within_mc_error",
  speedup_medium_min  = 5,
  speedup_large_min   = 10,
  measures            = c("wall_time_refit", "wall_time_core_update",
                           "decisions_agree", "mc_error_bound"),
  locked              = TRUE
)
