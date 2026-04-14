# S3: calibration studies for PCA / PLSC / CCA
#
# Paper Theorem 1 + Proposition 2 imply that under a true null the
# rejection rate of the collapsed P3 ladder should equal the nominal
# alpha up to Monte Carlo error. These tests generate many null
# datasets from the frozen Phase 0 bench generators, run the ladder
# at each, and check the empirical rate is within a binomial
# tolerance band around alpha.
#
# Skip on CRAN because they are expensive relative to unit tests.
# Local test run should finish in 1-3 minutes.

calibration_band <- function(alpha, n_sim, k = 3) {
  # k-SE binomial band around alpha. k = 3 makes the band loose enough
  # to survive CI noise; k = 2 is the formal 95% test but it is too
  # tight for n_sim as small as 200 and will flake on CI runners.
  se <- sqrt(alpha * (1 - alpha) / n_sim)
  c(lower = max(0, alpha - k * se), upper = min(1, alpha + k * se))
}

run_calibration <- function(alpha, n_sim, sim_fn, infer_fn) {
  rejections <- 0L
  for (i in seq_len(n_sim)) {
    dat <- sim_fn(i)
    res <- infer_fn(dat, alpha = alpha, seed = 1000L + i)
    if (sum(res$units$selected) >= 1L) {
      rejections <- rejections + 1L
    }
  }
  rejections / n_sim
}

test_that("oneblock_null calibration: empirical alpha within binomial band", {
  skip_on_cran()
  ensure_default_adapters()

  n_sim <- 200L
  alpha <- 0.05

  sim_fn <- function(i) {
    # Pure-noise oneblock via the frozen Phase 0 generator.
    bench_oneblock_null(n = 60, p = 12, seed = 5000L + i)$X
  }
  infer_fn <- function(dat, alpha, seed) {
    infer(
      adapter  = "svd_oneblock", data = dat,
      geometry = "oneblock", relation = "variance",
      B = 199L, R = 4L, alpha = alpha, seed = seed
    )
  }

  rate <- run_calibration(alpha, n_sim, sim_fn, infer_fn)
  band <- calibration_band(alpha, n_sim, k = 3)
  expect_gte(rate, band[["lower"]])
  expect_lte(rate, band[["upper"]])
})

test_that("cross_null covariance calibration: empirical alpha within binomial band", {
  skip_on_cran()
  ensure_default_adapters()

  n_sim <- 200L
  alpha <- 0.05

  sim_fn <- function(i) {
    dat <- bench_cross_null(
      n = 50, p_x = 6, p_y = 5,
      within_rank_x = 2, within_rank_y = 2,
      seed = 6000L + i
    )
    list(X = dat$X, Y = dat$Y)
  }
  infer_fn <- function(dat, alpha, seed) {
    infer(
      adapter  = "cross_svd", data = dat,
      geometry = "cross", relation = "covariance",
      B = 199L, R = 4L, alpha = alpha, seed = seed
    )
  }

  rate <- run_calibration(alpha, n_sim, sim_fn, infer_fn)
  band <- calibration_band(alpha, n_sim, k = 3)
  expect_gte(rate, band[["lower"]])
  expect_lte(rate, band[["upper"]])
})

test_that("cross_null correlation calibration: empirical alpha within binomial band", {
  skip_on_cran()
  ensure_default_adapters()

  n_sim <- 200L
  alpha <- 0.05

  sim_fn <- function(i) {
    dat <- bench_cross_null(
      n = 50, p_x = 5, p_y = 4,
      within_rank_x = 2, within_rank_y = 2,
      seed = 7000L + i
    )
    list(X = dat$X, Y = dat$Y)
  }
  infer_fn <- function(dat, alpha, seed) {
    infer(
      adapter  = "cross_svd", data = dat,
      geometry = "cross", relation = "correlation",
      B = 199L, R = 4L, alpha = alpha, seed = seed
    )
  }

  rate <- run_calibration(alpha, n_sim, sim_fn, infer_fn)
  band <- calibration_band(alpha, n_sim, k = 3)
  expect_gte(rate, band[["lower"]])
  expect_lte(rate, band[["upper"]])
})
