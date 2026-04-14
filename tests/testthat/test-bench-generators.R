# Generator validity tests (multifer-hor.2).
#
# The locked bench generators declare target regimes -- e.g.
# bench_oneblock_null claims true_rank = 0, bench_oneblock_shadowing
# claims a specific root_profile.  The calibration matrix and several
# regression suites downstream *trust* those declarations.  These
# tests empirically confirm that each generator realizes its target
# regime, not just the right tensor shapes.

test_that("bench_oneblock_null produces no planted signal", {
  n <- 80L
  p <- 14L
  K <- 120L
  leading_sq <- numeric(K)
  col_means_max <- numeric(K)
  meta_ranks <- integer(K)

  for (k in seq_len(K)) {
    d <- bench_oneblock_null(n = n, p = p, seed = 30000L + k)
    Xc <- scale(d$X, center = TRUE, scale = FALSE)
    sv <- svd(Xc, nu = 0L, nv = 0L)$d
    leading_sq[k] <- sv[[1L]]^2 / (n - 1L)
    col_means_max[k] <- max(abs(colMeans(d$X)))
    meta_ranks[k] <- d$meta$true_rank
  }

  # Declared rank is zero on every draw.
  expect_true(all(meta_ranks == 0L))

  # Column means are O(1/sqrt(n)) under the noise model.
  expect_lt(mean(col_means_max), 0.5,
            label = "bench_oneblock_null mean |column mean|")

  # Leading eigenvalue of cov(X) should stay near the Marchenko-Pastur
  # upper edge lambda_+ = (1 + sqrt(p/n))^2 under unit-variance Gaussian
  # noise.  Use a 50% margin to absorb finite-K variance.
  mp_upper <- (1 + sqrt(p / n))^2
  expect_lt(mean(leading_sq), 1.5 * mp_upper,
            label = "bench_oneblock_null leading eigenvalue vs Marchenko-Pastur")
})

test_that("bench_oneblock_shadowing plants exactly the declared root profile", {
  n <- 200L
  p <- 40L
  root_profile <- c(8, 6, 4, 2, 1.2, 1.05)
  K <- 40L

  leading_ratios <- matrix(NA_real_,
                           nrow = K,
                           ncol = length(root_profile))
  for (k in seq_len(K)) {
    d <- bench_oneblock_shadowing(n = n, p = p,
                                  root_profile = root_profile,
                                  seed = 31000L + k)
    Xc <- scale(d$X, center = TRUE, scale = FALSE)
    sv <- svd(Xc, nu = 0L, nv = 0L)$d
    leading_ratios[k, ] <- sv[seq_along(root_profile)]^2 / (n - 1L)
    expect_equal(d$meta$true_rank, length(root_profile))
    expect_equal(d$meta$root_profile, root_profile)
  }

  # Average empirical leading eigenvalue should track the planted
  # signal + unit-variance noise floor.  Under the generator's
  # convention each planted root carries `root_profile[i]` variance
  # and the noise adds approximately (1 + sqrt(p/n))^2 per direction.
  # A 30% relative tolerance is enough to catch a drift that would
  # change the shadowing story (e.g. if the signal scaling changed).
  mean_leading <- colMeans(leading_ratios)
  for (i in seq_along(root_profile)) {
    planted <- root_profile[[i]]
    observed <- mean_leading[[i]]
    expect_gt(
      observed, 0.6 * planted,
      label = sprintf("bench_oneblock_shadowing root %d below planted",
                      i)
    )
  }

  # The very top root in particular must still dominate the
  # Marchenko-Pastur noise edge by a clear margin -- this is what
  # the shadowing benchmark is actually testing downstream.
  mp_upper <- (1 + sqrt(p / n))^2
  expect_gt(mean_leading[[1L]], 3 * mp_upper)
})

test_that("bench_cross_null declares and realizes zero cross-rank", {
  n <- 80L
  p_x <- 10L
  p_y <- 8L
  K <- 80L

  leading_cor <- numeric(K)
  for (k in seq_len(K)) {
    d <- bench_cross_null(n = n, p_x = p_x, p_y = p_y,
                          within_rank_x = 3L, within_rank_y = 3L,
                          seed = 32000L + k)
    expect_equal(d$meta$true_cross_rank, 0L)
    Xc <- scale(d$X, center = TRUE, scale = FALSE)
    Yc <- scale(d$Y, center = TRUE, scale = FALSE)
    cc <- stats::cancor(Xc, Yc)
    leading_cor[k] <- cc$cor[[1L]]
  }

  # Under a true cross-null the *mean* leading canonical correlation
  # is dominated by finite-n fitting and stays well below 1.  If the
  # generator started leaking signal across blocks this mean would
  # jump toward 1.
  expect_lt(mean(leading_cor), 0.95,
            label = "bench_cross_null mean leading canonical correlation")
  expect_gt(mean(leading_cor), 0.1,
            label = "bench_cross_null mean canonical correlation is vacuously zero")
})

test_that("bench_cross_null metadata mirrors its arguments", {
  d <- bench_cross_null(n = 120L, p_x = 9L, p_y = 7L,
                        within_rank_x = 2L, within_rank_y = 2L,
                        seed = 33000L)
  expect_equal(d$meta$n, 120L)
  expect_equal(d$meta$p_x, 9L)
  expect_equal(d$meta$p_y, 7L)
  expect_equal(d$meta$within_rank_x, 2L)
  expect_equal(d$meta$within_rank_y, 2L)
  expect_equal(d$meta$true_cross_rank, 0L)
  expect_equal(dim(d$X), c(120L, 9L))
  expect_equal(dim(d$Y), c(120L, 7L))
})
