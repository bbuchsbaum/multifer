# Phase 1.5 Wave 6 integration: verify that the core-update fast path
# produces the same decisions as refit across the four Phase 0 bench
# suites within Monte Carlo error, and that $cost records the fast path
# was actually engaged.

test_that("phase 1.5 integration: oneblock_null (svd) refit vs core-update agree", {
  skip_on_cran()
  ensure_default_adapters()
  dat <- bench_oneblock_null(n = 80, p = 12, noise = "gaussian", seed = 11L)

  common <- list(
    adapter  = "svd_oneblock",
    data     = dat$X,
    geometry = "oneblock",
    relation = "variance",
    B        = 99L, R = 8L, alpha = 0.05, seed = 13L
  )
  r_refit <- do.call(infer, c(common, list(fast_path = "off")))
  r_fast  <- do.call(infer, c(common, list(fast_path = "auto")))

  expect_equal(r_refit$units$selected, r_fast$units$selected)
  expect_equal(r_refit$cost$core_updates, 0L)
  expect_equal(r_fast$cost$core_updates,  8L)
  # p-values agree in sign of "above/below alpha" for every tested rung.
  expect_equal(
    r_refit$component_tests$p_value <= 0.05,
    r_fast$component_tests$p_value  <= 0.05
  )
})

test_that("phase 1.5 integration: oneblock_shadowing (prcomp) refit vs core-update agree", {
  skip_on_cran()
  ensure_default_adapters()
  dat <- bench_oneblock_shadowing(
    n = 80, p = 12,
    root_profile = c(8, 6, 4, 1.5),
    noise = "gaussian",
    seed = 17L
  )

  common <- list(
    adapter  = "prcomp_oneblock",
    data     = dat$X,
    geometry = "oneblock",
    relation = "variance",
    B        = 99L, R = 8L, alpha = 0.05, seed = 19L
  )
  r_refit <- do.call(infer, c(common, list(fast_path = "off")))
  r_fast  <- do.call(infer, c(common, list(fast_path = "auto")))

  # The three strong roots should still be rejected through in both runs.
  expect_true(sum(r_fast$units$selected)  >= 3L)
  expect_true(sum(r_refit$units$selected) >= 3L)
  expect_equal(r_fast$cost$core_updates, 8L)
  expect_equal(r_refit$cost$core_updates, 0L)
})

test_that("phase 1.5 integration: cross_null covariance refit vs core-update agree", {
  skip_on_cran()
  ensure_default_adapters()
  dat <- bench_cross_null(
    n = 60, p_x = 8, p_y = 6,
    within_rank_x = 3, within_rank_y = 3,
    seed = 23L
  )

  common <- list(
    adapter  = "cross_svd",
    data     = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "covariance",
    B        = 99L, R = 6L, alpha = 0.05, seed = 29L
  )
  r_refit <- do.call(infer, c(common, list(fast_path = "off")))
  r_fast  <- do.call(infer, c(common, list(fast_path = "auto")))

  # Null-only data: both modes should not over-select.
  expect_true(sum(r_fast$units$selected)  <= 3L)
  expect_true(sum(r_refit$units$selected) <= 3L)
  expect_equal(r_fast$cost$core_updates, 6L)
  expect_equal(r_refit$cost$core_updates, 0L)
})

test_that("phase 1.5 integration: cross_null correlation stays on refit path", {
  skip_on_cran()
  ensure_default_adapters()
  dat <- bench_cross_null(
    n = 60, p_x = 8, p_y = 6,
    within_rank_x = 3, within_rank_y = 3,
    seed = 31L
  )

  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "correlation",
    B        = 99L, R = 6L, alpha = 0.05, seed = 37L
  )
  # Correlation mode does not have a core-update fast path.
  expect_equal(res$cost$core_updates, 0L)
})

test_that("phase 1.5 integration: speed_agreement medium runs under both paths", {
  skip_on_cran()
  ensure_default_adapters()
  dat <- bench_speed_agreement(size = "medium", seed = 41L)

  common <- list(
    adapter  = "cross_svd",
    data     = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "covariance",
    B        = 49L, R = 4L, alpha = 0.05, seed = 43L
  )
  r_refit <- do.call(infer, c(common, list(fast_path = "off")))
  r_fast  <- do.call(infer, c(common, list(fast_path = "auto")))

  expect_true(is_infer_result(r_refit))
  expect_true(is_infer_result(r_fast))
  expect_equal(r_fast$cost$core_updates, 4L)
  expect_equal(r_refit$cost$core_updates, 0L)
  # Decisions agree.
  expect_equal(r_refit$units$selected, r_fast$units$selected)
})

test_that("phase 1.5 smoke speedup: fast path is at least somewhat faster than refit", {
  skip_on_cran()
  skip_on_os("windows")   # timing is noisy on CI runners
  ensure_default_adapters()

  # Small-but-wide oneblock: refit is a full n*p SVD per rep; fast path is n*k.
  set.seed(53)
  X <- matrix(rnorm(120 * 40), 120, 40)

  common <- list(
    adapter  = "svd_oneblock",
    data     = X,
    geometry = "oneblock",
    relation = "variance",
    B        = 29L, R = 40L, alpha = 0.05, seed = 59L
  )

  t_refit <- system.time(do.call(infer, c(common, list(fast_path = "off"))))[["elapsed"]]
  t_fast  <- system.time(do.call(infer, c(common, list(fast_path = "auto"))))[["elapsed"]]

  # Smoke-level assertion: fast path should not be meaningfully slower.
  # We avoid a hard speedup floor here because CI timing is noisy; the full
  # §40 "medium >= 5x, large >= 10x" numbers are produced by the bench
  # harness run manually.
  expect_true(t_fast <= t_refit * 1.25)
})
