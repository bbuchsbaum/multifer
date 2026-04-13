# End-to-end Phase 1 integration tests: real adapters running through
# infer() against the four Phase 0 benchmark suites. These exercise
# the whole pipeline -- recipe compilation, sequential deflation
# engine, bootstrap perturbation, stability consumers, and result
# assembly -- on data drawn from each locked benchmark generator.
#
# B and R are kept modest so CI runs in seconds. The guarantees here
# are calibration / power on a single realization; they are not
# replacement for the larger-scale benchmark sweep that Phase 1.5+
# will run.

test_that("integration: bench_oneblock_null calibration via prcomp adapter", {
  skip_on_cran()
  ensure_default_adapters()

  dat <- bench_oneblock_null(n = 60, p = 12, seed = 11)

  res <- infer(
    adapter  = "prcomp_oneblock",
    data     = dat$X,
    geometry = "oneblock",
    relation = "variance",
    B        = 99L,
    R        = 8L,
    alpha    = 0.05,
    seed     = 17
  )

  expect_true(is_infer_result(res))
  # On a single null realization at alpha=0.05 with B=99 and p=12,
  # we tolerate up to ~3 false-positive component selections. Tight
  # bounds belong in the larger benchmark sweep.
  selected_n <- sum(res$units$selected)
  expect_true(selected_n <= 3L)
})

test_that("integration: bench_oneblock_shadowing recovers strong roots via prcomp", {
  skip_on_cran()
  ensure_default_adapters()

  dat <- bench_oneblock_shadowing(
    n = 80, p = 12,
    root_profile = c(8, 6, 4, 1.5),
    noise = "gaussian",
    seed = 23
  )

  res <- infer(
    adapter  = "prcomp_oneblock",
    data     = dat$X,
    geometry = "oneblock",
    relation = "variance",
    B        = 99L,
    R        = 8L,
    alpha    = 0.05,
    seed     = 29
  )

  # The three strong roots (8, 6, 4) should be recovered. The weakest
  # shadowed root (1.5) is harder; we don't assert on it under such
  # low B and R.
  expect_true(sum(res$units$selected) >= 3L)
})

test_that("integration: bench_cross_null calibration via cross_svd covariance", {
  skip_on_cran()
  ensure_default_adapters()

  dat <- bench_cross_null(
    n = 60, p_x = 8, p_y = 6,
    within_rank_x = 3, within_rank_y = 3,
    seed = 31
  )

  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "covariance",
    B        = 99L,
    R        = 6L,
    alpha    = 0.05,
    seed     = 37
  )

  selected_n <- sum(res$units$selected)
  expect_true(selected_n <= 3L)
})

test_that("integration: bench_cross_null calibration via cross_svd correlation", {
  skip_on_cran()
  ensure_default_adapters()

  dat <- bench_cross_null(
    n = 60, p_x = 8, p_y = 6,
    within_rank_x = 3, within_rank_y = 3,
    seed = 41
  )

  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "correlation",
    B        = 99L,
    R        = 6L,
    alpha    = 0.05,
    seed     = 43
  )

  selected_n <- sum(res$units$selected)
  expect_true(selected_n <= 3L)
})

test_that("integration: bench_speed_agreement medium runs and records cost", {
  skip_on_cran()
  ensure_default_adapters()

  dat <- bench_speed_agreement(size = "medium", seed = 47)

  # Phase 1 is refit-first only; we just check the engine completes
  # and records timing into the cost block. Speedup measurements vs
  # core-update are a Phase 1.5 deliverable.
  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "covariance",
    B        = 49L,
    R        = 4L,
    alpha    = 0.05,
    seed     = 53
  )

  expect_true(is_infer_result(res))
  expect_true(!is.null(res$cost$wall_time_phases))
  expect_true(any(names(res$cost$wall_time_phases) == "total"))
  expect_true(res$cost$wall_time_phases[["total"]] >= 0)
})

test_that("integration: multivarious::pca round-trips through infer()", {
  skip_on_cran()
  skip_if_not_installed("multivarious")
  ensure_default_adapters()

  set.seed(59)
  X <- matrix(rnorm(40 * 6), 40, 6)
  fit <- multivarious::pca(X, ncomp = 5)

  res <- infer(
    adapter  = "multivarious_pca",
    data     = X,
    geometry = "oneblock",
    relation = "variance",
    model    = fit,
    B        = 49L,
    R        = 4L,
    seed     = 61
  )
  expect_true(is_infer_result(res))
  expect_equal(res$provenance$adapter_id, "multivarious_pca")
})

test_that("integration: multivarious::plsc round-trips through infer()", {
  skip_on_cran()
  skip_if_not_installed("multivarious")
  ensure_default_adapters()

  set.seed(67)
  X <- matrix(rnorm(40 * 6), 40, 6)
  Y <- matrix(rnorm(40 * 5), 40, 5)
  fit <- multivarious::plsc(X, Y, ncomp = 4)

  res <- infer(
    adapter  = "multivarious_plsc",
    data     = list(X = X, Y = Y),
    geometry = "cross",
    model    = fit,
    B        = 49L,
    R        = 4L,
    seed     = 71
  )
  expect_true(is_infer_result(res))
  expect_equal(res$provenance$adapter_id, "multivarious_plsc")
  # multivarious_plsc declares only covariance; strict dispatch should
  # auto-resolve it without ambiguity.
  expect_equal(res$cost$wall_time_phases[["total"]] >= 0, TRUE)
})
