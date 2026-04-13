test_that("run_cross_ladder validates recipe geometry", {
  ensure_default_adapters()
  rec_one <- infer_recipe(geometry = "oneblock", relation = "variance",
                          adapter = "prcomp_oneblock")
  X <- matrix(rnorm(40), 10, 4)
  Y <- matrix(rnorm(40), 10, 4)
  expect_error(run_cross_ladder(rec_one, X, Y), "geometry must be")
})

test_that("run_cross_ladder validates input matrices", {
  ensure_default_adapters()
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  expect_error(run_cross_ladder(rec, "not a matrix", matrix(0, 2, 2)),
               "must be a numeric matrix")
  expect_error(
    run_cross_ladder(rec, matrix(0, 5, 3), matrix(0, 4, 3)),
    "same number of rows"
  )
})

test_that("run_cross_ladder runs on covariance recipe and returns expected schema", {
  ensure_default_adapters()
  set.seed(301)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, X, Y, B = 99L, alpha = 0.05, seed = 7)
  expect_true(is.list(res))
  for (nm in c("units", "component_tests", "roots_observed", "ladder_result")) {
    expect_true(nm %in% names(res))
  }
  expect_s3_class(res$units, "multifer_units")
  expect_true(nrow(res$component_tests) >= 1L)
})

test_that("run_cross_ladder runs on correlation recipe", {
  ensure_default_adapters()
  set.seed(302)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "correlation",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, X, Y, B = 99L, alpha = 0.05, seed = 7)
  expect_true(nrow(res$component_tests) >= 1L)
  expect_true(all(res$component_tests$p_value >= 0 &
                  res$component_tests$p_value <= 1))
})

test_that("run_cross_ladder is reproducible under fixed seed", {
  ensure_default_adapters()
  set.seed(303)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  r1 <- run_cross_ladder(rec, X, Y, B = 49L, seed = 11)
  r2 <- run_cross_ladder(rec, X, Y, B = 49L, seed = 11)
  expect_equal(r1$component_tests$p_value, r2$component_tests$p_value)
  expect_equal(r1$component_tests$observed_stat, r2$component_tests$observed_stat)
})

test_that("run_cross_ladder null calibration on bench_cross_null", {
  ensure_default_adapters()
  dat <- bench_cross_null(n = 60, p_x = 8, p_y = 6,
                          within_rank_x = 3, within_rank_y = 3,
                          seed = 777)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, dat$X, dat$Y, B = 99L, alpha = 0.05, seed = 7)
  # On a single null realization a few false positives are tolerable.
  expect_true(res$ladder_result$rejected_through <= 3L)
})
