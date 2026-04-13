test_that("infer() oneblock end-to-end produces a valid infer_result", {
  ensure_default_adapters()
  set.seed(701)
  signal <- matrix(rnorm(40 * 3), 40, 3) %*% diag(c(8, 6, 4)) %*%
    matrix(rnorm(3 * 6), 3, 6)
  X <- signal + 0.1 * matrix(rnorm(240), 40, 6)
  res <- infer(
    adapter  = "prcomp_oneblock",
    data     = X,
    geometry = "oneblock",
    relation = "variance",
    B        = 99L,
    R        = 10L,
    seed     = 7
  )
  expect_true(is_infer_result(res))
  expect_true(nrow(res$units) >= 1L)
  expect_true(nrow(res$component_tests) >= 1L)
  expect_true(nrow(res$variable_stability) >= 1L)
  expect_true(nrow(res$score_stability) >= 1L)
  expect_equal(nrow(res$subspace_stability), nrow(res$units))
  expect_equal(res$provenance$adapter_id, "prcomp_oneblock")
})

test_that("infer() cross end-to-end with explicit covariance", {
  ensure_default_adapters()
  set.seed(702)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = X, Y = Y),
    geometry = "cross",
    relation = "covariance",
    B        = 49L,
    R        = 6L,
    seed     = 9
  )
  expect_true(is_infer_result(res))
  expect_setequal(unique(res$variable_stability$domain), c("X", "Y"))
})

test_that("infer() reproduces under fixed seed", {
  ensure_default_adapters()
  set.seed(703)
  X <- matrix(rnorm(160), 40, 4)
  res1 <- infer(adapter = "prcomp_oneblock", data = X,
                geometry = "oneblock", relation = "variance",
                B = 49L, R = 5L, seed = 13)
  res2 <- infer(adapter = "prcomp_oneblock", data = X,
                geometry = "oneblock", relation = "variance",
                B = 49L, R = 5L, seed = 13)
  expect_equal(res1$component_tests$p_value, res2$component_tests$p_value)
  expect_equal(res1$component_tests$statistic, res2$component_tests$statistic)
})

test_that("infer() with strict dispatch errors on cross adapter without relation", {
  ensure_default_adapters()
  X <- matrix(rnorm(80), 20, 4)
  Y <- matrix(rnorm(60), 20, 3)
  expect_error(
    infer(adapter = "cross_svd", data = list(X = X, Y = Y),
          geometry = "cross", B = 9L, R = 2L, seed = 1),
    "covariance"
  )
})

test_that("infer() captures cost and mc reproducibility blocks", {
  ensure_default_adapters()
  set.seed(704)
  X <- matrix(rnorm(80), 20, 4)
  res <- infer(adapter = "prcomp_oneblock", data = X,
               geometry = "oneblock", relation = "variance",
               B = 49L, R = 5L, seed = 17)
  expect_true(!is.null(res$cost))
  expect_true(!is.null(res$mc))
  expect_equal(res$mc$rng_seed, 17L)
})

test_that("infer() rejects malformed adapter argument", {
  expect_error(infer(adapter = 123, data = matrix(0, 2, 2)),
               "registered adapter id")
})
