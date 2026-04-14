test_that("infer_pca matches direct infer() call", {
  ensure_default_adapters()
  set.seed(901)
  X <- matrix(rnorm(120), 30, 4)

  wrapped <- infer_pca(X, B = 29L, R = 4L, seed = 17L)
  direct <- infer(
    adapter = "prcomp_oneblock",
    data = X,
    geometry = "oneblock",
    relation = "variance",
    B = 29L,
    R = 4L,
    seed = 17L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
})

test_that("infer_plsc matches direct infer() call", {
  ensure_default_adapters()
  set.seed(902)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)

  wrapped <- infer_plsc(X, Y, B = 29L, R = 4L, seed = 19L)
  direct <- infer(
    adapter = "cross_svd",
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "covariance",
    B = 29L,
    R = 4L,
    seed = 19L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
})

test_that("infer_cca matches direct infer() call", {
  ensure_default_adapters()
  set.seed(903)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)

  wrapped <- infer_cca(X, Y, B = 29L, R = 4L, seed = 23L)
  direct <- infer(
    adapter = "cancor_cross",
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "correlation",
    B = 29L,
    R = 4L,
    seed = 23L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
  expect_equal(wrapped$provenance$adapter_id, "cancor_cross")
})

test_that("infer_cca forwards supported nuisance-adjusted designs through the default adapter", {
  ensure_default_adapters()
  set.seed(904)
  n <- 36
  Z <- cbind(1, scale(seq_len(n)), sin(seq_len(n) / 5))
  X <- matrix(rnorm(n * 5L), n, 5L) + 0.2 * Z %*% matrix(rnorm(ncol(Z) * 5L), ncol(Z), 5L)
  Y <- matrix(rnorm(n * 4L), n, 4L) + 0.2 * Z %*% matrix(rnorm(ncol(Z) * 4L), ncol(Z), 4L)

  wrapped <- infer_cca(
    X, Y,
    design = nuisance_adjusted(Z),
    B = 29L,
    R = 4L,
    seed = 31L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$provenance$adapter_id, "cancor_cross")
  expect_true(all(wrapped$component_tests$p_value >= 0 &
                  wrapped$component_tests$p_value <= 1))
})
