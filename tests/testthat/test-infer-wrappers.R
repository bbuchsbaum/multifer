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
    adapter = "cross_svd",
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
})
