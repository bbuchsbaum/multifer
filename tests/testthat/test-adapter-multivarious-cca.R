test_that("adapter_multivarious_cca constructs cleanly", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  a <- adapter_multivarious_cca()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "multivarious_cca")
  expect_equal(a$shape_kinds, "cross")
  expect_equal(a$validity_level, "conditional")
})

test_that("adapter_multivarious_cca declares ONLY correlation", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_cca()
  expect_true(adapter_supports(a, "cross", "correlation",
                               "component_significance"))
  expect_false(adapter_supports(a, "cross", "covariance",
                                "component_significance"))
})

test_that("adapter_multivarious_cca declares all four v1 targets for correlation", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_cca()
  for (tgt in c("component_significance", "variable_stability",
                "score_stability", "subspace_stability")) {
    expect_true(adapter_supports(a, "cross", "correlation", tgt))
  }
})

test_that("register_multivarious_cca_adapter registers when multivarious installed", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  register_multivarious_cca_adapter()
  expect_true("multivarious_cca" %in% list_infer_adapters())
})

test_that("strict-dispatch has NO ambiguity for multivarious_cca", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  register_multivarious_cca_adapter()
  rec <- infer_recipe(geometry = "cross", adapter = "multivarious_cca")
  expect_true(is_infer_recipe(rec))
  expect_equal(rec$shape$relation$kind, "correlation")
})

test_that("multivarious_cca matches cancor on full-rank lambda = 0", {
  skip_if_not_installed("multivarious")
  set.seed(301)
  X <- matrix(rnorm(160), 40, 4)
  Y <- matrix(rnorm(120), 40, 3)

  fit_mv <- multivarious::cca(X, Y, ncomp = 3L, lambda = 0)
  fit_cc <- stats::cancor(X, Y)
  a <- adapter_multivarious_cca(lambda = 0, ncomp = 3L)

  expect_equal(a$roots(fit_mv), fit_cc$cor[1:3], tolerance = 1e-8)
  expect_equal(dim(a$scores(fit_mv, "X")), c(40L, 3L))
  expect_equal(dim(a$scores(fit_mv, "Y")), c(40L, 3L))
  expect_equal(dim(a$loadings(fit_mv, "X")), c(4L, 3L))
  expect_equal(dim(a$loadings(fit_mv, "Y")), c(3L, 3L))
})

test_that("truncate hook keeps paired state consistent", {
  skip_if_not_installed("multivarious")
  set.seed(302)
  X <- matrix(rnorm(160), 40, 4)
  Y <- matrix(rnorm(120), 40, 3)
  fit <- multivarious::cca(X, Y, ncomp = 3L, lambda = 0)
  a <- adapter_multivarious_cca(lambda = 0, ncomp = 3L)
  trunc <- a$truncate(fit, 2L)

  expect_equal(length(a$roots(trunc)), 2L)
  expect_equal(dim(a$scores(trunc, "X")), c(40L, 2L))
  expect_equal(dim(a$scores(trunc, "Y")), c(40L, 2L))
  expect_equal(dim(a$loadings(trunc, "X")), c(4L, 2L))
  expect_equal(dim(a$loadings(trunc, "Y")), c(3L, 2L))
})

test_that("refit supports regularized p > n CCA", {
  skip_if_not_installed("multivarious")
  set.seed(303)
  X <- matrix(rnorm(60 * 14), 60, 14)
  Y <- matrix(rnorm(60 * 11), 60, 11)
  a <- adapter_multivarious_cca(lambda_x = 1e-2, lambda_y = 5e-3, ncomp = 4L)
  fit <- a$refit(NULL, list(X = X, Y = Y))

  expect_true(all(is.finite(a$roots(fit))))
  expect_equal(dim(a$scores(fit, "X")), c(60L, 4L))
  expect_equal(dim(a$scores(fit, "Y")), c(60L, 4L))
})
