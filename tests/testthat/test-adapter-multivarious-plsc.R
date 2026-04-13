test_that("adapter_multivarious_plsc constructs cleanly", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  a <- adapter_multivarious_plsc()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "multivarious_plsc")
  expect_equal(a$shape_kinds, "cross")
  expect_equal(a$validity_level, "conditional")
})

test_that("adapter_multivarious_plsc declares ONLY covariance", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_plsc()
  expect_true(adapter_supports(a, "cross", "covariance",
                               "component_significance"))
  expect_false(adapter_supports(a, "cross", "correlation",
                                "component_significance"))
})

test_that("adapter_multivarious_plsc declares all four v1 targets for covariance", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_plsc()
  for (tgt in c("component_significance", "variable_stability",
                "score_stability", "subspace_stability")) {
    expect_true(adapter_supports(a, "cross", "covariance", tgt))
  }
})

test_that("adapter_multivarious_plsc does NOT claim variable_significance", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_plsc()
  expect_false(adapter_supports(a, "cross", "covariance",
                                "variable_significance"))
})

test_that("register_multivarious_plsc_adapter registers when multivarious installed", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  register_multivarious_plsc_adapter()
  expect_true("multivarious_plsc" %in% list_infer_adapters())
})

test_that("strict-dispatch has NO ambiguity (single relation)", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  register_multivarious_plsc_adapter()
  rec <- infer_recipe(geometry = "cross", adapter = "multivarious_plsc")
  expect_true(is_infer_recipe(rec))
  expect_equal(rec$shape$relation$kind, "covariance")
})

test_that("roots hook returns squared singular values", {
  skip_if_not_installed("multivarious")
  set.seed(201)
  X <- matrix(rnorm(200), 20, 10)
  Y <- matrix(rnorm(160), 20, 8)
  fit <- multivarious::plsc(X, Y, ncomp = 5)
  a <- adapter_multivarious_plsc()
  rt <- a$roots(fit)
  expect_equal(length(rt), 5L)
  expect_equal(rt, fit$singvals^2, tolerance = 1e-10)
})

test_that("loadings hook returns p x k for X and q x k for Y", {
  skip_if_not_installed("multivarious")
  set.seed(202)
  X <- matrix(rnorm(200), 20, 10)
  Y <- matrix(rnorm(160), 20, 8)
  fit <- multivarious::plsc(X, Y, ncomp = 5)
  a <- adapter_multivarious_plsc()
  expect_equal(dim(a$loadings(fit, "X")), c(10L, 5L))
  expect_equal(dim(a$loadings(fit, "Y")), c(8L, 5L))
})

test_that("scores hook returns n x k for both domains", {
  skip_if_not_installed("multivarious")
  set.seed(203)
  X <- matrix(rnorm(200), 20, 10)
  Y <- matrix(rnorm(160), 20, 8)
  fit <- multivarious::plsc(X, Y, ncomp = 5)
  a <- adapter_multivarious_plsc()
  expect_equal(dim(a$scores(fit, "X")), c(20L, 5L))
  expect_equal(dim(a$scores(fit, "Y")), c(20L, 5L))
})

test_that("refit re-fits with the requested ncomp", {
  skip_if_not_installed("multivarious")
  set.seed(204)
  X <- matrix(rnorm(200), 20, 10)
  Y <- matrix(rnorm(160), 20, 8)
  a <- adapter_multivarious_plsc(ncomp = 3L)
  fit <- a$refit(NULL, list(X = X, Y = Y))
  expect_equal(length(fit$singvals), 3L)
})

test_that("component_stat returns a finite number for k = 1", {
  skip_if_not_installed("multivarious")
  set.seed(205)
  X <- matrix(rnorm(200), 20, 10)
  Y <- matrix(rnorm(160), 20, 8)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  a <- adapter_multivarious_plsc()
  fit <- multivarious::plsc(Xc, Yc, ncomp = 5)
  v <- a$component_stat(fit, list(X = Xc, Y = Yc), 1L)
  expect_true(is.finite(v))
  expect_true(v >= 0 && v <= 1)
})

test_that("null_action preserves shape and breaks pairing", {
  skip_if_not_installed("multivarious")
  set.seed(206)
  X <- matrix(rnorm(200), 20, 10)
  Y <- matrix(rnorm(160), 20, 8)
  a <- adapter_multivarious_plsc()
  perturbed <- a$null_action(NULL, list(X = X, Y = Y))
  expect_equal(dim(perturbed$X), dim(X))
  expect_equal(dim(perturbed$Y), dim(Y))
  expect_false(isTRUE(all.equal(crossprod(X, Y),
                                crossprod(perturbed$X, perturbed$Y))))
})
