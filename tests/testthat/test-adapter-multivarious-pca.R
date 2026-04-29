test_that("adapter_multivarious_pca constructs cleanly", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  a <- adapter_multivarious_pca()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "multivarious_pca")
  expect_equal(a$shape_kinds, "oneblock")
  expect_equal(a$validity_level, "conditional")
})

test_that("adapter_multivarious_pca declares all four v1 targets", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_pca()
  for (tgt in c("component_significance", "variable_stability",
                "score_stability", "subspace_stability")) {
    expect_true(adapter_supports(a, "oneblock", "variance", tgt))
  }
})

test_that("adapter_multivarious_pca does NOT claim variable_significance", {
  skip_if_not_installed("multivarious")
  a <- adapter_multivarious_pca()
  expect_false(adapter_supports(a, "oneblock", "variance",
                                "variable_significance"))
})

test_that("register_multivarious_pca_adapter registers when multivarious installed", {
  skip_if_not_installed("multivarious")
  clear_adapter_registry()
  register_multivarious_pca_adapter()
  expect_true("multivarious_pca" %in% list_infer_adapters())
})

test_that("roots hook returns sdev squared", {
  skip_if_not_installed("multivarious")
  set.seed(101)
  X <- matrix(rnorm(200), 20, 10)
  fit <- multivarious::pca(X, ncomp = 5)
  a <- adapter_multivarious_pca()
  rt <- a$roots(fit)
  expect_equal(length(rt), 5L)
  expect_equal(rt, multivarious::sdev(fit)^2, tolerance = 1e-10)
})

test_that("loadings hook returns p x k matrix", {
  skip_if_not_installed("multivarious")
  set.seed(102)
  X <- matrix(rnorm(200), 20, 10)
  fit <- multivarious::pca(X, ncomp = 5)
  a <- adapter_multivarious_pca()
  V <- a$loadings(fit)
  expect_equal(dim(V), c(10L, 5L))
})

test_that("scores hook returns n x k matrix", {
  skip_if_not_installed("multivarious")
  set.seed(103)
  X <- matrix(rnorm(200), 20, 10)
  fit <- multivarious::pca(X, ncomp = 5)
  a <- adapter_multivarious_pca()
  S <- a$scores(fit)
  expect_equal(dim(S), c(20L, 5L))
})

test_that("refit hook re-fits with the requested ncomp", {
  skip_if_not_installed("multivarious")
  set.seed(104)
  X <- matrix(rnorm(200), 20, 10)
  a <- adapter_multivarious_pca(ncomp = 3L)
  fit <- a$refit(NULL, X)
  expect_equal(length(multivarious::sdev(fit)), 3L)
})

test_that("component_stat returns a finite number on real data", {
  skip_if_not_installed("multivarious")
  set.seed(105)
  X <- matrix(rnorm(200), 20, 10)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  a <- adapter_multivarious_pca()
  fit <- multivarious::pca(X, ncomp = 5)
  v <- a$component_stat(fit, Xc, 1L)
  expect_true(is.finite(v))
  expect_true(v >= 0 && v <= 1)
})

test_that("null_action preserves dimensions and changes column covariance", {
  skip_if_not_installed("multivarious")
  set.seed(106)
  X <- matrix(rnorm(200), 20, 10) %*% matrix(rnorm(100), 10, 10)
  a <- adapter_multivarious_pca()
  Xp <- a$null_action(NULL, X)
  expect_equal(dim(Xp), dim(X))
  expect_false(isTRUE(all.equal(cov(X), cov(Xp))))
})

test_that("residualize removes the top component", {
  skip_if_not_installed("multivarious")
  set.seed(107)
  X <- matrix(rnorm(200), 20, 10)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  fit <- multivarious::pca(Xc, ncomp = 5)
  a <- adapter_multivarious_pca()
  resid <- a$residualize(fit, 1L, Xc)
  s_orig <- svd(Xc)$d[1]
  s_res  <- svd(resid)$d[1]
  expect_true(s_res < s_orig)
})

test_that("native PCA feature evidence action returns conditional bootstrap ratios", {
  skip_if_not_installed("multivarious")
  set.seed(108)
  X <- matrix(rnorm(160), 20, 8)
  a <- adapter_multivarious_pca(ncomp = 4L)
  fit <- a$refit(NULL, X)
  units <- form_units(a$roots(fit))

  out <- feature_evidence_from_adapter(
    adapter = a,
    fit = fit,
    data = X,
    units = units,
    statistic = "loading",
    orientation = "auto",
    R = 3L,
    seed = 108L
  )

  expect_s3_class(out, "multifer_feature_evidence")
  expect_true(nrow(out) > 0L)
  expect_equal(unique(out$method), "conditional_subspace_bootstrap")
  expect_equal(unique(out$calibration), "bootstrap")
  expect_true(all(is.na(out$p_value)))
  expect_true(all(out$validity_level == "conditional"))

  units_extra <- form_units(c(a$roots(fit), 0))
  out_extra <- feature_evidence_from_adapter(
    adapter = a,
    fit = fit,
    data = X,
    units = units_extra,
    statistic = "loading",
    orientation = "auto",
    R = 3L,
    seed = 109L
  )
  expect_false("u5" %in% out_extra$unit_id)
})
