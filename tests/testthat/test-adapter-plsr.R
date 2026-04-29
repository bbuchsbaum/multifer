make_plsr_fixture <- function(n = 48L, p = 6L, q = 2L, seed = 1L) {
  set.seed(seed)
  t1 <- scale(rnorm(n), center = TRUE, scale = FALSE)[, 1L]
  p_load <- c(1.3, -0.8, 0.6, 0, 0, 0)
  q_load <- c(1.2, -0.9)

  X <- matrix(rnorm(n * p, sd = 0.15), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * q, sd = 0.10), nrow = n, ncol = q)

  X <- X + t1 %o% p_load
  Y <- Y + t1 %o% q_load
  list(X = X, Y = Y)
}

test_that("adapter_plsr_refit constructs with predictive capabilities", {
  skip_if_not_installed("pls")

  a <- adapter_plsr_refit()

  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "plsr_refit")
  expect_equal(a$shape_kinds, "cross")
  expect_true(adapter_supports(a, "cross", "predictive", "component_significance"))
  expect_true(adapter_supports(a, "cross", "predictive", "variable_stability"))
  expect_true(adapter_supports(a, "cross", "predictive", "score_stability"))
  expect_true(adapter_supports(a, "cross", "predictive", "subspace_stability"))
  expect_null(a$core)
  expect_null(a$update_core)
  expect_true(all(c("paired_rows", "continuous_response") %in% a$declared_assumptions))
  expect_true(is.function(a$predict_response))
})

test_that("adapter_plsr_refit refit and predictive hooks round-trip", {
  skip_if_not_installed("pls")

  a <- adapter_plsr_refit()
  dat <- make_plsr_fixture(seed = 11L)
  fit <- a$refit(NULL, dat)

  expect_s3_class(fit, "mvr")
  expect_true(is.numeric(a$roots(fit)))
  expect_true(ncol(a$scores(fit, "X")) >= 1L)
  expect_true(ncol(a$loadings(fit, "X")) >= 1L)

  pred <- a$predict_response(
    fit,
    list(
      X = dat$X[1:7, , drop = FALSE],
      Y = dat$Y[1:7, , drop = FALSE]
    ),
    k = 1L
  )
  expect_true(is.matrix(pred))
  expect_identical(dim(pred), c(7L, ncol(dat$Y)))

  resid <- a$residualize(fit, 1L, dat)
  expect_true(is.list(resid))
  expect_identical(dim(resid$X), dim(dat$X))
  expect_identical(dim(resid$Y), dim(dat$Y))

  stat <- a$component_stat(
    fit,
    dat,
    k = 1L,
    split = list(train = 1:36, test = 37:48)
  )
  expect_true(is.numeric(stat))
  expect_length(stat, 1L)
  expect_gte(stat, 0)
})

test_that("plsr_refit null_action permutes Y relative to fixed X", {
  skip_if_not_installed("pls")

  a <- adapter_plsr_refit()
  dat <- make_plsr_fixture(seed = 13L)
  fit <- a$refit(NULL, dat)
  set.seed(19L)
  nul <- a$null_action(fit, dat)

  expect_equal(nul$X, dat$X)
  expect_false(identical(nul$Y, dat$Y))
  expect_equal(sort(nul$Y[, 1L]), sort(dat$Y[, 1L]))
  expect_equal(sort(nul$Y[, 2L]), sort(dat$Y[, 2L]))
})

test_that("plsr_refit component_stat is explicitly split-aware", {
  skip_if_not_installed("pls")

  a <- adapter_plsr_refit()
  dat <- make_plsr_fixture(seed = 14L)
  fit <- a$refit(NULL, dat)

  expect_true("split" %in% names(formals(a$component_stat)))
  expect_error(
    a$component_stat(fit, dat, k = 1L),
    "requires a split"
  )
})

test_that("register_plsr_refit_adapter registers the default id", {
  skip_if_not_installed("pls")
  clear_adapter_registry()

  register_plsr_refit_adapter()

  expect_true("plsr_refit" %in% list_infer_adapters())
  expect_true(is_infer_adapter(get_infer_adapter("plsr_refit")))
})

test_that("adapter_plsr_refit rejects unknown method via match.arg()", {
  skip_if_not_installed("pls")

  expect_error(
    adapter_plsr_refit(method = "kernelpls_typo"),
    "should be one of"
  )
  expect_error(
    adapter_plsr_refit(method = "model.frame"),
    "should be one of"
  )
})

for (m in c("simpls", "oscorespls", "kernelpls", "widekernelpls")) {
  local({
    method_name <- m
    test_that(sprintf("adapter_plsr_refit accepts method '%s'", method_name), {
      skip_if_not_installed("pls")
      dat <- make_plsr_fixture(seed = 17L)
      a <- adapter_plsr_refit(method = method_name)
      fit <- a$refit(NULL, dat)
      expect_s3_class(fit, "mvr")
      expect_true(is.numeric(a$roots(fit)))
      expect_gte(ncol(a$loadings(fit, "X")), 1L)
    })
  })
}
