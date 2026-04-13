# Tests for adapter_svd() and adapter_prcomp() -- structural and smoke tests only.
# No calls to infer() (owned by kxn.10). Engine round-trips belong to kxn.15.

# ---------------------------------------------------------------------------
# 1. adapter_svd() constructs cleanly
# ---------------------------------------------------------------------------

test_that("adapter_svd constructs a valid multifer_adapter", {
  clear_adapter_registry()
  a <- adapter_svd()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "svd_oneblock")
  expect_equal(a$shape_kinds, "oneblock")
  expect_true(is_capability_matrix(a$capabilities))
})

# ---------------------------------------------------------------------------
# 2. adapter_svd() declares all 4 v1 targets for (oneblock, variance)
# ---------------------------------------------------------------------------

test_that("adapter_svd supports all four v1 targets", {
  clear_adapter_registry()
  a <- adapter_svd()
  expect_true(adapter_supports(a, "oneblock", "variance", "component_significance"))
  expect_true(adapter_supports(a, "oneblock", "variance", "variable_stability"))
  expect_true(adapter_supports(a, "oneblock", "variance", "score_stability"))
  expect_true(adapter_supports(a, "oneblock", "variance", "subspace_stability"))
})

# ---------------------------------------------------------------------------
# 3. adapter_svd() does NOT claim variable_significance
# ---------------------------------------------------------------------------

test_that("adapter_svd does not claim variable_significance", {
  clear_adapter_registry()
  a <- adapter_svd()
  expect_false(adapter_supports(a, "oneblock", "variance", "variable_significance"))
})

# ---------------------------------------------------------------------------
# 4. adapter_prcomp() constructs cleanly
# ---------------------------------------------------------------------------

test_that("adapter_prcomp constructs a valid multifer_adapter", {
  clear_adapter_registry()
  a <- adapter_prcomp()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "prcomp_oneblock")
  expect_equal(a$shape_kinds, "oneblock")
  expect_true(is_capability_matrix(a$capabilities))
})

test_that("adapter_prcomp supports all four v1 targets", {
  clear_adapter_registry()
  a <- adapter_prcomp()
  expect_true(adapter_supports(a, "oneblock", "variance", "component_significance"))
  expect_true(adapter_supports(a, "oneblock", "variance", "variable_stability"))
  expect_true(adapter_supports(a, "oneblock", "variance", "score_stability"))
  expect_true(adapter_supports(a, "oneblock", "variance", "subspace_stability"))
})

test_that("adapter_prcomp does not claim variable_significance", {
  clear_adapter_registry()
  a <- adapter_prcomp()
  expect_false(adapter_supports(a, "oneblock", "variance", "variable_significance"))
})

# ---------------------------------------------------------------------------
# 5. Both adapters register via register_oneblock_baser_adapters()
# ---------------------------------------------------------------------------

test_that("register_oneblock_baser_adapters registers both adapters", {
  clear_adapter_registry()
  register_oneblock_baser_adapters()
  ids <- list_infer_adapters()
  expect_true("svd_oneblock" %in% ids)
  expect_true("prcomp_oneblock" %in% ids)
})

# ---------------------------------------------------------------------------
# 6. Fit smoke test for adapter_svd()
# ---------------------------------------------------------------------------

test_that("adapter_svd refit returns a list with u, d, v, center", {
  clear_adapter_registry()
  set.seed(42L)
  X   <- matrix(rnorm(50L), 10L, 5L)
  a   <- adapter_svd()
  fit <- a$refit(NULL, X)
  expect_true(is.list(fit))
  expect_true(!is.null(fit$u))
  expect_true(!is.null(fit$d))
  expect_true(!is.null(fit$v))
  expect_true(!is.null(fit$center))
  # u: n x min(n,p), d: min(n,p), v: p x min(n,p)
  expect_equal(nrow(fit$u), 10L)
  expect_equal(nrow(fit$v), 5L)
  expect_equal(length(fit$d), 5L)
  expect_equal(length(fit$center), 5L)
})

# ---------------------------------------------------------------------------
# 7. Fit smoke test for adapter_prcomp()
# ---------------------------------------------------------------------------

test_that("adapter_prcomp refit returns a prcomp-like object", {
  clear_adapter_registry()
  set.seed(42L)
  X   <- matrix(rnorm(50L), 10L, 5L)
  a   <- adapter_prcomp()
  fit <- a$refit(NULL, X)
  expect_true(is.list(fit))
  expect_true(!is.null(fit$sdev))
  expect_true(!is.null(fit$rotation))
  expect_true(!is.null(fit$x))
  expect_equal(length(fit$sdev), 5L)
  expect_equal(nrow(fit$rotation), 5L)
  expect_equal(nrow(fit$x), 10L)
})

# ---------------------------------------------------------------------------
# 8. component_stat returns a finite number in [0, 1]
# ---------------------------------------------------------------------------

test_that("adapter_svd component_stat returns finite value in [0, 1]", {
  clear_adapter_registry()
  set.seed(42L)
  X   <- matrix(rnorm(50L), 10L, 5L)
  a   <- adapter_svd()
  fit <- a$refit(NULL, X)
  val <- a$component_stat(x = fit, data = X, k = 1L)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1L)
  expect_true(is.finite(val))
  expect_true(val >= 0 && val <= 1)
})

test_that("adapter_prcomp component_stat returns finite value in [0, 1]", {
  clear_adapter_registry()
  set.seed(42L)
  X   <- matrix(rnorm(50L), 10L, 5L)
  a   <- adapter_prcomp()
  fit <- a$refit(NULL, X)
  val <- a$component_stat(x = fit, data = X, k = 1L)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1L)
  expect_true(is.finite(val))
  expect_true(val >= 0 && val <= 1)
})

# ---------------------------------------------------------------------------
# 9. residualize removes the first component (deflation reduces top root)
# ---------------------------------------------------------------------------

test_that("adapter_svd residualize deflates top component", {
  clear_adapter_registry()
  set.seed(42L)
  X    <- matrix(rnorm(50L), 10L, 5L)
  a    <- adapter_svd()
  fit  <- a$refit(NULL, X)

  resid <- a$residualize(x = fit, k = 1L, data = X)
  expect_equal(dim(resid), dim(X))

  fit2 <- a$refit(NULL, resid)

  # After deflating k=1, the absolute first root (eigenvalue) must shrink.
  # component_stat is a relative ratio and can increase after deflation, so
  # we compare the raw roots returned by the roots hook instead.
  root_before <- a$roots(fit)[1L]
  root_after  <- a$roots(fit2)[1L]
  expect_true(root_after < root_before)
})

test_that("adapter_prcomp residualize deflates top component", {
  clear_adapter_registry()
  set.seed(42L)
  X    <- matrix(rnorm(50L), 10L, 5L)
  a    <- adapter_prcomp()
  fit  <- a$refit(NULL, X)

  resid <- a$residualize(x = fit, k = 1L, data = X)
  expect_equal(dim(resid), dim(X))

  fit2 <- a$refit(NULL, resid)

  # After deflating k=1, the absolute first root (variance) must shrink.
  root_before <- a$roots(fit)[1L]
  root_after  <- a$roots(fit2)[1L]
  expect_true(root_after < root_before)
})

# ---------------------------------------------------------------------------
# 10. null_action preserves shape and destroys column covariance
# ---------------------------------------------------------------------------

test_that("adapter_svd null_action preserves dim and changes covariance", {
  clear_adapter_registry()
  set.seed(42L)
  # Structured matrix with clear covariance signal.
  signal  <- rnorm(10L)
  X       <- cbind(signal + rnorm(10L, sd = 0.1),
                   signal + rnorm(10L, sd = 0.1),
                   rnorm(10L),
                   rnorm(10L),
                   rnorm(10L))
  a       <- adapter_svd()
  fit     <- a$refit(NULL, X)
  Xperm   <- a$null_action(x = fit, data = X)

  expect_equal(dim(Xperm), dim(X))
  # Original cov between col1 and col2 should be clearly positive;
  # permuted columns should destroy it.
  cov_orig <- cov(X[, 1L], X[, 2L])
  cov_perm <- cov(Xperm[, 1L], Xperm[, 2L])
  expect_false(isTRUE(all.equal(cov_orig, cov_perm, tolerance = 0.2)))
})

test_that("adapter_prcomp null_action preserves dim and changes covariance", {
  clear_adapter_registry()
  set.seed(42L)
  signal  <- rnorm(10L)
  X       <- cbind(signal + rnorm(10L, sd = 0.1),
                   signal + rnorm(10L, sd = 0.1),
                   rnorm(10L),
                   rnorm(10L),
                   rnorm(10L))
  a       <- adapter_prcomp()
  fit     <- a$refit(NULL, X)
  Xperm   <- a$null_action(x = fit, data = X)

  expect_equal(dim(Xperm), dim(X))
  cov_orig <- cov(X[, 1L], X[, 2L])
  cov_perm <- cov(Xperm[, 1L], Xperm[, 2L])
  expect_false(isTRUE(all.equal(cov_orig, cov_perm, tolerance = 0.2)))
})

# ---------------------------------------------------------------------------
# 11. Registration-time validation still fires for under-provisioned adapters
# ---------------------------------------------------------------------------

test_that("dropping null_action when claiming component_significance errors", {
  clear_adapter_registry()
  expect_error(
    infer_adapter(
      adapter_id      = "bad_svd",
      adapter_version = "0.0.1",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = c("component_significance", "variable_stability",
                          "score_stability", "subspace_stability"))
      ),
      roots          = function(x, ...) x$d^2,
      scores         = function(x, domain = NULL, ...) x$u %*% diag(x$d, nrow = length(x$d)),
      loadings       = function(x, domain = NULL, ...) x$v,
      truncate       = function(x, k, ...) x,
      residualize    = function(x, k, data, ...) data,
      refit          = function(x, new_data, ...) svd(new_data),
      # null_action deliberately omitted
      component_stat = function(x, data, k, ...) 1.0,
      validity_level = "conditional"
    ),
    "component_significance"
  )
})

# Clean up after all tests in this file.
clear_adapter_registry()
