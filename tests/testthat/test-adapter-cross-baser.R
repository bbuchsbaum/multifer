# Tests for adapter_cross_svd() and adapter_cancor() -- structural and smoke
# tests only. No calls to infer() (engine not yet available).

# ---------------------------------------------------------------------------
# 1. adapter_cross_svd() constructs and declares BOTH relations
# ---------------------------------------------------------------------------

test_that("adapter_cross_svd constructs and declares both covariance and correlation", {
  clear_adapter_registry()
  a <- adapter_cross_svd()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "cross_svd")
  expect_equal(a$shape_kinds, "cross")
  expect_true(is_capability_matrix(a$capabilities))
  expect_true(adapter_supports(a, "cross", "covariance", "component_significance"))
  expect_true(adapter_supports(a, "cross", "correlation", "component_significance"))
})

# ---------------------------------------------------------------------------
# 2. adapter_cross_svd() does NOT claim variable_significance
# ---------------------------------------------------------------------------

test_that("adapter_cross_svd does not claim variable_significance", {
  clear_adapter_registry()
  a <- adapter_cross_svd()
  expect_false(adapter_supports(a, "cross", "covariance", "variable_significance"))
  expect_false(adapter_supports(a, "cross", "correlation", "variable_significance"))
})

# ---------------------------------------------------------------------------
# 3. adapter_cancor() constructs and declares ONLY correlation
# ---------------------------------------------------------------------------

test_that("adapter_cancor constructs and declares only correlation", {
  clear_adapter_registry()
  a <- adapter_cancor()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "cancor_cross")
  expect_equal(a$shape_kinds, "cross")
  expect_true(adapter_supports(a, "cross", "correlation", "component_significance"))
  expect_false(adapter_supports(a, "cross", "covariance", "component_significance"))
})

# ---------------------------------------------------------------------------
# 4. Both adapters register via register_cross_baser_adapters()
# ---------------------------------------------------------------------------

test_that("register_cross_baser_adapters registers both adapters", {
  clear_adapter_registry()
  register_cross_baser_adapters()
  ids <- list_infer_adapters()
  expect_true("cross_svd" %in% ids)
  expect_true("cancor_cross" %in% ids)
  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 5. Dual-relation adapter triggers strict-dispatch ambiguity error
# ---------------------------------------------------------------------------

test_that("adapter_cross_svd triggers strict-dispatch ambiguity error without relation", {
  clear_adapter_registry()
  register_infer_adapter("cross_svd", adapter_cross_svd())

  # Error must mention "covariance"
  expect_error(
    infer_recipe(geometry = "cross", adapter = "cross_svd"),
    regexp = "covariance"
  )

  # Same call must also mention "correlation"
  expect_error(
    infer_recipe(geometry = "cross", adapter = "cross_svd"),
    regexp = "correlation"
  )

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 6. Single-relation adapter_cancor() has no ambiguity without relation arg
# ---------------------------------------------------------------------------

test_that("adapter_cancor has no ambiguity under strict dispatch", {
  clear_adapter_registry()
  register_infer_adapter("cancor_cross", adapter_cancor())

  expect_no_error(
    infer_recipe(geometry = "cross", adapter = "cancor_cross")
  )

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 7. Explicit-relation path compiles cleanly for BOTH relations on cross_svd
# ---------------------------------------------------------------------------

test_that("infer_recipe compiles with explicit covariance on cross_svd", {
  clear_adapter_registry()
  register_infer_adapter("cross_svd", adapter_cross_svd())

  recipe <- infer_recipe(geometry = "cross", relation = "covariance",
                         adapter = "cross_svd")
  expect_true(is_infer_recipe(recipe))
  expect_equal(recipe$shape$relation$kind, "covariance")

  clear_adapter_registry()
})

test_that("infer_recipe compiles with explicit correlation on cross_svd", {
  clear_adapter_registry()
  register_infer_adapter("cross_svd", adapter_cross_svd())

  recipe <- infer_recipe(geometry = "cross", relation = "correlation",
                         adapter = "cross_svd")
  expect_true(is_infer_recipe(recipe))
  expect_equal(recipe$shape$relation$kind, "correlation")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 8. refit smoke test for cross_svd (covariance and correlation)
# ---------------------------------------------------------------------------

test_that("adapter_cross_svd refit returns correct slots for covariance", {
  clear_adapter_registry()
  set.seed(42L)
  X <- matrix(rnorm(50L), 10L, 5L)
  Y <- matrix(rnorm(40L), 10L, 4L)
  a <- adapter_cross_svd()

  fit <- a$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  expect_true(is.list(fit))
  expect_true(!is.null(fit$relation))
  expect_true(!is.null(fit$d))
  expect_true(!is.null(fit$Wx))
  expect_true(!is.null(fit$Wy))
  expect_true(!is.null(fit$Tx))
  expect_true(!is.null(fit$Ty))
  expect_equal(fit$relation, "covariance")
  # d should have min(p, q) = 4 singular values
  expect_equal(length(fit$d), 4L)
  # Wx: p x r, Wy: q x r
  expect_equal(nrow(fit$Wx), 5L)
  expect_equal(nrow(fit$Wy), 4L)
  # Tx: n x r, Ty: n x r
  expect_equal(nrow(fit$Tx), 10L)
  expect_equal(nrow(fit$Ty), 10L)
})

test_that("adapter_cross_svd refit returns correct slots for correlation", {
  clear_adapter_registry()
  set.seed(42L)
  X <- matrix(rnorm(50L), 10L, 5L)
  Y <- matrix(rnorm(40L), 10L, 4L)
  a <- adapter_cross_svd()

  fit <- a$refit(NULL, list(X = X, Y = Y, relation = "correlation"))
  expect_true(is.list(fit))
  expect_equal(fit$relation, "correlation")
  expect_true(!is.null(fit$d))
  expect_true(!is.null(fit$Wx))
  expect_true(!is.null(fit$Wy))
  expect_true(!is.null(fit$Tx))
  expect_true(!is.null(fit$Ty))
  # canonical correlations are in [0, 1]
  expect_true(all(fit$d >= -1e-10 & fit$d <= 1 + 1e-10))
})

# ---------------------------------------------------------------------------
# 9. refit smoke test for cancor
# ---------------------------------------------------------------------------

test_that("adapter_cancor refit returns correct slots", {
  clear_adapter_registry()
  set.seed(42L)
  X <- matrix(rnorm(50L), 10L, 5L)
  Y <- matrix(rnorm(40L), 10L, 4L)
  a <- adapter_cancor()

  fit <- a$refit(NULL, list(X = X, Y = Y))
  expect_true(is.list(fit))
  expect_true(!is.null(fit$cor))
  expect_true(!is.null(fit$xcoef))
  expect_true(!is.null(fit$ycoef))
  expect_true(!is.null(fit$xcoef_scores))
  expect_true(!is.null(fit$ycoef_scores))
  # cor: min(p, q) = 4 canonical correlations
  expect_equal(length(fit$cor), 4L)
  # canonical correlations in [0, 1]
  expect_true(all(fit$cor >= -1e-10 & fit$cor <= 1 + 1e-10))
  # score matrices: n x r
  expect_equal(nrow(fit$xcoef_scores), 10L)
  expect_equal(nrow(fit$ycoef_scores), 10L)
})

# ---------------------------------------------------------------------------
# 10. component_stat returns a finite number for both adapters
# ---------------------------------------------------------------------------

test_that("adapter_cross_svd component_stat returns finite value in [0, 1] for covariance", {
  clear_adapter_registry()
  set.seed(42L)
  X <- matrix(rnorm(50L), 10L, 5L)
  Y <- matrix(rnorm(40L), 10L, 4L)
  a   <- adapter_cross_svd()
  fit <- a$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  val <- a$component_stat(x = fit, data = list(X = X, Y = Y), k = 1L)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1L)
  expect_true(is.finite(val))
  expect_true(val >= 0 && val <= 1)
})

test_that("adapter_cross_svd component_stat returns finite value in [0, 1] for correlation", {
  clear_adapter_registry()
  set.seed(42L)
  X <- matrix(rnorm(50L), 10L, 5L)
  Y <- matrix(rnorm(40L), 10L, 4L)
  a   <- adapter_cross_svd()
  fit <- a$refit(NULL, list(X = X, Y = Y, relation = "correlation"))
  val <- a$component_stat(x = fit, data = list(X = X, Y = Y), k = 1L)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1L)
  expect_true(is.finite(val))
  expect_true(val >= 0 && val <= 1)
})

test_that("adapter_cancor component_stat returns finite value in [0, 1]", {
  clear_adapter_registry()
  set.seed(42L)
  X <- matrix(rnorm(50L), 10L, 5L)
  Y <- matrix(rnorm(40L), 10L, 4L)
  a   <- adapter_cancor()
  fit <- a$refit(NULL, list(X = X, Y = Y))
  val <- a$component_stat(x = fit, data = list(X = X, Y = Y), k = 1L)
  expect_true(is.numeric(val))
  expect_equal(length(val), 1L)
  expect_true(is.finite(val))
  expect_true(val >= 0 && val <= 1)
})

# ---------------------------------------------------------------------------
# 11. null_action preserves dimensions and breaks pairing
# ---------------------------------------------------------------------------

test_that("adapter_cross_svd null_action preserves dim and breaks pairing", {
  clear_adapter_registry()
  set.seed(42L)
  # Construct X, Y with a shared latent signal so crossprod differs under perm.
  signal <- rnorm(10L)
  X <- cbind(signal + rnorm(10L, sd = 0.1), rnorm(10L), rnorm(10L),
             rnorm(10L), rnorm(10L))
  Y <- cbind(signal + rnorm(10L, sd = 0.1), rnorm(10L), rnorm(10L),
             rnorm(10L))
  a   <- adapter_cross_svd()
  fit <- a$refit(NULL, list(X = X, Y = Y, relation = "covariance"))

  set.seed(1L)
  null_data <- a$null_action(x = fit, data = list(X = X, Y = Y))
  expect_equal(dim(null_data$Y), dim(Y))
  expect_equal(dim(null_data$X), dim(X))

  # Cross-product should differ after permutation
  cp_orig <- crossprod(X, Y)
  cp_perm <- crossprod(null_data$X, null_data$Y)
  expect_false(isTRUE(all.equal(cp_orig, cp_perm, tolerance = 1e-6)))
})

test_that("adapter_cancor null_action preserves dim and breaks pairing", {
  clear_adapter_registry()
  set.seed(42L)
  signal <- rnorm(10L)
  X <- cbind(signal + rnorm(10L, sd = 0.1), rnorm(10L), rnorm(10L),
             rnorm(10L), rnorm(10L))
  Y <- cbind(signal + rnorm(10L, sd = 0.1), rnorm(10L), rnorm(10L),
             rnorm(10L))
  a   <- adapter_cancor()
  fit <- a$refit(NULL, list(X = X, Y = Y))

  set.seed(1L)
  null_data <- a$null_action(x = fit, data = list(X = X, Y = Y))
  expect_equal(dim(null_data$Y), dim(Y))
  expect_equal(dim(null_data$X), dim(X))

  cp_orig <- crossprod(X, Y)
  cp_perm <- crossprod(null_data$X, null_data$Y)
  expect_false(isTRUE(all.equal(cp_orig, cp_perm, tolerance = 1e-6)))
})

# ---------------------------------------------------------------------------
# 12. Under-provisioned cross adapter error fires (dropping null_action)
# ---------------------------------------------------------------------------

test_that("dropping null_action when claiming component_significance errors", {
  clear_adapter_registry()
  expect_error(
    infer_adapter(
      adapter_id      = "bad_cross",
      adapter_version = "0.0.1",
      shape_kinds     = "cross",
      capabilities    = capability_matrix(
        list(geometry = "cross", relation = "covariance",
             targets  = c("component_significance", "variable_stability",
                          "score_stability", "subspace_stability"))
      ),
      roots          = function(x, ...) x$d,
      scores         = function(x, domain = c("X", "Y"), ...) x$Tx,
      loadings       = function(x, domain = c("X", "Y"), ...) x$Wx,
      truncate       = function(x, k, ...) x,
      residualize    = function(x, k, data, ...) data,
      refit          = function(x, new_data, ...) new_data,
      # null_action deliberately omitted
      component_stat = function(x, data, k, ...) 1.0,
      validity_level = "conditional"
    ),
    "component_significance"
  )
})

# Clean up after all tests in this file.
clear_adapter_registry()
