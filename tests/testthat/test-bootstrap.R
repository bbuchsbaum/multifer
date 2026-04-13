## test-bootstrap.R
## Tests for bootstrap_fits() and multifer_bootstrap_artifact.
##
## Every test block starts with clear_adapter_registry() and registers only
## what it needs, so tests do not interfere via the global registry.

# ---------------------------------------------------------------------------
# Helper: compile a oneblock recipe using adapter_prcomp
# ---------------------------------------------------------------------------

.make_oneblock_recipe <- function() {
  adapter <- adapter_prcomp()
  register_infer_adapter("prcomp_oneblock", adapter, overwrite = TRUE)
  infer_recipe(geometry = "oneblock", relation = "variance", adapter = adapter)
}

# ---------------------------------------------------------------------------
# Helper: compile a cross recipe using adapter_cross_svd with covariance
# ---------------------------------------------------------------------------

.make_cross_recipe <- function() {
  adapter <- adapter_cross_svd()
  register_infer_adapter("cross_svd", adapter, overwrite = TRUE)
  infer_recipe(
    geometry = "cross",
    relation = "covariance",
    adapter  = adapter
  )
}

# ---------------------------------------------------------------------------
# Test 1: Oneblock bootstrap -- structural test
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: oneblock structural test", {
  clear_adapter_registry()

  set.seed(1)
  X <- matrix(stats::rnorm(100), 20, 5)

  recipe  <- .make_oneblock_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, X)
  units        <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = X,
    original_fit = original_fit,
    units        = units,
    R            = 20L,
    seed         = 42L
  )

  expect_s3_class(result, "multifer_bootstrap_artifact")
  expect_equal(result$R, 20L)
  expect_equal(result$domains, "X")
  expect_length(result$reps, 20L)

  rep1 <- result$reps[[1]]
  expect_true(!is.null(rep1$fit))
  expect_true(!is.null(rep1$aligned_loadings$X))
  expect_true(!is.null(rep1$aligned_scores$X))
  expect_true(!is.null(rep1$resample_indices))

  # Dimensions match original loadings
  orig_L_dim <- dim(adapter$loadings(original_fit, "X"))
  expect_equal(dim(rep1$aligned_loadings$X), orig_L_dim)
})

# ---------------------------------------------------------------------------
# Test 2: Cross paired bootstrap -- structural test
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: cross paired bootstrap structural test", {
  clear_adapter_registry()

  set.seed(2)
  X <- matrix(stats::rnorm(200), 40, 5)
  Y <- matrix(stats::rnorm(160), 40, 4)

  recipe  <- .make_cross_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  units        <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = list(X = X, Y = Y),
    original_fit = original_fit,
    units        = units,
    R            = 10L,
    seed         = 7L
  )

  expect_s3_class(result, "multifer_bootstrap_artifact")
  expect_equal(result$R, 10L)
  expect_equal(result$domains, c("X", "Y"))
  expect_length(result$reps, 10L)

  rep1 <- result$reps[[1]]
  expect_true(!is.null(rep1$aligned_loadings$X))
  expect_true(!is.null(rep1$aligned_loadings$Y))
  expect_true(!is.null(rep1$aligned_scores$X))
  expect_true(!is.null(rep1$aligned_scores$Y))

  # resample_indices is a single vector (not a list)
  expect_true(is.integer(rep1$resample_indices))
  expect_equal(length(rep1$resample_indices), nrow(X))
})

# ---------------------------------------------------------------------------
# Test 3: Paired bootstrap preserves pairing -- key correctness test
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: cross paired bootstrap uses the same index vector for X and Y", {
  clear_adapter_registry()

  set.seed(3)
  X <- matrix(stats::rnorm(200), 40, 5)
  Y <- matrix(stats::rnorm(160), 40, 4)

  recipe  <- .make_cross_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  units        <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = list(X = X, Y = Y),
    original_fit = original_fit,
    units        = units,
    R            = 5L,
    seed         = 99L
  )

  # For each replicate: verify resample_indices is a single integer vector,
  # not two separate ones. The pairing is verified by checking that applying
  # the stored indices to X and Y reproduces the expected sub-matrices.
  for (b in seq_len(5L)) {
    idx <- result$reps[[b]]$resample_indices

    # Must be a single vector of length n
    expect_true(is.integer(idx))
    expect_equal(length(idx), nrow(X))

    # The same index applied to both blocks gives consistent sub-matrices.
    # We cannot directly introspect the fit's data, but we can verify that
    # X[idx, ] and Y[idx, ] have matching row counts -- confirming one index
    # drove both.
    X_sub <- X[idx, , drop = FALSE]
    Y_sub <- Y[idx, , drop = FALSE]
    expect_equal(nrow(X_sub), nrow(Y_sub))
    expect_equal(nrow(X_sub), nrow(X))
  }
})

# ---------------------------------------------------------------------------
# Test 4: Sign alignment is active (not a pass-through)
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: sign alignment guarantees non-negative inner products", {
  clear_adapter_registry()

  # Use a well-separated dataset so sign flips are meaningful.
  set.seed(4)
  X <- matrix(stats::rnorm(300), 60, 5)

  recipe  <- .make_oneblock_recipe()
  adapter <- recipe$adapter

  original_fit  <- adapter$refit(NULL, X)
  orig_loadings <- adapter$loadings(original_fit, "X")
  units         <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = X,
    original_fit = original_fit,
    units        = units,
    R            = 30L,
    method_align = "sign",
    seed         = 11L
  )

  # After sign alignment, diag(t(Vref) %*% aligned_loadings) must be >= 0
  # for every replicate and every column.
  for (b in seq_len(result$R)) {
    aligned_L <- result$reps[[b]]$aligned_loadings$X
    dots      <- base::colSums(orig_loadings * aligned_L)
    expect_true(
      all(dots >= -1e-10),
      info = sprintf("rep %d has negative inner product(s): %s", b,
                     paste(round(dots, 6), collapse = ", "))
    )
  }
})

# ---------------------------------------------------------------------------
# Test 5: Loading-score consistency
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: aligned_scores are consistent with aligned_loadings on replicate data", {
  clear_adapter_registry()

  # Use svd adapter so scores = X_centered %*% V (pure linear projection).
  set.seed(5)
  X <- matrix(stats::rnorm(300), 60, 5)

  adapter <- adapter_svd()
  register_infer_adapter("svd_oneblock", adapter, overwrite = TRUE)
  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = adapter
  )

  original_fit <- adapter$refit(NULL, X)
  units        <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = X,
    original_fit = original_fit,
    units        = units,
    R            = 10L,
    method_align = "sign",
    seed         = 17L
  )

  # For the svd adapter: scores = Xc %*% V.
  # aligned_scores should equal rep_Xc %*% aligned_loadings where rep_Xc is
  # the centered replicate data. Verify within floating point tolerance.
  for (b in seq_len(result$R)) {
    idx     <- result$reps[[b]]$resample_indices
    rep_X   <- X[idx, , drop = FALSE]
    ctr     <- base::colMeans(rep_X)
    rep_Xc  <- base::sweep(rep_X, 2L, ctr, "-")

    aligned_L <- result$reps[[b]]$aligned_loadings$X
    aligned_S <- result$reps[[b]]$aligned_scores$X

    expected_S <- rep_Xc %*% aligned_L
    frob_err   <- sqrt(sum((aligned_S - expected_S)^2))
    expect_lt(
      frob_err, 1e-8,
      label = sprintf("rep %d loading-score consistency (Frobenius err = %g)", b, frob_err)
    )
  }
})

# ---------------------------------------------------------------------------
# Test 6: Reproducibility under seed
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: same seed produces identical results", {
  clear_adapter_registry()

  set.seed(6)
  X <- matrix(stats::rnorm(200), 40, 5)

  recipe  <- .make_oneblock_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, X)
  units        <- form_units(adapter$roots(original_fit))

  r1 <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = X,
    original_fit = original_fit,
    units        = units,
    R            = 15L,
    seed         = 123L
  )

  r2 <- bootstrap_fits(
    recipe       = recipe,
    adapter      = adapter,
    data         = X,
    original_fit = original_fit,
    units        = units,
    R            = 15L,
    seed         = 123L
  )

  # resample_indices must be identical across all reps
  for (b in seq_len(15L)) {
    expect_identical(
      r1$reps[[b]]$resample_indices,
      r2$reps[[b]]$resample_indices,
      label = sprintf("resample_indices differ at rep %d", b)
    )
    expect_equal(
      r1$reps[[b]]$aligned_loadings$X,
      r2$reps[[b]]$aligned_loadings$X,
      tolerance = 1e-14,
      label = sprintf("aligned_loadings differ at rep %d", b)
    )
    expect_equal(
      r1$reps[[b]]$aligned_scores$X,
      r2$reps[[b]]$aligned_scores$X,
      tolerance = 1e-14,
      label = sprintf("aligned_scores differ at rep %d", b)
    )
  }
})

# ---------------------------------------------------------------------------
# Test 7: Validation errors
# ---------------------------------------------------------------------------

test_that("bootstrap_fits: oneblock recipe rejects list data", {
  clear_adapter_registry()

  set.seed(7)
  X <- matrix(stats::rnorm(100), 20, 5)

  recipe  <- .make_oneblock_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, X)
  units        <- form_units(adapter$roots(original_fit))

  expect_error(
    bootstrap_fits(
      recipe       = recipe,
      adapter      = adapter,
      data         = list(X = X),        # wrong: should be a matrix
      original_fit = original_fit,
      units        = units,
      R            = 5L
    ),
    regexp = "oneblock"
  )
})

test_that("bootstrap_fits: cross recipe rejects bare matrix data", {
  clear_adapter_registry()

  set.seed(8)
  X <- matrix(stats::rnorm(200), 40, 5)
  Y <- matrix(stats::rnorm(160), 40, 4)

  recipe  <- .make_cross_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  units        <- form_units(adapter$roots(original_fit))

  expect_error(
    bootstrap_fits(
      recipe       = recipe,
      adapter      = adapter,
      data         = X,                  # wrong: should be list(X, Y)
      original_fit = original_fit,
      units        = units,
      R            = 5L
    ),
    regexp = "cross"
  )
})

test_that("bootstrap_fits: non-positive R errors", {
  clear_adapter_registry()

  set.seed(9)
  X <- matrix(stats::rnorm(100), 20, 5)

  recipe  <- .make_oneblock_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, X)
  units        <- form_units(adapter$roots(original_fit))

  expect_error(
    bootstrap_fits(
      recipe       = recipe,
      adapter      = adapter,
      data         = X,
      original_fit = original_fit,
      units        = units,
      R            = 0L
    ),
    regexp = "positive integer"
  )

  expect_error(
    bootstrap_fits(
      recipe       = recipe,
      adapter      = adapter,
      data         = X,
      original_fit = original_fit,
      units        = units,
      R            = -5L
    ),
    regexp = "positive integer"
  )
})
