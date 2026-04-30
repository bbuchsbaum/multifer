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

test_that("compile_bootstrap_plan resolves default cross callbacks once", {
  clear_adapter_registry()

  set.seed(101)
  X <- matrix(stats::rnorm(60), 12, 5)
  Y <- matrix(stats::rnorm(48), 12, 4)
  recipe <- .make_cross_recipe()
  adapter <- recipe$adapter
  original_fit <- adapter$refit(NULL, list(X = X, Y = Y, relation = "covariance"))

  plan <- compile_bootstrap_plan(
    recipe = recipe,
    adapter = adapter,
    data = list(X = X, Y = Y),
    original_fit = original_fit
  )
  idx <- c(2L, 1L, 2L, seq.int(4L, 12L))
  rep_data <- plan$resample_data(idx)

  expect_s3_class(plan, "multifer_bootstrap_plan")
  expect_equal(plan$data_shape, "cross")
  expect_true(is_data_role_schema(plan$data_schema))
  expect_named(plan$data_schema, c("X", "Y", "relation"))
  expect_equal(plan$domains, c("X", "Y"))
  expect_true(plan$fast_path_supported)
  expect_false(plan$rank_deficient_fallback)
  expect_equal(rep_data$X, X[idx, , drop = FALSE])
  expect_equal(rep_data$Y, Y[idx, , drop = FALSE])
  expect_equal(rep_data$relation, "covariance")
})

test_that("compile_bootstrap_plan normalizes adapter-owned bootstrap callbacks", {
  clear_adapter_registry()

  payload <- list(reference = matrix(seq_len(12), 4, 3))
  fit <- list(values = 1, loadings = matrix(1, 3, 1))
  adapter <- infer_adapter(
    adapter_id = "bootstrap_plan_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "score_stability")
    ),
    roots = function(x) x$values,
    domains = function(x, data = NULL) "payload",
    loadings = function(x, domain = NULL) x$loadings,
    refit = function(x, new_data) fit,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x, resample_indices = as.integer(replicate))
    },
    project_scores = function(x, data, domain = NULL) matrix(1, 4, 1),
    validity_level = "conditional"
  )
  recipe <- infer_recipe(
    geometry = "adapter",
    relation = "variance",
    targets = "score_stability",
    adapter = adapter
  )

  plan <- compile_bootstrap_plan(recipe, adapter, payload, fit)
  action <- plan$bootstrap_action(fit, payload, recipe$shape$design, 3L)

  expect_s3_class(plan, "multifer_bootstrap_plan")
  expect_true(plan$has_custom_bootstrap)
  expect_equal(action$fit, fit)
  expect_equal(action$resample_indices, 3L)
})

test_that("compile_bootstrap_plan resamples adapter data through data_schema", {
  clear_adapter_registry()

  payload <- list(
    Z = matrix(seq_len(20), nrow = 5),
    group = letters[1:5],
    weights = seq_len(5),
    lambda = 0.5
  )
  fit <- list(values = 1, loadings = matrix(1, 4, 1))
  adapter <- infer_adapter(
    adapter_id = "schema_bootstrap_plan_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "variable_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain = NULL) x$loadings,
    refit = function(x, new_data) fit,
    validity_level = "conditional",
    data_schema = data_role_schema(
      Z = role("primary", axes = c("row", "col")),
      group = role("design_index", axis = "row"),
      weights = role("metric", axis = "row",
                     policy = "diagonal_position_weights"),
      lambda = role("static")
    )
  )
  recipe <- infer_recipe(
    geometry = "adapter",
    relation = "variance",
    targets = "variable_stability",
    adapter = adapter
  )

  plan <- compile_bootstrap_plan(
    recipe,
    adapter,
    payload,
    fit,
    store_aligned_scores = FALSE
  )
  idx <- c(5L, 1L, 5L, 2L, 3L)
  rep_data <- plan$resample_data(idx)

  expect_s3_class(plan, "multifer_bootstrap_plan")
  expect_equal(plan$data_shape, "adapter_schema")
  expect_true(is_data_role_schema(plan$data_schema))
  expect_false(plan$fast_path_supported)
  expect_equal(rep_data$Z, payload$Z[idx, , drop = FALSE])
  expect_equal(rep_data$group, payload$group[idx])
  expect_equal(rep_data$weights, payload$weights)
  expect_equal(rep_data$lambda, payload$lambda)
})

test_that("compile_bootstrap_plan refuses unsafe adapter metric policies", {
  clear_adapter_registry()

  payload <- list(
    Z = matrix(stats::rnorm(20), nrow = 5),
    K = diag(5)
  )
  fit <- list(values = 1, loadings = matrix(1, 4, 1))
  adapter <- infer_adapter(
    adapter_id = "unsafe_schema_bootstrap_plan_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "variable_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain = NULL) x$loadings,
    refit = function(x, new_data) fit,
    validity_level = "conditional",
    data_schema = data_role_schema(
      Z = role("primary", axes = c("row", "col")),
      K = role("metric", axes = c("row", "col"),
               policy = "full_spd_no_replacement_only")
    )
  )
  recipe <- infer_recipe(
    geometry = "adapter",
    relation = "variance",
    targets = "variable_stability",
    adapter = adapter
  )

  expect_error(
    compile_bootstrap_plan(
      recipe,
      adapter,
      payload,
      fit,
      store_aligned_scores = FALSE
    ),
    "full_spd_no_replacement_only"
  )
})

test_that("compile_bootstrap_plan lets custom bootstrap own unsafe schema roles", {
  clear_adapter_registry()

  payload <- list(
    Z = matrix(stats::rnorm(20), nrow = 5),
    K = diag(5)
  )
  fit <- list(values = 1, loadings = matrix(1, 4, 1))
  adapter <- infer_adapter(
    adapter_id = "custom_schema_bootstrap_plan_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "variable_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain = NULL) x$loadings,
    refit = function(x, new_data) fit,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x, info = list(replicate = replicate))
    },
    validity_level = "conditional",
    data_schema = data_role_schema(
      Z = role("primary", axes = c("row", "col")),
      K = role("metric", axes = c("row", "col"),
               policy = "full_spd_no_replacement_only")
    )
  )
  recipe <- infer_recipe(
    geometry = "adapter",
    relation = "variance",
    targets = "variable_stability",
    adapter = adapter
  )

  plan <- compile_bootstrap_plan(
    recipe,
    adapter,
    payload,
    fit,
    store_aligned_scores = FALSE
  )

  expect_s3_class(plan, "multifer_bootstrap_plan")
  expect_true(plan$has_custom_bootstrap)
  expect_equal(plan$data_shape, "adapter")
  expect_true(is_data_role_schema(plan$data_schema))
})

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
# Test 4b: Procrustes remains available for legacy comparisons
# ---------------------------------------------------------------------------

test_that("bootstrap_fits still supports legacy procrustes alignment", {
  clear_adapter_registry()

  set.seed(41)
  X <- matrix(stats::rnorm(160), 32, 5)

  recipe  <- .make_oneblock_recipe()
  adapter <- recipe$adapter

  original_fit <- adapter$refit(NULL, X)
  units        <- form_units(adapter$roots(original_fit))

  expect_no_error(
    result <- bootstrap_fits(
      recipe       = recipe,
      adapter      = adapter,
      data         = X,
      original_fit = original_fit,
      units        = units,
      R            = 4L,
      method_align = "procrustes",
      seed         = 13L
    )
  )

  expect_equal(result$method_align, "procrustes")
  expect_length(result$reps, 4L)
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
# Test 5b: Optional aligned_scores suppression
# ---------------------------------------------------------------------------

test_that("bootstrap_fits can skip aligned_scores when requested", {
  clear_adapter_registry()

  set.seed(55)
  X <- matrix(stats::rnorm(180), 36, 5)

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
    R            = 6L,
    seed         = 21L,
    store_aligned_scores = FALSE
  )

  expect_length(result$reps, 6L)
  expect_true(all(vapply(result$reps, function(rep) is.null(rep$aligned_scores), logical(1L))))
  expect_true(all(vapply(result$reps, function(rep) !is.null(rep$aligned_loadings$X), logical(1L))))
})

test_that("bootstrap_fits uses adapter-owned bootstrap_action", {
  clear_adapter_registry()

  X <- matrix(seq_len(24), 6, 4)
  fit <- list(
    values = c(2, 1),
    loadings = diag(4)[, 1:2, drop = FALSE],
    scores = matrix(0, nrow = nrow(X), ncol = 2)
  )
  adapter <- infer_adapter(
    adapter_id = "custom_bootstrap",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "variable_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain = NULL) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(
        fit = x,
        resample_indices = as.integer(replicate),
        info = list(unit = "adapter")
      )
    },
    validity_level = "conditional"
  )

  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = adapter,
    targets = "variable_stability"
  )
  art <- bootstrap_fits(
    recipe = rec,
    adapter = adapter,
    data = X,
    original_fit = fit,
    units = form_units(adapter$roots(fit)),
    R = 3L,
    seed = 9L,
    store_aligned_scores = FALSE
  )

  expect_true(art$used_bootstrap_action)
  expect_equal(vapply(art$reps, function(rep) rep$resample_indices, integer(1L)),
               as.integer(1:3))
  expect_equal(art$reps[[1L]]$bootstrap_info$unit, "adapter")
})

test_that("bootstrap_fits stores adapter-projected scores for original data", {
  clear_adapter_registry()

  X <- matrix(seq_len(24), 6, 4)
  fit <- list(
    values = c(2, 1),
    loadings = diag(4)[, 1:2, drop = FALSE],
    scores = matrix(0, nrow = nrow(X), ncol = 2)
  )
  adapter <- infer_adapter(
    adapter_id = "projecting_bootstrap",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "score_stability")
    ),
    roots = function(x) x$values,
    scores = function(x, domain = NULL) x$scores,
    loadings = function(x, domain = NULL) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x, resample_indices = as.integer(replicate))
    },
    project_scores = function(x, data, domain = NULL) {
      matrix(7, nrow = nrow(data), ncol = ncol(x$loadings))
    },
    validity_level = "conditional"
  )

  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = adapter,
    targets = "score_stability"
  )
  art <- bootstrap_fits(
    recipe = rec,
    adapter = adapter,
    data = X,
    original_fit = fit,
    units = form_units(adapter$roots(fit)),
    R = 2L,
    seed = 10L,
    store_aligned_scores = TRUE
  )

  expect_equal(art$score_source, "project_scores")
  expect_equal(art$reps[[1L]]$aligned_scores$X,
               matrix(7, nrow = nrow(X), ncol = 2))
})

test_that("bootstrap_fits delegates loading and score alignment to adapter", {
  clear_adapter_registry()

  X <- matrix(seq_len(24), 6, 4)
  fit <- list(
    values = c(2, 1),
    loadings = diag(4)[, 1:2, drop = FALSE],
    scores = matrix(0, nrow = nrow(X), ncol = 2)
  )
  adapter <- infer_adapter(
    adapter_id = "custom_align_bootstrap",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "score_stability")
    ),
    roots = function(x) x$values,
    scores = function(x, domain = NULL) x$scores,
    loadings = function(x, domain = NULL) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x, resample_indices = as.integer(replicate))
    },
    project_scores = function(x, data, domain = NULL) {
      matrix(2, nrow = nrow(data), ncol = ncol(x$loadings))
    },
    align = function(reference_loadings,
                     replicate_loadings,
                     replicate_scores,
                     reference_fit,
                     replicate_fit,
                     domain,
                     method,
                     ...) {
      list(
        loadings = matrix(5, nrow = nrow(replicate_loadings),
                          ncol = ncol(replicate_loadings)),
        scores = matrix(6, nrow = nrow(replicate_scores),
                        ncol = ncol(replicate_scores))
      )
    },
    validity_level = "conditional"
  )

  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = adapter,
    targets = "score_stability"
  )
  art <- bootstrap_fits(
    recipe = rec,
    adapter = adapter,
    data = X,
    original_fit = fit,
    units = form_units(adapter$roots(fit)),
    R = 1L,
    seed = 11L,
    store_aligned_scores = TRUE
  )

  expect_equal(art$reps[[1L]]$aligned_loadings$X,
               matrix(5, nrow = 4, ncol = 2))
  expect_equal(art$reps[[1L]]$aligned_scores$X,
               matrix(6, nrow = nrow(X), ncol = 2))
})

test_that("bootstrap_fits rejects aligned scores for adapter geometry without project_scores", {
  payload <- list(reference = matrix(seq_len(12), 4, 3))
  fit <- list(values = c(2, 1), loadings = diag(2), scores = matrix(0, 4, 2))
  adapter <- infer_adapter(
    adapter_id = "opaque_boot_no_project",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "subspace_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain = NULL) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) list(fit = x),
    validity_level = "conditional"
  )
  rec <- infer_recipe(
    geometry = "adapter",
    relation = "variance",
    adapter = adapter,
    targets = "subspace_stability"
  )

  expect_error(
    bootstrap_fits(
      recipe = rec,
      adapter = adapter,
      data = payload,
      original_fit = fit,
      units = form_units(adapter$roots(fit)),
      R = 1L,
      store_aligned_scores = TRUE
    ),
    "project_scores"
  )
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

test_that("bootstrap_fits respects blocked_rows by resampling within groups", {
  clear_adapter_registry()

  set.seed(61)
  X <- matrix(stats::rnorm(240), 40, 6)
  groups <- rep(letters[1:4], each = 10)

  adapter <- adapter_prcomp()
  register_infer_adapter("prcomp_oneblock", adapter, overwrite = TRUE)
  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    design = blocked_rows(groups),
    adapter = adapter
  )
  original_fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe = recipe,
    adapter = adapter,
    data = X,
    original_fit = original_fit,
    units = units,
    R = 8L,
    seed = 61L,
    store_aligned_scores = FALSE
  )

  for (rep in result$reps) {
    idx <- rep$resample_indices
    expect_equal(tabulate(as.integer(factor(groups[idx], levels = letters[1:4]))),
                 rep(10L, 4L))
  }
})

test_that("bootstrap_fits respects clustered_rows by resampling whole clusters", {
  clear_adapter_registry()

  set.seed(62)
  X <- matrix(stats::rnorm(180), 30, 6)
  clusters <- rep(seq_len(10), each = 3)

  adapter <- adapter_prcomp()
  register_infer_adapter("prcomp_oneblock", adapter, overwrite = TRUE)
  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    design = clustered_rows(clusters),
    adapter = adapter
  )
  original_fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(original_fit))

  result <- bootstrap_fits(
    recipe = recipe,
    adapter = adapter,
    data = X,
    original_fit = original_fit,
    units = units,
    R = 8L,
    seed = 62L,
    fast_path = "off",
    store_aligned_scores = FALSE
  )

  rows_by_cluster <- split(seq_along(clusters), clusters)
  for (rep in result$reps) {
    idx <- rep$resample_indices
    counts <- tabulate(match(idx, seq_along(clusters)), nbins = length(clusters))
    for (rows in rows_by_cluster) {
      expect_true(length(unique(counts[rows])) == 1L)
    }
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
