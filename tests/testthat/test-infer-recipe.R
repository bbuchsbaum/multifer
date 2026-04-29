# Clear the registry before every test in this file so tests are isolated.
clear_adapter_registry()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Minimal oneblock/variance adapter supporting component_significance only.
.make_pca_stub <- function(id = "stub_pca") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets  = "component_significance")
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1.0,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    validity_level = "conditional"
  )
}

# Oneblock adapter supporting all four v1 targets.
.make_pca_full <- function(id = "stub_pca_full") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1.0,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    loadings       = function(x, domain) numeric(0),
    scores         = function(x, domain) numeric(0),
    variable_stat  = function(x, data, domain, k) numeric(0),
    score_stat     = function(x, data, domain, k) numeric(0),
    validity_level = "conditional"
  )
}

# Cross adapter supporting BOTH (cross, covariance) and (cross, correlation),
# each with component_significance.
.make_cross_stub <- function(id = "stub_cross") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(geometry = "cross", relation = "covariance",
           targets  = "component_significance"),
      list(geometry = "cross", relation = "correlation",
           targets  = "component_significance")
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1.0,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    validity_level = "exact"
  )
}

.make_predictive_stub <- function(id = "stub_predictive") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(geometry = "cross", relation = "predictive",
           targets  = "component_significance")
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k, split = NULL) 1.0,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    predict_response = function(x, new_data, k = NULL) {
      if (is.null(new_data$Y)) {
        return(matrix(0, nrow = nrow(new_data$X), ncol = 1L))
      }
      matrix(0, nrow = nrow(new_data$Y), ncol = ncol(new_data$Y))
    },
    validity_level = "conditional"
  )
}

.make_geneig_stub <- function(id = "stub_geneig") {
  residualize_geneig <- function(x, k, data) data
  attr(residualize_geneig, "b_metric") <- TRUE

  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "geneig",
    capabilities    = capability_matrix(
      list(geometry = "geneig", relation = "generalized_eigen",
           targets  = "component_significance")
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1.0,
    residualize    = residualize_geneig,
    refit          = function(x, new_data) x,
    validity_level = "conditional"
  )
}

# ---------------------------------------------------------------------------
# 1. Happy path: oneblock PCA-style adapter
# ---------------------------------------------------------------------------

test_that("happy path oneblock: basic recipe compiles correctly", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  r <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca"
  )

  expect_true(is_infer_recipe(r))
  expect_true(is_infer_problem(r$problem))
  expect_equal(r$shape$geometry$kind, "oneblock")
  expect_equal(r$shape$relation$kind, "variance")
  expect_equal(r$adapter_id, "stub_pca")
  expect_true(is_infer_adapter(r$adapter))
  expect_equal(r$validity_level, "conditional")
  expect_true(r$strict)
  expect_true(is.na(r$downgrade_reason))
  expect_equal(r$problem$target_family, "latent_root")
  expect_equal(r$problem$engine_kind, "oneblock")
})

test_that("targets = 'default' expands to supported targets excluding variable_significance", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca_full", .make_pca_full())

  r <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca_full",
    targets  = "default"
  )

  expect_true("component_significance" %in% r$targets)
  expect_true("variable_stability"     %in% r$targets)
  expect_true("score_stability"        %in% r$targets)
  expect_true("subspace_stability"     %in% r$targets)
  expect_false("variable_significance" %in% r$targets)
})

# ---------------------------------------------------------------------------
# 2. Happy path cross -- explicit relation wins
# ---------------------------------------------------------------------------

test_that("cross adapter: explicit covariance relation compiles correctly", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r <- infer_recipe(
    geometry = "cross",
    relation = "covariance",
    adapter  = "stub_cross"
  )

  expect_true(is_infer_recipe(r))
  expect_equal(r$shape$relation$kind, "covariance")
})

test_that("cross adapter: explicit correlation relation compiles correctly", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    adapter  = "stub_cross"
  )

  expect_true(is_infer_recipe(r))
  expect_equal(r$shape$relation$kind, "correlation")
})

test_that("cross adapter: covariance and correlation recipes differ", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r_cov  <- infer_recipe(geometry = "cross", relation = "covariance",
                          adapter = "stub_cross")
  r_corr <- infer_recipe(geometry = "cross", relation = "correlation",
                          adapter = "stub_cross")

  expect_false(identical(r_cov$shape$relation$kind, r_corr$shape$relation$kind))
})

test_that("cross predictive adapter compiles when it declares the triple", {
  clear_adapter_registry()
  register_infer_adapter("stub_predictive", .make_predictive_stub())

  r <- infer_recipe(
    geometry = "cross",
    relation = "predictive",
    adapter  = "stub_predictive"
  )

  expect_true(is_infer_recipe(r))
  expect_equal(r$shape$geometry$kind, "cross")
  expect_equal(r$shape$relation$kind, "predictive")
  expect_true(is_infer_problem(r$problem))
  expect_equal(r$problem$target_family, "predictive_gain")
  expect_equal(r$problem$engine_kind, "predictive")
})

test_that("predictive relation is rejected for non-cross geometry at recipe compile time", {
  clear_adapter_registry()
  register_infer_adapter("stub_predictive", .make_predictive_stub())

  expect_error(
    infer_recipe(
      geometry = "oneblock",
      relation = "predictive",
      adapter  = "stub_predictive"
    ),
    regexp = "requires geometry 'cross'"
  )
})

test_that("predictive relation errors cleanly when adapter does not support the triple", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  expect_error(
    infer_recipe(
      geometry = "cross",
      relation = "predictive",
      adapter  = "stub_cross"
    ),
    regexp = "declares no supported targets|not supported"
  )
})

# ---------------------------------------------------------------------------
# 3. Strict mode ambiguity error
# ---------------------------------------------------------------------------

test_that("strict mode errors when relation is NULL and multiple are available", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  expect_error(
    infer_recipe(geometry = "cross", adapter = "stub_cross"),
    regexp = "covariance"   # error message must list at least one relation
  )
})

test_that("strict mode ambiguity error message lists both relations", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  err <- tryCatch(
    infer_recipe(geometry = "cross", adapter = "stub_cross"),
    error = function(e) conditionMessage(e)
  )
  expect_true(grepl("covariance", err))
  expect_true(grepl("correlation", err))
})

# ---------------------------------------------------------------------------
# 4. Non-strict mode: ambiguity warning + downgrade
# ---------------------------------------------------------------------------

test_that("non-strict mode emits warning and returns recipe on ambiguous relation", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  expect_warning(
    r <- infer_recipe(geometry = "cross", adapter = "stub_cross",
                      strict = FALSE),
    regexp = "Ambiguous relation"
  )
  expect_true(is_infer_recipe(r))
})

test_that("non-strict mode downgrades validity_level from exact", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r <- suppressWarnings(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = FALSE)
  )

  # adapter validity is "exact"; must be downgraded one step to "conditional"
  expect_equal(r$validity_level, "conditional")
})

test_that("non-strict mode populates downgrade_reason", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r <- suppressWarnings(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = FALSE)
  )

  expect_false(is.na(r$downgrade_reason))
  expect_true(nchar(r$downgrade_reason) > 0L)
})

# ---------------------------------------------------------------------------
# 5. Target not supported errors (regardless of strict)
# ---------------------------------------------------------------------------

test_that("unsupported target errors in strict mode", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  expect_error(
    infer_recipe(
      geometry = "oneblock",
      relation = "variance",
      targets  = c("component_significance", "variable_stability"),
      adapter  = "stub_pca"
    ),
    regexp = "variable_stability"
  )
})

test_that("unsupported target errors even with strict = FALSE", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  expect_error(
    infer_recipe(
      geometry = "oneblock",
      relation = "variance",
      targets  = c("component_significance", "variable_stability"),
      adapter  = "stub_pca",
      strict   = FALSE
    ),
    regexp = "variable_stability"
  )
})

# ---------------------------------------------------------------------------
# 6. variable_significance blocked even with strict = FALSE
# ---------------------------------------------------------------------------

test_that("variable_significance in targets errors with section 38 reference", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  expect_error(
    infer_recipe(
      geometry = "oneblock",
      relation = "variance",
      targets  = "variable_significance",
      adapter  = "stub_pca",
      strict   = FALSE
    ),
    regexp = "section 38"
  )
})

test_that("variable_significance errors even in strict mode", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  expect_error(
    infer_recipe(
      geometry = "oneblock",
      relation = "variance",
      targets  = c("component_significance", "variable_significance"),
      adapter  = "stub_pca",
      strict   = TRUE
    ),
    regexp = "section 38"
  )
})

# ---------------------------------------------------------------------------
# 7. Geometry not in adapter's shape_kinds
# ---------------------------------------------------------------------------

test_that("geometry not in adapter shape_kinds errors", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  expect_error(
    infer_recipe(
      geometry = "cross",
      relation = "covariance",
      adapter  = "stub_pca"
    ),
    regexp = "shape_kinds"
  )
})

# ---------------------------------------------------------------------------
# 8. Adapter lookup by string id
# ---------------------------------------------------------------------------

test_that("adapter can be specified as a string id and is looked up", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  r <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca"      # string -- not the object
  )
  expect_true(is_infer_recipe(r))
  expect_equal(r$adapter_id, "stub_pca")
})

# ---------------------------------------------------------------------------
# 9. Adapter lookup fails for unknown id
# ---------------------------------------------------------------------------

test_that("adapter string id not in registry errors clearly", {
  clear_adapter_registry()

  expect_error(
    infer_recipe(
      geometry = "oneblock",
      relation = "variance",
      adapter  = "does_not_exist"
    ),
    regexp = "does_not_exist"
  )
})

# ---------------------------------------------------------------------------
# 10. shape argument wins over component args
# ---------------------------------------------------------------------------

test_that("shape argument is used and individual components are ignored", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  sh <- typed_shape(
    geometry = geometry("oneblock"),
    relation = relation("variance"),
    design   = exchangeable_rows()
  )

  # Pass a conflicting geometry arg; shape should win
  r <- infer_recipe(
    shape    = sh,
    geometry = "cross",   # would fail on its own
    adapter  = .make_pca_stub()
  )

  expect_equal(r$shape$geometry$kind, "oneblock")
  expect_equal(r$shape$relation$kind, "variance")
})

# ---------------------------------------------------------------------------
# 11. Default design picks
# ---------------------------------------------------------------------------

test_that("NULL design for oneblock geometry defaults to exchangeable_rows", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  r <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    design   = NULL,
    adapter  = "stub_pca"
  )

  expect_equal(r$shape$design$kind, "exchangeable_rows")
})

test_that("NULL design for cross geometry defaults to paired_rows", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r <- infer_recipe(
    geometry = "cross",
    relation = "covariance",
    design   = NULL,
    adapter  = "stub_cross"
  )

  expect_equal(r$shape$design$kind, "paired_rows")
})

test_that("NULL design for multiblock geometry defaults to exchangeable_rows", {
  clear_adapter_registry()
  mb_adapter <- infer_adapter(
    adapter_id      = "stub_mb",
    adapter_version = "0.0.1",
    shape_kinds     = "multiblock",
    capabilities    = capability_matrix(
      list(geometry = "multiblock", relation = "covariance",
           targets  = "component_significance")
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1.0,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    validity_level = "heuristic"
  )

  r <- infer_recipe(
    geometry = "multiblock",
    relation = "covariance",
    design   = NULL,
    adapter  = mb_adapter
  )

  expect_equal(r$shape$design$kind, "exchangeable_rows")
})

test_that("geneig geometry defaults to exchangeable_rows and compiles", {
  clear_adapter_registry()
  register_infer_adapter("stub_geneig", .make_geneig_stub())

  r <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    design   = NULL,
    adapter  = "stub_geneig"
  )

  expect_true(is_infer_recipe(r))
  expect_equal(r$shape$design$kind, "exchangeable_rows")
  expect_equal(r$shape$relation$kind, "generalized_eigen")
})

test_that("adapter-owned geometry compiles with exchangeable_rows default", {
  a <- infer_adapter(
    adapter_id = "opaque_recipe_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "score_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain = NULL) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) list(fit = x),
    project_scores = function(x, data, domain = NULL) matrix(0, 2, 1),
    validity_level = "conditional"
  )

  rec <- infer_recipe(
    geometry = "adapter",
    relation = "variance",
    adapter = a,
    targets = "score_stability"
  )

  expect_true(is_infer_recipe(rec))
  expect_equal(rec$shape$geometry$kind, "adapter")
  expect_equal(rec$shape$design$kind, "exchangeable_rows")
})

test_that("infer_recipe refuses geneig with any non-generalized_eigen relation", {
  clear_adapter_registry()
  register_infer_adapter("stub_geneig", .make_geneig_stub())

  expect_error(
    infer_recipe(
      geometry = "geneig",
      relation = "variance",
      adapter  = "stub_geneig"
    ),
    regexp = "requires relation 'generalized_eigen'"
  )
})

# ---------------------------------------------------------------------------
# 12. Print method
# ---------------------------------------------------------------------------

test_that("print method runs without error and shows key fields", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  r <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca"
  )

  out <- capture.output(print(r))
  expect_true(any(grepl("multifer_infer_recipe", out)))
  expect_true(any(grepl("oneblock", out)))
  expect_true(any(grepl("variance", out)))
  expect_true(any(grepl("target_family", out)))
  expect_true(any(grepl("engine", out)))
  expect_true(any(grepl("stub_pca", out)))
  expect_true(any(grepl("component_significance", out)))
  expect_true(any(grepl("conditional", out)))
})

test_that("compiled infer_problem is printable and inspectable directly", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", .make_pca_stub())

  r <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca"
  )

  out <- capture.output(print(r$problem))
  expect_true(any(grepl("multifer_infer_problem", out)))
  expect_true(any(grepl("latent_root", out)))
  expect_true(any(grepl("oneblock", out)))
})

test_that("print method shows downgrade reason when strict = FALSE triggered", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross", .make_cross_stub())

  r <- suppressWarnings(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = FALSE)
  )

  out <- capture.output(print(r))
  expect_true(any(grepl("downgrade", out, ignore.case = TRUE)))
})

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------
clear_adapter_registry()
