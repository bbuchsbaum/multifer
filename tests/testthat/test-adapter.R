# Clear the registry before every test in this file so tests are isolated.
clear_adapter_registry()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Build a valid capability matrix for oneblock/variance supporting all four
# v1 targets.
.pca_caps <- function() {
  capability_matrix(
    list(
      geometry = "oneblock",
      relation = "variance",
      targets  = c("component_significance", "variable_stability",
                   "score_stability", "subspace_stability")
    )
  )
}

# Minimal PCA stub adapter: provides all hooks needed for all four v1 targets.
.make_pca_adapter <- function(id = "pca_stub") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = .pca_caps(),
    # hooks
    roots          = function(x) x$values,
    scores         = function(x, domain) x$scores,
    loadings       = function(x, domain) x$loadings,
    truncate       = function(x, k) x,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    align          = function(xb, xref) xb,
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1.0,
    variable_stat  = function(x, data, domain, k) rep(1.0, 3L),
    score_stat     = function(x, data, domain, k) rep(1.0, 3L),
    validity_level = "conditional"
  )
}

.make_geneig_residualize <- function(mode = c("b_metric", "euclidean", "delegate")) {
  mode <- match.arg(mode)
  fn <- function(x, k, data) data
  if (mode == "b_metric") {
    attr(fn, "b_metric") <- TRUE
  } else if (mode == "delegate") {
    attr(fn, "delegates_geneig_deflation") <- TRUE
  }
  fn
}

.predictive_caps <- function() {
  capability_matrix(
    list(
      geometry = "cross",
      relation = "predictive",
      targets  = "component_significance"
    )
  )
}

# ---------------------------------------------------------------------------
# 1. Valid minimal PCA stub constructs cleanly
# ---------------------------------------------------------------------------

test_that("valid PCA stub adapter constructs without error", {
  a <- .make_pca_adapter()
  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "pca_stub")
  expect_equal(a$adapter_version, "0.0.1")
  expect_equal(a$shape_kinds, "oneblock")
  expect_equal(a$validity_level, "conditional")
  expect_true(is_capability_matrix(a$capabilities))
  expect_equal(nrow(a$capabilities), 4L)
})

test_that("PCA stub adapter has all expected hook fields", {
  a <- .make_pca_adapter()
  for (hook in c("roots", "scores", "loadings", "truncate", "residualize",
                 "refit", "align", "null_action", "component_stat",
                 "variable_stat", "score_stat")) {
    expect_true(is.function(a[[hook]]),
                info = paste("hook", hook, "should be a function"))
  }
})

test_that("optional perturbation/projection hooks are accepted", {
  a <- infer_adapter(
    adapter_id      = "adapter_owned_bootstrap",
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = c("variable_stability", "score_stability"))
    ),
    roots = function(x) x$values,
    scores = function(x, domain) x$scores,
    loadings = function(x, domain) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x, info = list(replicate = replicate))
    },
    project_scores = function(x, data, domain = NULL) {
      matrix(0, nrow = nrow(data), ncol = ncol(x$loadings))
    },
    validity_level = "conditional"
  )

  expect_true(is_infer_adapter(a))
  expect_true(is.function(a$bootstrap_action))
  expect_true(is.function(a$project_scores))
})

test_that("project_scores can satisfy score_stability score extraction", {
  a <- infer_adapter(
    adapter_id      = "project_scores_only",
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "score_stability")
    ),
    roots = function(x) x$values,
    loadings = function(x, domain) x$loadings,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x)
    },
    project_scores = function(x, data, domain = NULL) {
      matrix(0, nrow = nrow(data), ncol = ncol(x$loadings))
    },
    validity_level = "conditional"
  )

  expect_true(is_infer_adapter(a))
  expect_true(is.function(a$project_scores))
})

# ---------------------------------------------------------------------------
# 2. Under-provisioned adapters error at construction
# ---------------------------------------------------------------------------

test_that("claiming component_significance without null_action errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "component_significance")
      ),
      component_stat = function(x, data, k) 1.0,
      residualize    = function(x, k, data) data,
      # null_action deliberately omitted
      refit          = function(x, new_data) x,
      validity_level = "conditional"
    ),
    "null_action"
  )
})

test_that("claiming component_significance without component_stat errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "component_significance")
      ),
      null_action    = function(x, data) data,
      residualize    = function(x, k, data) data,
      # component_stat deliberately omitted
      refit          = function(x, new_data) x,
      validity_level = "conditional"
    ),
    "component_stat"
  )
})

test_that("claiming component_significance without residualize errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "component_significance")
      ),
      null_action    = function(x, data) data,
      component_stat = function(x, data, k) 1.0,
      # residualize deliberately omitted
      refit          = function(x, new_data) x,
      validity_level = "conditional"
    ),
    "residualize"
  )
})

test_that("geneig component_significance refuses unmarked Euclidean residualize hooks", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad_geneig",
      shape_kinds     = "geneig",
      capabilities    = capability_matrix(
        list(
          geometry = "geneig",
          relation = "generalized_eigen",
          targets  = "component_significance"
        )
      ),
      null_action    = function(x, data) data,
      component_stat = function(x, data, k) 1,
      residualize    = .make_geneig_residualize("euclidean"),
      refit          = function(x, new_data) x,
      validity_level = "conditional"
    ),
    "B-metric deflation required - see notes/engine_geneig_spec.md"
  )
})

test_that("geneig component_significance accepts a marked B-metric residualize hook", {
  expect_no_error(
    infer_adapter(
      adapter_id      = "good_geneig",
      shape_kinds     = "geneig",
      capabilities    = capability_matrix(
        list(
          geometry = "geneig",
          relation = "generalized_eigen",
          targets  = "component_significance"
        )
      ),
      null_action    = function(x, data) data,
      component_stat = function(x, data, k) 1,
      residualize    = .make_geneig_residualize("b_metric"),
      refit          = function(x, new_data) x,
      validity_level = "conditional"
    )
  )
})

test_that("geneig component_significance accepts delegation to engine deflation", {
  expect_no_error(
    infer_adapter(
      adapter_id      = "delegate_geneig",
      shape_kinds     = "geneig",
      capabilities    = capability_matrix(
        list(
          geometry = "geneig",
          relation = "generalized_eigen",
          targets  = "component_significance"
        )
      ),
      null_action    = function(x, data) data,
      component_stat = function(x, data, k) 1,
      residualize    = .make_geneig_residualize("delegate"),
      refit          = function(x, new_data) x,
      validity_level = "conditional"
    )
  )
})

test_that("adapter_lda_refit is accepted by the geneig registration-time gate", {
  skip_if_not_installed("MASS")
  expect_no_error(adapter_lda_refit())
})

test_that("predictive component_significance refuses in-sample component_stat hooks", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad_predictive",
      shape_kinds     = "cross",
      capabilities    = .predictive_caps(),
      null_action     = function(x, data) data,
      component_stat  = function(x, data, k) 1,
      residualize     = function(x, k, data) data,
      refit           = function(x, new_data) x,
      predict_response = function(x, new_data, k = NULL) {
        matrix(0, nrow = nrow(new_data$X), ncol = ncol(new_data$Y))
      },
      validity_level  = "conditional"
    ),
    "Predictive cross-fit admissibility rule"
  )
})

test_that("predictive component_significance accepts split-aware held-out hooks", {
  expect_no_error(
    infer_adapter(
      adapter_id      = "good_predictive",
      shape_kinds     = "cross",
      capabilities    = .predictive_caps(),
      null_action     = function(x, data) data,
      component_stat  = function(x, data, k, split = NULL) 1,
      residualize     = function(x, k, data) data,
      refit           = function(x, new_data) x,
      predict_response = function(x, new_data, k = NULL) {
        matrix(0, nrow = nrow(new_data$X), ncol = ncol(new_data$Y))
      },
      validity_level  = "conditional"
    )
  )
})

test_that("claiming variable_stability without perturbation path errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "variable_stability")
      ),
      loadings       = function(x, domain) matrix(0),
      # neither refit nor core/update_core provided
      validity_level = "heuristic"
    ),
    "perturbation hook"
  )
})

test_that("claiming variable_stability without variable_stat or loadings errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "variable_stability")
      ),
      refit          = function(x, new_data) x,
      # neither variable_stat nor loadings provided
      validity_level = "heuristic"
    ),
    "variable_stat.*loadings|loadings.*variable_stat|variable_stat` or `loadings"
  )
})

test_that("claiming score_stability without loadings hook errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "score_stability")
      ),
      refit          = function(x, new_data) x,
      scores         = function(x, domain = NULL) matrix(0, 1, 1),
      # loadings deliberately omitted
      validity_level = "heuristic"
    ),
    "loadings"
  )
})

test_that("claiming score_stability without score extraction hook errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "score_stability")
      ),
      refit          = function(x, new_data) x,
      loadings       = function(x, domain = NULL) matrix(0, 1, 1),
      # scores and project_scores deliberately omitted
      validity_level = "heuristic"
    ),
    "scores.*project_scores|project_scores.*scores"
  )
})

test_that("claiming subspace_stability without loadings hook errors", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "subspace_stability")
      ),
      refit          = function(x, new_data) x,
      # loadings deliberately omitted
      validity_level = "heuristic"
    ),
    "loadings"
  )
})

# ---------------------------------------------------------------------------
# 3. Claiming variable_significance in v1 errors with Part 5 section 38 ref
# ---------------------------------------------------------------------------

test_that("claiming variable_significance errors with Part 5 section 38 reference", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "variable_significance")
      ),
      refit          = function(x, new_data) x,
      validity_level = "heuristic"
    ),
    "section 38"
  )
})

# ---------------------------------------------------------------------------
# 4. Unknown validity_level errors
# ---------------------------------------------------------------------------

test_that("unknown validity_level errors", {
  expect_error(
    infer_adapter(
      adapter_id   = "bad",
      shape_kinds  = "oneblock",
      capabilities = capability_matrix(),
      validity_level = "unknown_level"
    ),
    "validity_level"
  )
})

# ---------------------------------------------------------------------------
# 5. shape_kinds not a subset of valid geometries errors
# ---------------------------------------------------------------------------

test_that("shape_kinds with unknown geometry errors", {
  expect_error(
    infer_adapter(
      adapter_id   = "bad",
      shape_kinds  = c("oneblock", "pca"),
      capabilities = capability_matrix(),
      validity_level = "heuristic"
    ),
    "unknown geometry"
  )
})

# ---------------------------------------------------------------------------
# 6. Capability referencing geometry not in shape_kinds errors
# ---------------------------------------------------------------------------

test_that("capability geometry not in shape_kinds errors", {
  expect_error(
    infer_adapter(
      adapter_id   = "bad",
      shape_kinds  = "oneblock",
      capabilities = capability_matrix(
        list(geometry = "cross", relation = "covariance",
             targets  = "component_significance")
      ),
      null_action    = function(x, data) data,
      component_stat = function(x, data, k) 1.0,
      residualize    = function(x, k, data) data,
      validity_level = "heuristic"
    ),
    "shape_kinds"
  )
})

# ---------------------------------------------------------------------------
# 7. Registry: register, get, list, unregister
# ---------------------------------------------------------------------------

test_that("registry: register and get roundtrip", {
  clear_adapter_registry()
  a <- .make_pca_adapter("reg_pca")
  register_infer_adapter("reg_pca", a)
  a2 <- get_infer_adapter("reg_pca")
  expect_true(is_infer_adapter(a2))
  expect_equal(a2$adapter_id, "reg_pca")
})

test_that("registry: list_infer_adapters returns registered ids", {
  clear_adapter_registry()
  register_infer_adapter("reg_a", .make_pca_adapter("reg_a"))
  register_infer_adapter("reg_b", .make_pca_adapter("reg_b"))
  ids <- list_infer_adapters()
  expect_true("reg_a" %in% ids)
  expect_true("reg_b" %in% ids)
  expect_equal(length(ids), 2L)
})

test_that("registry: duplicate registration errors without overwrite", {
  clear_adapter_registry()
  register_infer_adapter("dup", .make_pca_adapter("dup"))
  expect_error(
    register_infer_adapter("dup", .make_pca_adapter("dup")),
    "already registered"
  )
})

test_that("registry: overwrite = TRUE replaces silently", {
  clear_adapter_registry()
  register_infer_adapter("ow", .make_pca_adapter("ow"))
  a2 <- .make_pca_adapter("ow")
  expect_no_error(register_infer_adapter("ow", a2, overwrite = TRUE))
  expect_equal(get_infer_adapter("ow")$adapter_id, "ow")
})

test_that("registry: unregister removes the entry", {
  clear_adapter_registry()
  register_infer_adapter("to_rm", .make_pca_adapter("to_rm"))
  unregister_infer_adapter("to_rm")
  expect_false("to_rm" %in% list_infer_adapters())
})

test_that("registry: get unknown id errors", {
  clear_adapter_registry()
  expect_error(get_infer_adapter("does_not_exist"), "No adapter registered")
})

test_that("registry: unregister unknown id errors", {
  clear_adapter_registry()
  expect_error(unregister_infer_adapter("does_not_exist"), "No adapter registered")
})

# ---------------------------------------------------------------------------
# 8. is_infer_adapter and adapter_supports
# ---------------------------------------------------------------------------

test_that("is_infer_adapter distinguishes adapters from other objects", {
  a <- .make_pca_adapter()
  expect_true(is_infer_adapter(a))
  expect_false(is_infer_adapter(list()))
  expect_false(is_infer_adapter("string"))
  expect_false(is_infer_adapter(NULL))
})

test_that("adapter_supports returns TRUE for declared triples", {
  a <- .make_pca_adapter()
  expect_true(adapter_supports(a, "oneblock", "variance",
                               "component_significance"))
  expect_true(adapter_supports(a, "oneblock", "variance",
                               "variable_stability"))
  expect_true(adapter_supports(a, "oneblock", "variance",
                               "score_stability"))
  expect_true(adapter_supports(a, "oneblock", "variance",
                               "subspace_stability"))
})

test_that("adapter_supports returns FALSE for undeclared triples", {
  a <- .make_pca_adapter()
  expect_false(adapter_supports(a, "cross", "covariance",
                                "component_significance"))
  expect_false(adapter_supports(a, "oneblock", "covariance",
                                "component_significance"))
})

test_that("adapter_capabilities returns the capability matrix", {
  a <- .make_pca_adapter()
  cm <- adapter_capabilities(a)
  expect_true(is_capability_matrix(cm))
  expect_equal(nrow(cm), 4L)
})

# ---------------------------------------------------------------------------
# 9. print method
# ---------------------------------------------------------------------------

test_that("print.multifer_adapter produces informative output", {
  a <- .make_pca_adapter()
  out <- capture.output(print(a))
  expect_true(any(grepl("multifer_adapter", out)))
  expect_true(any(grepl("pca_stub", out)))
  expect_true(any(grepl("oneblock", out)))
  expect_true(any(grepl("conditional", out)))
  expect_true(any(grepl("4", out)))
})

# ---------------------------------------------------------------------------
# 10. Extra: fast-path (core + update_core) satisfies perturbation requirement
# ---------------------------------------------------------------------------

test_that("core+update_core satisfies perturbation path for variable_stability", {
  expect_no_error(
    infer_adapter(
      adapter_id      = "fast_pca",
      shape_kinds     = "oneblock",
      capabilities    = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets  = "variable_stability")
      ),
      core            = function(x, data) list(),
      update_core     = function(core_obj, ...) core_obj,
      loadings        = function(x, domain) matrix(0),
      validity_level  = "conditional"
    )
  )
})

# Clean up after all tests in this file.
clear_adapter_registry()
