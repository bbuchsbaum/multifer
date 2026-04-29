make_fe_artifact <- function(loadings_by_rep, domains = names(loadings_by_rep[[1L]])) {
  structure(
    list(
      reps = lapply(loadings_by_rep, function(L) list(aligned_loadings = L)),
      R = length(loadings_by_rep),
      method_align = "sign",
      domains = domains,
      seed = NULL
    ),
    class = "multifer_bootstrap_artifact"
  )
}

test_that("infer_feature_evidence constructs and validates schema", {
  empty <- infer_feature_evidence()
  expect_s3_class(empty, "multifer_feature_evidence")
  expect_equal(nrow(empty), 0L)

  out <- infer_feature_evidence(
    scope = "unit",
    domain = "X",
    feature = "x1",
    unit_id = "u1",
    unit_set = "u1",
    statistic = "loading",
    estimate = 1,
    bootstrap_mean = 1.1,
    bias = 0.1,
    std_error = 0.2,
    lower = 0.7,
    upper = 1.3,
    z = 5,
    z_label = "estimate / bootstrap_se",
    p_value = NA_real_,
    p_adjusted = NA_real_,
    adjustment = "none",
    calibration = "bootstrap",
    interval_method = "percentile",
    null_label = NA_character_,
    method = "aligned_bootstrap_loading",
    orientation = "signed",
    identifiable = TRUE,
    validity_level = "conditional",
    warnings = NA_character_
  )
  expect_s3_class(out, "multifer_feature_evidence")
  expect_equal(out$feature, "x1")

  expect_error(
    infer_feature_evidence(scope = c("unit", "unit"), domain = "X"),
    "same length"
  )
  expect_error(
    infer_feature_evidence(
      scope = "unit", domain = "X", feature = "x1", unit_id = "u1",
      unit_set = "u1", statistic = "loading", estimate = 1,
      bootstrap_mean = 1, bias = 0, std_error = 1, lower = 0, upper = 2,
      z = 1, z_label = "z", p_value = NA_real_, p_adjusted = NA_real_,
      adjustment = "none", calibration = "guess", interval_method = "none",
      null_label = NA_character_, method = "m", orientation = "signed",
      identifiable = TRUE, validity_level = "conditional",
      warnings = NA_character_
    ),
    "unknown calibration"
  )
})

test_that("feature_evidence_from_bootstrap returns signed loading ratios", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = TRUE
  )
  L1 <- matrix(c(1, 2), ncol = 1, dimnames = list(c("a", "b"), "c1"))
  L2 <- matrix(c(3, 4), ncol = 1, dimnames = list(c("a", "b"), "c1"))
  art <- make_fe_artifact(list(list(X = L1), list(X = L2)))

  out <- feature_evidence_from_bootstrap(
    art, units, statistic = "loading", orientation = "auto",
    normalize = "none"
  )

  expect_s3_class(out, "multifer_feature_evidence")
  expect_equal(out$scope, c("unit", "unit"))
  expect_equal(out$feature, c("a", "b"))
  expect_equal(out$estimate, c(2, 3))
  expect_equal(out$bootstrap_mean, c(2, 3))
  expect_equal(out$std_error, c(sqrt(2), sqrt(2)))
  expect_equal(out$orientation, c("signed", "signed"))
  expect_true(all(is.na(out$p_value)))
  expect_equal(unique(out$calibration), "bootstrap")
})

test_that("feature_evidence_from_bootstrap guards signed non-identifiable units", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "subspace",
    members = list(c(1L, 2L)),
    identifiable = FALSE,
    selected = TRUE
  )
  L <- matrix(c(3, 4, 0, 0), nrow = 2, byrow = TRUE)
  rownames(L) <- c("a", "b")
  art <- make_fe_artifact(list(list(X = L), list(X = L)))

  expect_error(
    feature_evidence_from_bootstrap(
      art, units, statistic = "loading", orientation = "signed"
    ),
    "not defined"
  )

  out <- feature_evidence_from_bootstrap(
    art, units, statistic = "loading", orientation = "auto",
    normalize = "none"
  )
  expect_equal(out$statistic, c("subspace_norm", "subspace_norm"))
  expect_equal(out$estimate, c(5, 0))
  expect_equal(out$orientation, c("subspace_norm", "subspace_norm"))
})

test_that("squared and subspace feature evidence are sign and rotation invariant", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "subspace",
    members = list(c(1L, 2L)),
    identifiable = FALSE,
    selected = TRUE
  )
  L <- matrix(c(3, 4, 1, 2), nrow = 2, byrow = TRUE)
  rownames(L) <- c("a", "b")
  theta <- pi / 5
  Q <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  art <- make_fe_artifact(list(list(X = L), list(X = L %*% Q)))

  sq <- feature_evidence_from_bootstrap(
    art, units, statistic = "squared_loading", normalize = "none"
  )
  sn <- feature_evidence_from_bootstrap(
    art, units, statistic = "subspace_norm", normalize = "none"
  )

  expect_equal(sq$estimate, c(25, 5), tolerance = 1e-12)
  expect_equal(sn$estimate, c(5, sqrt(5)), tolerance = 1e-12)
  expect_equal(sq$std_error, c(0, 0), tolerance = 1e-12)
  expect_equal(sn$std_error, c(0, 0), tolerance = 1e-12)
})

test_that("aggregate feature evidence matches feature_importance_from_bootstrap", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )
  L1 <- matrix(c(1, 0, 0, 1), nrow = 2)
  L2 <- matrix(c(0, 1, 1, 0), nrow = 2)
  art <- make_fe_artifact(list(list(X = L1), list(X = L2)))

  ev <- feature_evidence_from_bootstrap(
    art, units, statistic = "squared_loading", scope = "aggregate",
    k = "all", weights = "equal", normalize = "none"
  )
  fi <- feature_importance_from_bootstrap(
    art, units, k = "all", weights = "equal", normalize = "none"
  )

  expect_equal(ev$estimate, fi$estimate, tolerance = 1e-14)
  expect_equal(ev$scope, rep("aggregate", nrow(ev)))
  expect_equal(ev$orientation, rep("unsigned", nrow(ev)))
})

test_that("feature_evidence_pvalues refuses signed loading statistics", {
  ensure_default_adapters()
  rec <- infer_recipe(
    geometry = "oneblock", relation = "variance",
    adapter = "prcomp_oneblock"
  )
  adapter <- get_infer_adapter("prcomp_oneblock")
  X <- matrix(rnorm(40), 10, 4)
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))

  expect_error(
    feature_evidence_pvalues(
      rec, adapter, X, units, original_fit = fit,
      statistic = "loading", B = 2
    ),
    "'arg' should be one of"
  )
})

test_that("feature_evidence_pvalues returns null-calibrated unsigned rows", {
  ensure_default_adapters()
  set.seed(2041)
  X <- matrix(rnorm(80), 20, 4)
  rec <- infer_recipe(
    geometry = "oneblock", relation = "variance",
    adapter = "prcomp_oneblock"
  )
  adapter <- get_infer_adapter("prcomp_oneblock")
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit))

  out <- feature_evidence_pvalues(
    rec, adapter, X, units,
    original_fit = fit,
    statistic = "squared_loading",
    scope = "aggregate",
    B = 5L,
    adjust = "BH",
    seed = 99L,
    normalize = "none"
  )

  expect_s3_class(out, "multifer_feature_evidence")
  expect_equal(unique(out$calibration), "null_refit")
  expect_equal(unique(out$adjustment), "BH")
  expect_true(all(out$p_value >= 1 / 6 & out$p_value <= 1))
  expect_true(all(out$p_adjusted >= out$p_value))
  expect_true(all(is.na(out$z)))
})

test_that("feature_evidence_pvalues supports adapter-owned nonnegative statistic", {
  rec <- infer_recipe(
    geometry = "oneblock", relation = "variance",
    adapter = infer_adapter(
      adapter_id = "feature_stat_pvalue_adapter",
      shape_kinds = "oneblock",
      capabilities = capability_matrix(
        list(geometry = "oneblock", relation = "variance",
             targets = "variable_stability")
      ),
      roots = function(x, ...) c(1, 0.5),
      loadings = function(x, domain = NULL, ...) x$loadings,
      null_action = function(x, data, ...) data + matrix(0.1, nrow(data), ncol(data)),
      refit = function(x, new_data, ...) {
        list(
          loadings = matrix(c(colMeans(new_data), colMeans(new_data)), ncol = 2)
        )
      },
      variable_stat = function(x, data, domain, k, ...) {
        stats::setNames(abs(x$loadings[, k[1L]]), paste0("v", seq_len(nrow(x$loadings))))
      },
      validity_level = "conditional"
    ),
    targets = "variable_stability",
    strict = FALSE
  )
  adapter <- rec$adapter
  X <- matrix(seq_len(12), nrow = 4)
  fit <- adapter$refit(NULL, X)
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )

  out <- feature_evidence_pvalues(
    rec, adapter, X, units,
    original_fit = fit,
    statistic = "adapter",
    scope = "unit",
    k = "all",
    B = 3L,
    adjust = "none",
    normalize = "none",
    seed = 12L
  )

  expect_s3_class(out, "multifer_feature_evidence")
  expect_equal(unique(out$statistic), "adapter")
  expect_equal(unique(out$orientation), "adapter_invariant")
  expect_equal(unique(out$method), "adapter_feature_stat")
  expect_true(all(out$p_value >= 1 / 4 & out$p_value <= 1))
})

test_that("feature_evidence_from_adapter consumes native evidence action", {
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = TRUE
  )
  native <- infer_feature_evidence(
    scope = "unit",
    domain = "X",
    feature = "x1",
    unit_id = "u1",
    unit_set = "u1",
    statistic = "native_ratio",
    estimate = 2,
    bootstrap_mean = 2,
    bias = 0,
    std_error = 0.5,
    lower = 1,
    upper = 3,
    z = 4,
    z_label = "native bootstrap ratio",
    p_value = NA_real_,
    p_adjusted = NA_real_,
    adjustment = "none",
    calibration = "bootstrap",
    interval_method = "percentile",
    null_label = NA_character_,
    method = "conditional_subspace_bootstrap",
    orientation = "adapter_invariant",
    identifiable = TRUE,
    validity_level = "conditional",
    warnings = NA_character_
  )
  adapter <- structure(
    list(
      feature_evidence_action = function(x, data, units, design, statistic,
                                         orientation, R = NULL, seed = NULL,
                                         ...) native
    ),
    class = "multifer_adapter"
  )

  out <- feature_evidence_from_adapter(
    adapter, fit = list(), data = matrix(0, 1, 1), units = units,
    statistic = "native_ratio"
  )

  expect_s3_class(out, "multifer_feature_evidence")
  expect_equal(out$method, "conditional_subspace_bootstrap")
  expect_true(is.na(out$p_value))
})
