test_that("feature_importance_from_bootstrap computes aggregate squared-loading importance", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, FALSE)
  )
  reps <- list(
    list(aligned_loadings = list(X = matrix(c(1, 0, 0, 1), nrow = 2))),
    list(aligned_loadings = list(X = matrix(c(1, 0, 0, 1), nrow = 2)))
  )
  art <- structure(
    list(reps = reps, R = 2L, method_align = "sign", domains = "X", seed = NULL),
    class = "multifer_bootstrap_artifact"
  )

  out <- feature_importance_from_bootstrap(
    art, units, k = "selected", normalize = "none"
  )

  expect_s3_class(out, "multifer_feature_importance")
  expect_equal(out$scope, c("aggregate", "aggregate"))
  expect_equal(out$estimate, c(1, 0))
  expect_false("p_value" %in% names(out))
})

test_that("feature_importance_from_bootstrap handles subspace units and selected semantics", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("subspace", "component"),
    members = list(c(1L, 2L), 3L),
    identifiable = c(FALSE, TRUE),
    selected = c(TRUE, FALSE)
  )
  L <- matrix(
    c(1, 2, 0,
      0, 0, 3),
    nrow = 2,
    byrow = TRUE
  )
  reps <- list(
    list(aligned_loadings = list(X = L)),
    list(aligned_loadings = list(X = L))
  )
  art <- structure(
    list(reps = reps, R = 2L, method_align = "sign", domains = "X", seed = NULL),
    class = "multifer_bootstrap_artifact"
  )

  out <- feature_importance_from_bootstrap(
    art, units, k = "selected", scope = "both", normalize = "none"
  )

  agg <- out[out$scope == "aggregate", ]
  unit <- out[out$scope == "unit", ]
  expect_equal(agg$estimate, c(5, 0))
  expect_equal(unit$unit_id, c("u1", "u1"))
  expect_equal(unit$estimate, c(5, 0))
})

test_that("feature_importance_from_bootstrap supports root weights, ranks, and top-m frequency", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )
  L1 <- matrix(c(1, 0, 0, 1), nrow = 2)
  L2 <- matrix(c(0, 1, 1, 0), nrow = 2)
  reps <- list(
    list(aligned_loadings = list(X = L1)),
    list(aligned_loadings = list(X = L2))
  )
  art <- structure(
    list(reps = reps, R = 2L, method_align = "sign", domains = "X", seed = NULL),
    class = "multifer_bootstrap_artifact"
  )

  out <- feature_importance_from_bootstrap(
    art, units, k = "all", weights = "root", roots = c(3, 1),
    normalize = "block", top_m = 1L
  )

  expect_equal(sum(out$estimate), 1)
  expect_true(all(out$rank_mean >= 1))
  expect_true(all(out$top_m_frequency >= 0 & out$top_m_frequency <= 1))
})

test_that("feature_importance_pvalues uses null_action without a bootstrap artifact", {
  clear_adapter_registry()

  set.seed(71)
  X <- matrix(stats::rnorm(80), 20, 4)
  adapter <- adapter_svd()
  register_infer_adapter("svd_oneblock", adapter, overwrite = TRUE)
  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = adapter
  )
  fit <- adapter$refit(NULL, X)
  units <- form_units(adapter$roots(fit), selected = c(TRUE, FALSE, FALSE, FALSE))

  out <- feature_importance_pvalues(
    recipe = recipe,
    adapter = adapter,
    data = X,
    units = units,
    original_fit = fit,
    B = 5L,
    k = "selected",
    adjust = "maxT",
    seed = 72L
  )

  expect_s3_class(out, "multifer_feature_importance_pvalues")
  expect_equal(nrow(out), ncol(X))
  expect_true(all(out$p_value >= 1 / 6 & out$p_value <= 1))
  expect_true(all(out$p_adjusted >= out$p_value))
})
