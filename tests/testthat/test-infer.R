test_that("infer() oneblock end-to-end produces a valid infer_result", {
  ensure_default_adapters()
  set.seed(701)
  signal <- matrix(rnorm(40 * 3), 40, 3) %*% diag(c(8, 6, 4)) %*%
    matrix(rnorm(3 * 6), 3, 6)
  X <- signal + 0.1 * matrix(rnorm(240), 40, 6)
  res <- infer(
    adapter  = "prcomp_oneblock",
    data     = X,
    geometry = "oneblock",
    relation = "variance",
    B        = 99L,
    R        = 10L,
    seed     = 7
  )
  expect_true(is_infer_result(res))
  expect_true(nrow(res$units) >= 1L)
  expect_true(nrow(res$component_tests) >= 1L)
  expect_true(nrow(res$variable_stability) >= 1L)
  expect_true(nrow(res$score_stability) >= 1L)
  expect_equal(nrow(res$subspace_stability), nrow(res$units))
  expect_equal(res$provenance$adapter_id, "prcomp_oneblock")
})

test_that("infer() can use an adapter-owned oneblock component ladder", {
  calls <- new.env(parent = emptyenv())
  calls$stat <- 0L
  X <- matrix(seq_len(24), 6, 4)
  adapter <- infer_adapter(
    adapter_id = "adapter_engine_oneblock",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "component_significance")
    ),
    roots = function(x) x$values,
    refit = function(x, new_data) list(values = 3),
    null_action = function(x, data) data,
    component_stat = function(x, data, k) {
      calls$stat <- calls$stat + 1L
      0.5
    },
    residualize = function(x, k, data) data * 0,
    component_engine = "adapter",
    validity_level = "conditional",
    checked_assumptions = .oneblock_baser_checks()
  )

  res <- infer(
    adapter = adapter,
    data = X,
    geometry = "oneblock",
    relation = "variance",
    targets = "component_significance",
    B = 3L,
    R = 0L,
    seed = 2L
  )

  expect_true(is_infer_result(res))
  expect_gt(calls$stat, 0L)
  expect_equal(res$component_tests$statistic, 0.5)
  expect_match(res$mc$statistic_label, "adapter component_stat on oneblock data")
})

test_that("infer() supports stability-only adapter-owned geometry", {
  payload <- list(
    reference = matrix(seq_len(12), 4, 3),
    blocks = list(a = matrix(seq_len(10), 5, 2)),
    map = list(a = c(1L, 1L, 2L, 3L, 4L))
  )
  fit <- list(
    values = c(4, 1),
    loadings = matrix(c(1, 0, 0, 1), 2, 2)
  )
  calls <- new.env(parent = emptyenv())
  calls$component_stat <- 0L
  adapter <- infer_adapter(
    adapter_id = "opaque_stability_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = c("score_stability", "subspace_stability"))
    ),
    roots = function(x) x$values,
    domains = function(x, data = NULL) "block",
    loadings = function(x, domain = NULL) x$loadings,
    refit = function(x, new_data) fit,
    bootstrap_action = function(x, data, design, replicate = NULL) {
      list(fit = x, resample_indices = replicate)
    },
    project_scores = function(x, data, domain = NULL) {
      matrix(1, nrow = nrow(data$reference), ncol = ncol(x$loadings))
    },
    null_action = function(x, data) data,
    component_stat = function(x, data, k) {
      calls$component_stat <- calls$component_stat + 1L
      1
    },
    residualize = function(x, k, data) data,
    validity_level = "conditional"
  )

  res <- infer(
    adapter = adapter,
    data = payload,
    geometry = "adapter",
    relation = "variance",
    targets = c("score_stability", "subspace_stability"),
    B = 3L,
    R = 2L,
    seed = 4L
  )

  expect_true(is_infer_result(res))
  expect_equal(nrow(res$component_tests), 0L)
  expect_equal(calls$component_stat, 0L)
  expect_equal(nrow(res$score_stability), nrow(res$units) * nrow(payload$reference))
  expect_equal(nrow(res$subspace_stability), nrow(res$units))
  expect_equal(res$provenance$capabilities,
               "adapter/variance:score_stability+subspace_stability")
})

test_that("infer() supports component-significance adapter-owned geometry", {
  payload <- list(reference = matrix(seq_len(12), 4, 3))
  adapter <- infer_adapter(
    adapter_id = "opaque_component_adapter",
    shape_kinds = "adapter",
    capabilities = capability_matrix(
      list(geometry = "adapter", relation = "variance",
           targets = "component_significance")
    ),
    roots = function(x) x$values,
    refit = function(x, new_data) list(values = 2),
    null_action = function(x, data) data,
    component_stat = function(x, data, k) 0.75,
    residualize = function(x, k, data) data,
    validity_level = "conditional"
  )

  res <- infer(
    adapter = adapter,
    data = payload,
    geometry = "adapter",
    relation = "variance",
    targets = "component_significance",
    B = 3L,
    R = 0L,
    seed = 5L
  )

  expect_true(is_infer_result(res))
  expect_equal(res$component_tests$statistic, 0.75)
  expect_match(res$mc$statistic_label, "adapter component_stat on adapter data")
})

test_that("infer() cross end-to-end with explicit covariance", {
  ensure_default_adapters()
  set.seed(702)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = X, Y = Y),
    geometry = "cross",
    relation = "covariance",
    B        = 49L,
    R        = 6L,
    seed     = 9
  )
  expect_true(is_infer_result(res))
  expect_setequal(unique(res$variable_stability$domain), c("X", "Y"))
})

test_that("infer() reproduces under fixed seed", {
  ensure_default_adapters()
  set.seed(703)
  X <- matrix(rnorm(160), 40, 4)
  res1 <- infer(adapter = "prcomp_oneblock", data = X,
                geometry = "oneblock", relation = "variance",
                B = 49L, R = 5L, seed = 13)
  res2 <- infer(adapter = "prcomp_oneblock", data = X,
                geometry = "oneblock", relation = "variance",
                B = 49L, R = 5L, seed = 13)
  expect_equal(res1$component_tests$p_value, res2$component_tests$p_value)
  expect_equal(res1$component_tests$statistic, res2$component_tests$statistic)
})

test_that("infer() with strict dispatch errors on cross adapter without relation", {
  ensure_default_adapters()
  X <- matrix(rnorm(80), 20, 4)
  Y <- matrix(rnorm(60), 20, 3)
  expect_error(
    infer(adapter = "cross_svd", data = list(X = X, Y = Y),
          geometry = "cross", B = 9L, R = 2L, seed = 1),
    "covariance"
  )
})

test_that("infer() captures cost and mc reproducibility blocks", {
  ensure_default_adapters()
  set.seed(704)
  X <- matrix(rnorm(80), 20, 4)
  res <- infer(adapter = "prcomp_oneblock", data = X,
               geometry = "oneblock", relation = "variance",
               B = 49L, R = 5L, seed = 17)
  expect_true(!is.null(res$cost))
  expect_true(!is.null(res$mc))
  expect_equal(res$mc$rng_seed, 17L)
})

test_that("infer() rejects malformed adapter argument", {
  expect_error(infer(adapter = 123, data = matrix(0, 2, 2)),
               "registered adapter id")
})
