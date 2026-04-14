# Sprint 1: executable validity contracts
#
# Every mature adapter ships concrete checked_assumptions that run
# against the raw `data` argument before the engine starts. Strict
# mode (the default) errors on any failure; permissive mode records
# violations in result$assumptions$checked.

test_that("run_adapter_checks returns empty list when adapter has no checks", {
  adapter <- infer_adapter(
    adapter_id     = "test_empty",
    shape_kinds    = "oneblock",
    capabilities   = capability_matrix(
      list(geometry = "oneblock", relation = "variance", targets = "component_significance")
    ),
    roots          = function(x, ...) x$d^2,
    residualize    = function(x, k, data, ...) data,
    null_action    = function(x, data, ...) data,
    component_stat = function(x, data, k, ...) 0,
    validity_level = "heuristic",
    checked_assumptions = list()
  )
  X <- matrix(1:12, 3, 4)
  expect_equal(run_adapter_checks(adapter, X, strict = TRUE), list())
})

test_that("oneblock adapter_svd checks pass on a clean matrix", {
  adapter <- adapter_svd()
  X <- matrix(rnorm(40), 10, 4)
  results <- run_adapter_checks(adapter, X, strict = TRUE)
  expect_true(length(results) >= 3L)
  expect_true(all(vapply(results, function(r) isTRUE(r$passed), logical(1L))))
})

test_that("oneblock adapter_svd strict check errors on non-matrix data", {
  adapter <- adapter_svd()
  expect_error(
    run_adapter_checks(adapter, data.frame(x = 1:5, y = 6:10), strict = TRUE),
    "must be a numeric matrix"
  )
})

test_that("oneblock adapter_svd strict check errors on NA values", {
  adapter <- adapter_svd()
  X <- matrix(rnorm(40), 10, 4); X[1L, 1L] <- NA_real_
  expect_error(
    run_adapter_checks(adapter, X, strict = TRUE),
    "must not contain NA"
  )
})

test_that("oneblock adapter_svd strict check errors on too-small matrix", {
  adapter <- adapter_svd()
  expect_error(
    run_adapter_checks(adapter, matrix(1:3, nrow = 1, ncol = 3), strict = TRUE),
    "at least 2 rows"
  )
})

test_that("permissive mode returns violation records without erroring", {
  adapter <- adapter_svd()
  X <- matrix(rnorm(40), 10, 4); X[1L, 1L] <- NA_real_
  results <- run_adapter_checks(adapter, X, strict = FALSE)
  failed <- Filter(function(r) !isTRUE(r$passed), results)
  expect_true(length(failed) >= 1L)
  expect_match(failed[[1L]]$detail, "NA")
})

test_that("cross adapter checks pass on a clean (X, Y) list", {
  adapter <- adapter_cross_svd()
  X <- matrix(rnorm(60), 15, 4)
  Y <- matrix(rnorm(45), 15, 3)
  results <- run_adapter_checks(adapter, list(X = X, Y = Y), strict = TRUE)
  expect_true(length(results) >= 4L)
  expect_true(all(vapply(results, function(r) isTRUE(r$passed), logical(1L))))
})

test_that("cross adapter errors when X and Y have mismatched row counts", {
  adapter <- adapter_cross_svd()
  X <- matrix(rnorm(60), 15, 4)
  Y <- matrix(rnorm(48), 12, 4)     # 12 rows, not 15
  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), strict = TRUE),
    "paired-row"
  )
})

test_that("cross adapter errors when data is not a list", {
  adapter <- adapter_cross_svd()
  expect_error(
    run_adapter_checks(adapter, matrix(1:12, 3, 4), strict = TRUE),
    "list containing numeric matrices"
  )
})

test_that("cross adapter errors on infinite values", {
  adapter <- adapter_cross_svd()
  X <- matrix(rnorm(60), 15, 4); Y <- matrix(rnorm(45), 15, 3)
  Y[3L, 2L] <- Inf
  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), strict = TRUE),
    "NA, NaN, or Inf"
  )
})

test_that("infer() runs checks and surfaces failures under strict mode", {
  ensure_default_adapters()
  X <- matrix(rnorm(40), 10, 4); X[2L, 3L] <- NaN
  expect_error(
    infer(adapter = "svd_oneblock", data = X,
          geometry = "oneblock", relation = "variance",
          B = 19L, R = 4L, seed = 1L),
    "must not contain NA, NaN, or Inf"
  )
})

test_that("infer() populates result$assumptions$checked with passed records on clean input", {
  ensure_default_adapters()
  set.seed(1)
  X <- matrix(rnorm(200), 40, 5)
  res <- infer(
    adapter = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 49L, R = 6L, alpha = 0.05, seed = 7L
  )
  expect_true(is.list(res$assumptions$checked))
  expect_true(length(res$assumptions$checked) >= 3L)
  expect_true(all(vapply(res$assumptions$checked,
                         function(r) isTRUE(r$passed),
                         logical(1L))))
})

test_that("oneblock adapter check fires when a column has zero variance", {
  adapter <- adapter_svd()
  X <- matrix(rnorm(40), 10, 4)
  X[, 2L] <- 1.0  # constant column, zero variance
  expect_error(
    run_adapter_checks(adapter, X, strict = TRUE),
    "non-zero sample variance"
  )
})

test_that("cross adapter check fires on zero-variance column", {
  adapter <- adapter_cross_svd()
  X <- matrix(rnorm(60), 15, 4)
  Y <- matrix(rnorm(45), 15, 3)
  Y[, 1L] <- 3.14  # constant column
  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), strict = TRUE),
    "non-zero sample variance"
  )
})

test_that("cross adapter check fires on collinear columns (rank deficiency)", {
  adapter <- adapter_cross_svd()
  set.seed(17L)
  X <- matrix(rnorm(60), 15, 4)
  X[, 4L] <- 2 * X[, 1L] + 0.5 * X[, 2L]  # exact linear combination
  Y <- matrix(rnorm(45), 15, 3)
  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), strict = TRUE),
    "full numerical column rank"
  )
})

test_that("cross adapter permissive mode records rank-deficiency violation without erroring", {
  adapter <- adapter_cross_svd()
  set.seed(19L)
  X <- matrix(rnorm(60), 15, 4)
  X[, 4L] <- 2 * X[, 1L] + 0.5 * X[, 2L]
  Y <- matrix(rnorm(45), 15, 3)
  results <- run_adapter_checks(adapter, list(X = X, Y = Y), strict = FALSE)
  rank_result <- results[["cross_blocks_full_column_rank"]]
  expect_false(isTRUE(rank_result$passed))
  expect_match(rank_result$detail, "full numerical column rank")
})

test_that("correlation design checks are present on both CCA adapters", {
  expected <- c(
    "cross_correlation_sample_size",
    "cross_nuisance_design_rank",
    "cross_grouped_design_consistency"
  )

  for (adapter in list(adapter_cross_svd(), adapter_cancor())) {
    got <- vapply(adapter$checked_assumptions, `[[`, character(1L), "name")
    expect_true(all(expected %in% got))
  }
})

test_that("correlation sample-size check passes on a sufficiently large paired design", {
  adapter <- adapter_cross_svd()
  set.seed(101)
  X <- matrix(rnorm(30 * 4), 30, 4)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  recipe <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = paired_rows(),
    adapter = adapter
  )

  results <- run_adapter_checks(adapter, list(X = X, Y = Y), recipe = recipe, strict = TRUE)
  expect_true(isTRUE(results[["cross_correlation_sample_size"]]$passed))
})

test_that("correlation sample-size check fails when n does not exceed p + q", {
  adapter <- adapter_cross_svd()
  set.seed(102)
  X <- matrix(rnorm(7 * 4), 7, 4)
  Y <- matrix(rnorm(7 * 3), 7, 3)
  recipe <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = paired_rows(),
    adapter = adapter
  )

  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), recipe = recipe, strict = TRUE),
    "cross_correlation_sample_size"
  )
})

test_that("nuisance-design rank check passes on full-rank Z with enough residual rows", {
  adapter <- adapter_cancor()
  set.seed(103)
  X <- matrix(rnorm(30 * 4), 30, 4)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  Z <- cbind(1, scale(seq_len(30)), sin(seq_len(30) / 5))
  recipe <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = nuisance_adjusted(Z),
    adapter = adapter
  )

  results <- run_adapter_checks(adapter, list(X = X, Y = Y), recipe = recipe, strict = TRUE)
  expect_true(isTRUE(results[["cross_nuisance_design_rank"]]$passed))
})

test_that("nuisance-design rank check fails on rank-deficient Z", {
  adapter <- adapter_cancor()
  set.seed(104)
  X <- matrix(rnorm(30 * 4), 30, 4)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  z <- scale(seq_len(30))
  Z <- cbind(1, z, 2 * z)
  recipe <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = nuisance_adjusted(Z),
    adapter = adapter
  )

  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), recipe = recipe, strict = TRUE),
    "cross_nuisance_design_rank"
  )
})

test_that("grouped-design consistency check passes when exchangeable blocks are real", {
  adapter <- adapter_cross_svd()
  set.seed(105)
  X <- matrix(rnorm(24 * 4), 24, 4)
  Y <- matrix(rnorm(24 * 3), 24, 3)
  recipe <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = blocked_rows(rep(1:8, each = 3)),
    adapter = adapter
  )

  results <- run_adapter_checks(adapter, list(X = X, Y = Y), recipe = recipe, strict = TRUE)
  expect_true(isTRUE(results[["cross_grouped_design_consistency"]]$passed))
})

test_that("grouped-design consistency check fails on singleton-only groups", {
  adapter <- adapter_cross_svd()
  set.seed(106)
  X <- matrix(rnorm(20 * 4), 20, 4)
  Y <- matrix(rnorm(20 * 3), 20, 3)
  recipe <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = blocked_rows(seq_len(20)),
    adapter = adapter
  )

  expect_error(
    run_adapter_checks(adapter, list(X = X, Y = Y), recipe = recipe, strict = TRUE),
    "cross_grouped_design_consistency"
  )
})

test_that("infer_cca strict mode surfaces grouped-design validity failures by name", {
  ensure_default_adapters()
  set.seed(107)
  X <- matrix(rnorm(20 * 4), 20, 4)
  Y <- matrix(rnorm(20 * 3), 20, 3)

  expect_error(
    infer_cca(
      X, Y,
      design = blocked_rows(seq_len(20)),
      B = 9L,
      R = 2L,
      seed = 1L
    ),
    "cross_grouped_design_consistency"
  )
})

test_that("singleton grouped permutations are warning-free in permissive mode", {
  ensure_default_adapters()
  set.seed(1071)
  n <- 12L
  X <- matrix(rnorm(n * 4), n, 4)
  Y <- matrix(rnorm(n * 3), n, 3)

  expect_no_warning(
    infer_cca(
      X, Y,
      design = blocked_rows(seq_len(n)),
      strict = FALSE,
      B = 9L,
      R = 0L,
      seed = 1L
    )
  )
})

test_that("infer_cca permissive mode records nuisance-design validity failures", {
  ensure_default_adapters()
  set.seed(108)
  X <- matrix(rnorm(30 * 4), 30, 4)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  z <- scale(seq_len(30))
  Z <- cbind(1, z, 2 * z)

  res <- infer_cca(
    X, Y,
    design = nuisance_adjusted(Z),
    strict = FALSE,
    B = 9L,
    R = 2L,
    seed = 2L
  )

  expect_true(is.list(res$assumptions$checked))
  expect_false(isTRUE(res$assumptions$checked[["cross_nuisance_design_rank"]]$passed))
})

test_that("infer() cross-covariance runs paired-row + finite + min-dim checks end-to-end", {
  ensure_default_adapters()
  set.seed(2)
  X <- matrix(rnorm(60), 15, 4)
  Y <- matrix(rnorm(45), 15, 3)
  res <- infer(
    adapter = "cross_svd", data = list(X = X, Y = Y),
    geometry = "cross", relation = "covariance",
    B = 49L, R = 4L, alpha = 0.05, seed = 3L
  )
  expect_true(is.list(res$assumptions$checked))
  checked_names <- vapply(res$assumptions$checked, function(r) r$name, character(1L))
  expect_true("cross_paired_rows" %in% checked_names)
  expect_true("cross_blocks_are_finite" %in% checked_names)
})
