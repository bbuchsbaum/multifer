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
