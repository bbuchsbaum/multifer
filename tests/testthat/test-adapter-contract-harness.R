test_that("check_infer_adapter smoke-checks all targets for oneblock adapters", {
  ensure_default_adapters()
  set.seed(7101)
  X <- matrix(rnorm(72), 18, 4)

  report <- check_infer_adapter(
    "svd_oneblock",
    X,
    B = 7L,
    R = 2L,
    seed = 11L
  )

  expect_s3_class(report, "multifer_adapter_check")
  expect_equal(
    sort(report$target),
    sort(c("component_significance", "variable_stability",
           "score_stability", "subspace_stability"))
  )
  expect_true(all(report$passed), info = paste(report$detail, collapse = "; "))
})

test_that("check_infer_adapter covers multi-relation cross adapters per target", {
  ensure_default_adapters()
  set.seed(7102)
  dat <- list(
    X = matrix(rnorm(120), 24, 5),
    Y = matrix(rnorm(96), 24, 4)
  )

  report <- check_infer_adapter(
    "cross_svd",
    dat,
    targets = "component_significance",
    B = 7L,
    R = 1L,
    seed = 12L
  )

  expect_equal(sort(report$relation), c("correlation", "covariance"))
  expect_true(all(report$passed), info = paste(report$detail, collapse = "; "))
})

test_that("check_infer_adapter can select fixtures by capability pair", {
  ensure_default_adapters()
  set.seed(7103)
  fixture <- list(
    default = "not valid data",
    "cross/correlation" = list(
      X = matrix(rnorm(100), 20, 5),
      Y = matrix(rnorm(80), 20, 4)
    )
  )

  report <- check_infer_adapter(
    "cancor_cross",
    fixture,
    targets = "component_significance",
    B = 7L,
    R = 1L,
    seed = 13L
  )

  expect_true(report$passed)
  expect_equal(report$geometry, "cross")
  expect_equal(report$relation, "correlation")
})

test_that("check_infer_adapter reports validity failures without stopping by default", {
  ensure_default_adapters()
  bad_X <- matrix(1, 10, 3)

  report <- check_infer_adapter(
    "svd_oneblock",
    bad_X,
    targets = "component_significance",
    run = "compile"
  )

  expect_false(report$passed)
  expect_equal(report$failed_stage, "checks")
  expect_match(report$detail, "failed validity check")
})

test_that("check_infer_adapter can fail hard for downstream contract tests", {
  ensure_default_adapters()
  bad_X <- matrix(1, 10, 3)

  expect_error(
    check_infer_adapter(
      "svd_oneblock",
      bad_X,
      targets = "component_significance",
      run = "compile",
      fail_on_error = TRUE
    ),
    "failed contract check"
  )
})

test_that("check_infer_adapter matches the shipped LDA capability boundary", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()
  set.seed(7104)
  y <- factor(rep(c("a", "b", "c"), each = 12L))
  X <- matrix(rnorm(length(y) * 4L), ncol = 4L)
  X[y == "a", 1L] <- X[y == "a", 1L] - 2
  X[y == "c", 1L] <- X[y == "c", 1L] + 2

  report <- check_infer_adapter(
    "lda_refit",
    list(X = X, y = y),
    B = 7L,
    R = 1L,
    seed = 14L
  )

  expect_equal(report$target, "component_significance")
  expect_true(report$passed, info = report$detail)
})

test_that("check_feature_evidence_adapter smoke-tests declared feature stats", {
  L <- matrix(c(1, 2, 3, 4), nrow = 2,
              dimnames = list(c("a", "b"), c("c1", "c2")))
  adapter <- infer_adapter(
    adapter_id = "feature_checker",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "variable_stability")
    ),
    roots = function(x, ...) c(1, 0.5),
    loadings = function(x, domain = NULL, ...) x$loadings,
    refit = function(x, new_data, ...) list(loadings = L),
    feature_stat_spec = function(x = NULL, data = NULL, ...) {
      data.frame(
        statistic = "vip",
        scope = "unit",
        orientation = "adapter_invariant",
        stringsAsFactors = FALSE
      )
    },
    feature_stat = function(x, data, domain, unit_id, members, statistic,
                            orientation, scope, ...) {
      stats::setNames(rowSums(abs(x$loadings[, members, drop = FALSE])),
                      rownames(x$loadings))
    },
    validity_level = "conditional"
  )
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, TRUE)
  )

  report <- check_feature_evidence_adapter(
    adapter = adapter,
    fit = list(loadings = L),
    data = matrix(0, 2, 2),
    units = units
  )

  expect_s3_class(report, "multifer_feature_evidence_adapter_check")
  expect_true(report$passed)
  expect_equal(report$statistic, "vip")
})

test_that("check_feature_evidence_adapter reports declared stat failures", {
  adapter <- infer_adapter(
    adapter_id = "feature_checker_bad",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "variable_stability")
    ),
    roots = function(x, ...) 1,
    loadings = function(x, domain = NULL, ...) matrix(1, 2, 1),
    refit = function(x, new_data, ...) list(),
    feature_stat_spec = function(x = NULL, data = NULL, ...) {
      data.frame(
        statistic = "bad",
        scope = "unit",
        orientation = "adapter_invariant",
        stringsAsFactors = FALSE
      )
    },
    feature_stat = function(...) c(NA_real_),
    validity_level = "conditional"
  )
  units <- infer_units(
    unit_id = "u1",
    unit_type = "component",
    members = list(1L),
    identifiable = TRUE,
    selected = TRUE
  )

  report <- check_feature_evidence_adapter(
    adapter = adapter,
    fit = list(),
    data = matrix(0, 2, 2),
    units = units
  )
  expect_false(report$passed)
  expect_match(report$detail, "finite numeric")
  expect_error(
    check_feature_evidence_adapter(
      adapter = adapter,
      fit = list(),
      data = matrix(0, 2, 2),
      units = units,
      fail_on_error = TRUE
    ),
    "failed feature evidence"
  )
})
