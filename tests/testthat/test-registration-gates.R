# Registration-time gate tests for the v1 boundary.
#
# These tests pin the "impossible to misuse" gates that live in the
# adapter contract. Each malformed adapter here must be refused at
# `infer_adapter()` time, with an error message that either points at
# the spec note by filename or uses the spec's exact wording, so a
# future contributor cannot silently ship an in-sample predictive
# significance path or a Euclidean-deflating geneig adapter.
#
# Spec references:
#   notes/engine_predictive_spec.md §5 -- predictive cross-fit gate
#   notes/engine_geneig_spec.md §5     -- b_metric residualize gate
#                                      -- SPD Cholesky gate

make_predictive_caps <- function() {
  capability_matrix(
    list(
      geometry = "cross",
      relation = "predictive",
      targets  = "component_significance"
    )
  )
}

make_geneig_caps <- function() {
  capability_matrix(
    list(
      geometry = "geneig",
      relation = "generalized_eigen",
      targets  = "component_significance"
    )
  )
}

test_that("predictive gate refuses a non-split-aware component_stat", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad_plsr_insample",
      shape_kinds     = "cross",
      capabilities    = make_predictive_caps(),
      roots           = function(x) x$values,
      loadings        = function(x, domain) matrix(0, 1, 1),
      refit           = function(x, new_data) x,
      predict_response = function(x, new_data, k = NULL) matrix(0, 1, 1),
      null_action     = function(x, data) data,
      residualize     = function(x, k, data) data,
      # component_stat has no `split` argument on purpose -- this is
      # the failure mode the gate exists to catch.
      component_stat  = function(x, data, k) 1.0,
      validity_level  = "conditional"
    ),
    regexp = "Predictive admissibility rule"
  )
})

test_that("predictive gate refusal points at engine_predictive_spec.md", {
  err <- tryCatch(
    infer_adapter(
      adapter_id      = "bad_plsr_insample",
      shape_kinds     = "cross",
      capabilities    = make_predictive_caps(),
      roots           = function(x) x$values,
      loadings        = function(x, domain) matrix(0, 1, 1),
      refit           = function(x, new_data) x,
      predict_response = function(x, new_data, k = NULL) matrix(0, 1, 1),
      null_action     = function(x, data) data,
      residualize     = function(x, k, data) data,
      component_stat  = function(x, data, k) 1.0,
      validity_level  = "conditional"
    ),
    error = identity
  )
  expect_true(inherits(err, "error"))
  expect_match(conditionMessage(err), "notes/engine_predictive_spec.md",
               fixed = TRUE)
  expect_match(conditionMessage(err), "split", fixed = TRUE)
})

test_that("predictive gate refuses a missing predict_response hook", {
  expect_error(
    infer_adapter(
      adapter_id      = "bad_plsr_no_predict",
      shape_kinds     = "cross",
      capabilities    = make_predictive_caps(),
      roots           = function(x) x$values,
      loadings        = function(x, domain) matrix(0, 1, 1),
      refit           = function(x, new_data) x,
      null_action     = function(x, data) data,
      residualize     = function(x, k, data) data,
      component_stat  = function(x, data, k, split = NULL) 1.0,
      validity_level  = "conditional"
      # predict_response deliberately omitted
    ),
    regexp = "predict_response"
  )
})

test_that("geneig gate refuses an adapter with a non-b_metric residualize", {
  bad_residualize <- function(x, k, data) data
  # deliberately do NOT set attr(bad_residualize, "b_metric") <- TRUE

  expect_error(
    infer_adapter(
      adapter_id      = "bad_geneig_euclidean",
      shape_kinds     = "geneig",
      capabilities    = make_geneig_caps(),
      roots           = function(x) x$values,
      loadings        = function(x, domain) matrix(0, 1, 1),
      refit           = function(x, new_data) x,
      null_action     = function(x, data) data,
      residualize     = bad_residualize,
      component_stat  = function(x, data, k) 1.0,
      validity_level  = "conditional"
    ),
    regexp = "b_metric"
  )
})

test_that("geneig gate refusal points at engine_geneig_spec.md", {
  bad_residualize <- function(x, k, data) data
  err <- tryCatch(
    infer_adapter(
      adapter_id      = "bad_geneig_euclidean",
      shape_kinds     = "geneig",
      capabilities    = make_geneig_caps(),
      roots           = function(x) x$values,
      loadings        = function(x, domain) matrix(0, 1, 1),
      refit           = function(x, new_data) x,
      null_action     = function(x, data) data,
      residualize     = bad_residualize,
      component_stat  = function(x, data, k) 1.0,
      validity_level  = "conditional"
    ),
    error = identity
  )
  expect_true(inherits(err, "error"))
  expect_match(conditionMessage(err), "notes/engine_geneig_spec.md",
               fixed = TRUE)
})

test_that("geneig gate accepts a residualize with b_metric = TRUE", {
  good_residualize <- function(x, k, data) data
  attr(good_residualize, "b_metric") <- TRUE

  a <- infer_adapter(
    adapter_id      = "good_geneig",
    shape_kinds     = "geneig",
    capabilities    = make_geneig_caps(),
    roots           = function(x) x$values,
    loadings        = function(x, domain) matrix(0, 1, 1),
    refit           = function(x, new_data) x,
    null_action     = function(x, data) data,
    residualize     = good_residualize,
    component_stat  = function(x, data, k) 1.0,
    validity_level  = "conditional"
  )
  expect_true(is_infer_adapter(a))
})

test_that("geneig gate accepts a residualize with delegates_geneig_deflation", {
  delegating <- function(x, k, data) data
  attr(delegating, "delegates_geneig_deflation") <- TRUE

  a <- infer_adapter(
    adapter_id      = "good_geneig_delegating",
    shape_kinds     = "geneig",
    capabilities    = make_geneig_caps(),
    roots           = function(x) x$values,
    loadings        = function(x, domain) matrix(0, 1, 1),
    refit           = function(x, new_data) x,
    null_action     = function(x, data) data,
    residualize     = delegating,
    component_stat  = function(x, data, k) 1.0,
    validity_level  = "conditional"
  )
  expect_true(is_infer_adapter(a))
})

test_that("geneig_operator refuses a non-SPD B with a Cholesky-failed message", {
  A <- diag(3)
  # B has a zero on the diagonal so it is positive semi-definite but
  # not strictly positive definite; chol() will fail.
  B <- diag(c(1, 0, 1))
  expect_error(
    geneig_operator(A, B),
    regexp = "positive definite"
  )
  err <- tryCatch(geneig_operator(A, B), error = identity)
  expect_match(conditionMessage(err), "Cholesky failed", fixed = TRUE)
})

test_that("geneig_operator refuses a clearly indefinite B", {
  A <- diag(3)
  B <- diag(c(1, -0.5, 1))
  err <- tryCatch(geneig_operator(A, B), error = identity)
  expect_true(inherits(err, "error"))
  expect_match(conditionMessage(err), "positive definite", fixed = TRUE)
})

test_that("plsr_refit still registers after the gate hardening", {
  skip_if_not_installed("pls")
  a <- adapter_plsr_refit()
  expect_true(is_infer_adapter(a))
  expect_true(adapter_supports(a, "cross", "predictive", "component_significance"))
})

test_that("lda_refit still registers after the gate hardening", {
  skip_if_not_installed("MASS")
  a <- adapter_lda_refit()
  expect_true(is_infer_adapter(a))
  expect_true(adapter_supports(a, "geneig", "generalized_eigen",
                               "component_significance"))
})
