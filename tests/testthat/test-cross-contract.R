# test-cross-contract.R
#
# Phase 0 cross-block contract-proof tests (issue 9qu.9).
# Mirrors test-oneblock-contract.R but for the cross geometry.
# The central value here is proving that strict dispatch actually catches
# the PLS-vs-CCA ambiguity at recipe compile time.
#
# Helper functions (make_stub_cross_adapter, make_stub_cross_infer_result,
# etc.) are provided by tests/testthat/helper-stub-adapters.R, which
# testthat auto-sources before any test runs.

# ---------------------------------------------------------------------------
# 1. Adapter registration lifecycle
# ---------------------------------------------------------------------------

test_that("stub_cross adapter registers and is retrievable", {
  clear_adapter_registry()

  adapter <- make_stub_cross_adapter("stub_cross")

  expect_true(is_infer_adapter(adapter))

  register_infer_adapter("stub_cross", adapter)

  expect_true("stub_cross" %in% list_infer_adapters())

  fetched <- get_infer_adapter("stub_cross")
  expect_true(is_infer_adapter(fetched))
  expect_equal(fetched$adapter_id, "stub_cross")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 2. Negative test: strict-mode ambiguity
#    A two-relation cross adapter + omitted relation + strict = TRUE must
#    error with a message that names BOTH candidate relations.
# ---------------------------------------------------------------------------

test_that("infer_recipe errors under strict = TRUE when cross relation is ambiguous", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross",
                         make_stub_cross_adapter("stub_cross",
                                                 relations = c("covariance",
                                                               "correlation")))

  # Error must mention "covariance"
  expect_error(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = TRUE),
    regexp = "covariance"
  )

  # Same call must also mention "correlation" (both candidates visible)
  expect_error(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = TRUE),
    regexp = "correlation"
  )

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 3. Positive test: covariance path resolves when relation is explicit
# ---------------------------------------------------------------------------

test_that("infer_recipe succeeds with explicit relation = covariance", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross",
                         make_stub_cross_adapter("stub_cross",
                                                 relations = c("covariance",
                                                               "correlation")))

  recipe <- infer_recipe(geometry = "cross",
                         relation  = "covariance",
                         adapter   = "stub_cross")

  expect_true(is_infer_recipe(recipe))
  expect_equal(recipe$shape$relation$kind, "covariance")
  expect_equal(recipe$shape$design$kind,   "paired_rows")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 4. Positive test: correlation path resolves and produces a different recipe
# ---------------------------------------------------------------------------

test_that("infer_recipe succeeds with explicit relation = correlation", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross",
                         make_stub_cross_adapter("stub_cross",
                                                 relations = c("covariance",
                                                               "correlation")))

  recipe_cov  <- infer_recipe(geometry = "cross",
                               relation  = "covariance",
                               adapter   = "stub_cross")

  recipe_corr <- infer_recipe(geometry = "cross",
                               relation  = "correlation",
                               adapter   = "stub_cross")

  expect_true(is_infer_recipe(recipe_corr))
  expect_equal(recipe_corr$shape$relation$kind, "correlation")

  # The two recipes must differ at the relation level
  expect_false(identical(recipe_cov$shape$relation$kind,
                         recipe_corr$shape$relation$kind))

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 5. Non-strict mode downgrade
#    strict = FALSE + omitted relation must warn, succeed, and downgrade
#    validity_level one step (conditional -> asymptotic).
# ---------------------------------------------------------------------------

test_that("infer_recipe warns and downgrades in non-strict mode with ambiguous relation", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross",
                         make_stub_cross_adapter("stub_cross",
                                                 relations     = c("covariance",
                                                                   "correlation"),
                                                 validity_level = "conditional"))

  expect_warning(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = FALSE),
    regexp = "covariance"
  )

  recipe <- withCallingHandlers(
    infer_recipe(geometry = "cross", adapter = "stub_cross", strict = FALSE),
    warning = function(w) invokeRestart("muffleWarning")
  )

  expect_true(is_infer_recipe(recipe))

  # Validity level must be one step worse: conditional -> asymptotic
  expect_equal(recipe$validity_level, "asymptotic")

  # downgrade_reason must be populated (non-NA)
  expect_false(is.na(recipe$downgrade_reason))

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 6. Single-relation cross adapter -- no ambiguity even without relation arg
# ---------------------------------------------------------------------------

test_that("infer_recipe strict mode succeeds when cross adapter has only one relation", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross_single",
                         make_stub_cross_adapter("stub_cross_single",
                                                 relations = "covariance"))

  # Must not error -- exactly one candidate relation, no ambiguity
  expect_no_error(
    infer_recipe(geometry = "cross",
                 adapter  = "stub_cross_single",
                 strict   = TRUE)
  )

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 7. Default design for cross geometry is paired_rows
# ---------------------------------------------------------------------------

test_that("infer_recipe default design for cross geometry is paired_rows", {
  clear_adapter_registry()
  register_infer_adapter("stub_cross_single",
                         make_stub_cross_adapter("stub_cross_single",
                                                 relations = "covariance"))

  recipe <- infer_recipe(geometry = "cross",
                         relation  = "covariance",
                         adapter   = "stub_cross_single")

  expect_equal(recipe$shape$design$kind, "paired_rows")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 8. Cross-block infer_result composes with two domains
# ---------------------------------------------------------------------------

test_that("make_stub_cross_infer_result produces a valid infer_result with X and Y domains", {
  res <- make_stub_cross_infer_result(
    unit_ids   = c("u1", "u2"),
    adapter_id = "stub_cross"
  )

  expect_true(is_infer_result(res))

  # Correct number of units
  expect_equal(nrow(res$units), 2L)

  # Foreign-key consistency across all downstream tables
  expect_true(all(res$component_tests$unit_id    %in% res$units$unit_id))
  expect_true(all(res$variable_stability$unit_id %in% res$units$unit_id))
  expect_true(all(res$score_stability$unit_id    %in% res$units$unit_id))
  expect_true(all(res$subspace_stability$unit_id %in% res$units$unit_id))

  # Both domains "X" and "Y" must appear in variable_stability
  expect_true("X" %in% res$variable_stability$domain)
  expect_true("Y" %in% res$variable_stability$domain)

  # Both domains must appear in score_stability
  expect_true("X" %in% res$score_stability$domain)
  expect_true("Y" %in% res$score_stability$domain)

  # Provenance carries the right adapter id
  expect_equal(res$provenance$adapter_id, "stub_cross")
})

# ---------------------------------------------------------------------------
# 9. Subspace unit in a cross context
# ---------------------------------------------------------------------------

test_that("infer_units accepts a subspace unit in a cross result", {
  units <- infer_units(
    unit_id      = c("u1", "u2"),
    unit_type    = c("component", "subspace"),
    members      = list(1L, c(1L, 2L)),
    identifiable = c(TRUE, FALSE),
    selected     = c(TRUE, TRUE)
  )

  expect_equal(units$unit_type[2L], "subspace")
  expect_false(units$identifiable[2L])
  expect_equal(units$members[[2L]], c(1L, 2L))

  # Full infer_result must accept this
  res <- infer_result(
    units = units,
    component_tests = infer_component_tests(
      unit_id        = c("u1", "u2"),
      statistic      = c(4.0, 2.0),
      p_value        = c(0.01, 0.04),
      mc_uncertainty = c(0.002, 0.006),
      stopped_at     = c(500L, 500L),
      null_label     = c("sign_flip", "sign_flip"),
      validity_level = c("conditional", "conditional")
    ),
    subspace_stability = infer_subspace_stability(
      unit_id              = c("u1", "u2"),
      principal_angle_mean = c(0.03, 0.20),
      principal_angle_max  = c(0.06, 0.40),
      alignment_method     = c("sign", "subspace"),
      stability_label      = c("stable", "tied")
    )
  )

  expect_true(is_infer_result(res))
  expect_equal(nrow(res$units), 2L)
})

# ---------------------------------------------------------------------------
# 10. Print output mentions STABILITY, not p-value
# ---------------------------------------------------------------------------

test_that("print.infer_result contains STABILITY measure warning for cross result", {
  res <- make_stub_cross_infer_result()
  out <- capture.output(print(res))
  expect_true(any(grepl("STABILITY measure, not a p-value", out)))
})

# ---------------------------------------------------------------------------
# 11. Cleanup guard: registry is empty at end of this file
# ---------------------------------------------------------------------------

test_that("adapter registry is clean after all cross-contract tests", {
  clear_adapter_registry()
  expect_equal(length(list_infer_adapters()), 0L)
})
