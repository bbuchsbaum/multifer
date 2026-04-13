# test-oneblock-contract.R
#
# Phase 0 canary: end-to-end contract-proof test for the oneblock PCA-like
# case using typed_shape + infer_recipe + stub adapter + frozen infer_result
# schema. No real inference engine runs -- all numbers are canned.

# ---------------------------------------------------------------------------
# 1. Adapter registration
# ---------------------------------------------------------------------------

test_that("stub_pca adapter registers and is retrievable", {
  clear_adapter_registry()

  adapter <- make_stub_oneblock_adapter("stub_pca")

  # Object passes is_infer_adapter before registration
  expect_true(is_infer_adapter(adapter))

  register_infer_adapter("stub_pca", adapter)

  # Listed after registration
  expect_true("stub_pca" %in% list_infer_adapters())

  # Retrievable by id
  fetched <- get_infer_adapter("stub_pca")
  expect_true(is_infer_adapter(fetched))
  expect_equal(fetched$adapter_id, "stub_pca")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 2. Typed shape composes
# ---------------------------------------------------------------------------

test_that("typed_shape composes and prints for oneblock/variance", {
  ts <- typed_shape(
    geometry("oneblock"),
    relation("variance"),
    exchangeable_rows()
  )

  expect_true(is_typed_shape(ts))
  expect_equal(ts$geometry$kind, "oneblock")
  expect_equal(ts$relation$kind, "variance")
  expect_equal(ts$design$kind, "exchangeable_rows")

  # print must not error
  expect_output(print(ts), "oneblock")
})

# ---------------------------------------------------------------------------
# 3. Recipe compiles against typed shape with strict = TRUE
# ---------------------------------------------------------------------------

test_that("infer_recipe compiles from typed_shape with default targets", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", make_stub_oneblock_adapter("stub_pca"))

  ts <- typed_shape(
    geometry("oneblock"),
    relation("variance"),
    exchangeable_rows()
  )

  recipe <- infer_recipe(shape = ts, targets = "default", adapter = "stub_pca")

  expect_true(is_infer_recipe(recipe))

  # Default targets expand to the four v1 targets; variable_significance
  # must NOT appear.
  expect_true("component_significance" %in% recipe$targets)
  expect_true("variable_stability"     %in% recipe$targets)
  expect_true("score_stability"        %in% recipe$targets)
  expect_true("subspace_stability"     %in% recipe$targets)
  expect_false("variable_significance" %in% recipe$targets)

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 4. Recipe compiles when shape is omitted -- geometry/relation passed directly
# ---------------------------------------------------------------------------

test_that("infer_recipe compiles from geometry + relation strings", {
  clear_adapter_registry()
  register_infer_adapter("stub_pca", make_stub_oneblock_adapter("stub_pca"))

  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca"
  )

  expect_true(is_infer_recipe(recipe))
  expect_equal(recipe$shape$geometry$kind, "oneblock")
  expect_equal(recipe$shape$relation$kind, "variance")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 5. infer_result composes and passes schema checks
# ---------------------------------------------------------------------------

test_that("make_stub_oneblock_infer_result produces a valid infer_result", {
  res <- make_stub_oneblock_infer_result(
    unit_ids   = c("u1", "u2"),
    adapter_id = "stub_pca"
  )

  expect_true(is_infer_result(res))

  # Correct number of units
  expect_equal(nrow(res$units), 2L)

  # Foreign-key consistency: all downstream unit_ids in the units table
  expect_true(all(res$component_tests$unit_id    %in% res$units$unit_id))
  expect_true(all(res$variable_stability$unit_id %in% res$units$unit_id))
  expect_true(all(res$score_stability$unit_id    %in% res$units$unit_id))
  expect_true(all(res$subspace_stability$unit_id %in% res$units$unit_id))

  # Provenance carries the right adapter id
  expect_equal(res$provenance$adapter_id, "stub_pca")
})

# ---------------------------------------------------------------------------
# 6. Unit schema is honest about subspace vs component
# ---------------------------------------------------------------------------

test_that("infer_units accepts a subspace unit with identifiable = FALSE", {
  units <- infer_units(
    unit_id      = c("u1", "u2", "u3"),
    unit_type    = c("component", "subspace", "component"),
    members      = list(1L, c(2L, 3L), 4L),
    identifiable = c(TRUE, FALSE, TRUE),
    selected     = c(TRUE, TRUE, FALSE)
  )

  expect_equal(units$unit_type[2L], "subspace")
  expect_false(units$identifiable[2L])
  expect_equal(units$members[[2L]], c(2L, 3L))

  # Constructor must accept this in a full infer_result
  res <- infer_result(
    units = units,
    component_tests = infer_component_tests(
      unit_id        = c("u1", "u2", "u3"),
      statistic      = c(4.0, 2.5, 1.0),
      p_value        = c(0.01, 0.04, 0.50),
      mc_uncertainty = c(0.002, 0.006, 0.01),
      stopped_at     = c(500L, 500L, 500L),
      null_label     = c("sign_flip", "sign_flip", "sign_flip"),
      validity_level = c("conditional", "conditional", "conditional")
    ),
    subspace_stability = infer_subspace_stability(
      unit_id              = c("u1", "u2", "u3"),
      principal_angle_mean = c(0.03, 0.25, 0.04),
      principal_angle_max  = c(0.06, 0.50, 0.08),
      alignment_method     = c("sign", "subspace", "sign"),
      stability_label      = c("stable", "tied", "stable")
    )
  )

  expect_true(is_infer_result(res))
  expect_equal(nrow(res$units), 3L)
})

# ---------------------------------------------------------------------------
# 7. Print output contains the stability-vs-p-value warning string
# ---------------------------------------------------------------------------

test_that("print.infer_result contains STABILITY measure warning", {
  res <- make_stub_oneblock_infer_result()
  out <- capture.output(print(res))
  expect_true(any(grepl("STABILITY measure, not a p-value", out)))
})

# ---------------------------------------------------------------------------
# 8. Recipe validity_level propagates from adapter
# ---------------------------------------------------------------------------

test_that("recipe validity_level matches adapter declaration", {
  clear_adapter_registry()

  adapter <- make_stub_oneblock_adapter("stub_pca", validity_level = "conditional")
  register_infer_adapter("stub_pca", adapter)

  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "stub_pca"
  )

  expect_equal(recipe$validity_level, "conditional")

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 9. Strict dispatch: single relation -- no error when relation omitted
# ---------------------------------------------------------------------------

test_that("infer_recipe strict mode succeeds when adapter has only one relation", {
  clear_adapter_registry()

  # stub_pca declares only (oneblock, variance) -- no ambiguity
  register_infer_adapter("stub_pca", make_stub_oneblock_adapter("stub_pca"))

  # Omitting relation under strict = TRUE must succeed because there is
  # exactly one candidate relation for this geometry.
  expect_no_error(
    infer_recipe(geometry = "oneblock", adapter = "stub_pca", strict = TRUE)
  )

  clear_adapter_registry()
})

# ---------------------------------------------------------------------------
# 10. Cleanup guard: registry is empty at end of this file
# ---------------------------------------------------------------------------

test_that("adapter registry is clean after all tests in this file", {
  clear_adapter_registry()
  expect_equal(length(list_infer_adapters()), 0L)
})
