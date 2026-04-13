test_that("empty_infer_result() produces a valid object", {
  r <- empty_infer_result()
  expect_true(is_infer_result(r))
  expect_equal(nrow(r$units), 0L)
  expect_equal(nrow(r$component_tests), 0L)
  expect_equal(nrow(r$variable_stability), 0L)
  expect_equal(nrow(r$score_stability), 0L)
  expect_equal(nrow(r$subspace_stability), 0L)
})

test_that("infer_units rejects length mismatches", {
  expect_error(
    infer_units(
      unit_id = c("u1", "u2"),
      unit_type = "component",
      members = list(1L),
      identifiable = TRUE,
      selected = TRUE
    ),
    "same length"
  )
})

test_that("infer_units rejects duplicate ids and unknown unit_type", {
  expect_error(
    infer_units(
      unit_id = c("u1", "u1"),
      unit_type = c("component", "component"),
      members = list(1L, 2L),
      identifiable = c(TRUE, TRUE),
      selected = c(TRUE, FALSE)
    ),
    "must be unique"
  )
  expect_error(
    infer_units(
      unit_id = "u1",
      unit_type = "bogus",
      members = list(1L),
      identifiable = TRUE,
      selected = TRUE
    ),
    "unknown unit_type"
  )
})

test_that("members list-column survives subsetting access", {
  u <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "subspace"),
    members = list(1L, c(2L, 3L)),
    identifiable = c(TRUE, FALSE),
    selected = c(TRUE, TRUE)
  )
  expect_equal(u$members[[1]], 1L)
  expect_equal(u$members[[2]], c(2L, 3L))
})

test_that("infer_component_tests rejects unknown validity_level", {
  expect_error(
    infer_component_tests(
      unit_id = "u1",
      statistic = 1.0,
      p_value = 0.01,
      mc_uncertainty = 0.005,
      stopped_at = 100L,
      null_label = "permute",
      validity_level = "believe-me"
    ),
    "unknown validity_level"
  )
})

test_that("infer_subspace_stability rejects unknown alignment_method", {
  expect_error(
    infer_subspace_stability(
      unit_id = "u1",
      principal_angle_mean = 0.1,
      principal_angle_max = 0.2,
      alignment_method = "guess",
      stability_label = "stable"
    ),
    "unknown alignment_method"
  )
})

test_that("infer_result enforces foreign keys on unit_id", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "component"),
    members = list(1L, 2L),
    identifiable = c(TRUE, TRUE),
    selected = c(TRUE, FALSE)
  )
  bad_tests <- infer_component_tests(
    unit_id = "u99",  # does not exist
    statistic = 1.0,
    p_value = 0.5,
    mc_uncertainty = 0.01,
    stopped_at = 100L,
    null_label = "permute",
    validity_level = "conditional"
  )
  expect_error(
    infer_result(units = units, component_tests = bad_tests),
    "references unknown unit_id"
  )
})

test_that("infer_result composes a full valid result", {
  units <- infer_units(
    unit_id = c("u1", "u2"),
    unit_type = c("component", "subspace"),
    members = list(1L, c(2L, 3L)),
    identifiable = c(TRUE, FALSE),
    selected = c(TRUE, TRUE)
  )
  tests <- infer_component_tests(
    unit_id = c("u1", "u2"),
    statistic = c(5.2, 3.1),
    p_value = c(0.001, 0.04),
    mc_uncertainty = c(0.001, 0.006),
    stopped_at = c(500L, 500L),
    null_label = c("sign_flip", "sign_flip"),
    validity_level = c("conditional", "conditional")
  )
  vs <- infer_variable_stability(
    unit_id = c("u1", "u1", "u2"),
    domain = c("X", "X", "X"),
    variable = c("v1", "v2", "v1"),
    estimate = c(0.7, 0.3, 0.5),
    stability = c(0.9, 0.6, 0.7),
    selection_freq = c(NA_real_, NA_real_, NA_real_),
    weight_sensitivity = c(NA_real_, NA_real_, NA_real_)
  )
  ss <- infer_score_stability(
    unit_id = "u1",
    domain = "X",
    observation = 1L,
    estimate = 0.1,
    lower = -0.1,
    upper = 0.3,
    leverage = NA_real_
  )
  sub <- infer_subspace_stability(
    unit_id = c("u1", "u2"),
    principal_angle_mean = c(0.05, 0.20),
    principal_angle_max = c(0.10, 0.45),
    alignment_method = c("sign", "subspace"),
    stability_label = c("stable", "tied")
  )

  r <- infer_result(
    units = units,
    component_tests = tests,
    variable_stability = vs,
    score_stability = ss,
    subspace_stability = sub,
    provenance = infer_provenance(
      adapter_id = "stub_pca",
      adapter_version = "0.0.1",
      capabilities = c("component_significance", "variable_stability")
    )
  )
  expect_true(is_infer_result(r))
  expect_equal(nrow(r$units), 2L)
  expect_equal(r$provenance$adapter_id, "stub_pca")
})

test_that("infer_result print method runs and mentions stability-vs-p-value", {
  r <- empty_infer_result()
  out <- capture.output(print(r))
  expect_true(any(grepl("infer_result", out)))
  expect_true(any(grepl("STABILITY measure, not a p-value", out)))
})

test_that("infer_result rejects non-schema sub-blocks", {
  expect_error(
    infer_result(units = data.frame(unit_id = "u1")),
    "must be built with infer_units"
  )
})

test_that("infer_assumptions, mc, cost, provenance construct with defaults", {
  a <- infer_assumptions(
    declared = "symmetric_noise",
    checked = list(rank = list(passed = TRUE, detail = "rank ok"))
  )
  m <- infer_mc(rng_seed = 42L, rng_kind = "Mersenne-Twister")
  ci <- infer_cost(full_data_ops = 1L, core_updates = 100L)
  p <- infer_provenance(adapter_id = "stub", adapter_version = "0.0.1")
  expect_s3_class(a, "multifer_assumptions")
  expect_s3_class(m, "multifer_mc")
  expect_s3_class(ci, "multifer_cost")
  expect_s3_class(p, "multifer_provenance")
})
