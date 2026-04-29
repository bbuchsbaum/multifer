test_that("empty capability_matrix is valid", {
  cm <- capability_matrix()
  expect_true(is_capability_matrix(cm))
  expect_equal(nrow(cm), 0L)
})

test_that("capability_matrix expands entries into one row per target", {
  cm <- capability_matrix(
    list(geometry = "oneblock", relation = "variance",
         targets = c("component_significance", "variable_stability"))
  )
  expect_equal(nrow(cm), 2L)
  expect_setequal(cm$target,
                  c("component_significance", "variable_stability"))
})

test_that("capability_matrix deduplicates repeated targets", {
  cm <- capability_matrix(
    list(geometry = "cross", relation = "covariance",
         targets = c("component_significance", "component_significance"))
  )
  expect_equal(nrow(cm), 1L)
})

test_that("capability_matrix rejects unknown geometry / relation / target", {
  expect_error(
    capability_matrix(list(geometry = "pca", relation = "variance",
                           targets = "component_significance")),
    "`geometry` must be"
  )
  expect_error(
    capability_matrix(list(geometry = "oneblock", relation = "cov",
                           targets = "component_significance")),
    "`relation` must be"
  )
  expect_error(
    capability_matrix(list(geometry = "oneblock", relation = "variance",
                           targets = "significance")),
    "unknown target"
  )
})

test_that("capability_matrix accepts predictive for cross and adapter", {
  cm <- capability_matrix(
    list(geometry = "cross", relation = "predictive",
         targets = "component_significance"),
    list(geometry = "adapter", relation = "predictive",
         targets = "component_significance")
  )

  expect_true(supports(cm, "cross", "predictive", "component_significance"))
  expect_true(supports(cm, "adapter", "predictive", "component_significance"))

  expect_error(
    capability_matrix(
      list(geometry = "oneblock", relation = "predictive",
           targets = "component_significance")
    ),
    "cross.*adapter"
  )

  expect_error(
    capability_matrix(
      list(geometry = "geneig", relation = "predictive",
           targets = "component_significance")
    ),
    "cross.*adapter"
  )
})

test_that("capability_matrix rejects malformed entries", {
  expect_error(
    capability_matrix(list(geometry = "oneblock")),
    "must be a list with fields"
  )
  expect_error(
    capability_matrix(list(geometry = "oneblock", relation = "variance",
                           targets = character(0))),
    "non-empty character vector"
  )
})

test_that("supports() reports truth for declared and absent triples", {
  cm <- capability_matrix(
    list(geometry = "cross", relation = "covariance",
         targets = c("component_significance", "variable_stability"))
  )
  expect_true(supports(cm, "cross", "covariance",
                       "component_significance"))
  expect_true(supports(cm, "cross", "covariance",
                       "variable_stability"))
  expect_false(supports(cm, "cross", "correlation",
                        "component_significance"))
  expect_false(supports(cm, "cross", "covariance",
                        "score_stability"))
})

test_that("capabilities_for() returns only declared targets", {
  cm <- capability_matrix(
    list(geometry = "cross", relation = "covariance",
         targets = c("component_significance", "variable_stability")),
    list(geometry = "cross", relation = "correlation",
         targets = "component_significance")
  )
  expect_setequal(
    capabilities_for(cm, "cross", "covariance"),
    c("component_significance", "variable_stability")
  )
  expect_equal(
    capabilities_for(cm, "cross", "correlation"),
    "component_significance"
  )
  expect_equal(
    capabilities_for(cm, "oneblock", "variance"),
    character(0)
  )
})

test_that("supports() rejects non-scalar input", {
  cm <- capability_matrix()
  expect_error(
    supports(cm, c("cross", "oneblock"), "variance",
             "component_significance"),
    "single string"
  )
})

test_that("valid_capability_targets has the Part 5 section 38 set", {
  expect_setequal(
    valid_capability_targets(),
    c("component_significance", "variable_stability", "score_stability",
      "subspace_stability", "variable_significance")
  )
})

test_that("print method runs", {
  cm <- capability_matrix(
    list(geometry = "oneblock", relation = "variance",
         targets = "component_significance")
  )
  out <- capture.output(print(cm))
  expect_true(any(grepl("oneblock", out)))
  expect_true(any(grepl("component_significance", out)))
})
