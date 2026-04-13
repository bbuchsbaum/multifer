test_that("default mode: clean roots produce singleton component units", {
  res <- form_units(c(5, 3, 1))

  expect_equal(nrow(res), 3L)
  expect_true(all(res$unit_type == "component"))
  expect_true(all(res$identifiable == TRUE))
  expect_true(all(res$selected == TRUE))

  m <- attr(res, "members")
  expect_equal(m[[1L]], 1L)
  expect_equal(m[[2L]], 2L)
  expect_equal(m[[3L]], 3L)
})

test_that("default mode: partial selection threads through correctly", {
  res <- form_units(c(5, 3, 1), selected = c(TRUE, TRUE, FALSE))

  expect_equal(nrow(res), 3L)
  expect_true(all(res$unit_type == "component"))
  expect_equal(res$selected, c(TRUE, TRUE, FALSE))
})

test_that("default mode: near-tied roots still produce singletons (key v1 test)", {
  res <- form_units(c(5.0001, 5, 1))

  expect_equal(nrow(res), 3L)
  expect_true(all(res$unit_type == "component"))
  expect_true(all(res$identifiable == TRUE))

  m <- attr(res, "members")
  expect_equal(m[[1L]], 1L)
  expect_equal(m[[2L]], 2L)
  expect_equal(m[[3L]], 3L)
})

test_that("opt-in grouping: no near-ties produces all singletons", {
  res <- form_units(c(5, 3, 1), group_near_ties = TRUE)

  expect_equal(nrow(res), 3L)
  expect_true(all(res$unit_type == "component"))
  expect_true(all(res$identifiable == TRUE))
})

test_that("opt-in grouping: near-tied pair merges into a subspace", {
  res <- form_units(c(5.0001, 5, 1),
                    group_near_ties = TRUE,
                    tie_threshold   = 0.01)

  expect_equal(nrow(res), 2L)

  expect_equal(res$unit_type[1L], "subspace")
  expect_equal(res$identifiable[1L], FALSE)
  expect_equal(attr(res, "members")[[1L]], c(1L, 2L))

  expect_equal(res$unit_type[2L], "component")
  expect_equal(res$identifiable[2L], TRUE)
  expect_equal(attr(res, "members")[[2L]], 3L)
})

test_that("opt-in grouping: three-way tie merges into a single subspace", {
  res <- form_units(c(5, 5, 5, 1), group_near_ties = TRUE)

  expect_equal(nrow(res), 2L)

  expect_equal(res$unit_type[1L], "subspace")
  expect_equal(attr(res, "members")[[1L]], c(1L, 2L, 3L))

  expect_equal(res$unit_type[2L], "component")
  expect_equal(attr(res, "members")[[2L]], 4L)
})

test_that("opt-in grouping: subspace selected only when all members selected", {
  # First two roots are tied and form a subspace; third is singleton.
  # Case 1: both tied roots selected -> subspace selected = TRUE
  res_all <- form_units(c(5.0001, 5, 1),
                        selected        = c(TRUE, TRUE, TRUE),
                        group_near_ties = TRUE,
                        tie_threshold   = 0.01)
  expect_equal(res_all$selected[1L], TRUE)
  expect_equal(res_all$selected[2L], TRUE)

  # Case 2: first root selected, second not -> subspace selected = FALSE
  res_partial <- form_units(c(5.0001, 5, 1),
                            selected        = c(TRUE, FALSE, TRUE),
                            group_near_ties = TRUE,
                            tie_threshold   = 0.01)
  expect_equal(res_partial$selected[1L], FALSE)
  expect_equal(res_partial$selected[2L], TRUE)

  # Case 3: neither tied root selected -> subspace selected = FALSE
  res_none <- form_units(c(5.0001, 5, 1),
                         selected        = c(FALSE, FALSE, TRUE),
                         group_near_ties = TRUE,
                         tie_threshold   = 0.01)
  expect_equal(res_none$selected[1L], FALSE)
})

test_that("validation: non-sorted roots error", {
  expect_error(form_units(c(1, 3, 5)),
               regexp = "sorted in descending order",
               fixed  = FALSE)
})

test_that("validation: negative roots error", {
  expect_error(form_units(c(5, 3, -1)),
               regexp = "non-negative",
               fixed  = FALSE)
})

test_that("validation: mismatched selected length errors", {
  expect_error(form_units(c(5, 3, 1), selected = c(TRUE, FALSE)),
               regexp = "same length",
               fixed  = FALSE)
})

test_that("validation: negative tie_threshold errors", {
  expect_error(form_units(c(5, 3, 1), group_near_ties = TRUE,
                          tie_threshold = -0.01),
               regexp = "non-negative",
               fixed  = FALSE)
})

test_that("output is a valid multifer_units and composes into infer_result", {
  res <- form_units(c(5, 3, 1))

  expect_true(inherits(res, "multifer_units"))

  # Frozen schema round-trip: should not throw.
  ir <- infer_result(units = res)
  expect_true(inherits(ir, "infer_result"))
})
