test_that("typed_shape composes from a valid triple", {
  ts <- typed_shape(
    geometry = geometry("cross"),
    relation = relation("correlation"),
    design   = paired_rows()
  )
  expect_true(is_typed_shape(ts))
  expect_equal(ts$geometry$kind, "cross")
  expect_equal(ts$relation$kind, "correlation")
  expect_equal(ts$design$kind,   "paired_rows")
})

test_that("typed_shape rejects wrong slot types", {
  expect_error(
    typed_shape("oneblock", relation("variance"), exchangeable_rows()),
    "must be a multifer_geometry"
  )
  expect_error(
    typed_shape(geometry("oneblock"), "variance", exchangeable_rows()),
    "must be a multifer_relation"
  )
  expect_error(
    typed_shape(geometry("oneblock"), relation("variance"), "rows"),
    "must be a multifer_design"
  )
})

test_that("typed_shape does not enforce compatibility at construction", {
  # Compatibility checks belong to recipe compilation (issue 9qu.6),
  # not to typed_shape() itself. Part 5 §34 is explicit about this.
  ts <- typed_shape(
    geometry = geometry("cross"),
    relation = relation("variance"),  # unusual but not structurally invalid
    design   = paired_rows()
  )
  expect_true(is_typed_shape(ts))
})

test_that("typed_shape prints all three slots", {
  ts <- typed_shape(geometry("oneblock"), relation("variance"),
                    exchangeable_rows())
  out <- capture.output(print(ts))
  expect_true(any(grepl("geometry: oneblock", out)))
  expect_true(any(grepl("relation: variance", out)))
  expect_true(any(grepl("design:   exchangeable_rows", out)))
})
