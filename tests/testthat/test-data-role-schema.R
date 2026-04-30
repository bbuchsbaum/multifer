test_that("role validates kind, axes, and metric policy", {
  r <- role("primary", axes = c("row", "col"))
  expect_true(is_data_role(r))
  expect_equal(r$kind, "primary")
  expect_equal(r$axes, c("row", "col"))
  expect_null(r$policy)

  m <- role("metric", axis = "row",
            policy = "diagonal_observation_weights")
  expect_true(is_data_role(m))
  expect_equal(m$axes, "row")
  expect_equal(m$policy, "diagonal_observation_weights")

  expect_error(role("bogus"), "kind")
  expect_error(role("primary", axes = "time"), "unknown role axis")
  expect_error(role("static", axes = "row"), "Static roles")
  expect_error(role("primary", policy = "identity"), "only valid for metric")
  expect_error(role("metric", policy = "full_spd"), "metric `policy`")
})

test_that("data_role_schema requires named role objects", {
  schema <- data_role_schema(
    X = role("primary", axis = "row"),
    weights = role("metric", axis = "row",
                   policy = "diagonal_position_weights"),
    family = role("static")
  )

  expect_true(is_data_role_schema(schema))
  expect_named(schema, c("X", "weights", "family"))
  expect_equal(schema$weights$policy, "diagonal_position_weights")

  expect_error(data_role_schema(), "requires at least one")
  expect_error(data_role_schema(role("primary")), "must be named")
  expect_error(data_role_schema(X = list(kind = "primary")),
               "built with `role\\(\\)`")
})

test_that("null_spec records declarative null semantics", {
  spec <- null_spec(
    kind = "between_role_decoupling",
    randomized_role = "Y",
    reference_role = "X",
    axis = "row",
    preserves = c("X_rows", "Y_marginals")
  )

  expect_true(is_null_spec(spec))
  expect_equal(spec$kind, "between_role_decoupling")
  expect_equal(spec$axis, "row")
  expect_equal(spec$preserves, c("X_rows", "Y_marginals"))

  expect_error(null_spec("shuffle"), "kind")
  expect_error(null_spec("within_role_randomization", axis = "time"),
               "axis")
  expect_error(null_spec("parametric_draw", preserves = NA_character_),
               "preserves")
})
