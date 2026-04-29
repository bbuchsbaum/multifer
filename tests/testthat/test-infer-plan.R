test_that("compile_infer_plan resolves adapter ladder callbacks", {
  a <- infer_adapter(
    adapter_id = "plan_oneblock",
    shape_kinds = "oneblock",
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets = "component_significance")
    ),
    roots = function(x, ...) 1,
    refit = function(x, new_data, ...) list(fit = TRUE),
    null_action = function(x, data, ...) data,
    component_stat = function(x, data, k, ...) 2,
    residualize = function(x, k, data, ...) data,
    validity_level = "conditional"
  )
  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    targets = "component_significance",
    adapter = a
  )

  seen <- new.env(parent = emptyenv())
  plan <- compile_infer_plan(
    rec,
    adapter = a,
    validate_data = function(data) {
      seen$validated <- TRUE
      invisible(data)
    },
    labels = list(statistic = "custom statistic")
  )

  expect_true(is_infer_plan(plan))
  expect_equal(plan$component_stat(list(), matrix(0, 2, 2), 1L), 2)
  expect_equal(plan$labels$statistic, "custom statistic")
  plan$validate_data(matrix(0, 2, 2))
  expect_true(seen$validated)
})

test_that("compile_infer_plan curries predictive component_stat signature", {
  seen <- new.env(parent = emptyenv())
  a <- infer_adapter(
    adapter_id = "plan_predictive",
    shape_kinds = "cross",
    capabilities = capability_matrix(
      list(geometry = "cross", relation = "predictive",
           targets = "component_significance")
    ),
    roots = function(x, ...) 1,
    refit = function(x, new_data, ...) list(fit = TRUE),
    null_action = function(x, data, ...) data,
    component_stat = function(x, data, k, split = NULL, ...) {
      seen$split <- split
      3
    },
    residualize = function(x, k, data, ...) data,
    predict_response = function(x, new_data, k = NULL, ...) {
      matrix(0, nrow = nrow(new_data$Y), ncol = ncol(new_data$Y))
    },
    validity_level = "conditional"
  )
  rec <- infer_recipe(
    geometry = "cross",
    relation = "predictive",
    targets = "component_significance",
    adapter = a
  )

  plan <- compile_infer_plan(rec, adapter = a)
  dat <- list(X = matrix(0, 4, 2), Y = matrix(0, 4, 1))

  expect_equal(plan$component_stat(list(), dat, 1L), 3)
  expect_null(seen$split)
  expect_match(plan$labels$statistic, "predictive")
})

test_that("compile_oneblock_ladder_plan reproduces first-rung P3 statistic", {
  set.seed(7201)
  X <- matrix(rnorm(60), 12, 5)
  ensure_default_adapters()
  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "prcomp_oneblock"
  )

  plan <- compile_oneblock_ladder_plan(rec, X)
  s2 <- svd(sweep(X, 2L, colMeans(X), "-"))$d^2

  expect_true(is_ladder_plan(plan))
  expect_equal(plan$roots_observed, s2, tolerance = 1e-12)
  expect_equal(plan$observed_stat_fn(1L, plan$initial_data),
               s2[1L] / sum(s2),
               tolerance = 1e-12)
  expect_equal(plan$labels$null, "column permutation of residual matrix")
})

test_that("run_oneblock_ladder delegates its callbacks through a ladder plan", {
  set.seed(7202)
  X <- matrix(rnorm(96), 16, 6)
  ensure_default_adapters()
  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "prcomp_oneblock"
  )

  plan <- compile_oneblock_ladder_plan(rec, X)
  ladder <- ladder_driver(
    observed_stat_fn = plan$observed_stat_fn,
    null_stat_fn = plan$null_stat_fn,
    deflate_fn = plan$deflate_fn,
    initial_data = plan$initial_data,
    max_steps = plan$max_steps,
    B = 29L,
    alpha = 0.05,
    seed = 3L
  )
  expected <- .ladder_plan_result(plan, ladder)
  observed <- run_oneblock_ladder(rec, X, B = 29L, alpha = 0.05, seed = 3L)

  expect_equal(observed$component_tests, expected$component_tests)
  expect_equal(observed$roots_observed, expected$roots_observed)
  expect_equal(observed$labels, expected$labels)
})
