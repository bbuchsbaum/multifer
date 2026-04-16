test_that("multifer_task_seeds is deterministic given a master seed", {
  s1 <- multifer_task_seeds(10L, master_seed = 42L)
  s2 <- multifer_task_seeds(10L, master_seed = 42L)
  expect_equal(s1, s2)
  s3 <- multifer_task_seeds(10L, master_seed = 43L)
  expect_false(identical(s1, s3))
})

test_that("multifer_task_seeds preserves caller's RNG stream", {
  set.seed(100)
  before <- stats::runif(1)
  set.seed(100)
  .ignore <- multifer_task_seeds(50L, master_seed = 7L)
  after <- stats::runif(1)
  expect_equal(before, after)
})

test_that("multifer_parallel_lapply sequential backend applies seeds", {
  seeds <- c(1L, 2L, 3L)
  result <- multifer_parallel_lapply(
    X = list("a", "b", "c"),
    FUN = function(x) stats::runif(1),
    seeds = seeds,
    backend = "sequential"
  )
  expect_length(result, 3L)
  # Same seeds => same draws.
  result2 <- multifer_parallel_lapply(
    X = list("a", "b", "c"),
    FUN = function(x) stats::runif(1),
    seeds = seeds,
    backend = "sequential"
  )
  expect_equal(result, result2)
})

test_that("multifer_parallel_lapply validation", {
  expect_error(
    multifer_parallel_lapply(X = 1:3, FUN = identity,
                             seeds = 1:2, backend = "sequential"),
    "same length"
  )
})

test_that("multifer_parallel_lapply is bit-identical across backends with fixed seeds", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  on.exit(try(multifer_parallel_shutdown(), silent = TRUE), add = TRUE)

  seeds <- c(11L, 12L, 13L)
  seq_out <- multifer_parallel_lapply(
    X = list("a", "b", "c"),
    FUN = function(x) stats::runif(1),
    seeds = seeds,
    backend = "sequential"
  )
  mir_out <- multifer_parallel_lapply(
    X = list("a", "b", "c"),
    FUN = function(x) stats::runif(1),
    seeds = seeds,
    backend = "mirai"
  )

  expect_equal(seq_out, mir_out)
})

test_that("bootstrap_fits parallel=mirai produces valid artifact", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  ensure_default_adapters()
  on.exit(try(multifer_parallel_shutdown(), silent = TRUE), add = TRUE)

  set.seed(1)
  X <- matrix(rnorm(40 * 5), 40, 5)
  adapter <- adapter_svd()
  recipe  <- infer_recipe(
    geometry = "oneblock", relation = "variance",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(NULL, X)
  units <- form_units(original_fit$d^2)

  artifact <- bootstrap_fits(
    recipe = recipe, adapter = adapter, data = X,
    original_fit = original_fit, units = units,
    R = 6L, method_align = "sign",
    seed = 123L, parallel = "mirai"
  )

  expect_equal(artifact$parallel, "mirai")
  expect_equal(length(artifact$reps), 6L)
  # Each rep has the four canonical fields.
  keys <- c("fit", "aligned_loadings", "aligned_scores", "resample_indices")
  expect_true(all(vapply(
    artifact$reps,
    function(r) all(keys %in% names(r)),
    logical(1L)
  )))
})

test_that("bootstrap_fits parallel=mirai is reproducible under fixed seed", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  ensure_default_adapters()
  on.exit(try(multifer_parallel_shutdown(), silent = TRUE), add = TRUE)

  set.seed(2)
  X <- matrix(rnorm(30 * 4), 30, 4)
  adapter <- adapter_svd()
  recipe  <- infer_recipe(
    geometry = "oneblock", relation = "variance",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(NULL, X)
  units <- form_units(original_fit$d^2)

  a1 <- bootstrap_fits(
    recipe = recipe, adapter = adapter, data = X,
    original_fit = original_fit, units = units,
    R = 4L, method_align = "sign", seed = 77L, parallel = "mirai"
  )
  a2 <- bootstrap_fits(
    recipe = recipe, adapter = adapter, data = X,
    original_fit = original_fit, units = units,
    R = 4L, method_align = "sign", seed = 77L, parallel = "mirai"
  )
  # Same resample indices.
  expect_equal(
    lapply(a1$reps, `[[`, "resample_indices"),
    lapply(a2$reps, `[[`, "resample_indices")
  )
  # Same aligned loadings.
  expect_equal(
    lapply(a1$reps, `[[`, "aligned_loadings"),
    lapply(a2$reps, `[[`, "aligned_loadings")
  )
})

test_that("infer() ladder $mc output is bit-identical across parallel backends", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  ensure_default_adapters()
  on.exit(try(multifer_parallel_shutdown(), silent = TRUE), add = TRUE)

  set.seed(3)
  X <- matrix(rnorm(30 * 4), 30, 4)

  r_seq <- infer(
    adapter  = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 49L, R = 4L, alpha = 0.05, seed = 99L,
    parallel = "sequential"
  )
  r_mir <- infer(
    adapter  = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 49L, R = 4L, alpha = 0.05, seed = 99L,
    parallel = "mirai"
  )
  # Ladder-path outputs are single-threaded in both runs.
  expect_equal(r_seq$mc$exceedance_counts, r_mir$mc$exceedance_counts)
  expect_equal(r_seq$mc$total_draws_used,  r_mir$mc$total_draws_used)
  expect_equal(r_seq$mc$stopping_boundary, r_mir$mc$stopping_boundary)
  expect_equal(r_seq$component_tests$p_value, r_mir$component_tests$p_value)
})
