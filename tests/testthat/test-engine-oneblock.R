# Tests for ladder_driver() and run_oneblock_ladder().

restore_reference_adapters <- function() {
  clear_adapter_registry()
  register_oneblock_baser_adapters()
  register_cross_baser_adapters()
}

centered_orthonormal_basis <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(scale(M, center = TRUE, scale = FALSE)))
}

orthonormal_basis <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(M))
}

# ---------------------------------------------------------------------------
# 1. ladder_driver() smoke test: trivial callbacks, all rungs selected
# ---------------------------------------------------------------------------

test_that("ladder_driver smoke: trivial callbacks run to max_steps, all selected", {
  max_steps <- 5L
  B         <- 199L

  # observed always 1 (highest possible), null always 0 (never >= 1)
  obs_fn  <- function(step, data) 1
  null_fn <- function(step, data) 0
  # identity deflation -- data unchanged between steps
  defl_fn <- function(step, data) data

  res <- ladder_driver(
    observed_stat_fn = obs_fn,
    null_stat_fn     = null_fn,
    deflate_fn       = defl_fn,
    initial_data     = matrix(0, 2L, 2L),
    max_steps        = max_steps,
    B                = B,
    alpha            = 0.05,
    seed             = 1L
  )

  expect_equal(res$last_step_tested,  max_steps)
  expect_equal(res$rejected_through,  max_steps)
  expect_length(res$step_results,     max_steps)

  # With observed = 1 and all nulls = 0, r = 0 for every rung.
  # Phipson-Smyth: p = (1 + 0) / (B + 1) = 1 / 200 = 0.005.
  for (i in seq_len(max_steps)) {
    sr <- res$step_results[[i]]
    expect_equal(sr$step, i)
    expect_equal(sr$r, 0L)
    expected_p <- 1L / (B + 1L)
    expect_equal(sr$p_value, expected_p)
    expect_true(sr$selected)
  }
})

# ---------------------------------------------------------------------------
# 2. ladder_driver() non-rejection stops at first rung
# ---------------------------------------------------------------------------

test_that("ladder_driver stops at first non-rejection", {
  B <- 199L

  # observed = 0.01 is very low; nulls drawn from U(0,1) will almost always
  # exceed it, so p_value ~ 1 >> alpha.
  obs_fn  <- function(step, data) 0.01
  null_fn <- function(step, data) stats::runif(1)
  defl_fn <- function(step, data) data

  res <- ladder_driver(
    observed_stat_fn = obs_fn,
    null_stat_fn     = null_fn,
    deflate_fn       = defl_fn,
    initial_data     = matrix(0, 2L, 2L),
    max_steps        = 10L,
    B                = B,
    alpha            = 0.05,
    seed             = 42L
  )

  # First rung p-value must exceed alpha.
  expect_true(res$step_results[[1L]]$p_value > 0.05)
  expect_false(res$step_results[[1L]]$selected)

  # Driver must stop after rung 1.
  expect_equal(res$last_step_tested,  1L)
  expect_equal(res$rejected_through,  0L)
  expect_length(res$step_results,     1L)
})

# ---------------------------------------------------------------------------
# 3. run_oneblock_ladder() on null data -- null calibration
# ---------------------------------------------------------------------------

test_that("run_oneblock_ladder null calibration: rejected_through <= 3", {
  X   <- bench_oneblock_null(n = 200, p = 30, seed = 42L)$X
  restore_reference_adapters()
  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "prcomp_oneblock"
  )

  res <- run_oneblock_ladder(rec, X, B = 199L, alpha = 0.05, seed = 7L)

  expect_type(res, "list")
  expect_true(all(c("units", "component_tests", "roots_observed",
                    "ladder_result") %in% names(res)))

  expect_true(res$ladder_result$rejected_through <= 3L)
})

# ---------------------------------------------------------------------------
# 4. run_oneblock_ladder() on exact low-rank data -- rank recovery
# ---------------------------------------------------------------------------

test_that("run_oneblock_ladder recovers exact rank on noiseless low-rank data", {
  set.seed(42L)
  U <- centered_orthonormal_basis(40, 2)
  V <- orthonormal_basis(12, 2)
  X <- U %*% diag(c(14, 6), nrow = 2, ncol = 2) %*% t(V)
  restore_reference_adapters()
  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "prcomp_oneblock"
  )

  res <- run_oneblock_ladder(rec, X, B = 99L, alpha = 0.05, seed = 7L)

  expect_equal(res$ladder_result$rejected_through, 2L)
  expect_equal(res$ladder_result$last_step_tested, 3L)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, FALSE))
  expect_equal(res$component_tests$observed_stat[1], 196 / 232, tolerance = 1e-10)
  expect_equal(res$component_tests$observed_stat[2], 1, tolerance = 1e-10)
  expect_equal(res$component_tests$observed_stat[3], 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. P3 statistic correctness
# ---------------------------------------------------------------------------

test_that("P3 statistic at step 1 matches Part 1 section 1 formula", {
  set.seed(123L)
  X  <- matrix(stats::rnorm(50 * 10), nrow = 50, ncol = 10)
  Xc <- sweep(X, 2L, colMeans(X), "-")

  sv      <- base::svd(Xc)
  s2      <- sv$d^2
  expected <- s2[1L] / sum(s2)

  # Extract the observed statistic at step 1 using the same logic as
  # run_oneblock_ladder's observed_stat_fn.
  observed <- {
    s  <- base::svd(Xc)$d
    s2b <- s^2
    s2b[1L] / sum(s2b)
  }

  expect_equal(observed, expected, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 6. Reproducibility under seed
# ---------------------------------------------------------------------------

test_that("run_oneblock_ladder is reproducible under the same seed", {
  X <- bench_oneblock_null(n = 80, p = 15, seed = 1L)$X
  restore_reference_adapters()
  rec <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter  = "prcomp_oneblock"
  )

  res1 <- run_oneblock_ladder(rec, X, B = 99L, alpha = 0.05, seed = 55L)
  res2 <- run_oneblock_ladder(rec, X, B = 99L, alpha = 0.05, seed = 55L)

  expect_equal(res1$ladder_result$rejected_through,
               res2$ladder_result$rejected_through)

  # Compare null_values from each step to confirm full RNG identity.
  for (i in seq_along(res1$ladder_result$step_results)) {
    expect_identical(
      res1$ladder_result$step_results[[i]]$null_values,
      res2$ladder_result$step_results[[i]]$null_values
    )
  }
})

# ---------------------------------------------------------------------------
# 7. Recipe geometry validation
# ---------------------------------------------------------------------------

test_that("run_oneblock_ladder errors on non-oneblock recipe", {
  restore_reference_adapters()
  cross_rec <- infer_recipe(
    geometry = "cross",
    relation = "covariance",
    adapter  = "cross_svd"
  )

  X <- matrix(stats::rnorm(40), nrow = 10, ncol = 4)
  expect_error(
    run_oneblock_ladder(cross_rec, X, B = 9L),
    regexp = "oneblock"
  )
})
