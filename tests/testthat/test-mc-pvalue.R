test_that("known p-value: observed far above null, alternative = greater", {
  result <- mc_p_value(
    observed_stat = 5,
    null_gen_fn   = function() rnorm(1),
    B             = 999L,
    alternative   = "greater",
    seed          = 1L
  )
  expect_true(result$p_value < 0.01)
  expect_equal(result$B, 999L)
  expect_equal(result$observed_stat, 5)
  expect_equal(result$alternative, "greater")
  expect_equal(result$seed, 1L)
  expect_length(result$null_values, 999L)
})

test_that("exact formula: all null > observed, r = B, p_value = 1", {
  B <- 50L
  result <- mc_p_value(
    observed_stat = 0,
    null_gen_fn   = function() 1,
    B             = B,
    alternative   = "greater"
  )
  expect_equal(result$r, B)
  expect_equal(result$p_value, (1 + B) / (B + 1))
  expect_equal(result$p_value, 1)
})

test_that("exact formula: no null >= observed, r = 0, p_value = 1/(B+1)", {
  B <- 99L
  result <- mc_p_value(
    observed_stat = 100,
    null_gen_fn   = function() 1,
    B             = B,
    alternative   = "greater"
  )
  expect_equal(result$r, 0L)
  expect_equal(result$p_value, 1 / (B + 1))
})

test_that("alternative = less works correctly", {
  result <- mc_p_value(
    observed_stat = -5,
    null_gen_fn   = function() rnorm(1),
    B             = 999L,
    alternative   = "less",
    seed          = 2L
  )
  expect_true(result$p_value < 0.01)
  expect_equal(result$alternative, "less")
})

test_that("alternative = two_sided works correctly (observed ~ 2 sigma)", {
  result <- mc_p_value(
    observed_stat = 2,
    null_gen_fn   = function() rnorm(1),
    B             = 999L,
    alternative   = "two_sided",
    seed          = 3L
  )
  # two-sided p-value for |z| >= 2 from N(0,1) is ~0.046
  expect_true(result$p_value > 0.01)
  expect_true(result$p_value < 0.15)
  expect_equal(result$alternative, "two_sided")
})

test_that("deterministic under same seed", {
  gen <- function() rnorm(1)
  r1 <- mc_p_value(observed_stat = 1.5, null_gen_fn = gen, B = 200L,
                   seed = 99L)
  r2 <- mc_p_value(observed_stat = 1.5, null_gen_fn = gen, B = 200L,
                   seed = 99L)
  expect_identical(r1$null_values, r2$null_values)
  expect_identical(r1$p_value, r2$p_value)
})

test_that("RNG state is restored after call with seed", {
  set.seed(42L)
  before <- .Random.seed

  mc_p_value(observed_stat = 0, null_gen_fn = function() rnorm(1),
             B = 10L, seed = 7L)

  after <- .Random.seed
  expect_identical(before, after)
})

test_that("mc_se formula is correct", {
  B <- 199L
  result <- mc_p_value(
    observed_stat = 100,
    null_gen_fn   = function() 1,
    B             = B,
    alternative   = "greater"
  )
  p <- result$p_value
  expected_se <- sqrt(p * (1 - p) / B)
  expect_equal(result$mc_se, expected_se)
})

test_that("B = 1 edge case: does not crash, returns valid p_value", {
  result <- mc_p_value(
    observed_stat = 0,
    null_gen_fn   = function() 1,
    B             = 1L,
    alternative   = "greater"
  )
  # null (1) >= observed (0): r=1, p = (1+1)/(1+1) = 1
  expect_equal(result$p_value, 1)
  expect_length(result$null_values, 1L)

  result2 <- mc_p_value(
    observed_stat = 100,
    null_gen_fn   = function() 1,
    B             = 1L,
    alternative   = "greater"
  )
  # null (1) < observed (100): r=0, p = 1/2
  expect_equal(result2$p_value, 0.5)
})

test_that("validation: non-integer B errors", {
  expect_error(
    mc_p_value(0, function() 1, B = 1.5),
    regexp = "positive integer",
    fixed  = FALSE
  )
})

test_that("validation: B <= 0 errors", {
  expect_error(
    mc_p_value(0, function() 1, B = 0L),
    regexp = "positive integer"
  )
  expect_error(
    mc_p_value(0, function() 1, B = -5L),
    regexp = "positive integer"
  )
})

test_that("validation: non-function null_gen_fn errors", {
  expect_error(
    mc_p_value(0, null_gen_fn = 42, B = 10L),
    regexp = "function"
  )
})

test_that("validation: NaN observed_stat errors", {
  expect_error(
    mc_p_value(NaN, function() 1, B = 10L),
    regexp = "finite"
  )
})

test_that("validation: null_gen_fn returning a vector errors", {
  expect_error(
    mc_p_value(0, null_gen_fn = function() c(1, 2), B = 5L),
    regexp = "single finite numeric"
  )
})

test_that("validation: null_gen_fn returning non-finite errors", {
  expect_error(
    mc_p_value(0, null_gen_fn = function() Inf, B = 5L),
    regexp = "single finite numeric"
  )
})
