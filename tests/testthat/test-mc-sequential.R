test_that("mc_sequential_bc matches mc_p_value when early stop never fires", {
  set.seed(1)
  # Observed 1e6, nulls are standard normal: r is always 0, no early stop.
  obs <- 1e6
  gen <- function() stats::rnorm(1)

  seq_res <- mc_sequential_bc(
    observed_stat = obs, null_gen_fn = gen, B_max = 199L, alpha = 0.05,
    batch_size = 32L, seed = 2L
  )
  expect_equal(seq_res$r, 0L)
  expect_equal(seq_res$drawn, 199L)
  expect_equal(seq_res$stop_reason, "exhausted")
  expect_equal(seq_res$p_value, 1 / (199L + 1L))
})

test_that("mc_sequential_bc stops early when r reaches h (non-rejection)", {
  # Observed = 0, nulls = 1: every draw exceeds, r grows as fast as drawn.
  obs <- 0
  gen <- function() 1
  # h = ceiling(0.05 * 201) = 11
  res <- mc_sequential_bc(
    observed_stat = obs, null_gen_fn = gen, B_max = 200L, alpha = 0.05,
    batch_size = 4L, seed = 3L
  )
  expect_equal(res$stop_reason, "non_reject_early")
  expect_true(res$drawn < 200L)
  # When stopped early, p = (1 + r) / (drawn + 1) should exceed alpha.
  expect_true(res$p_value > 0.05)
  expect_true(res$r >= res$h)
})

test_that("mc_sequential_bc default h agrees with Phipson-Smyth alpha boundary", {
  res <- mc_sequential_bc(
    observed_stat = 1, null_gen_fn = function() 0,
    B_max = 999L, alpha = 0.05, batch_size = 128L
  )
  expect_equal(res$h, as.integer(ceiling(0.05 * 1000L)))
})

test_that("mc_sequential_bc validation", {
  expect_error(mc_sequential_bc("bad", function() 1, 10L, 0.05),
               "single finite numeric")
  expect_error(mc_sequential_bc(1, "bad", 10L, 0.05),
               "must be a function")
  expect_error(mc_sequential_bc(1, function() 1, 0L, 0.05),
               "positive integer")
  expect_error(mc_sequential_bc(1, function() 1, 10L, 1.5),
               "in \\(0, 1\\)")
  expect_error(mc_sequential_bc(1, function() 1, 10L, 0.05, h = -1),
               "positive integer")
  expect_error(mc_sequential_bc(1, function() 1, 10L, 0.05, batch_size = 0),
               "positive integer")
})

test_that("mc_sequential_bc preserves reproducibility under fixed seed", {
  obs <- 0.5
  gen <- function() stats::runif(1)
  r1 <- mc_sequential_bc(obs, gen, 200L, 0.05, seed = 99L)
  r2 <- mc_sequential_bc(obs, gen, 200L, 0.05, seed = 99L)
  expect_equal(r1$p_value,  r2$p_value)
  expect_equal(r1$r,        r2$r)
  expect_equal(r1$drawn,    r2$drawn)
  expect_equal(r1$null_values, r2$null_values)
})

test_that("mc_sequential_bc batch_schedule sums to drawn", {
  res <- mc_sequential_bc(
    observed_stat = 0.2, null_gen_fn = function() stats::runif(1),
    B_max = 500L, alpha = 0.05, batch_size = 64L, seed = 7L
  )
  expect_equal(sum(res$batch_schedule), res$drawn)
})

test_that("mc_sequential_bc rejects non-finite null draws", {
  gen_bad <- function() NaN
  expect_error(mc_sequential_bc(1, gen_bad, 10L, 0.05),
               "finite numeric")
})
