test_that("mc_budget_allocator tracks pool balance", {
  alloc <- mc_budget_allocator(B_total = 1000L, per_rung_cap = 300L)
  expect_equal(alloc$remaining_fn(), 1000L)
  expect_equal(alloc$used_fn(), 0L)

  c1 <- alloc$checkout()
  expect_equal(c1, 300L)   # capped by per_rung_cap
  alloc$record_use(c1, used = 250L)
  expect_equal(alloc$used_fn(), 250L)
  expect_equal(alloc$remaining_fn(), 750L)

  c2 <- alloc$checkout()
  expect_equal(c2, 300L)
  alloc$record_use(c2, used = 300L)
  expect_equal(alloc$used_fn(), 550L)
})

test_that("mc_budget_allocator caps by remaining pool when pool runs low", {
  alloc <- mc_budget_allocator(B_total = 100L, per_rung_cap = 80L)
  c1 <- alloc$checkout()
  alloc$record_use(c1, 60L)
  # remaining = 40 now, so next checkout is capped at 40 not 80.
  c2 <- alloc$checkout()
  expect_equal(c2, 40L)
})

test_that("mc_budget_allocator honors request argument", {
  alloc <- mc_budget_allocator(B_total = 500L, per_rung_cap = 100L)
  c1 <- alloc$checkout(request = 50L)
  expect_equal(c1, 50L)
  alloc$record_use(c1, 50L)
})

test_that("mc_budget_allocator schedule records actual draws per rung", {
  alloc <- mc_budget_allocator(B_total = 500L, per_rung_cap = 100L)
  for (i in 1:3) {
    cap <- alloc$checkout()
    alloc$record_use(cap, used = cap - 10L)
  }
  expect_equal(alloc$schedule, c(90L, 90L, 90L))
  expect_equal(alloc$used_fn(), 270L)
})

test_that("mc_budget_allocator validation", {
  expect_error(mc_budget_allocator(0),    "positive integer")
  expect_error(mc_budget_allocator(-5),   "positive integer")
  expect_error(mc_budget_allocator(100, per_rung_cap = 0), "positive integer")
  alloc <- mc_budget_allocator(100)
  expect_error(alloc$checkout(request = 0), "positive integer")
  expect_error(alloc$record_use(10, -1),    "non-negative")
  expect_error(alloc$record_use(10, 20),    "<= `allocated`")
})

test_that("print method does not error", {
  alloc <- mc_budget_allocator(100, per_rung_cap = 50)
  expect_output(print(alloc), "multifer_budget_allocator")
  expect_output(print(alloc), "total:")
  expect_output(print(alloc), "remaining:")
})

test_that("ladder_driver integrates allocator correctly", {
  set.seed(99)
  obs_fn <- function(step, data) 1
  null_fn <- function(step, data) stats::rnorm(1)  # almost always < 1
  defl_fn <- function(step, data) data

  res <- ladder_driver(
    observed_stat_fn = obs_fn,
    null_stat_fn     = null_fn,
    deflate_fn       = defl_fn,
    initial_data     = matrix(0, 2L, 2L),
    max_steps        = 4L,
    B                = 49L,
    B_total          = 100L,
    batch_size       = 8L,
    alpha            = 0.05,
    seed             = 11L
  )
  expect_true(inherits(res$allocator, "multifer_budget_allocator"))
  expect_true(res$allocator$used_fn() <= 100L)
  expect_true(res$total_draws_used == sum(res$allocator$schedule))
  expect_match(res$stopping_boundary, "besag_clifford")
})
