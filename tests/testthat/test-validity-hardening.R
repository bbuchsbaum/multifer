# High-discrimination validity tests for numerical helper paths and
# sequential Monte Carlo control flow. These tests target invariants that
# are easy to break with seemingly harmless optimizations.

make_stream_gen <- function(values) {
  i <- 0L
  force(values)
  function() {
    i <<- i + 1L
    if (i > length(values)) {
      stop("stream exhausted", call. = FALSE)
    }
    values[[i]]
  }
}

test_that("top_singular_values matches full SVD on canonical matrices", {
  set.seed(7001)

  tall <- matrix(stats::rnorm(80 * 12), nrow = 80, ncol = 12)
  wide <- matrix(stats::rnorm(12 * 80), nrow = 12, ncol = 80)

  U <- qr.Q(qr(matrix(stats::rnorm(20 * 3), nrow = 20, ncol = 3)))
  V <- qr.Q(qr(matrix(stats::rnorm(10 * 3), nrow = 10, ncol = 3)))
  lowrank <- U %*% diag(c(9, 1e-6, 1e-10), nrow = 3, ncol = 3) %*% t(V)

  zero <- matrix(0, nrow = 15, ncol = 10)

  cases <- list(
    list(name = "tall",    X = tall,    ks = c(1L, 2L, 3L), tol = 1e-10),
    list(name = "wide",    X = wide,    ks = c(1L, 2L, 3L), tol = 1e-10),
    list(name = "lowrank", X = lowrank, ks = c(1L, 2L, 3L), tol = 1e-6),
    list(name = "zero",    X = zero,    ks = c(1L, 2L, 3L), tol = 0)
  )

  for (case in cases) {
    X <- case$X
    ref <- base::svd(X)$d
    for (k in case$ks) {
      got <- top_singular_values(X, k)
      expect_true(
        max(abs(got - ref[seq_len(k)])) <= case$tol,
        info = sprintf("%s k=%d", case$name, k)
      )
    }
  }
})

test_that("mc_sequential_bc respects consumed-prefix semantics under early stop", {
  observed <- 0.5
  values <- c(0.6, 0.7, 0.1, 0.8, 0.2, 0.4, 0.9, 0.3, 0.1, 0.2)

  res <- mc_sequential_bc(
    observed_stat = observed,
    null_gen_fn   = make_stream_gen(values),
    B_max         = 10L,
    alpha         = 0.2,
    batch_size    = 2L
  )

  expect_equal(res$h, 3L)
  expect_equal(res$stop_reason, "non_reject_early")
  expect_equal(res$drawn, 4L)
  expect_equal(res$r, 3L)
  expect_equal(res$null_values, values[1:4])
  expect_equal(res$p_value, (1 + 3) / (4 + 1), tolerance = 1e-12)
  expect_equal(sum(res$batch_schedule), res$drawn)
})

test_that("mc_sequential_bc batch granularity only delays stopping, not the boundary logic", {
  observed <- 0.5
  values <- c(0.6, 0.7, 0.1, 0.8, 0.2, 0.4, 0.9, 0.3, 0.1, 0.2)

  fine <- mc_sequential_bc(
    observed_stat = observed,
    null_gen_fn   = make_stream_gen(values),
    B_max         = 10L,
    alpha         = 0.2,
    batch_size    = 1L
  )
  coarse <- mc_sequential_bc(
    observed_stat = observed,
    null_gen_fn   = make_stream_gen(values),
    B_max         = 10L,
    alpha         = 0.2,
    batch_size    = 5L
  )

  expect_equal(fine$stop_reason, "non_reject_early")
  expect_equal(coarse$stop_reason, "non_reject_early")
  expect_true(fine$drawn <= coarse$drawn)
  expect_true(fine$r >= fine$h)
  expect_true(coarse$r >= coarse$h)
  expect_true(all(fine$null_values %in% values))
  expect_true(all(coarse$null_values %in% values))
})

test_that("ladder_driver records explicit budget exhaustion and exact pool usage", {
  res <- ladder_driver(
    observed_stat_fn = function(step, data) 1,
    null_stat_fn     = function(step, data) 0,
    deflate_fn       = function(step, data) data,
    initial_data     = matrix(0, nrow = 2, ncol = 2),
    max_steps        = 4L,
    B                = 3L,
    B_total          = 5L,
    batch_size       = 1L,
    alpha            = 0.5,
    seed             = 99L
  )

  expect_equal(res$last_step_tested, 3L)
  expect_equal(res$rejected_through, 2L)
  expect_equal(length(res$step_results), 3L)

  expect_equal(res$step_results[[1]]$drawn, 3L)
  expect_equal(res$step_results[[2]]$drawn, 2L)
  expect_equal(res$step_results[[3]]$stop_reason, "budget_exhausted")
  expect_equal(res$step_results[[3]]$B, 0L)
  expect_equal(res$step_results[[3]]$drawn, 0L)
  expect_true(is.na(res$step_results[[3]]$p_value))

  expect_equal(res$total_draws_used, 5L)
  expect_equal(res$allocator$used_fn(), 5L)
  expect_equal(res$allocator$remaining_fn(), 0L)
  expect_equal(res$allocator$schedule, c(3L, 2L))
})
