test_that("list_benchmarks() returns exactly 4 entries all locked", {
  df <- list_benchmarks()
  expect_equal(nrow(df), 4L)
  expect_true(all(df$locked))
  expect_true(all(c("name", "geometry", "purpose", "locked") %in% names(df)))
})

test_that("list_benchmarks() contains the four expected suite names", {
  df <- list_benchmarks()
  expect_true("oneblock_null"       %in% df$name)
  expect_true("oneblock_shadowing"  %in% df$name)
  expect_true("cross_null"          %in% df$name)
  expect_true("speed_agreement"     %in% df$name)
})

test_that("bench_oneblock_null runs with defaults and returns correct fields", {
  res <- bench_oneblock_null()
  expect_type(res, "list")
  expect_true(!is.null(res$X))
  expect_true(!is.null(res$meta))
  expect_true(is.matrix(res$X))
  expect_equal(dim(res$X), c(200L, 50L))
  expect_equal(res$meta$n, 200L)
  expect_equal(res$meta$p, 50L)
  expect_equal(res$meta$true_rank, 0L)
  expect_equal(res$meta$noise, "gaussian")
})

test_that("bench_oneblock_shadowing runs with defaults and returns correct fields", {
  res <- bench_oneblock_shadowing()
  expect_type(res, "list")
  expect_true(is.matrix(res$X))
  expect_equal(dim(res$X), c(200L, 50L))
  expect_equal(res$meta$n, 200L)
  expect_equal(res$meta$p, 50L)
  expect_equal(res$meta$true_rank, 6L)
  expect_equal(res$meta$root_profile, c(8, 6, 4, 2, 1.2, 1.05))
})

test_that("bench_cross_null runs with defaults and returns correct fields", {
  res <- bench_cross_null()
  expect_type(res, "list")
  expect_true(is.matrix(res$X))
  expect_true(is.matrix(res$Y))
  expect_equal(dim(res$X), c(200L, 40L))
  expect_equal(dim(res$Y), c(200L, 30L))
  expect_equal(res$meta$n, 200L)
  expect_equal(res$meta$p_x, 40L)
  expect_equal(res$meta$p_y, 30L)
  expect_equal(res$meta$true_cross_rank, 0L)
})

test_that("bench_speed_agreement runs with medium and returns correct fields", {
  res <- bench_speed_agreement("medium")
  expect_type(res, "list")
  expect_true(is.matrix(res$X))
  expect_true(is.matrix(res$Y))
  expect_equal(res$meta$size, "medium")
  expect_equal(res$meta$n, 500L)
  expect_equal(res$meta$p_x, 100L)
  expect_equal(res$meta$p_y, 80L)
  expect_equal(res$meta$true_rank, 4L)
})

test_that("bench_speed_agreement runs with large and returns correct fields", {
  res <- bench_speed_agreement("large")
  expect_equal(res$meta$size, "large")
  expect_equal(res$meta$n, 2000L)
  expect_equal(res$meta$p_x, 500L)
  expect_equal(res$meta$p_y, 400L)
  expect_equal(res$meta$true_rank, 6L)
})

test_that("generators respect seed: same seed gives identical output", {
  r1 <- bench_oneblock_null(seed = 42L)
  r2 <- bench_oneblock_null(seed = 42L)
  expect_identical(r1$X, r2$X)

  r3 <- bench_oneblock_shadowing(seed = 7L)
  r4 <- bench_oneblock_shadowing(seed = 7L)
  expect_identical(r3$X, r4$X)

  r5 <- bench_cross_null(seed = 99L)
  r6 <- bench_cross_null(seed = 99L)
  expect_identical(r5$X, r6$X)
  expect_identical(r5$Y, r6$Y)

  r7 <- bench_speed_agreement("medium", seed = 1L)
  r8 <- bench_speed_agreement("medium", seed = 1L)
  expect_identical(r7$X, r8$X)
  expect_identical(r7$Y, r8$Y)
})

test_that("bench_oneblock_null_targets has correct over_selection_rate", {
  expect_equal(bench_oneblock_null_targets$over_selection_rate, 0.05)
  expect_true(bench_oneblock_null_targets$locked)
})

test_that("get_benchmark returns generator and targets for oneblock_null", {
  b <- get_benchmark("oneblock_null")
  expect_type(b, "list")
  expect_true(!is.null(b$generator))
  expect_true(!is.null(b$targets))
  expect_equal(b$targets$over_selection_rate, 0.05)
  expect_true(is.function(b$generator))
})

test_that("run_benchmark_generator works as thin wrapper", {
  res <- run_benchmark_generator("oneblock_null", n = 50, p = 10, seed = 5L)
  expect_equal(dim(res$X), c(50L, 10L))
})

test_that("get_benchmark errors clearly for unknown name", {
  expect_error(get_benchmark("nonexistent_suite"), "Unknown benchmark")
  expect_error(get_benchmark("nonexistent_suite"), "nonexistent_suite")
})

test_that("get_benchmark errors on non-character input", {
  expect_error(get_benchmark(123), "suite_name")
})

test_that("bench_oneblock_null noise variants all work", {
  r_ht <- bench_oneblock_null(n = 50, p = 10, noise = "heavy_tailed", seed = 1L)
  expect_equal(r_ht$meta$noise, "heavy_tailed")
  expect_equal(dim(r_ht$X), c(50L, 10L))

  r_hs <- bench_oneblock_null(n = 50, p = 10, noise = "heteroscedastic", seed = 1L)
  expect_equal(r_hs$meta$noise, "heteroscedastic")
  expect_equal(dim(r_hs$X), c(50L, 10L))
})

test_that("bench_oneblock_shadowing heteroscedastic noise works", {
  res <- bench_oneblock_shadowing(n = 50, p = 20,
                                   root_profile = c(4, 2),
                                   noise = "heteroscedastic",
                                   seed = 3L)
  expect_equal(res$meta$true_rank, 2L)
  expect_equal(dim(res$X), c(50L, 20L))
})
