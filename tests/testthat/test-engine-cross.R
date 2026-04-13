centered_orthonormal_cross <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(scale(M, center = TRUE, scale = FALSE)))
}

orthonormal_basis_cross <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(M))
}

make_exact_cross_covariance_engine <- function(n, p_x, p_y, scale_x, scale_y) {
  k <- length(scale_x)
  stopifnot(length(scale_y) == k)
  T  <- centered_orthonormal_cross(n, k)
  Wx <- orthonormal_basis_cross(p_x, k)
  Wy <- orthonormal_basis_cross(p_y, k)
  list(
    X = T %*% diag(scale_x, nrow = k, ncol = k) %*% t(Wx),
    Y = T %*% diag(scale_y, nrow = k, ncol = k) %*% t(Wy)
  )
}

make_exact_cross_correlation_engine <- function(n, scales_x, scales_y) {
  k <- length(scales_x)
  stopifnot(length(scales_y) == k)
  T  <- centered_orthonormal_cross(n, k)
  Rx <- orthonormal_basis_cross(k, k)
  Ry <- orthonormal_basis_cross(k, k)
  list(
    X = T %*% diag(scales_x, nrow = k, ncol = k) %*% Rx,
    Y = T %*% diag(scales_y, nrow = k, ncol = k) %*% Ry
  )
}

test_that("run_cross_ladder validates recipe geometry", {
  ensure_default_adapters()
  rec_one <- infer_recipe(geometry = "oneblock", relation = "variance",
                          adapter = "prcomp_oneblock")
  X <- matrix(rnorm(40), 10, 4)
  Y <- matrix(rnorm(40), 10, 4)
  expect_error(run_cross_ladder(rec_one, X, Y), "geometry must be")
})

test_that("run_cross_ladder validates input matrices", {
  ensure_default_adapters()
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  expect_error(run_cross_ladder(rec, "not a matrix", matrix(0, 2, 2)),
               "must be a numeric matrix")
  expect_error(
    run_cross_ladder(rec, matrix(0, 5, 3), matrix(0, 4, 3)),
    "same number of rows"
  )
})

test_that("run_cross_ladder runs on covariance recipe and returns expected schema", {
  ensure_default_adapters()
  set.seed(301)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, X, Y, B = 99L, alpha = 0.05, seed = 7)
  expect_true(is.list(res))
  for (nm in c("units", "component_tests", "roots_observed", "ladder_result")) {
    expect_true(nm %in% names(res))
  }
  expect_s3_class(res$units, "multifer_units")
  expect_true(nrow(res$component_tests) >= 1L)
})

test_that("run_cross_ladder runs on correlation recipe", {
  ensure_default_adapters()
  set.seed(302)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "correlation",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, X, Y, B = 99L, alpha = 0.05, seed = 7)
  expect_true(nrow(res$component_tests) >= 1L)
  expect_true(all(res$component_tests$p_value >= 0 &
                  res$component_tests$p_value <= 1))
})

test_that("run_cross_ladder is reproducible under fixed seed", {
  ensure_default_adapters()
  set.seed(303)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  r1 <- run_cross_ladder(rec, X, Y, B = 49L, seed = 11)
  r2 <- run_cross_ladder(rec, X, Y, B = 49L, seed = 11)
  expect_equal(r1$component_tests$p_value, r2$component_tests$p_value)
  expect_equal(r1$component_tests$observed_stat, r2$component_tests$observed_stat)
})

test_that("run_cross_ladder null calibration on bench_cross_null", {
  ensure_default_adapters()
  dat <- bench_cross_null(n = 60, p_x = 8, p_y = 6,
                          within_rank_x = 3, within_rank_y = 3,
                          seed = 777)
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, dat$X, dat$Y, B = 99L, alpha = 0.05, seed = 7)
  # On a single null realization a few false positives are tolerable.
  expect_true(res$ladder_result$rejected_through <= 3L)
})

test_that("run_cross_ladder covariance recovers exact noiseless shared cross rank", {
  ensure_default_adapters()
  set.seed(401)
  dat <- make_exact_cross_covariance_engine(
    n = 28,
    p_x = 6,
    p_y = 5,
    scale_x = c(7, 3, 1.5),
    scale_y = c(5, 2, 0.5)
  )
  rec <- infer_recipe(geometry = "cross", relation = "covariance",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, dat$X, dat$Y, B = 39L, alpha = 0.05, seed = 7L)

  expect_equal(res$ladder_result$rejected_through, 3L)
  expect_equal(res$ladder_result$last_step_tested, 4L)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("run_cross_ladder correlation is validity-capped at the first root", {
  ensure_default_adapters()
  set.seed(402)
  dat <- make_exact_cross_correlation_engine(
    n = 30,
    scales_x = c(5, 2, 1),
    scales_y = c(4, 3, 0.5)
  )
  rec <- infer_recipe(geometry = "cross", relation = "correlation",
                      adapter = "cross_svd")
  res <- run_cross_ladder(rec, dat$X, dat$Y, B = 39L, alpha = 0.05, seed = 7L)

  expect_equal(res$ladder_result$last_step_tested, 1L)
  expect_true(res$ladder_result$rejected_through %in% c(0L, 1L))
  expect_equal(nrow(res$component_tests), 1L)
  expect_true(isTRUE(res$component_tests$selected[1L]))
})
