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

make_exact_canonical_correlation_engine <- function(n, canonical_corrs,
                                                    scale_x = NULL,
                                                    scale_y = NULL,
                                                    extra_scale_x = 1,
                                                    extra_scale_y = 1) {
  k <- length(canonical_corrs)
  if (is.null(scale_x)) scale_x <- rep(1, k)
  if (is.null(scale_y)) scale_y <- rep(1, k)
  stopifnot(length(scale_x) == k, length(scale_y) == k)

  B <- centered_orthonormal_cross(n, 2L * k + 2L)
  A <- B[, seq_len(k), drop = FALSE]
  C <- B[, k + seq_len(k), drop = FALSE]
  extra_x <- B[, 2L * k + 1L, drop = FALSE]
  extra_y <- B[, 2L * k + 2L, drop = FALSE]

  Qx_shared <- A
  Qy_shared <- A %*% diag(canonical_corrs, nrow = k, ncol = k) +
    C %*% diag(sqrt(1 - canonical_corrs^2), nrow = k, ncol = k)

  X <- cbind(Qx_shared %*% diag(scale_x, nrow = k, ncol = k),
             extra_scale_x * extra_x)
  Y <- cbind(Qy_shared %*% diag(scale_y, nrow = k, ncol = k),
             extra_scale_y * extra_y)

  list(X = X, Y = Y, canonical_corrs = c(canonical_corrs, 0))
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

test_that("run_cross_ladder correlation recovers exact shared canonical rank under paired rows", {
  ensure_default_adapters()
  set.seed(402)
  dat <- make_exact_canonical_correlation_engine(
    n = 36,
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5),
    extra_scale_x = 1.5,
    extra_scale_y = 1.2
  )
  rec <- infer_recipe(geometry = "cross", relation = "correlation",
                      adapter = "cross_svd")
  res <- run_cross_ladder(
    rec, dat$X, dat$Y,
    B = 39L, alpha = 0.05, seed = 7L, max_steps = 4L
  )

  expect_equal(res$ladder_result$rejected_through, 3L)
  expect_equal(res$ladder_result$last_step_tested, 4L)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("run_cross_ladder correlation remains conservative for nuisance-adjusted designs", {
  ensure_default_adapters()
  set.seed(403)
  dat <- make_exact_canonical_correlation_engine(
    n = 36,
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5),
    extra_scale_x = 1.5,
    extra_scale_y = 1.2
  )
  Z <- cbind(1, scale(seq_len(nrow(dat$X))))
  rec <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = nuisance_adjusted(Z),
    adapter = "cross_svd"
  )
  res <- run_cross_ladder(
    rec, dat$X, dat$Y,
    B = 39L, alpha = 0.05, seed = 7L, max_steps = 4L
  )

  expect_equal(res$ladder_result$last_step_tested, 1L)
  expect_true(res$ladder_result$rejected_through %in% c(0L, 1L))
  expect_equal(nrow(res$component_tests), 1L)
  expect_true(isTRUE(res$component_tests$selected[1L]))
})

test_that("run_cross_ladder correlation controls false positives on replicated null draws", {
  ensure_default_adapters()
  rec <- infer_recipe(geometry = "cross", relation = "correlation",
                      adapter = "cross_svd")

  n_rep <- 40L
  hits <- integer(n_rep)
  for (i in seq_len(n_rep)) {
    dat <- bench_cross_null(
      n = 50, p_x = 6, p_y = 6,
      within_rank_x = 3, within_rank_y = 3,
      seed = 800 + i
    )
    res <- run_cross_ladder(rec, dat$X, dat$Y, B = 39L, alpha = 0.05, seed = i)
    hits[i] <- as.integer(res$ladder_result$rejected_through > 0L)
  }

  fpr <- mean(hits)
  expect_true(fpr >= 0 && fpr <= 0.15, info = sprintf("empirical FPR = %.3f", fpr))
})
