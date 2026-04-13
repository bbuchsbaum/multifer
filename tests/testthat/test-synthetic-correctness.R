# Synthetic-data correctness tests for the computational core.
#
# These tests lock down method behavior with analytic or differential
# oracles, metamorphic invariants, and degenerate edge cases. They are
# intentionally separate from the existing smoke/contract tests.

centered_orthonormal <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(scale(M, center = TRUE, scale = FALSE)))
}

orthonormal_basis <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(M))
}

make_exact_oneblock <- function(n, p, singular_values) {
  k <- length(singular_values)
  U <- centered_orthonormal(n, k)
  V <- orthonormal_basis(p, k)
  X <- U %*% diag(singular_values, nrow = k, ncol = k) %*% t(V)
  list(X = X, U = U, V = V)
}

make_exact_cross_covariance <- function(n, p_x, p_y, scale_x, scale_y) {
  k <- length(scale_x)
  stopifnot(length(scale_y) == k)

  T  <- centered_orthonormal(n, k)
  Wx <- orthonormal_basis(p_x, k)
  Wy <- orthonormal_basis(p_y, k)

  X <- T %*% diag(scale_x, nrow = k, ncol = k) %*% t(Wx)
  Y <- T %*% diag(scale_y, nrow = k, ncol = k) %*% t(Wy)

  list(X = X, Y = Y, T = T, Wx = Wx, Wy = Wy)
}

make_exact_cross_correlation <- function(n, scales_x, scales_y) {
  k <- length(scales_x)
  stopifnot(length(scales_y) == k)

  T  <- centered_orthonormal(n, k)
  Rx <- orthonormal_basis(k, k)
  Ry <- orthonormal_basis(k, k)

  X <- T %*% diag(scales_x, nrow = k, ncol = k) %*% Rx
  Y <- T %*% diag(scales_y, nrow = k, ncol = k) %*% Ry

  list(X = X, Y = Y, T = T)
}

make_noisy_cross_pair <- function(n, p_x, p_y, latent_scales, noise_sd) {
  k <- length(latent_scales)
  T <- centered_orthonormal(n, k)
  Ax <- matrix(stats::rnorm(k * p_x), nrow = k, ncol = p_x)
  Ay <- matrix(stats::rnorm(k * p_y), nrow = k, ncol = p_y)

  signal_x <- T %*% diag(latent_scales, nrow = k, ncol = k) %*% Ax
  signal_y <- T %*% diag(latent_scales, nrow = k, ncol = k) %*% Ay

  X <- signal_x + matrix(stats::rnorm(n * p_x, sd = noise_sd), nrow = n)
  Y <- signal_y + matrix(stats::rnorm(n * p_y, sd = noise_sd), nrow = n)

  list(
    X = scale(X, center = TRUE, scale = FALSE),
    Y = scale(Y, center = TRUE, scale = FALSE)
  )
}

expect_binomial_fpr <- function(n_false, n_rep, alpha, conf = 0.95) {
  lower <- stats::qbinom((1 - conf) / 2, n_rep, alpha) / n_rep
  upper <- stats::qbinom(1 - (1 - conf) / 2, n_rep, alpha) / n_rep
  observed <- n_false / n_rep

  expect_true(
    observed >= lower && observed <= upper,
    info = sprintf(
      "Observed false-positive rate %.3f fell outside binomial %.0f%% band [%.3f, %.3f] around alpha = %.3f.",
      observed, conf * 100, lower, upper, alpha
    )
  )
}

test_that("oneblock adapters match exact tail-ratio oracles on noiseless low-rank data", {
  set.seed(1001)

  singular_values <- c(9, 4, 1.5)
  dat <- make_exact_oneblock(n = 18, p = 7, singular_values = singular_values)
  X <- dat$X

  expected_roots_svd    <- singular_values^2
  expected_roots_prcomp <- singular_values^2 / (nrow(X) - 1)
  expected_tail_2       <- expected_roots_svd[2] / sum(expected_roots_svd[2:3])

  svd_adapter <- adapter_svd()
  prcomp_adapter <- adapter_prcomp()

  fit_svd <- svd_adapter$refit(NULL, X)
  fit_prcomp <- prcomp_adapter$refit(NULL, X)

  expect_equal(svd_adapter$roots(fit_svd)[1:3], expected_roots_svd, tolerance = 1e-10)
  expect_equal(prcomp_adapter$roots(fit_prcomp)[1:3], expected_roots_prcomp, tolerance = 1e-10)

  expect_equal(
    svd_adapter$component_stat(fit_svd, X, k = 2L),
    expected_tail_2,
    tolerance = 1e-10
  )
  expect_equal(
    prcomp_adapter$component_stat(fit_prcomp, X, k = 2L),
    expected_tail_2,
    tolerance = 1e-10
  )
})

test_that("oneblock residualization and ladder recover exact rank on noiseless data", {
  set.seed(1002)

  singular_values <- c(12, 5)
  dat <- make_exact_oneblock(n = 24, p = 8, singular_values = singular_values)
  X <- dat$X

  svd_adapter <- adapter_svd()
  prcomp_adapter <- adapter_prcomp()

  fit_svd <- svd_adapter$refit(NULL, X)
  fit_prcomp <- prcomp_adapter$refit(NULL, X)

  resid_svd <- svd_adapter$residualize(fit_svd, k = 1L, data = X)
  resid_prcomp <- prcomp_adapter$residualize(fit_prcomp, k = 1L, data = X)

  fit_svd_resid <- svd_adapter$refit(NULL, resid_svd)
  fit_prcomp_resid <- prcomp_adapter$refit(NULL, resid_prcomp)

  expect_equal(svd_adapter$roots(fit_svd_resid)[1L], singular_values[2]^2, tolerance = 1e-10)
  expect_equal(
    prcomp_adapter$roots(fit_prcomp_resid)[1L],
    singular_values[2]^2 / (nrow(X) - 1),
    tolerance = 1e-10
  )
  expect_true(max(abs(svd_adapter$roots(fit_svd_resid)[-1])) < 1e-10)
  expect_true(max(abs(prcomp_adapter$roots(fit_prcomp_resid)[-1])) < 1e-10)

  clear_adapter_registry()
  register_oneblock_baser_adapters()

  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "prcomp_oneblock"
  )

  result <- run_oneblock_ladder(recipe, X, B = 39L, alpha = 0.05, seed = 77L)

  expect_equal(result$ladder_result$rejected_through, 2L)
  expect_equal(result$ladder_result$last_step_tested, 3L)
  expect_identical(result$component_tests$selected, c(TRUE, TRUE, FALSE))
  expect_equal(result$component_tests$observed_stat[1], 144 / 169, tolerance = 1e-10)
  expect_equal(result$component_tests$observed_stat[2], 1, tolerance = 1e-10)
  expect_equal(result$component_tests$observed_stat[3], 0, tolerance = 1e-10)
  expect_equal(result$component_tests$p_value[1:2], rep(1 / 40, 2), tolerance = 1e-12)
  expect_equal(result$component_tests$p_value[3], 1)
})

test_that("oneblock ladder is invariant to column permutation and handles the zero matrix", {
  set.seed(1003)

  dat <- make_exact_oneblock(n = 20, p = 9, singular_values = c(10, 3))
  X <- dat$X
  X_perm <- X[, c(4, 1, 8, 2, 9, 3, 5, 7, 6), drop = FALSE]

  clear_adapter_registry()
  register_oneblock_baser_adapters()

  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "svd_oneblock"
  )

  res_a <- run_oneblock_ladder(recipe, X, B = 29L, alpha = 0.05, seed = 404L)
  res_b <- run_oneblock_ladder(recipe, X_perm, B = 29L, alpha = 0.05, seed = 404L)

  expect_equal(res_a$ladder_result$rejected_through, 2L)
  expect_equal(res_b$ladder_result$rejected_through, 2L)
  expect_equal(res_a$roots_observed, res_b$roots_observed, tolerance = 1e-10)
  expect_equal(res_a$component_tests$observed_stat, res_b$component_tests$observed_stat, tolerance = 1e-12)
  expect_identical(res_a$component_tests$selected, c(TRUE, TRUE, FALSE))
  expect_identical(res_b$component_tests$selected, c(TRUE, TRUE, FALSE))

  zero_res <- run_oneblock_ladder(
    recipe,
    matrix(0, nrow = 12, ncol = 5),
    B = 19L,
    alpha = 0.05,
    seed = 11L
  )

  expect_equal(zero_res$ladder_result$rejected_through, 0L)
  expect_equal(zero_res$ladder_result$last_step_tested, 1L)
  expect_equal(zero_res$component_tests$observed_stat, 0)
  expect_equal(zero_res$component_tests$p_value, 1)
  expect_true(all(zero_res$ladder_result$step_results[[1]]$null_values == 0))
})

test_that("oneblock adapters obey scaling laws and preserve column marginals under null action", {
  set.seed(10035)

  dat <- make_exact_oneblock(n = 22, p = 6, singular_values = c(8, 3, 1))
  X <- dat$X
  X_scaled <- 2.5 * X

  for (adapter in list(adapter_svd(), adapter_prcomp())) {
    fit <- adapter$refit(NULL, X)
    fit_scaled <- adapter$refit(NULL, X_scaled)

    expect_equal(adapter$roots(fit_scaled), 2.5^2 * adapter$roots(fit), tolerance = 1e-10)
    expect_equal(
      adapter$component_stat(fit, X, k = 2L),
      adapter$component_stat(fit_scaled, X_scaled, k = 2L),
      tolerance = 1e-12
    )
    expect_true(is.na(adapter$component_stat(fit, X, k = 7L)))

    set.seed(91L)
    X_null <- adapter$null_action(fit, X)

    expect_equal(dim(X_null), dim(X))
    expect_equal(colSums(X_null), colSums(X), tolerance = 1e-12)
    for (j in seq_len(ncol(X))) {
      expect_equal(sort(X_null[, j]), sort(X[, j]), tolerance = 1e-12)
    }
  }
})

test_that("cross covariance adapter matches exact singular-value and deflation oracles", {
  set.seed(1004)

  scale_x <- c(7, 3, 1.5)
  scale_y <- c(5, 2, 0.5)
  dat <- make_exact_cross_covariance(
    n = 28,
    p_x = 6,
    p_y = 5,
    scale_x = scale_x,
    scale_y = scale_y
  )

  expected_singular_values <- scale_x * scale_y
  expected_tail_2 <- expected_singular_values[2]^2 / sum(expected_singular_values[2:3]^2)

  adapter <- adapter_cross_svd()
  fit <- adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))

  expect_equal(fit$d[1:3], expected_singular_values, tolerance = 1e-10)
  expect_equal(
    adapter$component_stat(fit, list(X = dat$X, Y = dat$Y), k = 2L),
    expected_tail_2,
    tolerance = 1e-10
  )

  resid <- adapter$residualize(fit, k = 1L, data = list(X = dat$X, Y = dat$Y))
  fit_resid <- adapter$refit(fit, resid)

  expect_equal(fit_resid$d[1:2], expected_singular_values[2:3], tolerance = 1e-10)
  expect_true(max(abs(fit_resid$d[-(1:2)])) < 1e-10)
})

test_that("cross correlation adapters recover exact unit canonical correlations", {
  set.seed(1005)

  dat <- make_exact_cross_correlation(
    n = 30,
    scales_x = c(5, 2, 1),
    scales_y = c(4, 3, 0.5)
  )

  cross_adapter <- adapter_cross_svd()
  cancor_adapter <- adapter_cancor()

  fit_cross <- cross_adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "correlation"))
  fit_cancor <- cancor_adapter$refit(NULL, list(X = dat$X, Y = dat$Y))

  expect_equal(fit_cross$d[1:3], c(1, 1, 1), tolerance = 1e-10)
  expect_equal(fit_cancor$cor[1:3], c(1, 1, 1), tolerance = 1e-10)
  expect_equal(fit_cross$d, fit_cancor$cor, tolerance = 1e-10)

  expect_equal(
    cross_adapter$component_stat(fit_cross, list(X = dat$X, Y = dat$Y), k = 2L),
    0.5,
    tolerance = 1e-10
  )
  expect_equal(
    cancor_adapter$component_stat(fit_cancor, list(X = dat$X, Y = dat$Y), k = 2L),
    0.5,
    tolerance = 1e-10
  )
})

test_that("cross covariance and correlation paths obey distinct scaling laws", {
  set.seed(10055)

  dat <- make_exact_cross_correlation(
    n = 32,
    scales_x = c(6, 3, 1),
    scales_y = c(5, 2, 0.5)
  )

  x_scale <- c(2, 0.5, 3)
  y_scale <- c(4, 1.5, 0.25)

  X_scaled <- sweep(dat$X, 2L, x_scale, "*")
  Y_scaled <- sweep(dat$Y, 2L, y_scale, "*")

  cross_adapter <- adapter_cross_svd()
  cancor_adapter <- adapter_cancor()

  fit_cov <- cross_adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))
  fit_cov_global <- cross_adapter$refit(NULL, list(X = 3 * dat$X, Y = 0.2 * dat$Y, relation = "covariance"))
  fit_corr <- cross_adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "correlation"))
  fit_corr_scaled <- cross_adapter$refit(NULL, list(X = X_scaled, Y = Y_scaled, relation = "correlation"))
  fit_cancor <- cancor_adapter$refit(NULL, list(X = dat$X, Y = dat$Y))
  fit_cancor_scaled <- cancor_adapter$refit(NULL, list(X = X_scaled, Y = Y_scaled))

  expect_equal(fit_cov_global$d, 0.6 * fit_cov$d, tolerance = 1e-10)
  expect_equal(
    cross_adapter$component_stat(fit_cov, list(X = dat$X, Y = dat$Y), k = 2L),
    cross_adapter$component_stat(fit_cov_global, list(X = 3 * dat$X, Y = 0.2 * dat$Y), k = 2L),
    tolerance = 1e-12
  )

  expect_equal(fit_corr_scaled$d, fit_corr$d, tolerance = 1e-10)
  expect_equal(fit_cancor_scaled$cor, fit_cancor$cor, tolerance = 1e-10)
  expect_equal(fit_corr_scaled$d, fit_cancor_scaled$cor, tolerance = 1e-10)
})

test_that("cross null action preserves within-block geometry while breaking cross pairing", {
  set.seed(10065)

  dat <- make_exact_cross_covariance(
    n = 26,
    p_x = 5,
    p_y = 4,
    scale_x = c(9, 2),
    scale_y = c(7, 1)
  )

  adapter <- adapter_cross_svd()
  fit <- adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))

  set.seed(13L)
  permuted <- adapter$null_action(fit, list(X = dat$X, Y = dat$Y))
  fit_permuted <- adapter$refit(fit, c(permuted, list(relation = "covariance")))

  expect_equal(permuted$X, dat$X, tolerance = 0)
  expect_equal(crossprod(permuted$Y), crossprod(dat$Y), tolerance = 1e-12)
  expect_equal(colSums(permuted$Y), colSums(dat$Y), tolerance = 1e-12)
  expect_true(fit_permuted$d[1L] < fit$d[1L])
  expect_equal(
    sort(rowSums(permuted$Y^2)),
    sort(rowSums(dat$Y^2)),
    tolerance = 1e-12
  )
})

test_that("cross fits are invariant to common row permutation and correlation path matches cancor", {
  set.seed(1006)

  dat <- make_noisy_cross_pair(
    n = 40,
    p_x = 6,
    p_y = 5,
    latent_scales = c(4, 2, 1),
    noise_sd = 0.15
  )
  perm <- c(11, 5, 27, 2, 19, 33, 8, 24, 14, 1,
            36, 6, 29, 18, 40, 9, 3, 21, 13, 31,
            7, 28, 16, 38, 10, 25, 4, 34, 17, 12,
            35, 26, 20, 39, 15, 23, 32, 30, 22, 37)

  cross_adapter <- adapter_cross_svd()
  cancor_adapter <- adapter_cancor()

  fit_cov <- cross_adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))
  fit_cov_perm <- cross_adapter$refit(
    NULL,
    list(X = dat$X[perm, , drop = FALSE], Y = dat$Y[perm, , drop = FALSE], relation = "covariance")
  )

  fit_corr <- cross_adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "correlation"))
  fit_corr_perm <- cross_adapter$refit(
    NULL,
    list(X = dat$X[perm, , drop = FALSE], Y = dat$Y[perm, , drop = FALSE], relation = "correlation")
  )
  fit_cancor <- cancor_adapter$refit(NULL, list(X = dat$X, Y = dat$Y))

  expect_equal(fit_cov$d, fit_cov_perm$d, tolerance = 1e-10)
  expect_equal(fit_corr$d, fit_corr_perm$d, tolerance = 1e-10)
  expect_equal(fit_corr$d, fit_cancor$cor, tolerance = 1e-8)

  expect_equal(
    cross_adapter$component_stat(fit_cov, list(X = dat$X, Y = dat$Y), k = 1L),
    cross_adapter$component_stat(
      fit_cov_perm,
      list(X = dat$X[perm, , drop = FALSE], Y = dat$Y[perm, , drop = FALSE]),
      k = 1L
    ),
    tolerance = 1e-10
  )
  expect_equal(
    cross_adapter$component_stat(fit_corr, list(X = dat$X, Y = dat$Y), k = 2L),
    cancor_adapter$component_stat(fit_cancor, list(X = dat$X, Y = dat$Y), k = 2L),
    tolerance = 1e-8
  )
})

test_that("oneblock ladder empirical false-positive rate is controlled under Gaussian null", {
  clear_adapter_registry()
  register_oneblock_baser_adapters()

  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "prcomp_oneblock"
  )

  alpha <- 0.05
  n_rep <- 80L
  false_positives <- integer(n_rep)

  for (i in seq_len(n_rep)) {
    X <- bench_oneblock_null(n = 60, p = 12, noise = "gaussian", seed = 5000L + i)$X
    res <- run_oneblock_ladder(recipe, X, B = 39L, alpha = alpha, seed = 9000L + i)
    false_positives[i] <- as.integer(res$ladder_result$rejected_through > 0L)
  }

  expect_binomial_fpr(sum(false_positives), n_rep, alpha)
})

test_that("oneblock ladder empirical false-positive rate is controlled under heteroscedastic null", {
  clear_adapter_registry()
  register_oneblock_baser_adapters()

  recipe <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "prcomp_oneblock"
  )

  alpha <- 0.05
  n_rep <- 80L
  false_positives <- integer(n_rep)

  for (i in seq_len(n_rep)) {
    X <- bench_oneblock_null(n = 60, p = 12, noise = "heteroscedastic", seed = 6000L + i)$X
    res <- run_oneblock_ladder(recipe, X, B = 39L, alpha = alpha, seed = 10000L + i)
    false_positives[i] <- as.integer(res$ladder_result$rejected_through > 0L)
  }

  expect_binomial_fpr(sum(false_positives), n_rep, alpha)
})

test_that("cross covariance ladder empirical false-positive rate is controlled under null", {
  clear_adapter_registry()
  register_cross_baser_adapters()

  recipe <- infer_recipe(
    geometry = "cross",
    relation = "covariance",
    adapter = "cross_svd"
  )

  alpha <- 0.05
  n_rep <- 60L
  false_positives <- integer(n_rep)

  for (i in seq_len(n_rep)) {
    dat <- bench_cross_null(
      n = 50,
      p_x = 8,
      p_y = 6,
      within_rank_x = 3,
      within_rank_y = 3,
      seed = 7000L + i
    )
    res <- run_cross_ladder(recipe, dat$X, dat$Y, B = 39L, alpha = alpha, seed = 11000L + i)
    false_positives[i] <- as.integer(res$ladder_result$rejected_through > 0L)
  }

  expect_binomial_fpr(sum(false_positives), n_rep, alpha)
})
