make_cross_bootstrap_fixture <- function(n, p_x, p_y, signal, noise_x, noise_y) {
  k <- length(signal)
  latent <- qr.Q(qr(scale(matrix(rnorm(n * k), nrow = n, ncol = k),
                            center = TRUE, scale = FALSE)))
  Wx <- qr.Q(qr(matrix(rnorm(p_x * k), nrow = p_x, ncol = k)))
  Wy <- qr.Q(qr(matrix(rnorm(p_y * k), nrow = p_y, ncol = k)))

  X <- latent %*% diag(signal, nrow = k, ncol = k) %*% t(Wx) +
    matrix(rnorm(n * p_x, sd = noise_x), nrow = n, ncol = p_x)
  Y <- latent %*% diag(signal, nrow = k, ncol = k) %*% t(Wy) +
    matrix(rnorm(n * p_y, sd = noise_y), nrow = n, ncol = p_y)

  list(
    X = scale(X, center = TRUE, scale = FALSE),
    Y = scale(Y, center = TRUE, scale = FALSE)
  )
}

test_that("adapter_cross_svd$core returns block SVDs for covariance mode", {
  set.seed(1)
  X <- matrix(rnorm(30 * 5), 30, 5)
  Y <- matrix(rnorm(30 * 4), 30, 4)
  adapter <- adapter_cross_svd()
  fit <- adapter$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  core <- adapter$core(fit, list(X = X, Y = Y, relation = "covariance"))

  expect_equal(core$relation, "covariance")
  expect_equal(dim(core$Ux), c(30L, 5L))
  expect_equal(dim(core$Uy), c(30L, 4L))
  expect_equal(dim(core$Vx), c(5L, 5L))
  expect_equal(dim(core$Vy), c(4L, 4L))
  expect_equal(length(core$dx), 5L)
  expect_equal(length(core$dy), 4L)
})

test_that("cross covariance update_core reproduces refit bootstrap to 1e-10", {
  set.seed(2)
  n <- 40; px <- 6; py <- 5
  X <- matrix(rnorm(n * px), n, px)
  Y <- matrix(rnorm(n * py), n, py)
  adapter <- adapter_cross_svd()
  fit <- adapter$refit(NULL, list(X = X, Y = Y, relation = "covariance"))
  core <- adapter$core(fit, list(X = X, Y = Y, relation = "covariance"))

  set.seed(111)
  idx <- sample.int(n, n, replace = TRUE)
  X_boot <- X[idx, , drop = FALSE]
  Y_boot <- Y[idx, , drop = FALSE]

  ref <- adapter$refit(
    NULL, list(X = X_boot, Y = Y_boot, relation = "covariance")
  )
  upd <- adapter$update_core(core, indices = idx)

  # Singular values agree.
  expect_equal(sort(upd$d, decreasing = TRUE),
               sort(ref$d, decreasing = TRUE),
               tolerance = 1e-10)

  # Outer product of the cross-product SVD is sign-invariant.
  recon_ref <- ref$Wx %*% diag(ref$d) %*% t(ref$Wy)
  recon_upd <- upd$Wx %*% diag(upd$d) %*% t(upd$Wy)
  expect_equal(recon_upd, recon_ref, tolerance = 1e-10)

  # Centers agree.
  expect_equal(upd$center_x, ref$center_x, tolerance = 1e-10)
  expect_equal(upd$center_y, ref$center_y, tolerance = 1e-10)

  # Scores: Tx_new = X_boot_c %*% Wx_new. Verify against fresh computation.
  X_boot_c <- sweep(X_boot, 2L, ref$center_x, "-")
  Y_boot_c <- sweep(Y_boot, 2L, ref$center_y, "-")
  Tx_direct <- X_boot_c %*% upd$Wx
  Ty_direct <- Y_boot_c %*% upd$Wy
  expect_equal(upd$Tx, Tx_direct, tolerance = 1e-10)
  expect_equal(upd$Ty, Ty_direct, tolerance = 1e-10)
})

test_that("cross covariance refit trims exact zero tail to active numerical rank", {
  set.seed(21)
  dat <- make_cross_bootstrap_fixture(
    n = 60L, p_x = 80L, p_y = 80L,
    signal = c(4.0, 3.6), noise_x = 1.5, noise_y = 1.5
  )
  adapter <- adapter_cross_svd()
  fit <- adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))

  expect_equal(length(fit$d), nrow(dat$X) - 1L)
  expect_true(all(fit$d > 0))
  expect_equal(ncol(fit$Wx), length(fit$d))
  expect_equal(ncol(fit$Wy), length(fit$d))
  expect_equal(ncol(fit$Tx), length(fit$d))
  expect_equal(ncol(fit$Ty), length(fit$d))
})

test_that("bootstrap_fits uses fast path for cross covariance", {
  set.seed(3)
  X <- matrix(rnorm(40 * 6), 40, 6)
  Y <- matrix(rnorm(40 * 5), 40, 5)
  adapter <- adapter_cross_svd()
  recipe <- infer_recipe(
    geometry = "cross", relation = "covariance",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(
    NULL, list(X = X, Y = Y, relation = "covariance")
  )
  units <- form_units(original_fit$d^2)

  artifact <- bootstrap_fits(
    recipe = recipe, adapter = adapter,
    data = list(X = X, Y = Y),
    original_fit = original_fit, units = units,
    R = 6L, method_align = "sign", seed = 31
  )
  expect_true(isTRUE(artifact$used_fast_path))
  expect_equal(length(artifact$reps), 6L)

  # Each rep's fit should have the covariance fit shape.
  fit_keys <- c("relation", "d", "Wx", "Wy", "Tx", "Ty", "center_x", "center_y")
  expect_true(all(vapply(
    artifact$reps,
    function(r) all(fit_keys %in% names(r$fit)),
    logical(1L)
  )))
})

test_that("cross covariance fast-path bootstrap matches refit artifact and stability summaries", {
  scenarios <- list(
    list(name = "balanced", n = 48L, p_x = 7L, p_y = 6L,
         signal = c(5, 3, 1.5), noise_x = 0.8, noise_y = 0.9),
    list(name = "wide_y_noisy_y", n = 60L, p_x = 5L, p_y = 11L,
         signal = c(4, 2.7), noise_x = 0.7, noise_y = 1.8)
  )

  adapter <- adapter_cross_svd()
  recipe <- infer_recipe(
    geometry = "cross", relation = "covariance",
    adapter = adapter, strict = TRUE
  )

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    set.seed(100 + i)
    dat <- make_cross_bootstrap_fixture(
      n = sc$n, p_x = sc$p_x, p_y = sc$p_y,
      signal = sc$signal, noise_x = sc$noise_x, noise_y = sc$noise_y
    )

    original_fit <- adapter$refit(
      NULL, list(X = dat$X, Y = dat$Y, relation = "covariance")
    )
    units <- form_units(adapter$roots(original_fit))

    fast <- bootstrap_fits(
      recipe = recipe, adapter = adapter,
      data = list(X = dat$X, Y = dat$Y),
      original_fit = original_fit, units = units,
      R = 8L, method_align = "sign", seed = 700L + i,
      fast_path = "auto", core_rank = NULL
    )
    slow <- bootstrap_fits(
      recipe = recipe, adapter = adapter,
      data = list(X = dat$X, Y = dat$Y),
      original_fit = original_fit, units = units,
      R = 8L, method_align = "sign", seed = 700L + i,
      fast_path = "off", core_rank = NULL
    )

    expect_true(isTRUE(fast$used_fast_path), info = sc$name)
    expect_false(isTRUE(slow$used_fast_path), info = sc$name)

    for (b in seq_len(fast$R)) {
      rep_fast <- fast$reps[[b]]
      rep_slow <- slow$reps[[b]]

      expect_equal(rep_fast$resample_indices, rep_slow$resample_indices,
                   info = sprintf("%s rep %d indices", sc$name, b))
      expect_equal(rep_fast$fit$d, rep_slow$fit$d, tolerance = 1e-10,
                   info = sprintf("%s rep %d singular values", sc$name, b))
      expect_equal(rep_fast$fit$center_x, rep_slow$fit$center_x, tolerance = 1e-10,
                   info = sprintf("%s rep %d center_x", sc$name, b))
      expect_equal(rep_fast$fit$center_y, rep_slow$fit$center_y, tolerance = 1e-10,
                   info = sprintf("%s rep %d center_y", sc$name, b))

      recon_fast <- rep_fast$fit$Wx %*% diag(rep_fast$fit$d, nrow = length(rep_fast$fit$d)) %*%
        t(rep_fast$fit$Wy)
      recon_slow <- rep_slow$fit$Wx %*% diag(rep_slow$fit$d, nrow = length(rep_slow$fit$d)) %*%
        t(rep_slow$fit$Wy)
      expect_equal(recon_fast, recon_slow, tolerance = 1e-10,
                   info = sprintf("%s rep %d reconstructed cross operator", sc$name, b))

      expect_equal(rep_fast$aligned_loadings$X, rep_slow$aligned_loadings$X, tolerance = 1e-10,
                   info = sprintf("%s rep %d aligned X loadings", sc$name, b))
      expect_equal(rep_fast$aligned_loadings$Y, rep_slow$aligned_loadings$Y, tolerance = 1e-10,
                   info = sprintf("%s rep %d aligned Y loadings", sc$name, b))
      expect_equal(rep_fast$aligned_scores$X, rep_slow$aligned_scores$X, tolerance = 1e-10,
                   info = sprintf("%s rep %d aligned X scores", sc$name, b))
      expect_equal(rep_fast$aligned_scores$Y, rep_slow$aligned_scores$Y, tolerance = 1e-10,
                   info = sprintf("%s rep %d aligned Y scores", sc$name, b))
    }

    var_fast <- variable_stability_from_bootstrap(fast, units)
    var_slow <- variable_stability_from_bootstrap(slow, units)
    expect_equal(var_fast, var_slow, tolerance = 1e-10,
                 info = sprintf("%s variable stability", sc$name))

    score_fast <- score_stability_from_bootstrap(fast, dat, units)
    score_slow <- score_stability_from_bootstrap(slow, dat, units)
    expect_equal(score_fast, score_slow, tolerance = 1e-10,
                 info = sprintf("%s score stability", sc$name))

    sub_fast <- subspace_stability_from_bootstrap(fast, original_fit, adapter, units)
    sub_slow <- subspace_stability_from_bootstrap(slow, original_fit, adapter, units)
    expect_equal(sub_fast, sub_slow, tolerance = 1e-10,
                 info = sprintf("%s subspace stability", sc$name))
  }
})

test_that("cross covariance fast-path matches refit on leading signal units in small-n wide near-tied regime", {
  set.seed(212)
  dat <- make_cross_bootstrap_fixture(
    n = 60L, p_x = 80L, p_y = 80L,
    signal = c(4.0, 3.6), noise_x = 1.5, noise_y = 1.5
  )

  adapter <- adapter_cross_svd()
  recipe <- infer_recipe(
    geometry = "cross", relation = "covariance",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(
    NULL, list(X = dat$X, Y = dat$Y, relation = "covariance")
  )
  signal_units <- form_units(adapter$roots(original_fit)[1:2])

  fast <- bootstrap_fits(
    recipe = recipe, adapter = adapter,
    data = list(X = dat$X, Y = dat$Y),
    original_fit = original_fit, units = signal_units,
    R = 8L, method_align = "sign", seed = 812L,
    fast_path = "auto", core_rank = NULL
  )
  slow <- bootstrap_fits(
    recipe = recipe, adapter = adapter,
    data = list(X = dat$X, Y = dat$Y),
    original_fit = original_fit, units = signal_units,
    R = 8L, method_align = "sign", seed = 812L,
    fast_path = "off", core_rank = NULL
  )

  for (b in seq_len(fast$R)) {
    expect_equal(fast$reps[[b]]$resample_indices, slow$reps[[b]]$resample_indices)
    expect_equal(length(fast$reps[[b]]$fit$d), length(slow$reps[[b]]$fit$d))
    expect_equal(fast$reps[[b]]$fit$d[1:2], slow$reps[[b]]$fit$d[1:2], tolerance = 1e-10)
    expect_equal(
      fast$reps[[b]]$aligned_loadings$X[, 1:2, drop = FALSE],
      slow$reps[[b]]$aligned_loadings$X[, 1:2, drop = FALSE],
      tolerance = 1e-10
    )
    expect_equal(
      fast$reps[[b]]$aligned_loadings$Y[, 1:2, drop = FALSE],
      slow$reps[[b]]$aligned_loadings$Y[, 1:2, drop = FALSE],
      tolerance = 1e-10
    )
    expect_equal(
      fast$reps[[b]]$aligned_scores$X[, 1:2, drop = FALSE],
      slow$reps[[b]]$aligned_scores$X[, 1:2, drop = FALSE],
      tolerance = 1e-10
    )
    expect_equal(
      fast$reps[[b]]$aligned_scores$Y[, 1:2, drop = FALSE],
      slow$reps[[b]]$aligned_scores$Y[, 1:2, drop = FALSE],
      tolerance = 1e-10
    )
  }

  expect_equal(
    variable_stability_from_bootstrap(fast, signal_units),
    variable_stability_from_bootstrap(slow, signal_units),
    tolerance = 1e-10
  )
  expect_equal(
    score_stability_from_bootstrap(fast, dat, signal_units),
    score_stability_from_bootstrap(slow, dat, signal_units),
    tolerance = 1e-10
  )
  expect_equal(
    subspace_stability_from_bootstrap(fast, original_fit, adapter, signal_units),
    subspace_stability_from_bootstrap(slow, original_fit, adapter, signal_units),
    tolerance = 1e-10
  )
})

test_that("bootstrap_fits falls back to refit for cross correlation", {
  set.seed(4)
  X <- matrix(rnorm(40 * 6), 40, 6)
  Y <- matrix(rnorm(40 * 5), 40, 5)
  adapter <- adapter_cross_svd()
  recipe <- infer_recipe(
    geometry = "cross", relation = "correlation",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(
    NULL, list(X = X, Y = Y, relation = "correlation")
  )
  units <- form_units(original_fit$d^2)

  artifact <- bootstrap_fits(
    recipe = recipe, adapter = adapter,
    data = list(X = X, Y = Y),
    original_fit = original_fit, units = units,
    R = 4L, method_align = "sign", seed = 37
  )
  expect_false(isTRUE(artifact$used_fast_path))
})

test_that("infer(cross, covariance) reports core_updates in $cost", {
  skip_on_cran()
  ensure_default_adapters()
  set.seed(5)
  X <- matrix(rnorm(40 * 6), 40, 6)
  Y <- matrix(rnorm(40 * 5), 40, 5)
  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = X, Y = Y),
    geometry = "cross", relation = "covariance",
    B = 49L, R = 5L, alpha = 0.05, seed = 41
  )
  expect_true(is_infer_result(res))
  expect_equal(res$cost$core_updates, 5L)
})

test_that("infer(cross, correlation) does NOT report core_updates", {
  skip_on_cran()
  ensure_default_adapters()
  set.seed(6)
  X <- matrix(rnorm(40 * 6), 40, 6)
  Y <- matrix(rnorm(40 * 5), 40, 5)
  res <- infer(
    adapter  = "cross_svd",
    data     = list(X = X, Y = Y),
    geometry = "cross", relation = "correlation",
    B = 49L, R = 5L, alpha = 0.05, seed = 43
  )
  expect_equal(res$cost$core_updates, 0L)
})
