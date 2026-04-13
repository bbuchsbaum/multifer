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
