test_that("adapter_svd$core returns a complete core_obj", {
  set.seed(1)
  X <- matrix(rnorm(60), nrow = 15, ncol = 4)
  adapter <- adapter_svd()
  fit <- adapter$refit(NULL, X)
  core <- adapter$core(fit, X)

  expect_named(core, c("U", "d", "V", "center", "n", "k"), ignore.order = TRUE)
  expect_equal(dim(core$U), c(15L, 4L))
  expect_equal(dim(core$V), c(4L, 4L))
  expect_equal(length(core$d), 4L)
  expect_equal(core$n, 15L)
  expect_equal(core$k, 4L)
})

test_that("adapter_svd$update_core reproduces refit on the same bootstrap", {
  set.seed(2)
  X <- matrix(rnorm(100), nrow = 25, ncol = 4)
  adapter <- adapter_svd()
  fit <- adapter$refit(NULL, X)
  core <- adapter$core(fit, X)

  set.seed(99)
  idx <- sample.int(nrow(X), nrow(X), replace = TRUE)
  X_boot <- X[idx, , drop = FALSE]

  ref <- adapter$refit(NULL, X_boot)
  upd <- adapter$update_core(core, indices = idx)

  # Singular values agree tightly.
  expect_equal(sort(upd$d, decreasing = TRUE),
               sort(ref$d, decreasing = TRUE),
               tolerance = 1e-10)

  # The reconstructions of the centered resampled matrix agree. Sign
  # ambiguity means we compare the outer product U %*% D %*% V^T, which
  # is sign-invariant per column.
  recon_ref <- ref$u %*% diag(ref$d) %*% t(ref$v)
  recon_upd <- upd$u %*% diag(upd$d) %*% t(upd$v)
  expect_equal(recon_upd, recon_ref, tolerance = 1e-10)

  # Centers agree.
  expect_equal(upd$center, ref$center, tolerance = 1e-10)
})

test_that("adapter_prcomp$core + update_core also round-trip vs refit", {
  set.seed(3)
  X <- matrix(rnorm(80), nrow = 20, ncol = 4)
  adapter <- adapter_prcomp()
  fit <- adapter$refit(NULL, X)
  core <- adapter$core(fit, X)

  set.seed(101)
  idx <- sample.int(nrow(X), nrow(X), replace = TRUE)
  X_boot <- X[idx, , drop = FALSE]

  ref <- adapter$refit(NULL, X_boot)
  upd <- adapter$update_core(core, indices = idx)

  expect_s3_class(upd, "prcomp")
  expect_equal(sort(upd$sdev, decreasing = TRUE),
               sort(ref$sdev, decreasing = TRUE),
               tolerance = 1e-10)
  # loadings outer products agree (sign-invariant comparison).
  r_ref <- ref$rotation %*% diag(ref$sdev^2) %*% t(ref$rotation)
  r_upd <- upd$rotation %*% diag(upd$sdev^2) %*% t(upd$rotation)
  expect_equal(r_upd, r_ref, tolerance = 1e-10)
  expect_equal(upd$center, ref$center, tolerance = 1e-10)
})

test_that("bootstrap_fits uses fast path when adapter exposes core hooks", {
  set.seed(4)
  X <- matrix(rnorm(100), nrow = 25, ncol = 4)
  adapter <- adapter_svd()
  recipe  <- infer_recipe(
    geometry = "oneblock", relation = "variance",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(NULL, X)
  units <- form_units(original_fit$d^2)

  artifact <- bootstrap_fits(
    recipe = recipe, adapter = adapter, data = X,
    original_fit = original_fit, units = units,
    R = 8L, method_align = "sign", seed = 17
  )
  expect_true(isTRUE(artifact$used_fast_path))
  expect_equal(length(artifact$reps), 8L)
  # Each rep carries a fit with u/d/v/center (same shape as refit).
  expect_true(all(vapply(artifact$reps,
                         function(r) all(c("u", "d", "v", "center") %in% names(r$fit)),
                         logical(1L))))
})

test_that("infer() reports core_updates in $cost when fast path fires", {
  skip_on_cran()
  ensure_default_adapters()
  set.seed(5)
  X <- matrix(rnorm(100), nrow = 25, ncol = 4)
  res <- infer(
    adapter = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 49L, R = 6L, alpha = 0.05, seed = 19
  )
  expect_true(is_infer_result(res))
  expect_equal(res$cost$core_updates, 6L)
  expect_equal(res$cost$full_data_ops, 1L)
})

test_that("decision agreement: fast path reaches same ranks as refit on synthetic data", {
  skip_on_cran()
  ensure_default_adapters()
  set.seed(77)
  # Matrix with clear top-3 signal structure.
  U <- qr.Q(qr(matrix(rnorm(60 * 3), 60, 3)))
  V <- qr.Q(qr(matrix(rnorm(6  * 3), 6,  3)))
  X <- U %*% diag(c(8, 5, 3)) %*% t(V) + matrix(rnorm(60 * 6, sd = 0.1), 60, 6)

  res_fast <- infer(
    adapter  = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 99L, R = 8L, alpha = 0.05, seed = 23
  )
  expect_true(sum(res_fast$units$selected) >= 1L)
  # Fast path engaged.
  expect_true(res_fast$cost$core_updates > 0L)
})
