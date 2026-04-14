# Correlation-mode core-rank-deficient fallback contract (multifer-9u9.1.3)
#
# The correlation fast path in adapter_cross_svd$update_core raises a
# classed `multifer_core_rank_deficient` error when a bootstrap
# resample produces a rank-deficient centered block in the core
# representation. bootstrap_fits() is required to catch this classed
# error and fall back to adapter$refit for that replicate, without
# aborting the entire call. This file pins both halves of that
# contract.

test_that("update_core raises classed error on an all-duplicate resample", {
  set.seed(1)
  n  <- 24L
  px <- 5L
  py <- 4L
  X <- matrix(rnorm(n * px), n, px)
  Y <- matrix(rnorm(n * py), n, py)

  adapter <- adapter_cross_svd()
  original_fit <- adapter$refit(
    NULL, list(X = X, Y = Y, relation = "correlation")
  )
  core <- adapter$core(original_fit, list(X = X, Y = Y))

  # All-duplicate resample: every bootstrap index points at row 1.
  # Ux[rep(1, n), ] then centers to the zero matrix, so Bx has all
  # zero eigenvalues and the rank-deficiency guard must fire.
  expect_error(
    adapter$update_core(core, indices = rep(1L, n)),
    class = "multifer_core_rank_deficient"
  )
})

test_that("bootstrap_fits catches the classed error and falls back to refit", {
  # Inject a synthetic adapter whose update_core always throws the
  # classed condition. Its other hooks delegate to adapter_cross_svd
  # so bootstrap_fits can still align and assemble reps. Every
  # replicate must take the refit fallback path, the artifact must
  # finish without erroring, and the fallback counter must reach R.
  set.seed(2)
  n  <- 30L
  px <- 5L
  py <- 4L
  X <- matrix(rnorm(n * px), n, px)
  Y <- matrix(rnorm(n * py), n, py)

  base_adapter <- adapter_cross_svd()
  injected <- base_adapter
  injected$update_core <- function(core_obj, indices = NULL, ...) {
    stop(structure(
      class = c("multifer_core_rank_deficient", "error", "condition"),
      list(
        message = "synthetic rank-deficient fault for fallback testing",
        call = NULL
      )
    ))
  }

  recipe <- infer_recipe(
    geometry = "cross", relation = "correlation",
    adapter = injected, strict = TRUE
  )
  original_fit <- base_adapter$refit(
    NULL, list(X = X, Y = Y, relation = "correlation")
  )
  units <- form_units(original_fit$d[1:3]^2)

  artifact <- bootstrap_fits(
    recipe = recipe, adapter = injected,
    data = list(X = X, Y = Y),
    original_fit = original_fit, units = units,
    R = 6L, method_align = "sign", seed = 701L,
    fast_path = "auto", core_rank = NULL
  )

  # The fast-path metadata should still report TRUE because the
  # adapter exposed core()/update_core() hooks, even though every
  # replicate ended up on the refit fallback.
  expect_true(isTRUE(artifact$used_fast_path))
  expect_equal(artifact$cost$core_rank_deficient_fallbacks, 6L)
  expect_equal(artifact$R, 6L)

  # Each replicate's fit must be a valid correlation fit from refit.
  for (b in seq_len(artifact$R)) {
    f <- artifact$reps[[b]]$fit
    expect_equal(f$relation, "correlation")
    expect_true(all(is.finite(f$d)))
    expect_true(all(is.finite(f$Wx)))
    expect_true(all(is.finite(f$Wy)))
  }
})

test_that("bootstrap_fits reports zero fallbacks on a non-degenerate fixture", {
  # Sanity pin: the fallback counter must stay at 0 when nothing
  # actually goes wrong, so a regression that accidentally routes
  # everything through refit would be caught here.
  set.seed(3)
  n  <- 60L
  px <- 8L
  py <- 6L
  X <- matrix(rnorm(n * px), n, px)
  Y <- matrix(rnorm(n * py), n, py)

  adapter <- adapter_cross_svd()
  recipe <- infer_recipe(
    geometry = "cross", relation = "correlation",
    adapter = adapter, strict = TRUE
  )
  original_fit <- adapter$refit(
    NULL, list(X = X, Y = Y, relation = "correlation")
  )
  units <- form_units(original_fit$d[1:3]^2)

  artifact <- bootstrap_fits(
    recipe = recipe, adapter = adapter,
    data = list(X = X, Y = Y),
    original_fit = original_fit, units = units,
    R = 8L, method_align = "sign", seed = 809L,
    fast_path = "auto", core_rank = NULL
  )

  expect_true(isTRUE(artifact$used_fast_path))
  expect_equal(artifact$cost$core_rank_deficient_fallbacks, 0L)
})
