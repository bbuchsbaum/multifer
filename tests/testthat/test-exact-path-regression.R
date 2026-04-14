# S5: exact-path correctness regression pin for cross covariance
#
# With all Phase 1.5 + solidification commits landed, the cross-
# covariance fast-path infer() must produce the same ladder decisions
# and the same unit selections as the refit-only path across a
# representative grid of problem shapes. This is the exact-only
# correctness pin.
#
# The §40 speed target (≥10× on cross cov large) requires structural
# approximation and is explicitly out of scope for the "solidify"
# track -- the fast path may be as fast or slightly slower than
# refit at small square-ish shapes where k_x * k_y equals p * q and
# the core gate (use_core = k_x*k_y < p*q) keeps us on the direct
# crossprod path. The test therefore does NOT assert a wall-clock
# speedup; it only pins that the two paths agree on decisions so
# regressions in the exact core-space math surface immediately.

test_that("cross covariance fast path agrees with refit path on a grid of shapes", {
  skip_on_cran()
  ensure_default_adapters()

  shapes <- list(
    list(name = "tall_narrow",  n = 200L, p = 30L, q = 25L),
    list(name = "square_mid",   n = 200L, p = 60L, q = 50L),
    list(name = "p_dominant",   n = 200L, p = 100L, q = 40L),
    list(name = "q_dominant",   n = 200L, p = 40L, q = 100L)
  )

  for (sc in shapes) {
    set.seed(2026L + nchar(sc$name))
    X <- matrix(stats::rnorm(sc$n * sc$p), sc$n, sc$p)
    Y <- matrix(stats::rnorm(sc$n * sc$q), sc$n, sc$q)

    common <- list(
      adapter  = "cross_svd", data = list(X = X, Y = Y),
      geometry = "cross", relation = "covariance",
      B = 49L, R = 6L, alpha = 0.05, seed = 41L
    )

    r_off  <- do.call(infer, c(common, list(fast_path = "off")))
    r_auto <- do.call(infer, c(common, list(fast_path = "auto")))

    expect_equal(
      r_off$units$selected,
      r_auto$units$selected,
      info = sprintf("%s selection agreement", sc$name)
    )
    expect_equal(
      r_off$component_tests$p_value <= 0.05,
      r_auto$component_tests$p_value <= 0.05,
      info = sprintf("%s rejection-decision agreement", sc$name)
    )
  }
})

test_that("blocked_rows correlation fast path agrees with refit path on a shape grid (multifer-9u9.1.1)", {
  skip_on_cran()
  ensure_default_adapters()

  shapes <- list(
    list(name = "blocked_6x6",  n = 36L, p = 8L, q = 6L, groups = rep(1:6, each = 6)),
    list(name = "blocked_4x9",  n = 36L, p = 8L, q = 6L, groups = rep(1:4, each = 9)),
    list(name = "blocked_3x12", n = 36L, p = 8L, q = 6L, groups = rep(1:3, each = 12))
  )

  for (sc in shapes) {
    set.seed(2026L + nchar(sc$name))
    X <- matrix(stats::rnorm(sc$n * sc$p), sc$n, sc$p)
    Y <- matrix(stats::rnorm(sc$n * sc$q), sc$n, sc$q)

    common <- list(
      adapter  = "cross_svd",
      data     = list(X = X, Y = Y),
      geometry = "cross",
      relation = "correlation",
      design   = blocked_rows(sc$groups),
      B = 49L, R = 6L, alpha = 0.05, seed = 41L
    )

    r_off  <- do.call(infer, c(common, list(fast_path = "off")))
    r_auto <- do.call(infer, c(common, list(fast_path = "auto")))

    expect_equal(
      r_off$units$selected,
      r_auto$units$selected,
      info = sprintf("%s selection agreement", sc$name)
    )
    expect_equal(
      r_off$component_tests$p_value <= 0.05,
      r_auto$component_tests$p_value <= 0.05,
      info = sprintf("%s rejection-decision agreement", sc$name)
    )
  }
})
