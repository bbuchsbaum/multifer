# Cross-adapter parity across the CCA support matrix (multifer-aai.2).
#
# cancor_cross and cross_svd (correlation mode) compute the same
# canonical-correlation roots.  This test confirms they also produce
# the same ladder output when given the same data, the same design,
# and the same Monte Carlo seed -- across every supported row of the
# CCA support matrix AND the stress fixtures from multifer-aai.1.
# Failure here means the two shipped CCA adapters have started to
# diverge, and the error message points at notes/cca_support_and_calibration.md
# so a future contributor knows where to start.

fit_cca_backend <- function(adapter, data, design, max_steps = 3L,
                            seed = 31L) {
  infer(
    adapter = adapter,
    data = list(X = data$X, Y = data$Y),
    geometry = "cross",
    relation = "correlation",
    design = design,
    B = 29L,
    R = 0L,
    seed = seed,
    max_steps = max_steps
  )
}

expect_cca_adapter_parity <- function(data, design, label,
                                      adapter = "cross_svd",
                                      max_steps = 3L,
                                      tol = 1e-8) {
  ref <- fit_cca_backend("cancor_cross", data, design, max_steps)
  alt <- fit_cca_backend(adapter, data, design, max_steps)

  expect_true(is_infer_result(ref), info = label)
  expect_true(is_infer_result(alt), info = label)
  expect_equal(
    ref$component_tests$statistic,
    alt$component_tests$statistic,
    tolerance = 1e-10,
    info = paste0(
      label,
      ": canonical-correlation statistics differ for adapter ",
      adapter,
      " -- see notes/cca_multivarious_parity.md"
    )
  )
  expect_equal(
    ref$component_tests$p_value,
    alt$component_tests$p_value,
    tolerance = tol,
    info = paste0(
      label,
      ": p-values differ for adapter ",
      adapter,
      " -- see notes/cca_multivarious_parity.md"
    )
  )
  expect_equal(ref$units, alt$units, info = label)
  expect_equal(
    ref$units$selected,
    alt$units$selected,
    info = label
  )
  expect_equal(
    ref$component_tests$stopped_at,
    alt$component_tests$stopped_at,
    info = label
  )
  expect_equal(
    nrow(ref$component_tests),
    nrow(alt$component_tests),
    info = label
  )
}

expect_cca_backends_parity <- function(data, design, label,
                                       max_steps = 3L,
                                       tol = 1e-8) {
  for (adapter in c("cross_svd", "multivarious_cca")) {
    expect_cca_adapter_parity(
      data = data,
      design = design,
      label = paste0(label, " / ", adapter),
      adapter = adapter,
      max_steps = max_steps,
      tol = tol
    )
  }
}

add_cca_noise_columns <- function(data, extra_x = 0L, extra_y = 0L,
                                  seed = 1L, scale = 0.35) {
  set.seed(seed)
  if (extra_x > 0L) {
    data$X <- cbind(
      data$X,
      matrix(rnorm(nrow(data$X) * extra_x, sd = scale),
             nrow(data$X), extra_x)
    )
  }
  if (extra_y > 0L) {
    data$Y <- cbind(
      data$Y,
      matrix(rnorm(nrow(data$Y) * extra_y, sd = scale),
             nrow(data$Y), extra_y)
    )
  }
  data
}

rank_deficient_cca_fixture <- function() {
  set.seed(530)
  dat <- make_exact_canonical_correlation_support(
    n = 48,
    canonical_corrs = c(0.92, 0.75),
    scale_x = c(4, 2),
    scale_y = c(3, 2)
  )
  dat$X <- cbind(dat$X, dat$X[, 1L] + dat$X[, 2L])
  dat
}

expect_rank_deficient_boundary <- function(adapter) {
  dat <- rank_deficient_cca_fixture()
  expect_error(
    fit_cca_backend(adapter, dat, paired_rows(), max_steps = 2L),
    "full numerical column rank",
    info = adapter
  )
}

stress_fixtures_for_parity <- function() {
  fx <- cca_stress_fixtures()
  fx
}

test_that("CCA parity: paired_rows on all stress fixtures", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()
  fixtures <- stress_fixtures_for_parity()

  for (nm in names(fixtures)) {
    expect_cca_backends_parity(
      data = fixtures[[nm]],
      design = paired_rows(),
      label = paste0("paired_rows / ", nm)
    )
  }
})

test_that("CCA parity: blocked_rows on all stress fixtures", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()
  fixtures <- stress_fixtures_for_parity()

  block_size <- 6L
  for (nm in names(fixtures)) {
    dat <- fixtures[[nm]]
    n_rows <- nrow(dat$X)
    n_groups <- n_rows %/% block_size
    groups <- rep(seq_len(n_groups), each = block_size, length.out = n_rows)
    expect_cca_backends_parity(
      data = dat,
      design = blocked_rows(groups),
      label = paste0("blocked_rows / ", nm)
    )
  }
})

test_that("CCA parity: nuisance_adjusted(Z) on the tight-df fixture", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()
  set.seed(520)
  dat <- make_nuisance_adjusted_correlation_support(
    n = 48,
    canonical_corrs = c(0.95, 0.80, 0.55),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.8)
  )
  expect_cca_backends_parity(
    data = dat,
    design = nuisance_adjusted(dat$Z),
    label = "nuisance_adjusted"
  )
})

test_that("CCA parity: nuisance_adjusted(Z, groups) on a grouped fixture", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()
  set.seed(521)
  groups <- rep(1:8, each = 6)
  dat <- make_nuisance_adjusted_correlation_support(
    n = length(groups),
    canonical_corrs = c(0.95, 0.80, 0.55),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.8),
    groups = groups
  )
  expect_cca_backends_parity(
    data = dat,
    design = nuisance_adjusted(dat$Z, groups = dat$groups),
    label = "nuisance_adjusted(groups)"
  )
})

test_that("CCA parity: exchangeable_rows first-root cap matches across both adapters", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()
  fixtures <- stress_fixtures_for_parity()
  dat <- fixtures$near_tied

  expect_cca_backends_parity(
    data = dat,
    design = exchangeable_rows(),
    label = "exchangeable_rows (first-root cap)",
    max_steps = 1L
  )
})

test_that("CCA parity: aspect-ratio stress fixtures match across backends", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()
  base <- make_exact_canonical_correlation_support(
    n = 80,
    canonical_corrs = c(0.94, 0.81, 0.47),
    scale_x = c(4, 2, 1),
    scale_y = c(3, 2, 1)
  )
  fixtures <- list(
    wide_x_narrow_y = add_cca_noise_columns(base, extra_x = 8L,
                                            extra_y = 1L, seed = 531),
    narrow_x_wide_y = add_cca_noise_columns(base, extra_x = 1L,
                                            extra_y = 8L, seed = 532)
  )

  for (nm in names(fixtures)) {
    expect_cca_backends_parity(
      data = fixtures[[nm]],
      design = paired_rows(),
      label = paste0("aspect ratio / ", nm),
      max_steps = 3L
    )
  }
})

test_that("CCA parity: rank-deficient data fail at the shared validity gate", {
  skip_if_not_installed("multivarious")
  ensure_default_adapters()

  for (adapter in c("cancor_cross", "cross_svd", "multivarious_cca")) {
    expect_rank_deficient_boundary(adapter)
  }
})
