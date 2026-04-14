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

expect_cca_adapter_parity <- function(data, design, label,
                                      max_steps = 3L,
                                      tol = 1e-8) {
  cc <- infer(
    adapter = "cancor_cross",
    data = list(X = data$X, Y = data$Y),
    geometry = "cross",
    relation = "correlation",
    design = design,
    B = 29L,
    R = 0L,
    seed = 31L,
    max_steps = max_steps
  )
  cs <- infer(
    adapter = "cross_svd",
    data = list(X = data$X, Y = data$Y),
    geometry = "cross",
    relation = "correlation",
    design = design,
    B = 29L,
    R = 0L,
    seed = 31L,
    max_steps = max_steps
  )

  expect_true(is_infer_result(cc), info = label)
  expect_true(is_infer_result(cs), info = label)
  expect_equal(
    cc$component_tests$statistic,
    cs$component_tests$statistic,
    tolerance = 1e-10,
    info = paste0(
      label,
      ": canonical-correlation statistics differ -- see ",
      "notes/cca_support_and_calibration.md"
    )
  )
  expect_equal(
    cc$component_tests$p_value,
    cs$component_tests$p_value,
    tolerance = tol,
    info = paste0(
      label,
      ": p-values differ -- see notes/cca_support_and_calibration.md"
    )
  )
  expect_equal(cc$units, cs$units, info = label)
}

stress_fixtures_for_parity <- function() {
  fx <- cca_stress_fixtures()
  fx
}

test_that("CCA parity: paired_rows on all stress fixtures", {
  ensure_default_adapters()
  fixtures <- stress_fixtures_for_parity()

  for (nm in names(fixtures)) {
    expect_cca_adapter_parity(
      data = fixtures[[nm]],
      design = paired_rows(),
      label = paste0("paired_rows / ", nm)
    )
  }
})

test_that("CCA parity: blocked_rows on all stress fixtures", {
  ensure_default_adapters()
  fixtures <- stress_fixtures_for_parity()

  block_size <- 6L
  for (nm in names(fixtures)) {
    dat <- fixtures[[nm]]
    n_rows <- nrow(dat$X)
    n_groups <- n_rows %/% block_size
    groups <- rep(seq_len(n_groups), each = block_size, length.out = n_rows)
    expect_cca_adapter_parity(
      data = dat,
      design = blocked_rows(groups),
      label = paste0("blocked_rows / ", nm)
    )
  }
})

test_that("CCA parity: nuisance_adjusted(Z) on the tight-df fixture", {
  ensure_default_adapters()
  set.seed(520)
  dat <- make_nuisance_adjusted_correlation_support(
    n = 48,
    canonical_corrs = c(0.95, 0.80, 0.55),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.8)
  )
  expect_cca_adapter_parity(
    data = dat,
    design = nuisance_adjusted(dat$Z),
    label = "nuisance_adjusted"
  )
})

test_that("CCA parity: nuisance_adjusted(Z, groups) on a grouped fixture", {
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
  expect_cca_adapter_parity(
    data = dat,
    design = nuisance_adjusted(dat$Z, groups = dat$groups),
    label = "nuisance_adjusted(groups)"
  )
})

test_that("CCA parity: exchangeable_rows first-root cap matches across both adapters", {
  ensure_default_adapters()
  fixtures <- stress_fixtures_for_parity()
  dat <- fixtures$near_tied

  expect_cca_adapter_parity(
    data = dat,
    design = exchangeable_rows(),
    label = "exchangeable_rows (first-root cap)",
    max_steps = 1L
  )
})
