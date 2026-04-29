test_that("infer dispatches multiblock data through adapter hooks", {
  set.seed(2001)
  adapter <- make_stub_multiblock_adapter()
  dat <- make_stub_multiblock_data()

  res <- infer(
    adapter = adapter,
    data = dat,
    geometry = "multiblock",
    relation = "variance",
    B = 9L,
    R = 3L,
    seed = 19L
  )

  expect_s3_class(res, "infer_result")
  expect_equal(res$provenance$adapter_id, "stub_multiblock")
  expect_true(all(names(dat) %in% unique(res$variable_stability$domain)))
  expect_true(all(names(dat) %in% unique(res$score_stability$domain)))
  expect_gt(nrow(res$component_tests), 0L)
})

test_that("bootstrap_fits uses arbitrary multiblock domains", {
  set.seed(2002)
  adapter <- make_stub_multiblock_adapter()
  dat <- make_stub_multiblock_data()
  rec <- infer_recipe(
    geometry = "multiblock",
    relation = "variance",
    adapter = adapter,
    targets = c("variable_stability", "score_stability", "subspace_stability")
  )
  fit <- adapter$refit(NULL, dat)
  units <- form_units(adapter$roots(fit))

  art <- bootstrap_fits(
    recipe = rec,
    adapter = adapter,
    data = dat,
    original_fit = fit,
    units = units,
    R = 2L,
    seed = 33L
  )

  expect_s3_class(art, "multifer_bootstrap_artifact")
  expect_equal(art$domains, names(dat))
  expect_equal(names(art$reps[[1L]]$aligned_loadings), names(dat))
})

test_that("infer passes supplied multiblock model into the ladder", {
  set.seed(2003)
  adapter <- make_stub_multiblock_adapter()
  dat <- make_stub_multiblock_data()
  fit <- adapter$refit(NULL, dat)
  base_refit <- adapter$refit
  adapter$refit <- function(x, new_data, ...) {
    if (is.null(x)) {
      stop("unexpected original refit", call. = FALSE)
    }
    base_refit(x, new_data, ...)
  }

  expect_no_error(res <- infer(
    adapter = adapter,
    data = dat,
    model = fit,
    geometry = "multiblock",
    relation = "variance",
    targets = "component_significance",
    B = 5L,
    R = 0L,
    seed = 41L
  ))

  expect_s3_class(res, "infer_result")
})

test_that("multiblock geometry requires aligned matrix blocks", {
  adapter <- make_stub_multiblock_adapter()
  bad <- list(alpha = matrix(1, 3, 2), beta = matrix(1, 4, 2))

  expect_error(
    infer(
      adapter = adapter,
      data = bad,
      geometry = "multiblock",
      relation = "variance",
      targets = "component_significance",
      B = 3L,
      R = 0L
    ),
    "same row count"
  )
})
