test_that("infer_bundle returns result plus reusable bootstrap artifact", {
  ensure_default_adapters()
  set.seed(1901)
  X <- matrix(rnorm(120), 30, 4)

  bundle <- infer_bundle(
    adapter = "prcomp_oneblock",
    data = X,
    geometry = "oneblock",
    relation = "variance",
    B = 9L,
    R = 3L,
    seed = 19L
  )

  expect_s3_class(bundle, "multifer_infer_bundle")
  expect_true(is_infer_result(bundle$result))
  expect_s3_class(bundle$artifacts$bootstrap, "multifer_bootstrap_artifact")
  expect_false(is.null(bundle$artifacts$original_fit))
  expect_s3_class(bundle$artifacts$recipe, "multifer_infer_recipe")
  expect_s3_class(bundle$artifacts$adapter, "multifer_adapter")

  ev <- feature_evidence_from_bootstrap(
    bundle$artifacts$bootstrap,
    bundle$result$units,
    original_fit = bundle$artifacts$original_fit,
    adapter = bundle$artifacts$adapter,
    statistic = "squared_loading",
    scope = "unit",
    k = "all",
    normalize = "none"
  )
  expect_s3_class(ev, "multifer_feature_evidence")
  expect_true(nrow(ev) > 0L)
})

test_that("infer default result shape remains unchanged", {
  ensure_default_adapters()
  set.seed(1902)
  X <- matrix(rnorm(100), 25, 4)

  result <- infer(
    adapter = "prcomp_oneblock",
    data = X,
    geometry = "oneblock",
    relation = "variance",
    B = 9L,
    R = 2L,
    seed = 23L
  )

  expect_true(is_infer_result(result))
  expect_false(inherits(result, "multifer_infer_bundle"))
  expect_false("artifacts" %in% names(result))
})
