test_that("multifer_maturity_levels() returns the canonical three-label vocabulary", {
  expect_equal(multifer_maturity_levels(),
               c("mature", "narrow", "planned"))
})

test_that("multifer_adapter_maturity() returns the known label for each v1 adapter", {
  expect_equal(multifer_adapter_maturity("svd_oneblock"), "mature")
  expect_equal(multifer_adapter_maturity("prcomp_oneblock"), "mature")
  expect_equal(multifer_adapter_maturity("cross_svd"), "mature")
  expect_equal(multifer_adapter_maturity("cancor_cross"), "mature")
  expect_equal(multifer_adapter_maturity("multivarious_pca"), "mature")
  expect_equal(multifer_adapter_maturity("multivarious_plsc"), "mature")
  expect_equal(multifer_adapter_maturity("multivarious_cca"), "narrow")
  expect_equal(multifer_adapter_maturity("lda_refit"), "narrow")
  expect_equal(multifer_adapter_maturity("plsr_refit"), "narrow")
})

test_that("multifer_adapter_maturity() returns 'unknown' for unrecognized ids", {
  expect_equal(multifer_adapter_maturity("imaginary_adapter"), "unknown")
})

test_that("multifer_adapter_maturity() without args returns the full v1 map", {
  tbl <- multifer_adapter_maturity()
  expect_type(tbl, "character")
  expect_true("cross_svd" %in% names(tbl))
  expect_true(all(tbl %in% multifer_maturity_levels()))
})

test_that("list_infer_adapters() defaults to a character vector (backcompat)", {
  ensure_default_adapters()
  v <- list_infer_adapters()
  expect_type(v, "character")
  expect_true(length(v) > 0L)
})

test_that("list_infer_adapters(details = TRUE) returns adapter_id/maturity/geometry/relation", {
  ensure_default_adapters()
  df <- list_infer_adapters(details = TRUE)
  expect_s3_class(df, "data.frame")
  expect_setequal(
    colnames(df),
    c("adapter_id", "maturity", "geometry", "relation")
  )
  expect_true(all(df$maturity %in% c(multifer_maturity_levels(), "unknown")))

  # Every registered v1 adapter should have a canonical maturity
  # (no 'unknown' rows) -- that is what 384.2 pins.
  expect_false(any(df$maturity == "unknown"))

  # Mature tier must include at least the oneblock and cross
  # latent-root backbone.
  mature_ids <- df$adapter_id[df$maturity == "mature"]
  expect_true("cross_svd" %in% mature_ids)
  expect_true("prcomp_oneblock" %in% mature_ids ||
                "svd_oneblock" %in% mature_ids)

  # Narrow tier must include the geneig and predictive wrappers.
  narrow_ids <- df$adapter_id[df$maturity == "narrow"]
  expect_true("lda_refit" %in% narrow_ids)
  expect_true("plsr_refit" %in% narrow_ids)
})
