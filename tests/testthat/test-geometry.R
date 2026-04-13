test_that("geometry accepts the four valid kinds", {
  for (kind in c("oneblock", "cross", "multiblock", "geneig")) {
    g <- geometry(kind)
    expect_true(is_geometry(g))
    expect_equal(g$kind, kind)
    expect_s3_class(g, paste0("multifer_geometry_", kind))
  }
})

test_that("geometry rejects unknown kinds", {
  expect_error(geometry("pca"), "Unknown geometry")
  expect_error(geometry(""), "Unknown geometry")
})

test_that("geometry rejects non-scalar / non-character input", {
  expect_error(geometry(c("oneblock", "cross")), "single non-NA")
  expect_error(geometry(NA_character_), "single non-NA")
  expect_error(geometry(1), "single non-NA")
  expect_error(geometry(NULL), "single non-NA")
})

test_that("geometry prints", {
  g <- geometry("cross")
  expect_output(print(g), "multifer_geometry: cross")
})
