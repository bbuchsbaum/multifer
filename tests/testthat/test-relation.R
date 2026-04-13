test_that("relation accepts the four valid kinds", {
  for (kind in c("variance", "covariance", "correlation", "generalized_eigen")) {
    r <- relation(kind)
    expect_true(is_relation(r))
    expect_equal(r$kind, kind)
    expect_s3_class(r, paste0("multifer_relation_", kind))
  }
})

test_that("relation rejects unknown kinds", {
  expect_error(relation("pls"), "Unknown relation")
  expect_error(relation("cov"), "Unknown relation")
})

test_that("relation rejects non-scalar / non-character input", {
  expect_error(relation(c("variance", "covariance")), "single non-NA")
  expect_error(relation(NA_character_), "single non-NA")
  expect_error(relation(NULL), "single non-NA")
})

test_that("relation prints", {
  r <- relation("correlation")
  expect_output(print(r), "multifer_relation: correlation")
})
