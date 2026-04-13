test_that("exchangeable_rows and paired_rows construct", {
  d1 <- exchangeable_rows()
  d2 <- paired_rows()
  expect_true(is_design(d1))
  expect_true(is_design(d2))
  expect_equal(d1$kind, "exchangeable_rows")
  expect_equal(d2$kind, "paired_rows")
})

test_that("blocked_rows accepts valid group vectors", {
  d <- blocked_rows(c(1, 1, 2, 2, 3, 3))
  expect_true(is_design(d))
  expect_s3_class(d$groups, "factor")
  expect_equal(length(levels(d$groups)), 3L)

  d2 <- blocked_rows(factor(c("a", "b", "a", "b")))
  expect_equal(length(levels(d2$groups)), 2L)
})

test_that("blocked_rows rejects invalid groups", {
  expect_error(blocked_rows(NULL), "must not be NULL")
  expect_error(blocked_rows(integer(0)), "length >= 1")
  expect_error(blocked_rows(c(1L, NA_integer_)), "must not contain NA")
  expect_error(blocked_rows(list("a", "b")), "must be a factor")
})

test_that("nuisance_adjusted accepts a numeric matrix", {
  Z <- matrix(rnorm(12), nrow = 4, ncol = 3)
  d <- nuisance_adjusted(Z)
  expect_true(is_design(d))
  expect_equal(dim(d$Z), c(4L, 3L))
  expect_null(d$groups)

  d2 <- nuisance_adjusted(Z, groups = c("a", "a", "b", "b"))
  expect_s3_class(d2$groups, "factor")
  expect_equal(length(levels(d2$groups)), 2L)
})

test_that("nuisance_adjusted rejects bad Z", {
  expect_error(nuisance_adjusted(NULL), "must not be NULL")
  expect_error(nuisance_adjusted(1:4), "numeric matrix")
  Z <- matrix(c(1, NA, 3, 4), nrow = 2)
  expect_error(nuisance_adjusted(Z), "must not contain NA")
  expect_error(nuisance_adjusted(matrix(1:9, 3, 3), groups = 1:2),
               "length nrow\\(Z\\)")
  expect_error(nuisance_adjusted(matrix(1:9, 3, 3), groups = c(1, NA, 2)),
               "must not contain NA")
  expect_error(nuisance_adjusted(matrix(1:9, 3, 3), groups = list("a", "b", "c")),
               "must be a factor")
})

test_that("design print methods run", {
  expect_output(print(exchangeable_rows()), "exchangeable_rows")
  expect_output(print(paired_rows()), "paired_rows")
  expect_output(print(blocked_rows(c(1, 1, 2, 2))), "blocks")
  expect_output(print(nuisance_adjusted(matrix(1:4, 2, 2))), "Z:")
  expect_output(print(nuisance_adjusted(matrix(1:8, 4, 2), groups = c(1, 1, 2, 2))),
                "blocks")
})
