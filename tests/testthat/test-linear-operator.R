test_that("linear_operator rejects bad inputs", {
  expect_error(linear_operator("not a fn", identity, c(2L, 2L)), "must be a function")
  expect_error(linear_operator(identity, "not a fn", c(2L, 2L)), "must be a function")
  expect_error(linear_operator(identity, identity, c(2L)),         "length 2")
  expect_error(linear_operator(identity, identity, c(2L, -1L)),    "non-negative")
  expect_error(linear_operator(identity, identity, c(2L, 2L), gram = 1), "NULL or a function")
})

test_that("as_linear_operator.matrix wraps a dense matrix", {
  set.seed(1)
  A <- matrix(rnorm(20), nrow = 5, ncol = 4)
  op <- as_linear_operator(A)

  expect_true(is_linear_operator(op))
  expect_equal(dim(op), c(5L, 4L))
  expect_equal(nrow(op), 5L)
  expect_equal(ncol(op), 4L)
})

test_that("lo_apply matches A %*% v", {
  set.seed(2)
  A <- matrix(rnorm(20), nrow = 5, ncol = 4)
  op <- as_linear_operator(A)
  v  <- matrix(rnorm(8), nrow = 4, ncol = 2)
  expect_equal(lo_apply(op, v),  A %*% v)
  expect_equal(lo_apply(op, v[, 1L]), A %*% v[, 1L, drop = FALSE])
})

test_that("lo_apply_t matches t(A) %*% v", {
  set.seed(3)
  A <- matrix(rnorm(20), nrow = 5, ncol = 4)
  op <- as_linear_operator(A)
  v  <- matrix(rnorm(10), nrow = 5, ncol = 2)
  expect_equal(lo_apply_t(op, v), crossprod(A, v))
})

test_that("lo_gram matches crossprod(A) and uses hook when present", {
  set.seed(4)
  A <- matrix(rnorm(20), nrow = 5, ncol = 4)
  op <- as_linear_operator(A)
  expect_equal(lo_gram(op), crossprod(A))
})

test_that("lo_gram falls back to apply_t(apply(I)) when hook absent", {
  set.seed(5)
  A <- matrix(rnorm(20), nrow = 5, ncol = 4)
  op <- linear_operator(
    apply   = function(v) A %*% v,
    apply_t = function(v) crossprod(A, v),
    dim     = c(5L, 4L)
  )
  expect_null(op$gram)
  expect_equal(lo_gram(op), crossprod(A))
})

test_that("as.matrix materializes the operator", {
  set.seed(6)
  A <- matrix(rnorm(12), nrow = 3, ncol = 4)
  op <- as_linear_operator(A)
  expect_equal(as.matrix(op), A)
})

test_that("apply rejects shape mismatches", {
  A <- matrix(1:12, nrow = 3, ncol = 4)
  op <- as_linear_operator(A)
  bad <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(lo_apply(op, bad), "apply: input has 2 rows")
  expect_error(lo_apply_t(op, bad), "apply_t: input has 2 rows")
})

test_that("is_linear_operator returns FALSE for non-operators", {
  expect_false(is_linear_operator(matrix(1)))
  expect_false(is_linear_operator(NULL))
  expect_false(is_linear_operator(list(apply = identity)))
})

test_that("print method does not error", {
  op <- as_linear_operator(matrix(1:4, 2, 2))
  expect_output(print(op), "multifer_linear_operator")
  expect_output(print(op), "2 x 2")
})
