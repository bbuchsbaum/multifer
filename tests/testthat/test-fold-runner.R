test_that("fold ids resolve deterministically without disturbing caller RNG", {
  ids1 <- .fold_resolve_ids(n = 12L, n_folds = 4L, seed = 91L)
  ids2 <- .fold_resolve_ids(n = 12L, n_folds = 4L, seed = 91L)

  expect_identical(ids1, ids2)
  expect_lte(max(tabulate(ids1)) - min(tabulate(ids1)), 1L)

  set.seed(501L)
  .fold_resolve_ids(n = 12L, n_folds = 4L, seed = 91L)
  after <- sample.int(1000L, 5L)

  set.seed(501L)
  expected <- sample.int(1000L, 5L)

  expect_identical(after, expected)
})

test_that("row-aligned subset handles structured payloads and rejects misaligned data", {
  payload <- list(
    X = matrix(seq_len(12L), nrow = 6L, ncol = 2L),
    Y = matrix(seq_len(18L), nrow = 6L, ncol = 3L),
    group = factor(letters[seq_len(6L)])
  )

  sub <- .fold_subset_row_aligned(payload, c(2L, 4L, 6L))

  expect_equal(dim(sub$X), c(3L, 2L))
  expect_equal(dim(sub$Y), c(3L, 3L))
  expect_identical(as.character(sub$group), c("b", "d", "f"))

  bad_payload <- list(
    X = matrix(seq_len(12L), nrow = 6L, ncol = 2L),
    Y = matrix(seq_len(15L), nrow = 5L, ncol = 3L)
  )

  expect_error(
    .fold_subset_row_aligned(bad_payload, 1:3),
    "row-aligned"
  )
})

test_that("fold map threads structured train/test payloads through the worker", {
  payload <- list(
    X = matrix(seq_len(12L), nrow = 6L, ncol = 2L),
    Y = matrix(seq_len(18L), nrow = 6L, ncol = 3L)
  )
  fold_ids <- c(1L, 1L, 2L, 2L, 3L, 3L)

  out <- .fold_map(
    data = payload,
    fold_ids = fold_ids,
    fold_fun = function(train_data, test_data, split) {
      c(
        fold = split$fold,
        train_n = nrow(train_data$X),
        test_n = nrow(test_data$X),
        y_sum = sum(test_data$Y)
      )
    }
  )

  expect_length(out, 3L)
  expect_equal(vapply(out, `[[`, numeric(1L), "train_n"), c(4, 4, 4))
  expect_equal(vapply(out, `[[`, numeric(1L), "test_n"), c(2, 2, 2))
  expect_equal(vapply(out, `[[`, numeric(1L), "fold"), c(1, 2, 3))
})
