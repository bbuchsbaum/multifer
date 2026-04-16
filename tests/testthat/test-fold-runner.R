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

test_that("fold ids reject malformed specifications", {
  expect_error(
    .fold_resolve_ids(n = 6L, folds = c(1, 1, 2, 2, NA, 3)),
    "must not contain missing values"
  )

  expect_error(
    .fold_resolve_ids(n = 6L, folds = rep(1L, 6L)),
    "at least 2 non-empty folds"
  )

  expect_error(
    .fold_make_splits(c(1L, 1L, 2L, 2L, NA, 3L)),
    "must not contain missing values"
  )
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

  expect_error(
    .fold_subset_row_aligned(list(), 1:3),
    "must not be an empty list"
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

test_that("fold map matches sequential and mirai backends under fixed seeds", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  on.exit(try(multifer_parallel_shutdown(), silent = TRUE), add = TRUE)

  payload <- list(
    X = matrix(seq_len(12L), nrow = 6L, ncol = 2L),
    Y = matrix(seq_len(18L), nrow = 6L, ncol = 3L)
  )
  fold_ids <- c(1L, 1L, 2L, 2L, 3L, 3L)
  seeds <- c(11L, 12L, 13L)

  fold_worker <- function(train_data, test_data, split) {
    c(
      fold = split$fold,
      draw = stats::runif(1L),
      train_sum = sum(train_data$X),
      test_sum = sum(test_data$Y)
    )
  }

  seq_out <- .fold_map(
    data = payload,
    fold_ids = fold_ids,
    fold_fun = fold_worker,
    backend = "sequential",
    seeds = seeds
  )

  mir_out <- .fold_map(
    data = payload,
    fold_ids = fold_ids,
    fold_fun = fold_worker,
    backend = "mirai",
    seeds = seeds
  )

  expect_equal(seq_out, mir_out)
})
