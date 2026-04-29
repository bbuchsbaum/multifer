make_lda_signal_fixture <- function(n_per_class = 24L, seed = 1L) {
  set.seed(seed)
  y <- factor(rep(c("a", "b", "c"), each = n_per_class))
  means <- rbind(
    c(-4, 0, 0, 0),
    c( 0, 0, 0, 0),
    c( 4, 0, 0, 0)
  )
  X <- do.call(rbind, lapply(seq_len(nrow(means)), function(i) {
    matrix(stats::rnorm(n_per_class * ncol(means), sd = 0.4), ncol = ncol(means)) +
      matrix(means[i, ], nrow = n_per_class, ncol = ncol(means), byrow = TRUE)
  }))
  list(X = X, y = y)
}

make_lda_null_fixture <- function(n_per_class = 24L, seed = 1L) {
  set.seed(seed)
  y <- factor(rep(c("a", "b", "c"), each = n_per_class))
  X <- matrix(stats::rnorm(length(y) * 4L, sd = 1), ncol = 4L)
  list(X = X, y = y)
}

test_that("adapter_lda_refit constructs and exposes only shipped geneig significance", {
  skip_if_not_installed("MASS")
  a <- adapter_lda_refit()

  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "lda_refit")
  expect_equal(a$shape_kinds, "geneig")
  expect_equal(unique(a$capabilities$target), "component_significance")
  expect_equal(a$geneig_deflation, "b_metric")
})

test_that("lda_refit null_action permutes labels rather than rows", {
  skip_if_not_installed("MASS")
  dat <- make_lda_signal_fixture(seed = 2L)
  a <- adapter_lda_refit()
  op <- .lda_geneig_operator(dat$X, dat$y)

  permuted <- a$null_action(NULL, list(
    A = op$A,
    B = op$B,
    metric = op$metric,
    state = dat
  ))

  expect_equal(permuted$state$X, dat$X)
  expect_false(identical(as.character(permuted$state$y), as.character(dat$y)))
})

test_that("lda_refit component_stat returns the generalized-root tail ratio", {
  skip_if_not_installed("MASS")
  dat <- make_lda_signal_fixture(seed = 3L)
  fit <- .fit_lda_refit(dat$X, dat$y)
  a <- adapter_lda_refit()

  expected_1 <- fit$roots[1L] / sum(fit$roots)
  expected_2 <- fit$roots[2L] / sum(fit$roots[2:length(fit$roots)])

  expect_equal(a$component_stat(fit, dat, 1L), expected_1, tolerance = 1e-12)
  expect_equal(a$component_stat(fit, dat, 2L), expected_2, tolerance = 1e-12)
})

test_that("lda_refit fails early for unsupported within-class rank regimes", {
  skip_if_not_installed("MASS")
  a <- adapter_lda_refit()
  y <- factor(rep(c("a", "b", "c"), each = 4L))

  set.seed(31L)
  high_dimensional <- list(
    X = matrix(stats::rnorm(length(y) * 20L), nrow = length(y)),
    y = y
  )
  high_dim_checks <- run_adapter_checks(a, high_dimensional, strict = FALSE)
  expect_false(high_dim_checks$lda_within_class_full_rank$passed)
  expect_match(
    high_dim_checks$lda_within_class_full_rank$detail,
    "within-class centered rank"
  )
  expect_error(
    infer_lda(high_dimensional$X, y, adapter = a, B = 9L, seed = 31L),
    "lda_within_class_full_rank"
  )

  set.seed(32L)
  collinear_x <- matrix(stats::rnorm(length(y) * 4L), nrow = length(y))
  collinear_x[, 4L] <- collinear_x[, 1L]
  collinear_checks <- run_adapter_checks(
    a,
    list(X = collinear_x, y = y),
    strict = FALSE
  )
  expect_false(collinear_checks$lda_within_class_full_rank$passed)
})

test_that("lda_refit distinguishes class imbalance from too-small classes", {
  skip_if_not_installed("MASS")
  a <- adapter_lda_refit()

  set.seed(33L)
  imbalanced_y <- factor(rep(c("a", "b", "c"), times = c(2L, 5L, 9L)))
  imbalanced_x <- matrix(stats::rnorm(length(imbalanced_y) * 4L),
                         nrow = length(imbalanced_y))
  ok <- run_adapter_checks(
    a,
    list(X = imbalanced_x, y = imbalanced_y),
    strict = FALSE
  )
  expect_true(ok$lda_within_class_sample_size$passed)
  expect_true(ok$lda_within_class_full_rank$passed)

  res <- infer_lda(imbalanced_x, imbalanced_y, adapter = a, B = 9L, seed = 33L)
  expect_true(is_infer_result(res))
  expect_true(all(is.finite(res$component_tests$p_value)))

  too_small_y <- factor(rep(c("a", "b", "c"), times = c(1L, 5L, 5L)))
  too_small_x <- matrix(stats::rnorm(length(too_small_y) * 3L),
                        nrow = length(too_small_y))
  fail <- run_adapter_checks(
    a,
    list(X = too_small_x, y = too_small_y),
    strict = FALSE
  )
  expect_false(fail$lda_within_class_sample_size$passed)
  expect_error(
    infer_lda(too_small_x, too_small_y, adapter = a, B = 9L, seed = 34L),
    "lda_within_class_sample_size"
  )
})

test_that("lda_refit + geneig engine recovers rank bounded by K - 1", {
  skip_if_not_installed("MASS")
  dat <- make_lda_signal_fixture(seed = 4L)
  ensure_default_adapters()
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = "lda_refit"
  )

  res <- run_geneig_ladder(
    rec,
    .lda_geneig_operator(dat$X, dat$y),
    state = dat,
    B = 39L,
    alpha = 0.05,
    seed = 7L
  )

  expect_true(res$ladder_result$rejected_through >= 1L)
  expect_true(res$ladder_result$rejected_through <= (nlevels(dat$y) - 1L))
  expect_true(isTRUE(res$component_tests$selected[1L]))
  expect_false(tail(res$component_tests$selected, 1L))
})

test_that("lda_refit label-permutation null is approximately calibrated", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = "lda_refit"
  )

  pvals <- vapply(seq_len(20L), function(i) {
    dat <- make_lda_null_fixture(seed = 100L + i)
    res <- run_geneig_ladder(
      rec,
      .lda_geneig_operator(dat$X, dat$y),
      state = dat,
      B = 39L,
      alpha = 0.05,
      seed = 200L + i
    )
    res$component_tests$p_value[1L]
  }, double(1L))

  expect_true(all(pvals >= 0 & pvals <= 1))
  expect_true(mean(pvals) > 0.3 && mean(pvals) < 0.7)
})
