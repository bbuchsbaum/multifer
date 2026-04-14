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

test_that("adapter_lda_refit constructs and exposes the geneig capability quartet", {
  skip_if_not_installed("MASS")
  a <- adapter_lda_refit()

  expect_true(is_infer_adapter(a))
  expect_equal(a$adapter_id, "lda_refit")
  expect_equal(a$shape_kinds, "geneig")
  expect_true(all(
    c("component_significance", "variable_stability",
      "score_stability", "subspace_stability") %in% a$capabilities$target
  ))
  expect_true(isTRUE(attr(a$residualize, "b_metric")))
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
