make_lda_wrapper_fixture <- function(n_per_class = 24L, seed = 1L) {
  set.seed(seed)
  labels <- factor(rep(c("a", "b", "c"), each = n_per_class))
  means <- rbind(
    c(-4, 0, 0, 0),
    c( 0, 0, 0, 0),
    c( 4, 0, 0, 0)
  )
  X <- do.call(rbind, lapply(seq_len(nrow(means)), function(i) {
    matrix(stats::rnorm(n_per_class * ncol(means), sd = 0.4), ncol = ncol(means)) +
      matrix(means[i, ], nrow = n_per_class, ncol = ncol(means), byrow = TRUE)
  }))
  list(X = X, labels = labels)
}

make_plsr_wrapper_fixture <- function(n = 54L, p = 6L, q = 2L, seed = 1L) {
  set.seed(seed)
  t1 <- scale(rnorm(n), center = TRUE, scale = FALSE)[, 1L]
  X <- matrix(rnorm(n * p, sd = 0.14), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * q, sd = 0.10), nrow = n, ncol = q)
  X <- X + t1 %o% c(1.4, -0.9, 0.7, 0, 0, 0)
  Y <- Y + t1 %o% c(1.2, -1.0)
  list(X = X, Y = Y)
}

test_that("infer_pca matches direct infer() call", {
  ensure_default_adapters()
  set.seed(901)
  X <- matrix(rnorm(120), 30, 4)

  wrapped <- infer_pca(X, B = 29L, R = 4L, seed = 17L)
  direct <- infer(
    adapter = "prcomp_oneblock",
    data = X,
    geometry = "oneblock",
    relation = "variance",
    B = 29L,
    R = 4L,
    seed = 17L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
})

test_that("infer_plsc matches direct infer() call", {
  ensure_default_adapters()
  set.seed(902)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)

  wrapped <- infer_plsc(X, Y, B = 29L, R = 4L, seed = 19L)
  direct <- infer(
    adapter = "cross_svd",
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "covariance",
    B = 29L,
    R = 4L,
    seed = 19L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
})

test_that("infer_cca matches direct infer() call", {
  ensure_default_adapters()
  set.seed(903)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)

  wrapped <- infer_cca(X, Y, B = 29L, R = 4L, seed = 23L)
  direct <- infer(
    adapter = "cancor_cross",
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "correlation",
    B = 29L,
    R = 4L,
    seed = 23L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
  expect_equal(wrapped$provenance$adapter_id, "cancor_cross")
})

test_that("infer_cca forwards supported nuisance-adjusted designs through the default adapter", {
  ensure_default_adapters()
  set.seed(904)
  n <- 36
  Z <- cbind(1, scale(seq_len(n)), sin(seq_len(n) / 5))
  X <- matrix(rnorm(n * 5L), n, 5L) + 0.2 * Z %*% matrix(rnorm(ncol(Z) * 5L), ncol(Z), 5L)
  Y <- matrix(rnorm(n * 4L), n, 4L) + 0.2 * Z %*% matrix(rnorm(ncol(Z) * 4L), ncol(Z), 4L)

  wrapped <- infer_cca(
    X, Y,
    design = nuisance_adjusted(Z),
    B = 29L,
    R = 4L,
    seed = 31L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$provenance$adapter_id, "cancor_cross")
  expect_true(all(wrapped$component_tests$p_value >= 0 &
                  wrapped$component_tests$p_value <= 1))
})

test_that("infer_plsr matches direct infer() call on a synthetic fixture", {
  skip_if_not_installed("pls")
  ensure_default_adapters()
  dat <- make_plsr_wrapper_fixture(seed = 904)

  wrapped <- infer_plsr(
    dat$X,
    dat$Y,
    targets = "component_significance",
    B = 29L,
    seed = 37L
  )
  direct <- infer(
    adapter = "plsr_refit",
    data = list(X = dat$X, Y = dat$Y),
    geometry = "cross",
    relation = "predictive",
    targets = "component_significance",
    B = 29L,
    seed = 37L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
  expect_equal(wrapped$provenance$adapter_id, "plsr_refit")
})

test_that("infer_plsr refuses covariance requests through the predictive wrapper", {
  skip_if_not_installed("pls")
  dat <- make_plsr_wrapper_fixture(seed = 905)

  expect_error(
    infer_plsr(dat$X, dat$Y, relation = "covariance"),
    "relation = \"predictive\""
  )
})

test_that("infer_lda recovers planted discriminant rank through the geneig wrapper", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()
  dat <- make_lda_wrapper_fixture(seed = 905)

  wrapped <- infer_lda(
    dat$X,
    dat$labels,
    B = 39L,
    seed = 41L
  )
  direct <- infer(
    adapter = "lda_refit",
    data = list(X = dat$X, y = dat$labels),
    geometry = "geneig",
    relation = "generalized_eigen",
    targets = "component_significance",
    B = 39L,
    seed = 41L
  )

  expect_true(is_infer_result(wrapped))
  expect_equal(wrapped$component_tests, direct$component_tests)
  expect_equal(wrapped$units, direct$units)
  expect_true(wrapped$units$selected[1L])
  expect_false(tail(wrapped$units$selected, 1L))
  expect_equal(wrapped$provenance$adapter_id, "lda_refit")
})

test_that("infer_lda refuses invalid labels", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()
  X <- matrix(rnorm(60), 15, 4)

  expect_error(
    infer_lda(X, labels = rep("a", nrow(X))),
    "`labels` must be a factor"
  )
  expect_error(
    infer_lda(X, labels = factor(rep("a", nrow(X) - 1L))),
    "one value per row"
  )
})
