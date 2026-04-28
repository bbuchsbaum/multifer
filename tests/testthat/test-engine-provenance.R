test_that("infer_mc carries engine-provenance labels for oneblock results", {
  ensure_default_adapters()

  set.seed(7001)
  X <- matrix(rnorm(120), 30, 4)

  res <- infer_pca(X, B = 29L, R = 4L, seed = 11L)

  expect_true(is_infer_result(res))
  expect_equal(res$mc$estimand_label, "latent variance roots")
  expect_match(res$mc$statistic_label, "Vitale P3 tail-ratio")
  expect_match(res$mc$null_label, "column permutation")
})

test_that("infer_mc carries engine-provenance labels for PLSC covariance results", {
  ensure_default_adapters()

  set.seed(7002)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)

  res <- infer_plsc(X, Y, B = 29L, R = 4L, seed = 12L)

  expect_equal(res$mc$estimand_label, "cross-covariance singular roots")
  expect_match(res$mc$statistic_label, "cross-covariance roots")
  expect_match(res$mc$null_label, "row permutation of Y")
})

test_that("infer_mc carries engine-provenance labels for CCA correlation results", {
  ensure_default_adapters()

  set.seed(7003)
  X <- matrix(rnorm(200), 40, 5)
  Y <- matrix(rnorm(160), 40, 4)

  res <- infer_cca(X, Y, B = 29L, R = 4L, seed = 13L)

  expect_equal(res$mc$estimand_label, "canonical correlations")
  expect_match(res$mc$statistic_label, "canonical correlations")
  expect_match(res$mc$null_label, "row permutation of Y")
})

test_that("infer_mc carries engine-provenance labels for LDA geneig results", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()

  set.seed(7004)
  labels <- factor(rep(c("a", "b", "c"), each = 15L))
  X <- matrix(rnorm(45 * 4), 45, 4)
  X[labels == "b", 1L] <- X[labels == "b", 1L] + 1.5
  X[labels == "c", 2L] <- X[labels == "c", 2L] - 1.5

  res <- infer_lda(X, labels, B = 29L)

  expect_equal(res$mc$estimand_label, "generalized eigenvalues (B-metric)")
  expect_match(res$mc$statistic_label, "generalized eigenvalues")
  expect_equal(res$mc$null_label, "label permutation")
})

test_that("infer_mc carries engine-provenance labels for PLSR predictive results", {
  skip_if_not_installed("pls")
  ensure_default_adapters()

  set.seed(7005)
  n <- 60L
  X <- matrix(rnorm(n * 5), n, 5)
  beta <- matrix(c(1, 0, -1, 0, 0,
                   0, 1, 0, -1, 0), ncol = 2L)
  Y <- X %*% beta + matrix(rnorm(n * 2, sd = 0.5), n, 2)

  res <- infer_plsr(X, Y, B = 29L, R = 0L, seed = 15L)

  expect_equal(res$mc$estimand_label, "held-out predictive gain")
  expect_match(res$mc$statistic_label, "incremental held-out predictive gain")
  expect_match(res$mc$null_label, "row permutation of Y")
})

test_that("component_tests$null_label is supplied by the engine, not dispatch", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()

  set.seed(7100)
  labels <- factor(rep(c("a", "b", "c"), each = 15L))
  X <- matrix(rnorm(45 * 4), 45, 4)
  X[labels == "b", 1L] <- X[labels == "b", 1L] + 1.5

  res <- infer_lda(X, labels, B = 29L)

  # The per-rung null label in component_tests is now the engine's own
  # string, not the legacy "permute_labels" dispatch constant.
  expect_true(all(res$component_tests$null_label == "label permutation"))
  expect_false(any(res$component_tests$null_label == "permute_labels"))
})

test_that("R/infer.R dispatch contains no geometry-specific null-label branches", {
  infer_path <- test_path("..", "..", "R", "infer.R")
  skip_if_not(file.exists(infer_path), "R/infer.R is not available in installed-package checks")

  infer_src <- readLines(infer_path)
  expect_false(any(grepl("geom_kind == \"geneig\"", infer_src, fixed = TRUE) &
                   grepl("permute_labels", infer_src, fixed = FALSE)))
  expect_false(any(grepl("\"permute_labels\"", infer_src, fixed = TRUE)))
  expect_false(any(grepl("\"row_permute_y_resid_basis\"", infer_src, fixed = TRUE)))
  expect_false(any(grepl("\"column_permute\"", infer_src, fixed = TRUE)))
})

test_that("print.infer_result header includes the engine estimand label", {
  ensure_default_adapters()

  set.seed(7006)
  X <- matrix(rnorm(120), 30, 4)
  res <- infer_pca(X, B = 29L, R = 4L, seed = 17L)

  output <- capture.output(print(res))
  expect_true(any(grepl("estimand:", output, fixed = TRUE)))
  expect_true(any(grepl("latent variance roots", output, fixed = TRUE)))
})
