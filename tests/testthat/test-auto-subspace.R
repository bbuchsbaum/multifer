# S4: subspace inference as the default higher-level output
#
# With auto_subspace = TRUE (the default), infer() routes through
# form_units(group_near_ties = TRUE) so near-tied roots bundle into
# subspace units automatically, and print.infer_result leads with
# their principal-angle stability.

test_that("auto_subspace = TRUE groups exactly-tied roots into a subspace unit", {
  ensure_default_adapters()

  # Construct an oneblock matrix with two exactly-tied signal strengths
  # (8, 8) followed by a well-separated smaller root. Even with light
  # noise, the observed top-two singular values should have a relative
  # gap well below the default 1 percent tie_threshold.
  set.seed(101L)
  n <- 200L; p <- 8L
  k <- 3L
  U <- qr.Q(qr(matrix(stats::rnorm(n * k), nrow = n)))
  V <- qr.Q(qr(matrix(stats::rnorm(p * k), nrow = p)))
  sig <- c(8.0, 8.0, 3.0)
  X <- U %*% diag(sig * sqrt(n - 1), nrow = k) %*% t(V) +
    matrix(stats::rnorm(n * p, sd = 0.01), nrow = n, ncol = p)

  res <- infer(
    adapter  = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 99L, R = 6L, alpha = 0.05, seed = 19L,
    tie_threshold = 0.05   # 5 percent band to survive noise robustly
  )

  expect_true(sum(res$units$unit_type == "subspace") >= 1L,
              info = "at least one subspace unit should form from the tied pair")
  subspace_row <- which(res$units$unit_type == "subspace")[1L]
  expect_true(length(res$units$members[[subspace_row]]) >= 2L,
              info = "the subspace unit should bundle >= 2 tied roots")
})

test_that("auto_subspace = FALSE preserves Phase 1 component-only behavior", {
  ensure_default_adapters()
  set.seed(103L)
  n <- 200L; p <- 8L
  k <- 3L
  U <- qr.Q(qr(matrix(stats::rnorm(n * k), nrow = n)))
  V <- qr.Q(qr(matrix(stats::rnorm(p * k), nrow = p)))
  sig <- c(8.0, 8.0, 3.0)
  X <- U %*% diag(sig * sqrt(n - 1), nrow = k) %*% t(V) +
    matrix(stats::rnorm(n * p, sd = 0.01), nrow = n, ncol = p)

  res <- infer(
    adapter  = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 99L, R = 6L, alpha = 0.05, seed = 23L,
    auto_subspace = FALSE
  )

  expect_equal(sum(res$units$unit_type == "subspace"), 0L)
})

test_that("print.infer_result leads with subspace inference when subspaces are present", {
  ensure_default_adapters()
  set.seed(107L)
  n <- 200L; p <- 8L
  k <- 3L
  U <- qr.Q(qr(matrix(stats::rnorm(n * k), nrow = n)))
  V <- qr.Q(qr(matrix(stats::rnorm(p * k), nrow = p)))
  sig <- c(6.0, 6.0, 2.5)
  X <- U %*% diag(sig * sqrt(n - 1), nrow = k) %*% t(V) +
    matrix(stats::rnorm(n * p, sd = 0.01), nrow = n, ncol = p)

  res <- infer(
    adapter  = "svd_oneblock", data = X,
    geometry = "oneblock", relation = "variance",
    B = 99L, R = 6L, alpha = 0.05, seed = 29L,
    tie_threshold = 0.05
  )

  output <- capture.output(print(res))
  joined <- paste(output, collapse = "\n")
  expect_match(joined, "Subspace inference")
  expect_match(joined, "mean angle")
})

test_that("auto_subspace passes through cross-covariance paired-root case", {
  ensure_default_adapters()
  set.seed(109L)
  n <- 80L; p_x <- 6L; p_y <- 5L
  k <- 2L
  U <- qr.Q(qr(matrix(stats::rnorm(n * k), nrow = n)))
  Wx <- qr.Q(qr(matrix(stats::rnorm(p_x * k), nrow = p_x)))
  Wy <- qr.Q(qr(matrix(stats::rnorm(p_y * k), nrow = p_y)))
  sig <- c(5.0, 4.99)
  X <- U %*% diag(sig, nrow = k) %*% t(Wx) +
    matrix(stats::rnorm(n * p_x, sd = 0.1), n, p_x)
  Y <- U %*% diag(sig, nrow = k) %*% t(Wy) +
    matrix(stats::rnorm(n * p_y, sd = 0.1), n, p_y)

  res <- infer(
    adapter  = "cross_svd", data = list(X = X, Y = Y),
    geometry = "cross", relation = "covariance",
    B = 99L, R = 6L, alpha = 0.05, seed = 31L
  )
  # Not asserting subspace presence here because cross generators may not
  # produce perfectly tied roots at this seed, but verify infer() still
  # runs end-to-end with auto_subspace as the default.
  expect_true(is_infer_result(res))
  expect_equal(nrow(res$units), length(unique(c(
    unlist(res$units$members)
  ))))
})
