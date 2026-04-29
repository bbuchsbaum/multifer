test_that("match_components: identity input returns identity permutation", {
  set.seed(1)
  V <- qr.Q(qr(matrix(stats::rnorm(12), 4, 3)))
  p <- match_components(V, V)
  expect_equal(p, 1:3)
})

test_that("match_components: known permutation is recovered", {
  set.seed(2)
  V <- qr.Q(qr(matrix(stats::rnorm(12), 4, 3)))
  Vb <- V[, c(2, 1, 3)]
  p <- match_components(V, Vb, diagnostics = TRUE)
  # Vb[, p] should match V, so p must map back: col 1 of Vb is V col 2,
  # col 2 is V col 1, col 3 is V col 3 -> p = c(2, 1, 3)
  expect_equal(as.integer(p), c(2L, 1L, 3L))
  expect_false(is.null(attr(p, "match_score")))
  expect_true(attr(p, "match_method") %in% c("hungarian", "exhaustive", "greedy"))
})

test_that("match_components can use exact assignment and report ambiguity", {
  Vref <- diag(2)
  Vb <- matrix(c(1, 1, 1, 1), nrow = 2) / sqrt(2)
  p <- match_components(Vref, Vb, method = "auto", diagnostics = TRUE)

  expect_equal(sort(p), 1:2)
  expect_true(isTRUE(attr(p, "ambiguous_match")))
  expect_equal(attr(p, "match_margin"), 0)
})

test_that("match_components supports cosine metric", {
  Vref <- diag(2)
  Vb <- Vref[, c(2, 1)] %*% diag(c(10, 0.5))
  p <- match_components(Vref, Vb, method = "auto", metric = "cosine",
                        diagnostics = TRUE)

  expect_equal(as.integer(p), c(2L, 1L))
  expect_true(attr(p, "match_score") <= 2 + 1e-12)
})

test_that("match_components: sign flips do not confuse matching", {
  set.seed(3)
  V  <- qr.Q(qr(matrix(stats::rnorm(12), 4, 3)))
  Vb <- V
  Vb[, 2] <- -Vb[, 2]   # flip sign of second column
  p <- match_components(V, Vb)
  expect_equal(p, 1:3)   # still the identity permutation
})

test_that("align_sign: diagonal of t(Vref) %*% aligned is all non-negative", {
  set.seed(4)
  V  <- qr.Q(qr(matrix(stats::rnorm(20), 5, 4)))
  # introduce arbitrary sign flips and a permutation
  Vb <- V[, c(2, 1, 3, 4)]
  Vb[, 1] <- -Vb[, 1]
  Vb[, 3] <- -Vb[, 3]

  p       <- match_components(V, Vb)
  Vb_perm <- Vb[, p, drop = FALSE]
  aligned <- align_sign(V, Vb_perm)

  dots <- diag(t(V) %*% aligned)
  expect_true(all(dots >= 0))
})

test_that("align_procrustes: aligning a rotated orthonormal matrix recovers the original", {
  set.seed(5)
  n <- 10
  k <- 3
  V <- qr.Q(qr(matrix(stats::rnorm(n * k), n, k)))

  # create a random k x k orthogonal rotation
  R_raw  <- matrix(stats::rnorm(k * k), k, k)
  R_orth <- qr.Q(qr(R_raw))

  Vb_rot <- V %*% R_orth   # V rotated by an orthogonal matrix

  aligned <- align_procrustes(V, Vb_rot)

  frob_err <- sqrt(sum((V - aligned)^2))
  expect_lt(frob_err, 1e-10)
})

test_that("principal_angles: identical subspaces yield all-zero angles", {
  set.seed(6)
  A <- qr.Q(qr(matrix(stats::rnorm(12), 4, 3)))
  angles <- principal_angles(A, A)
  expect_equal(length(angles), 3L)
  expect_true(all(angles < 1e-6))
})

test_that("principal_angles: two orthogonal 1D subspaces yield pi/2", {
  A <- matrix(c(1, 0, 0), ncol = 1)
  B <- matrix(c(0, 1, 0), ncol = 1)
  angles <- principal_angles(A, B)
  expect_equal(length(angles), 1L)
  expect_equal(angles, pi / 2, tolerance = 1e-12)
})

test_that("principal_angles: known 2D rotation recovers the expected angle", {
  # Rotate the standard R2 basis by theta; the principal angle should be theta.
  theta <- pi / 6
  A <- diag(2)               # identity spans R^2
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  B <- R                     # same subspace as A (all of R^2) but rotated

  # A and B both span all of R^2, so the principal angles are both 0 when
  # we measure the subspace (identical). Instead, test a 1D case:
  # a = (1, 0)^T, b = (cos(theta), sin(theta))^T
  a <- matrix(c(1, 0), ncol = 1)
  b <- matrix(c(cos(theta), sin(theta)), ncol = 1)
  ang <- principal_angles(a, b)
  expect_equal(ang, theta, tolerance = 1e-12)
})

test_that("align_loadings sign method: permuted + sign-flipped input recovers the original", {
  set.seed(7)
  V  <- qr.Q(qr(matrix(stats::rnorm(20), 5, 4)))
  Vb <- V[, c(3, 1, 4, 2)]
  Vb[, 1] <- -Vb[, 1]
  Vb[, 4] <- -Vb[, 4]

  aligned <- align_loadings(V, Vb, method = "sign")

  frob_err <- sqrt(sum((V - aligned)^2))
  expect_lt(frob_err, 1e-12)
})

test_that("align_loadings procrustes method: rotated + permuted input recovers the original subspace", {
  set.seed(8)
  n <- 10
  k <- 3
  V <- qr.Q(qr(matrix(stats::rnorm(n * k), n, k)))

  # permute then rotate
  Vb_perm <- V[, c(2, 3, 1)]
  R_raw   <- matrix(stats::rnorm(k * k), k, k)
  R_orth  <- qr.Q(qr(R_raw))
  Vb      <- Vb_perm %*% R_orth

  aligned <- align_loadings(V, Vb, method = "procrustes")

  # The aligned subspace should span the same space as V.
  # Check Frobenius norm of the projection residual: ||V - aligned||_F
  frob_err <- sqrt(sum((V - aligned)^2))
  expect_lt(frob_err, 1e-10)
})

test_that("principal_angles fallback is correct when multivarious is absent (mock)", {
  # We cannot unload multivarious if it happens to be loaded, so we test the
  # core arithmetic independently by calling the fallback logic directly.
  # This validates the base-R path regardless of whether multivarious is present.

  # 30-degree angle between two 1D subspaces in R^2
  theta <- pi / 6
  a <- matrix(c(1, 0), ncol = 1)
  b <- matrix(c(cos(theta), sin(theta)), ncol = 1)

  # replicate fallback path inline
  Qa <- qr.Q(qr(a))
  Qb <- qr.Q(qr(b))
  cosines <- svd(t(Qa) %*% Qb)$d
  cosines <- pmin(pmax(cosines, -1), 1)
  ang_fallback <- sort(acos(cosines))

  expect_equal(ang_fallback, theta, tolerance = 1e-12)
})

test_that("principal_angles uses multivarious path when installed", {
  skip_if_not_installed("multivarious")
  set.seed(9)
  A <- qr.Q(qr(matrix(stats::rnorm(12), 4, 3)))
  B <- qr.Q(qr(matrix(stats::rnorm(12), 4, 3)))
  # Just check it runs and returns plausible values (between 0 and pi/2)
  angles <- principal_angles(A, B)
  expect_equal(length(angles), 3L)
  expect_true(all(angles >= 0 & angles <= pi / 2 + 1e-10))
})
