center_cols <- function(X) {
  sweep(X, 2L, colMeans(X), "-")
}

.random_orthonormal <- function(nrow, ncol) {
  qr.Q(qr(matrix(stats::rnorm(nrow * ncol), nrow = nrow, ncol = ncol)))
}

.make_shadowing_matrix <- function(n, p, root_profile, heteroscedastic = FALSE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  k <- length(root_profile)
  U <- .random_orthonormal(n, k)
  V <- .random_orthonormal(p, k)
  signal <- U %*% diag(sqrt((n - 1) * root_profile), nrow = k) %*% t(V)

  noise <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  if (heteroscedastic) {
    noise <- sweep(noise, 2L, sqrt(stats::rexp(p, rate = 1)), "*")
  }

  signal + noise
}

deflate_rank1_steps <- function(X, steps) {
  resid <- center_cols(X)
  if (steps <= 0L) {
    return(resid)
  }
  for (i in seq_len(steps)) {
    sv <- svd(resid)
    resid <- resid - sv$d[1L] * (sv$u[, 1L, drop = FALSE] %*% t(sv$v[, 1L, drop = FALSE]))
  }
  resid
}

tail_ratio_from_s2 <- function(s2, k) {
  if (length(s2) < k) {
    return(0)
  }
  denom <- sum(s2[k:length(s2)])
  if (denom <= 0) {
    return(0)
  }
  s2[k] / denom
}

vitale_p3_original <- function(permuted_residual, step) {
  sv <- svd(permuted_residual)
  if (step == 1L) {
    return(tail_ratio_from_s2(sv$d^2, 1L))
  }
  U_drop <- sv$u[, seq_len(step - 1L), drop = FALSE]
  projected <- permuted_residual - U_drop %*% crossprod(U_drop, permuted_residual)
  tail_ratio_from_s2(svd(projected)$d^2, 1L)
}

vitale_p3_collapsed <- function(permuted_residual, step) {
  tail_ratio_from_s2(svd(permuted_residual)$d^2, step)
}

observed_original <- function(X, step) {
  tail_ratio_from_s2(svd(center_cols(X))$d^2, step)
}

observed_collapsed <- function(X, step) {
  tail_ratio_from_s2(svd(deflate_rank1_steps(X, step - 1L))$d^2, 1L)
}

test_that("observed Vitale ladder statistic matches deflated-residual formulation across steps", {
  set.seed(101L)
  X <- .make_shadowing_matrix(
    n = 120,
    p = 30,
    root_profile = c(9, 5, 2.5, 1.25),
    heteroscedastic = FALSE,
    seed = 101L
  )

  for (step in 1:4) {
    expect_equal(
      observed_collapsed(X, step),
      observed_original(X, step),
      tolerance = 1e-10
    )
  }
})

test_that("original projected-null P3 equals collapsed tail-ratio on the same permuted residual", {
  set.seed(202L)
  X <- .make_shadowing_matrix(
    n = 100,
    p = 25,
    root_profile = c(8, 4, 2, 1.2),
    heteroscedastic = TRUE,
    seed = 202L
  )

  for (step in 1:4) {
    residual <- deflate_rank1_steps(X, step - 1L)
    for (draw in 1:8) {
      permuted <- apply(residual, 2L, sample)
      expect_equal(
        vitale_p3_collapsed(permuted, step),
        vitale_p3_original(permuted, step),
        tolerance = 1e-10
      )
    }
  }
})
