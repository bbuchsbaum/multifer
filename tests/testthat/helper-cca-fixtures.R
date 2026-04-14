# Shared CCA fixture builders and stress-matrix helpers.
#
# Used by test-cca-support-matrix.R (the primary support-boundary
# regression) and test-cca-parity.R (the cross-adapter parity sweep
# in multifer-aai.2).  Keeping these in a testthat helper file avoids
# sourcing one test file from another.

make_exact_canonical_correlation_support <- function(n, canonical_corrs,
                                                     scale_x = NULL,
                                                     scale_y = NULL,
                                                     extra_scale_x = 1,
                                                     extra_scale_y = 1) {
  k <- length(canonical_corrs)
  if (is.null(scale_x)) scale_x <- rep(1, k)
  if (is.null(scale_y)) scale_y <- rep(1, k)
  stopifnot(length(scale_x) == k, length(scale_y) == k)

  centered_orthonormal <- function(n, k) {
    M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
    qr.Q(qr(scale(M, center = TRUE, scale = FALSE)))
  }

  B <- centered_orthonormal(n, 2L * k + 2L)
  A <- B[, seq_len(k), drop = FALSE]
  C <- B[, k + seq_len(k), drop = FALSE]
  extra_x <- B[, 2L * k + 1L, drop = FALSE]
  extra_y <- B[, 2L * k + 2L, drop = FALSE]

  Qx_shared <- A
  Qy_shared <- A %*% diag(canonical_corrs, nrow = k, ncol = k) +
    C %*% diag(sqrt(1 - canonical_corrs^2), nrow = k, ncol = k)

  X <- cbind(Qx_shared %*% diag(scale_x, nrow = k, ncol = k),
             extra_scale_x * extra_x)
  Y <- cbind(Qy_shared %*% diag(scale_y, nrow = k, ncol = k),
             extra_scale_y * extra_y)

  list(X = X, Y = Y, canonical_corrs = c(canonical_corrs, 0))
}

make_nuisance_adjusted_correlation_support <- function(n, canonical_corrs,
                                                       scale_x = NULL,
                                                       scale_y = NULL,
                                                       groups = NULL) {
  z1 <- scale(seq_len(n))
  z2 <- cos(seq_len(n) / 4)
  Z <- cbind(1, z1, z2)

  if (is.null(groups)) {
    qz_full <- qr.Q(qr(Z), complete = TRUE)
    Qz <- qz_full[, seq.int(ncol(Z) + 1L, n), drop = FALSE]
    base <- make_exact_canonical_correlation_support(
      n = ncol(Qz),
      canonical_corrs = canonical_corrs,
      scale_x = scale_x,
      scale_y = scale_y,
      extra_scale_x = 0.8,
      extra_scale_y = 1.1
    )
    basis <- Qz
    reduced_groups <- NULL
  } else {
    resid_basis <- .nuisance_residual_basis(Z, n = n, groups = groups)
    base <- make_exact_canonical_correlation_support(
      n = ncol(resid_basis$Q),
      canonical_corrs = canonical_corrs,
      scale_x = scale_x,
      scale_y = scale_y,
      extra_scale_x = 0.8,
      extra_scale_y = 1.1
    )
    basis <- resid_basis$Q
    reduced_groups <- resid_basis$groups
  }

  bx <- matrix(c(2.0, -1.5, 1.0), nrow = ncol(Z), ncol = 1L)
  by <- matrix(c(-1.0, 1.2, -0.8), nrow = ncol(Z), ncol = 1L)

  list(
    X = basis %*% base$X + Z %*% bx %*% matrix(1, nrow = 1L, ncol = ncol(base$X)),
    Y = basis %*% base$Y + Z %*% by %*% matrix(1, nrow = 1L, ncol = ncol(base$Y)),
    Z = Z,
    groups = if (is.null(groups)) NULL else as.factor(groups),
    reduced_groups = reduced_groups
  )
}

expect_multiroot_support <- function(design, X, Y, label) {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = design,
    adapter = "cross_svd"
  )

  res <- run_cross_ladder(
    rec, X, Y,
    B = 39L, alpha = 0.05, seed = 7L, max_steps = 4L
  )

  expect_equal(res$ladder_result$rejected_through, 3L, info = label)
  expect_equal(res$ladder_result$last_step_tested, 4L, info = label)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, TRUE, FALSE), info = label)
}

expect_multiroot_walks_ladder <- function(design, X, Y, label,
                                          expected_through = 1L,
                                          max_steps = 3L) {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = design,
    adapter = "cross_svd"
  )
  res <- run_cross_ladder(
    rec, X, Y,
    B = 39L, alpha = 0.05, seed = 11L, max_steps = max_steps
  )
  expect_gt(res$ladder_result$last_step_tested, 1L)
  expect_gte(res$ladder_result$rejected_through, expected_through)
}

expect_first_root_cap <- function(X, Y, label) {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = exchangeable_rows(),
    adapter = "cross_svd"
  )
  res <- run_cross_ladder(
    rec, X, Y,
    B = 39L, alpha = 0.05, seed = 11L, max_steps = 3L
  )
  expect_equal(res$ladder_result$last_step_tested, 1L, info = label)
  expect_equal(nrow(res$component_tests), 1L, info = label)
}

cca_stress_fixtures <- function() {
  list(
    near_tied = local({
      set.seed(510)
      make_exact_canonical_correlation_support(
        n = 48,
        canonical_corrs = c(0.92, 0.91, 0.40),
        scale_x = c(4, 4, 1),
        scale_y = c(3, 3, 0.5),
        extra_scale_x = 1.4,
        extra_scale_y = 1.1
      )
    }),
    weak_trailing = local({
      set.seed(511)
      make_exact_canonical_correlation_support(
        n = 60,
        canonical_corrs = c(0.95, 0.70, 0.18),
        scale_x = c(5, 2, 0.6),
        scale_y = c(4, 3, 0.3),
        extra_scale_x = 1.2,
        extra_scale_y = 1.0
      )
    }),
    wide_x_narrow_y = local({
      set.seed(512)
      make_exact_canonical_correlation_support(
        n = 60,
        canonical_corrs = c(0.95, 0.85, 0.55),
        scale_x = c(5, 3, 1.5),
        scale_y = c(3, 2, 1),
        extra_scale_x = 2.0,
        extra_scale_y = 0.8
      )
    }),
    wide_y_narrow_x = local({
      set.seed(513)
      make_exact_canonical_correlation_support(
        n = 60,
        canonical_corrs = c(0.95, 0.85, 0.55),
        scale_x = c(3, 2, 1),
        scale_y = c(5, 3, 1.5),
        extra_scale_x = 0.8,
        extra_scale_y = 2.0
      )
    })
  )
}
