# CCA support-matrix regression tests
#
# These tests pin the shipped support boundary for the `(cross,
# correlation)` family. Supported designs should permit multi-root
# inference. Designs outside that matrix must conservatively cap the
# ladder to the first root.

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

test_that("CCA support matrix: shipped supported designs permit multi-root recovery", {
  ensure_default_adapters()

  set.seed(402)
  paired <- make_exact_canonical_correlation_support(
    n = 36,
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5),
    extra_scale_x = 1.5,
    extra_scale_y = 1.2
  )
  expect_multiroot_support(paired_rows(), paired$X, paired$Y, "paired_rows")

  set.seed(403)
  nuisance <- make_nuisance_adjusted_correlation_support(
    n = 36,
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5)
  )
  expect_multiroot_support(
    nuisance_adjusted(nuisance$Z),
    nuisance$X, nuisance$Y,
    "nuisance_adjusted"
  )

  set.seed(406)
  groups <- rep(1:8, each = 6)
  nuisance_grouped <- make_nuisance_adjusted_correlation_support(
    n = length(groups),
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5),
    groups = groups
  )
  expect_multiroot_support(
    nuisance_adjusted(nuisance_grouped$Z, groups = nuisance_grouped$groups),
    nuisance_grouped$X, nuisance_grouped$Y,
    "nuisance_adjusted(groups)"
  )

  set.seed(402)
  blocked <- make_exact_canonical_correlation_support(
    n = 36,
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5),
    extra_scale_x = 1.5,
    extra_scale_y = 1.2
  )
  expect_multiroot_support(
    blocked_rows(rep(1:6, each = 6)),
    blocked$X, blocked$Y,
    "blocked_rows"
  )
})

test_that("CCA support matrix: unsupported designs cap the ladder to the first root", {
  ensure_default_adapters()
  set.seed(403)

  dat <- make_exact_canonical_correlation_support(
    n = 36,
    canonical_corrs = c(0.95, 0.7, 0.4),
    scale_x = c(5, 2, 1),
    scale_y = c(4, 3, 0.5),
    extra_scale_x = 1.5,
    extra_scale_y = 1.2
  )

  rec <- infer_recipe(
    geometry = "cross",
    relation = "correlation",
    design = exchangeable_rows(),
    adapter = "cross_svd"
  )
  res <- run_cross_ladder(
    rec, dat$X, dat$Y,
    B = 39L, alpha = 0.05, seed = 7L, max_steps = 4L
  )

  expect_equal(res$ladder_result$last_step_tested, 1L)
  expect_equal(nrow(res$component_tests), 1L)
  expect_true(res$component_tests$selected[[1L]])
  expect_true(sum(res$units$selected) <= 1L)
})

test_that("CCA support matrix: nuisance-adjusted grouped null stays near nominal alpha", {
  skip_on_cran()
  ensure_default_adapters()

  n_rep <- 40L
  hits <- integer(n_rep)
  groups <- rep(1:8, each = 6)

  for (i in seq_len(n_rep)) {
    dat <- bench_cross_null(
      n = length(groups), p_x = 6, p_y = 5,
      within_rank_x = 2, within_rank_y = 2,
      seed = 1200L + i
    )
    z1 <- scale(seq_len(length(groups)))
    z2 <- sin(seq_len(length(groups)) / 5)
    Z <- cbind(1, z1, z2)

    rec <- infer_recipe(
      geometry = "cross",
      relation = "correlation",
      design = nuisance_adjusted(Z, groups = groups),
      adapter = "cross_svd"
    )
    res <- run_cross_ladder(
      rec, dat$X, dat$Y,
      B = 39L, alpha = 0.05, seed = 2000L + i
    )
    hits[i] <- as.integer(res$ladder_result$rejected_through > 0L)
  }

  fpr <- mean(hits)
  expect_true(fpr >= 0 && fpr <= 0.2,
              info = sprintf("grouped nuisance-adjusted empirical FPR = %.3f", fpr))
})
