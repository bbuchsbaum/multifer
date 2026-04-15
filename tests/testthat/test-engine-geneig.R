centered_orthonormal_basis_geneig <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(scale(M, center = TRUE, scale = FALSE)))
}

orthonormal_basis_geneig <- function(n, k) {
  M <- matrix(stats::rnorm(n * k), nrow = n, ncol = k)
  qr.Q(qr(M))
}

make_rank2_lda_fixture <- function(n_per_class = 28L, seed = 1L) {
  set.seed(seed)
  means <- rbind(
    c(-3.5, -2.5, 0, 0, 0),
    c(-3.5,  2.5, 0, 0, 0),
    c( 3.5, -2.5, 0, 0, 0),
    c( 3.5,  2.5, 0, 0, 0)
  )
  y <- factor(rep(seq_len(nrow(means)), each = n_per_class))
  X <- do.call(rbind, lapply(seq_len(nrow(means)), function(i) {
    noise <- matrix(stats::rnorm(n_per_class * ncol(means), sd = 0.35),
                    nrow = n_per_class, ncol = ncol(means))
    sweep(noise, 2L, means[i, ], `+`)
  }))
  list(X = X, y = y)
}

make_exact_lda_pencil <- function(class_means, within_metric) {
  grand_mean <- colMeans(class_means)
  between <- matrix(0, nrow = ncol(class_means), ncol = ncol(class_means))
  for (i in seq_len(nrow(class_means))) {
    shift <- matrix(class_means[i, ] - grand_mean, ncol = 1L)
    between <- between + shift %*% t(shift)
  }
  geneig_operator(
    A = (between + t(between)) / 2,
    B = within_metric,
    metric = "within_class"
  )
}

make_geneig_null_fixture <- function(n_per_class = 12L, p = 5L, seed = 1L) {
  set.seed(seed)
  y <- factor(rep(seq_len(3L), each = n_per_class))
  X <- matrix(stats::rnorm(length(y) * p), nrow = length(y), ncol = p)
  list(X = X, y = y)
}

make_exact_geneig_fixture <- function(lambda, B_metric, seed = 1L) {
  set.seed(seed)
  p <- length(lambda)
  U <- qr.Q(qr(matrix(stats::rnorm(p * p), nrow = p, ncol = p)))
  eig_B <- eigen((B_metric + t(B_metric)) / 2, symmetric = TRUE)
  B_half <- eig_B$vectors %*%
    diag(sqrt(eig_B$values), nrow = length(eig_B$values)) %*%
    t(eig_B$vectors)

  A <- B_half %*% U %*% diag(lambda, nrow = p, ncol = p) %*% t(U) %*% B_half
  geneig_operator((A + t(A)) / 2, B_metric, metric = "within_class")
}

make_geneig_zero_null_adapter <- function(id = "stub_geneig_zero_null") {
  residualize_geneig <- function(x, k, data) data
  attr(residualize_geneig, "b_metric") <- TRUE

  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "geneig",
    capabilities    = capability_matrix(
      list(
        geometry = "geneig",
        relation = "generalized_eigen",
        targets = "component_significance"
      )
    ),
    null_action    = function(x, data) {
      list(
        A = matrix(0, nrow = nrow(data$A), ncol = ncol(data$A)),
        B = data$B,
        metric = data$metric,
        state = data$state
      )
    },
    component_stat = function(x, data, k) 1,
    residualize    = residualize_geneig,
    refit          = function(x, new_data) x,
    validity_level = "exact"
  )
}

euclidean_deflate_geneig_state <- function(data) {
  C <- .geneig_whiten(data)
  eig <- eigen(C, symmetric = TRUE)
  z1 <- eig$vectors[, 1L, drop = FALSE]
  v1 <- data$B_inv_half %*% z1
  scale <- as.numeric(sqrt(t(v1) %*% data$B %*% v1))
  v1 <- v1 / scale
  P_e <- diag(nrow(data$A)) - v1 %*% t(v1)
  .geneig_normalize_state(list(
    A = (P_e %*% data$A %*% P_e + t(P_e %*% data$A %*% P_e)) / 2,
    B = data$B,
    metric = data$metric,
    state = data$state
  ))
}

make_geneig_oneblock_bridge_adapter <- function(id = "stub_geneig_oneblock_bridge") {
  residualize_geneig <- function(x, k, data) data
  attr(residualize_geneig, "b_metric") <- TRUE

  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "geneig",
    capabilities    = capability_matrix(
      list(
        geometry = "geneig",
        relation = "generalized_eigen",
        targets = "component_significance"
      )
    ),
    null_action    = function(x, data) {
      X <- data$state$factor_matrix
      X_perm <- base::apply(X, 2L, base::sample)
      list(
        A = crossprod(X_perm),
        B = data$B,
        metric = data$metric,
        state = list(factor_matrix = X_perm)
      )
    },
    component_stat = function(x, data, k) 1,
    residualize    = residualize_geneig,
    refit          = function(x, new_data) x,
    validity_level = "exact"
  )
}

test_that("run_geneig_ladder validates recipe and operator inputs", {
  skip_if_not_installed("MASS")
  adapter <- make_geneig_zero_null_adapter()
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter
  )
  oneblock_adapter <- infer_adapter(
    adapter_id      = "stub_pca_geneig_guard",
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(
        geometry = "oneblock",
        relation = "variance",
        targets = "component_significance"
      )
    ),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1,
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) x,
    validity_level = "conditional"
  )
  expect_error(
    run_geneig_ladder(
      infer_recipe(geometry = "oneblock", relation = "variance",
                   adapter = oneblock_adapter),
      geneig_operator(diag(2), diag(2))
    ),
    "geometry must be"
  )
  expect_error(
    run_geneig_ladder(rec, diag(2)),
    "must be a multifer_geneig_operator"
  )
})

test_that("run_geneig_ladder returns expected schema on a planted rank-k pencil", {
  skip_if_not_installed("MASS")
  adapter <- make_geneig_zero_null_adapter()
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter
  )
  op <- make_exact_geneig_fixture(
    lambda = c(9, 4, 1, 0, 0),
    B_metric = matrix(c(
      2.0, 0.3, 0.0, 0.0, 0.0,
      0.3, 1.6, 0.2, 0.0, 0.0,
      0.0, 0.2, 1.4, 0.1, 0.0,
      0.0, 0.0, 0.1, 1.3, 0.2,
      0.0, 0.0, 0.0, 0.2, 1.2
    ), nrow = 5, byrow = TRUE),
    seed = 11L
  )

  res <- run_geneig_ladder(rec, op, B = 39L, alpha = 0.05, seed = 7L)

  expect_true(is.list(res))
  expect_true(all(c("units", "component_tests", "roots_observed",
                    "ladder_result") %in% names(res)))
  expect_equal(res$ladder_result$rejected_through, 3L)
  expect_equal(res$ladder_result$last_step_tested, 4L)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(res$roots_observed[1:5], c(9, 4, 1, 0, 0), tolerance = 1e-10)
})

test_that("run_geneig_ladder matches run_oneblock_ladder when B = I", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()

  set.seed(42L)
  U <- centered_orthonormal_basis_geneig(40, 2)
  V <- orthonormal_basis_geneig(12, 2)
  X_raw <- U %*% diag(c(14, 6), nrow = 2, ncol = 2) %*% t(V)
  X_centered <- sweep(X_raw, 2L, colMeans(X_raw), "-")

  rec_one <- infer_recipe(
    geometry = "oneblock",
    relation = "variance",
    adapter = "prcomp_oneblock"
  )
  rec_geneig <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = make_geneig_oneblock_bridge_adapter()
  )

  oneblock <- run_oneblock_ladder(rec_one, X_raw, B = 99L, alpha = 0.05, seed = 7L)
  geneig <- run_geneig_ladder(
    rec_geneig,
    geneig_operator(crossprod(X_centered), diag(ncol(X_centered)), metric = "identity"),
    state = list(factor_matrix = X_centered),
    B = 99L,
    alpha = 0.05,
    seed = 7L
  )

  expect_identical(geneig$roots_observed, oneblock$roots_observed)
  expect_identical(geneig$component_tests$observed_stat,
                   oneblock$component_tests$observed_stat)
  expect_identical(geneig$component_tests$p_value,
                   oneblock$component_tests$p_value)
  expect_identical(geneig$component_tests$selected,
                   oneblock$component_tests$selected)
  expect_equal(geneig$ladder_result$rejected_through,
               oneblock$ladder_result$rejected_through)
})

test_that("run_geneig_ladder recovers a planted rank-2 four-class LDA pencil", {
  skip_if_not_installed("MASS")
  adapter <- make_geneig_zero_null_adapter(id = "stub_geneig_zero_null_rank2")
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter
  )
  op <- make_exact_lda_pencil(
    class_means = rbind(
      c(-3, -2, 0, 0, 0),
      c(-3,  2, 0, 0, 0),
      c( 3, -2, 0, 0, 0),
      c( 3,  2, 0, 0, 0)
    ),
    within_metric = matrix(c(
      1.7, 0.2, 0.0, 0.0, 0.0,
      0.2, 1.4, 0.1, 0.0, 0.0,
      0.0, 0.1, 1.3, 0.2, 0.0,
      0.0, 0.0, 0.2, 1.2, 0.1,
      0.0, 0.0, 0.0, 0.1, 1.1
    ), nrow = 5, byrow = TRUE)
  )

  res <- run_geneig_ladder(
    rec,
    op,
    B = 79L,
    alpha = 0.05,
    seed = 9L
  )

  expect_equal(res$ladder_result$rejected_through, 2L)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, FALSE))
  expect_lte(max(res$roots_observed[-(1:2)]), 1e-8)
  expect_gte(3L, sum(res$roots_observed > 1e-8))
})

test_that("run_geneig_ladder uses B-metric deflation rather than Euclidean deflation", {
  skip_if_not_installed("MASS")

  op <- make_exact_geneig_fixture(
    lambda = c(9, 4, 1, 0, 0),
    B_metric = matrix(c(
      2.0, 0.4, 0.0, 0.0, 0.0,
      0.4, 1.8, 0.2, 0.0, 0.0,
      0.0, 0.2, 1.5, 0.3, 0.0,
      0.0, 0.0, 0.3, 1.4, 0.2,
      0.0, 0.0, 0.0, 0.2, 1.2
    ), nrow = 5, byrow = TRUE),
    seed = 33L
  )
  state0 <- .geneig_normalize_state(op)
  state_b <- .geneig_deflate_state(state0, zero_tol = 1e-12)
  state_e <- euclidean_deflate_geneig_state(state0)

  expect_equal(.geneig_roots(state_b), c(4, 1, 0, 0, 0), tolerance = 1e-8)
  expect_gt(max(abs(.geneig_roots(state_e) - c(4, 1, 0, 0, 0))), 1e-3)
  expect_gt(abs(.geneig_observed_stat(state_e, zero_tol = 1e-12) -
                  .geneig_observed_stat(state_b, zero_tol = 1e-12)),
            1e-3)
})

test_that("run_geneig_ladder is invariant under orthogonal congruence transforms", {
  skip_if_not_installed("MASS")
  adapter <- make_geneig_zero_null_adapter(id = "stub_geneig_zero_null_congruence")
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter
  )

  op <- make_exact_geneig_fixture(
    lambda = c(9, 4, 1, 0, 0),
    B_metric = matrix(c(
      2.0, 0.3, 0.0, 0.0, 0.0,
      0.3, 1.6, 0.2, 0.0, 0.0,
      0.0, 0.2, 1.4, 0.1, 0.0,
      0.0, 0.0, 0.1, 1.3, 0.2,
      0.0, 0.0, 0.0, 0.2, 1.2
    ), nrow = 5, byrow = TRUE),
    seed = 17L
  )

  set.seed(91L)
  Q <- qr.Q(qr(matrix(stats::rnorm(25L), nrow = 5L, ncol = 5L)))
  op_rot <- geneig_operator(
    A = t(Q) %*% op$A %*% Q,
    B = t(Q) %*% op$B %*% Q,
    metric = op$metric
  )

  res_base <- run_geneig_ladder(rec, op, B = 59L, alpha = 0.05, seed = 12L)
  res_rot <- run_geneig_ladder(rec, op_rot, B = 59L, alpha = 0.05, seed = 12L)

  expect_equal(res_rot$roots_observed, res_base$roots_observed, tolerance = 1e-10)
  expect_equal(res_rot$component_tests$observed_stat,
               res_base$component_tests$observed_stat,
               tolerance = 1e-12)
  expect_identical(res_rot$component_tests$p_value, res_base$component_tests$p_value)
  expect_identical(res_rot$component_tests$selected, res_base$component_tests$selected)
})

test_that("run_geneig_ladder is invariant under common positive scaling of A and B", {
  skip_if_not_installed("MASS")
  adapter <- make_geneig_zero_null_adapter(id = "stub_geneig_zero_null_scaled")
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter
  )

  op <- make_exact_geneig_fixture(
    lambda = c(9, 4, 1, 0, 0),
    B_metric = matrix(c(
      2.0, 0.3, 0.0, 0.0, 0.0,
      0.3, 1.6, 0.2, 0.0, 0.0,
      0.0, 0.2, 1.4, 0.1, 0.0,
      0.0, 0.0, 0.1, 1.3, 0.2,
      0.0, 0.0, 0.0, 0.2, 1.2
    ), nrow = 5, byrow = TRUE),
    seed = 27L
  )
  scale_factor <- 13.5
  op_scaled <- geneig_operator(
    A = scale_factor * op$A,
    B = scale_factor * op$B,
    metric = op$metric
  )

  res_base <- run_geneig_ladder(rec, op, B = 59L, alpha = 0.05, seed = 13L)
  res_scaled <- run_geneig_ladder(rec, op_scaled, B = 59L, alpha = 0.05, seed = 13L)

  expect_equal(res_scaled$roots_observed, res_base$roots_observed, tolerance = 1e-10)
  expect_equal(res_scaled$component_tests$observed_stat,
               res_base$component_tests$observed_stat,
               tolerance = 1e-12)
  expect_identical(res_scaled$component_tests$p_value, res_base$component_tests$p_value)
  expect_identical(res_scaled$component_tests$selected, res_base$component_tests$selected)
})

test_that("run_geneig_ladder stays finite for a near-singular but SPD metric", {
  skip_if_not_installed("MASS")
  adapter <- make_geneig_zero_null_adapter(id = "stub_geneig_zero_null_near_singular")
  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = adapter
  )

  near_singular_B <- diag(c(1, 1e-5, 5e-6, 2, 3))
  op <- make_exact_geneig_fixture(
    lambda = c(5, 2, 0.5, 0, 0),
    B_metric = near_singular_B,
    seed = 41L
  )

  res <- run_geneig_ladder(rec, op, B = 39L, alpha = 0.05, seed = 5L)

  expect_true(all(is.finite(res$roots_observed)))
  expect_true(all(is.finite(res$component_tests$observed_stat)))
  expect_true(all(is.finite(res$component_tests$p_value)))
  expect_true(all(res$component_tests$p_value >= 0 & res$component_tests$p_value <= 1))
})

test_that("label-permutation null is approximately calibrated at rung 1", {
  skip_if_not_installed("MASS")
  ensure_default_adapters()

  rec <- infer_recipe(
    geometry = "geneig",
    relation = "generalized_eigen",
    adapter = "lda_refit"
  )
  alpha <- 0.05
  n_rep <- 32L

  rejected <- vapply(seq_len(n_rep), function(i) {
    dat <- make_geneig_null_fixture(seed = 1000L + i)
    res <- run_geneig_ladder(
      rec,
      .lda_geneig_operator(dat$X, dat$y),
      state = dat,
      B = 511L,
      alpha = alpha,
      seed = 2000L + i
    )
    isTRUE(res$component_tests$selected[1L])
  }, logical(1L))

  fpr <- mean(rejected)
  se <- sqrt(alpha * (1 - alpha) / n_rep)
  expect_lte(abs(fpr - alpha), 2 * se)
})
