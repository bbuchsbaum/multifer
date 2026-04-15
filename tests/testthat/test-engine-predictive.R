make_predictive_fixture <- function(n = 72L,
                                    p = 10L,
                                    q = 4L,
                                    rank = 2L,
                                    x_strength = c(3.6, 3.1),
                                    y_strength = c(3.8, 3.2),
                                    noise_x = 0.10,
                                    noise_y = 0.08,
                                    seed = 1L) {
  set.seed(seed)
  rank <- as.integer(rank)
  Tx <- if (rank > 0L) qr.Q(qr(matrix(rnorm(n * rank), nrow = n, ncol = rank))) else NULL
  Px <- if (rank > 0L) qr.Q(qr(matrix(rnorm(p * rank), nrow = p, ncol = rank))) else NULL
  Cy <- if (rank > 0L) qr.Q(qr(matrix(rnorm(q * rank), nrow = q, ncol = rank))) else NULL

  X <- matrix(rnorm(n * p, sd = noise_x), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * q, sd = noise_y), nrow = n, ncol = q)

  if (rank > 0L) {
    X <- X + Tx %*% diag(x_strength[seq_len(rank)], nrow = rank) %*% t(Px)
    Y <- Y + Tx %*% diag(y_strength[seq_len(rank)], nrow = rank) %*% t(Cy)
  }

  list(X = X, Y = Y)
}

make_plsr_predictive_fixture <- function(n = 72L,
                                         p = 8L,
                                         q = 4L,
                                         seed = 1L) {
  set.seed(seed)
  scores <- qr.Q(qr(scale(matrix(rnorm(n * 2L), nrow = n, ncol = 2L),
                           center = TRUE, scale = FALSE)))

  x_loadings <- rbind(
    c(1.3, -0.9, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0),
    c(0.0,  0.0, 0.0, 1.1, -0.7, 0.5, 0.0, 0.0)
  )
  y_loadings <- rbind(
    c(1.2, -1.0, 0.0, 0.0),
    c(0.0,  0.0, 0.9, -0.7)
  )

  X <- matrix(rnorm(n * p, sd = 0.10), nrow = n, ncol = p) +
    scores %*% diag(c(3.5, 1.8), nrow = 2L, ncol = 2L) %*% x_loadings
  Y <- matrix(rnorm(n * q, sd = 0.08), nrow = n, ncol = q) +
    scores %*% diag(c(3.2, 1.4), nrow = 2L, ncol = 2L) %*% y_loadings

  list(X = X, Y = Y)
}

fit_predictive_rrr <- function(data, ridge = 1e-8) {
  Xc <- sweep(data$X, 2L, colMeans(data$X), "-")
  Yc <- sweep(data$Y, 2L, colMeans(data$Y), "-")
  xtx <- crossprod(Xc) + diag(ridge, ncol(Xc))
  B_ols <- solve(xtx, crossprod(Xc, Yc))
  Yhat <- Xc %*% B_ols
  sv <- svd(Yhat)
  keep <- which(sv$d > 1e-10)

  if (length(keep) == 0L) {
    V <- matrix(0, nrow = ncol(Yc), ncol = 0L)
    scores <- matrix(0, nrow = nrow(Xc), ncol = 0L)
  } else {
    V <- sv$v[, keep, drop = FALSE]
    scores <- sv$u[, keep, drop = FALSE] %*%
      diag(sv$d[keep], nrow = length(keep), ncol = length(keep))
  }

  list(
    B_ols = B_ols,
    V = V,
    scores = scores,
    Xc = Xc,
    Yc = Yc,
    x_center = colMeans(data$X),
    y_center = colMeans(data$Y)
  )
}

make_predictive_rrr_adapter <- function(id = "predictive_rrr_stub") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(
        geometry = "cross",
        relation = "predictive",
        targets = "component_significance"
      )
    ),
    residualize = function(x, k, data) {
      if (ncol(x$V) == 0L || k <= 0L) {
        return(list(X = x$Xc, Y = x$Yc))
      }
      kk <- min(k, ncol(x$V))
      Vk <- x$V[, seq_len(kk), drop = FALSE]
      y_hat <- x$Xc %*% x$B_ols %*% Vk %*% t(Vk)
      T_scores <- x$scores[, seq_len(kk), drop = FALSE]
      if (ncol(T_scores) == 0L) {
        x_res <- x$Xc
      } else {
        P_load <- qr.solve(crossprod(T_scores), crossprod(T_scores, x$Xc))
        x_res <- x$Xc - T_scores %*% P_load
      }
      list(X = x_res, Y = x$Yc - y_hat)
    },
    refit = function(x, new_data) fit_predictive_rrr(new_data),
    predict_response = function(x, new_data, k = NULL) {
      Xc_new <- sweep(new_data$X, 2L, x$x_center, "-")
      if (is.null(k)) {
        k <- ncol(x$V)
      }
      if (k <= 0L || ncol(x$V) == 0L) {
        return(matrix(0, nrow = nrow(Xc_new), ncol = length(x$y_center)))
      }
      kk <- min(k, ncol(x$V))
      Vk <- x$V[, seq_len(kk), drop = FALSE]
      Xc_new %*% x$B_ols %*% Vk %*% t(Vk)
    },
    null_action = function(x, data) {
      perm <- sample.int(nrow(data$Y))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },
    component_stat = function(x, data, k, split = NULL) 1,
    validity_level = "conditional"
  )
}

make_nonpredictive_cross_adapter <- function(id = "predictive_guard_covariance") {
  infer_adapter(
    adapter_id      = id,
    adapter_version = "0.0.1",
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(
        geometry = "cross",
        relation = "covariance",
        targets = "component_significance"
      )
    ),
    residualize    = function(x, k, data) data,
    refit          = function(x, new_data) list(),
    null_action    = function(x, data) data,
    component_stat = function(x, data, k) 1,
    validity_level = "conditional"
  )
}

test_that("run_predictive_ladder refuses non-predictive recipes", {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "covariance",
    adapter = make_nonpredictive_cross_adapter()
  )
  dat <- make_predictive_fixture(rank = 1L, seed = 10L)

  expect_error(
    run_predictive_ladder(rec, dat$X, dat$Y),
    "relation must be"
  )
})

test_that("run_predictive_ladder recovers planted supervised rank", {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "predictive",
    adapter = make_predictive_rrr_adapter()
  )
  dat <- make_predictive_fixture(rank = 2L, seed = 101L)

  res <- run_predictive_ladder(
    rec,
    dat$X,
    dat$Y,
    n_folds = 4L,
    B = 95L,
    alpha = 0.05,
    seed = 77L
  )

  expect_true(is.list(res))
  expect_true(all(c("units", "component_tests", "roots_observed",
                    "ladder_result") %in% names(res)))
  expect_equal(res$ladder_result$rejected_through, 2L)
  expect_identical(res$component_tests$selected, c(TRUE, TRUE, FALSE))
  expect_equal(res$roots_observed, res$component_tests$observed_stat)
})

test_that("run_predictive_ladder matches explicit folds to generated fold ids", {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "predictive",
    adapter = make_predictive_rrr_adapter(id = "predictive_rrr_explicit")
  )
  dat <- make_predictive_fixture(rank = 2L, seed = 131L)
  fold_ids <- .fold_resolve_ids(n = nrow(dat$X), n_folds = 4L, seed = 77L)

  res_explicit <- run_predictive_ladder(
    rec,
    dat$X,
    dat$Y,
    folds = fold_ids,
    B = 95L,
    alpha = 0.05,
    seed = 77L
  )

  res_generated <- run_predictive_ladder(
    rec,
    dat$X,
    dat$Y,
    n_folds = 4L,
    B = 95L,
    alpha = 0.05,
    seed = 77L
  )

  expect_equal(res_explicit$roots_observed, res_generated$roots_observed)
  expect_equal(res_explicit$component_tests, res_generated$component_tests)
  expect_equal(res_explicit$ladder_result$rejected_through,
               res_generated$ladder_result$rejected_through)
})

test_that("predictive ladder first-rung rejection rate is calibrated under shuffled Y", {
  rec <- infer_recipe(
    geometry = "cross",
    relation = "predictive",
    adapter = make_predictive_rrr_adapter(id = "predictive_rrr_null")
  )

  alpha <- 0.05
  n_rep <- 24L

  rejected <- vapply(seq_len(n_rep), function(i) {
    dat <- make_predictive_fixture(rank = 2L, seed = 3000L + i)
    set.seed(4000L + i)
    dat$Y <- dat$Y[sample.int(nrow(dat$Y)), , drop = FALSE]

    res <- run_predictive_ladder(
      rec,
      dat$X,
      dat$Y,
      n_folds = 4L,
      B = 127L,
      alpha = alpha,
      seed = 5000L + i
    )

    isTRUE(res$component_tests$selected[1L])
  }, logical(1L))

  fpr <- mean(rejected)
  se <- sqrt(alpha * (1 - alpha) / n_rep)
  expect_lte(abs(fpr - alpha), 2 * se)
})

test_that("plsr_refit recovers planted predictive rank on a supervised fixture", {
  skip_if_not_installed("pls")
  ensure_default_adapters()

  dat <- make_plsr_predictive_fixture(seed = 901L)
  res <- infer_plsr(
    dat$X,
    dat$Y,
    targets = "component_significance",
    B = 39L,
    seed = 19L
  )

  expect_true(is_infer_result(res))
  expect_identical(res$units$selected, c(TRUE, TRUE, FALSE))
  expect_equal(res$provenance$adapter_id, "plsr_refit")
})

test_that("plsr_refit rung-1 null calibration is within 2 SE at B = 1000", {
  skip_on_cran()
  skip_if_not_installed("pls")
  ensure_default_adapters()

  alpha <- 0.05
  n_rep <- 24L

  rejected <- vapply(seq_len(n_rep), function(i) {
    dat <- make_plsr_predictive_fixture(seed = 1200L + i)
    set.seed(2200L + i)
    dat$Y <- dat$Y[sample.int(nrow(dat$Y)), , drop = FALSE]

    res <- infer_plsr(
      dat$X,
      dat$Y,
      targets = "component_significance",
      B = 1000L,
      seed = 3200L + i
    )

    isTRUE(res$component_tests$p_value[1L] <= alpha)
  }, logical(1L))

  fpr <- mean(rejected)
  se <- sqrt(alpha * (1 - alpha) / n_rep)
  expect_lte(abs(fpr - alpha), 2 * se)
})

test_that("plsr_refit stability outputs separate selected signal units from the noise unit", {
  skip_on_cran()
  skip_if_not_installed("pls")
  ensure_default_adapters()

  dat <- make_plsr_predictive_fixture(seed = 1501L)
  res <- infer_plsr(
    dat$X,
    dat$Y,
    B = 29L,
    R = 12L,
    seed = 29L
  )

  expect_identical(res$units$selected, c(TRUE, TRUE, FALSE))

  var_summary <- stats::aggregate(
    stability ~ unit_id + domain,
    data = res$variable_stability,
    FUN = mean
  )
  score_tbl <- res$score_stability
  score_tbl$rel_width <- (score_tbl$upper - score_tbl$lower) /
    (abs(score_tbl$estimate) + 1e-8)
  score_summary <- stats::aggregate(
    rel_width ~ unit_id + domain,
    data = score_tbl,
    FUN = stats::median
  )

  for (domain in c("X", "Y")) {
    dom_var <- var_summary[var_summary$domain == domain, ]
    dom_score <- score_summary[score_summary$domain == domain, ]

    expect_gt(
      dom_var$stability[dom_var$unit_id == "u1"],
      dom_var$stability[dom_var$unit_id == "u3"]
    )
    expect_gt(
      dom_var$stability[dom_var$unit_id == "u2"],
      dom_var$stability[dom_var$unit_id == "u3"]
    )
    expect_lt(
      dom_score$rel_width[dom_score$unit_id == "u1"],
      dom_score$rel_width[dom_score$unit_id == "u3"]
    )
    expect_lt(
      dom_score$rel_width[dom_score$unit_id == "u2"],
      dom_score$rel_width[dom_score$unit_id == "u3"]
    )
  }
})
