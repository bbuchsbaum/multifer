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

# ---------------------------------------------------------------------------
# Paper Theorem 1 parity calibration suite
# ---------------------------------------------------------------------------
# The two tests above prove the P3 collapse on one synthetic matrix each.
# The two tests below promote that to a full calibration suite:
#   (a) a grid of shape / SNR / noise configurations
#   (b) end-to-end ladder-decision equivalence under a shared seed
# This is what the paper's Theorem 1 claim rests on, so it must hold
# across a range of input regimes, not just a single test matrix.

test_that("Vitale P3 parity holds across a grid of synthetic scenarios", {
  scenarios <- list(
    list(name = "strong_signal",          n = 100, p = 20,
         root_profile = c(10, 5, 2, 1),     hetero = FALSE),
    list(name = "near_ties",               n = 150, p = 30,
         root_profile = c(5, 4.9, 4.8, 1),  hetero = FALSE),
    list(name = "shadowing_strong_tail",   n = 120, p = 25,
         root_profile = c(12, 8, 0.5, 0.4, 0.3), hetero = FALSE),
    list(name = "heteroscedastic",         n = 100, p = 20,
         root_profile = c(6, 3, 1.2),       hetero = TRUE),
    list(name = "tall_n_small_p",          n = 400, p = 8,
         root_profile = c(4, 2, 1),         hetero = FALSE),
    list(name = "pure_noise",              n = 80,  p = 15,
         root_profile = NULL,               hetero = FALSE)
  )

  for (sc in scenarios) {
    set.seed(7L + nchar(sc$name))
    X <- if (is.null(sc$root_profile)) {
      matrix(stats::rnorm(sc$n * sc$p), sc$n, sc$p)
    } else {
      .make_shadowing_matrix(
        n = sc$n, p = sc$p, root_profile = sc$root_profile,
        heteroscedastic = sc$hetero, seed = 13L
      )
    }

    max_step <- min(
      max(1L, length(sc$root_profile %||% integer(0))) + 1L,
      min(dim(X)) - 1L,
      5L
    )

    for (step in seq_len(max_step)) {
      expect_equal(
        observed_collapsed(X, step),
        observed_original(X, step),
        tolerance = 1e-9,
        info = sprintf("%s observed step=%d", sc$name, step)
      )

      residual <- deflate_rank1_steps(X, step - 1L)
      for (draw in seq_len(6L)) {
        set.seed(1000L + draw)
        permuted <- base::apply(residual, 2L, base::sample)
        expect_equal(
          vitale_p3_collapsed(permuted, step),
          vitale_p3_original(permuted, step),
          tolerance = 1e-9,
          info = sprintf("%s null step=%d draw=%d", sc$name, step, draw)
        )
      }
    }
  }
})

test_that("P3 parity implies equal ladder decisions under a shared seed", {
  # Full ladder run in both formulations with identical RNG state; the
  # collapsed null and the projected-null P3 must yield identical
  # observed statistics, identical null trajectories, and identical
  # rejection counts across a grid of shapes. This is the end-to-end
  # decision-equivalence claim backing Theorem 1.
  make_callbacks <- function(formulation) {
    list(
      # Observed at rung a: lambda_a(X) / sum_{q >= a} lambda_q(X).
      # By the collapse identity this equals lambda_1(E_{a-1}) /
      # ||E_{a-1}||_F^2 on the current rung's deflated data. The two
      # formulations agree on the observed side; only the null
      # differs.
      observed_stat_fn = function(step, data) {
        tail_ratio_from_s2(svd(data)$d^2, 1L)
      },
      null_stat_fn = function(step, data) {
        permuted <- base::apply(data, 2L, base::sample)
        if (formulation == "collapsed") {
          # Collapsed P3: lambda_a(M_b) / sum_{q >= a} lambda_q(M_b)
          # where M_b is the permuted current-rung data. Uses step-
          # indexing into the permuted SVD spectrum directly.
          return(tail_ratio_from_s2(svd(permuted)$d^2, step))
        }
        # Original projected-null P3: project out the leading step-1
        # left singular vectors of the permuted residual, then take
        # the top-1 ratio of the projected matrix. By the rank-(step-1)
        # removal identity this equals the collapsed form above.
        if (step == 1L) {
          return(tail_ratio_from_s2(svd(permuted)$d^2, 1L))
        }
        sv <- svd(permuted)
        U_drop <- sv$u[, seq_len(step - 1L), drop = FALSE]
        proj <- permuted - U_drop %*% base::crossprod(U_drop, permuted)
        tail_ratio_from_s2(svd(proj)$d^2, 1L)
      },
      deflate_fn = function(step, data) {
        sv <- svd(data)
        data - sv$d[1L] *
          (sv$u[, 1L, drop = FALSE] %*% t(sv$v[, 1L, drop = FALSE]))
      }
    )
  }

  scenarios <- list(
    list(name = "strong_signal",      n = 100, p = 20,
         root_profile = c(9, 5, 2, 1)),
    list(name = "shadowing",          n = 120, p = 25,
         root_profile = c(12, 8, 0.5, 0.4)),
    list(name = "near_ties",          n = 100, p = 18,
         root_profile = c(5, 4.9, 4.8, 1))
  )

  for (sc in scenarios) {
    X <- .make_shadowing_matrix(
      n = sc$n, p = sc$p, root_profile = sc$root_profile,
      heteroscedastic = FALSE, seed = 31L
    )
    Xc <- center_cols(X)

    cb_coll <- make_callbacks("collapsed")
    cb_orig <- make_callbacks("original")

    res_coll <- ladder_driver(
      observed_stat_fn = cb_coll$observed_stat_fn,
      null_stat_fn     = cb_coll$null_stat_fn,
      deflate_fn       = cb_coll$deflate_fn,
      initial_data     = Xc,
      max_steps        = 4L,
      B                = 59L,
      alpha            = 0.05,
      seed             = 97L
    )
    res_orig <- ladder_driver(
      observed_stat_fn = cb_orig$observed_stat_fn,
      null_stat_fn     = cb_orig$null_stat_fn,
      deflate_fn       = cb_orig$deflate_fn,
      initial_data     = Xc,
      max_steps        = 4L,
      B                = 59L,
      alpha            = 0.05,
      seed             = 97L
    )

    expect_equal(res_coll$rejected_through, res_orig$rejected_through,
                 info = sprintf("%s rejected_through", sc$name))
    expect_equal(length(res_coll$step_results),
                 length(res_orig$step_results),
                 info = sprintf("%s steps_tested", sc$name))

    for (i in seq_along(res_coll$step_results)) {
      expect_equal(
        res_coll$step_results[[i]]$observed_stat,
        res_orig$step_results[[i]]$observed_stat,
        tolerance = 1e-12,
        info = sprintf("%s observed_stat step=%d", sc$name, i)
      )
      expect_equal(
        res_coll$step_results[[i]]$r,
        res_orig$step_results[[i]]$r,
        info = sprintf("%s r step=%d", sc$name, i)
      )
      expect_equal(
        res_coll$step_results[[i]]$p_value,
        res_orig$step_results[[i]]$p_value,
        tolerance = 1e-12,
        info = sprintf("%s p_value step=%d", sc$name, i)
      )
    }
  }
})
