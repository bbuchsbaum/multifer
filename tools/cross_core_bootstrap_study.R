#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- if (length(file_arg) > 0L) normalizePath(file_arg[1L]) else NA_character_
root <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..")) else getwd()

user_args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- user_args[startsWith(user_args, prefix)]
  if (length(hit) == 0L) return(default)
  sub(prefix, "", hit[[1L]])
}

n_boot <- as.integer(arg_value("n-boot", "80"))
outdir <- arg_value("outdir", file.path(root, "tools", "results", "cross_core_bootstrap"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  pkgload::load_all(root, quiet = TRUE)
})

simulate_cross_covariance <- function(n, p_x, p_y, signal, noise_x, noise_y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- length(signal)
  latent <- qr.Q(qr(scale(matrix(stats::rnorm(n * k), nrow = n, ncol = k),
                            center = TRUE, scale = FALSE)))
  Wx <- qr.Q(qr(matrix(stats::rnorm(p_x * k), nrow = p_x, ncol = k)))
  Wy <- qr.Q(qr(matrix(stats::rnorm(p_y * k), nrow = p_y, ncol = k)))

  X <- latent %*% diag(signal, nrow = k, ncol = k) %*% t(Wx) +
    matrix(stats::rnorm(n * p_x, sd = noise_x), nrow = n, ncol = p_x)
  Y <- latent %*% diag(signal, nrow = k, ncol = k) %*% t(Wy) +
    matrix(stats::rnorm(n * p_y, sd = noise_y), nrow = n, ncol = p_y)

  list(
    X = scale(X, center = TRUE, scale = FALSE),
    Y = scale(Y, center = TRUE, scale = FALSE)
  )
}

scenario_grid <- function() {
  list(
    list(name = "balanced_separated", n = 150L, p_x = 40L, p_y = 40L,
         signal = c(6, 2), noise_x = 1.0, noise_y = 1.0),
    list(name = "balanced_near_tied", n = 150L, p_x = 40L, p_y = 40L,
         signal = c(4.0, 3.6), noise_x = 1.0, noise_y = 1.0),
    list(name = "wide_y_noisy_y", n = 150L, p_x = 20L, p_y = 120L,
         signal = c(4.0, 3.6), noise_x = 1.0, noise_y = 2.5),
    list(name = "wide_x_noisy_x", n = 150L, p_x = 120L, p_y = 20L,
         signal = c(4.0, 3.6), noise_x = 2.5, noise_y = 1.0),
    list(name = "small_n_wide_near_tied", n = 60L, p_x = 80L, p_y = 80L,
         signal = c(4.0, 3.6), noise_x = 1.5, noise_y = 1.5),
    list(name = "large_n_balanced", n = 300L, p_x = 60L, p_y = 60L,
         signal = c(6.0, 2.0), noise_x = 1.0, noise_y = 1.0)
  )
}

fit_reconstruction <- function(fit) {
  fit$Wx %*% diag(fit$d, nrow = length(fit$d), ncol = length(fit$d)) %*% t(fit$Wy)
}

active_root_count <- function(d, tol = 1e-10) {
  if (length(d) == 0L) return(0L)
  thresh <- max(1, d[[1L]]) * tol
  sum(d > thresh)
}

take_cols <- function(M, k) {
  if (is.null(M)) return(NULL)
  if (k <= 0L) {
    return(M[, 0L, drop = FALSE])
  }
  M[, seq_len(k), drop = FALSE]
}

max_abs_diff <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(0)
  if (length(a) != length(b)) return(Inf)
  max(abs(as.numeric(a) - as.numeric(b)))
}

artifact_diffs <- function(fast, slow, k_eval) {
  per_rep <- vector("list", fast$R)
  for (b in seq_len(fast$R)) {
    rf <- fast$reps[[b]]
    rs <- slow$reps[[b]]
    k_fast <- active_root_count(rf$fit$d)
    k_slow <- active_root_count(rs$fit$d)
    k_cmp <- min(k_eval, k_fast, k_slow)
    per_rep[[b]] <- data.frame(
      rep = b,
      active_rank_fast = k_fast,
      active_rank_slow = k_slow,
      singular_value_max_abs = max_abs_diff(rf$fit$d[seq_len(k_cmp)], rs$fit$d[seq_len(k_cmp)]),
      slow_trailing_sv_max_abs = if (k_slow > k_cmp) max(abs(rs$fit$d[seq.int(k_cmp + 1L, k_slow)])) else 0,
      center_x_max_abs = max_abs_diff(rf$fit$center_x, rs$fit$center_x),
      center_y_max_abs = max_abs_diff(rf$fit$center_y, rs$fit$center_y),
      cross_recon_max_abs = max_abs_diff(fit_reconstruction(rf$fit), fit_reconstruction(rs$fit)),
      aligned_loading_x_max_abs = max_abs_diff(take_cols(rf$aligned_loadings$X, k_cmp),
                                               take_cols(rs$aligned_loadings$X, k_cmp)),
      aligned_loading_y_max_abs = max_abs_diff(take_cols(rf$aligned_loadings$Y, k_cmp),
                                               take_cols(rs$aligned_loadings$Y, k_cmp)),
      aligned_score_x_max_abs = max_abs_diff(take_cols(rf$aligned_scores$X, k_cmp),
                                             take_cols(rs$aligned_scores$X, k_cmp)),
      aligned_score_y_max_abs = max_abs_diff(take_cols(rf$aligned_scores$Y, k_cmp),
                                             take_cols(rs$aligned_scores$Y, k_cmp)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, per_rep)
}

stability_diffs <- function(fast, slow, data, units, original_fit, adapter) {
  variable_fast <- variable_stability_from_bootstrap(fast, units)
  variable_slow <- variable_stability_from_bootstrap(slow, units)
  score_fast <- score_stability_from_bootstrap(fast, data, units)
  score_slow <- score_stability_from_bootstrap(slow, data, units)
  subspace_fast <- subspace_stability_from_bootstrap(fast, original_fit, adapter, units)
  subspace_slow <- subspace_stability_from_bootstrap(slow, original_fit, adapter, units)

  data.frame(
    variable_estimate_max_abs = max_abs_diff(variable_fast$estimate, variable_slow$estimate),
    variable_stability_max_abs = max_abs_diff(variable_fast$stability, variable_slow$stability),
    score_estimate_max_abs = max_abs_diff(score_fast$estimate, score_slow$estimate),
    score_lower_max_abs = max_abs_diff(score_fast$lower, score_slow$lower),
    score_upper_max_abs = max_abs_diff(score_fast$upper, score_slow$upper),
    subspace_mean_max_abs = max_abs_diff(subspace_fast$principal_angle_mean,
                                         subspace_slow$principal_angle_mean),
    subspace_max_max_abs = max_abs_diff(subspace_fast$principal_angle_max,
                                        subspace_slow$principal_angle_max),
    stringsAsFactors = FALSE
  )
}

run_scenario <- function(sc, n_boot, seed) {
  dat <- simulate_cross_covariance(
    n = sc$n, p_x = sc$p_x, p_y = sc$p_y,
    signal = sc$signal, noise_x = sc$noise_x, noise_y = sc$noise_y,
    seed = seed
  )
  adapter <- adapter_cross_svd()
  recipe <- infer_recipe(geometry = "cross", relation = "covariance", adapter = adapter)
  original_fit <- adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))
  k_units <- min(length(sc$signal), active_root_count(adapter$roots(original_fit)))
  units <- form_units(adapter$roots(original_fit)[seq_len(k_units)])

  fast_time <- system.time({
    fast <- bootstrap_fits(
      recipe = recipe, adapter = adapter,
      data = list(X = dat$X, Y = dat$Y),
      original_fit = original_fit, units = units,
      R = n_boot, method_align = "sign", seed = seed + 1000L,
      fast_path = "auto", core_rank = NULL
    )
  })[["elapsed"]]
  slow_time <- system.time({
    slow <- bootstrap_fits(
      recipe = recipe, adapter = adapter,
      data = list(X = dat$X, Y = dat$Y),
      original_fit = original_fit, units = units,
      R = n_boot, method_align = "sign", seed = seed + 1000L,
      fast_path = "off", core_rank = NULL
    )
  })[["elapsed"]]

  rep_df <- artifact_diffs(fast, slow, k_eval = k_units)
  stability_df <- stability_diffs(fast, slow, dat, units, original_fit, adapter)

  summary_df <- data.frame(
    scenario = sc$name,
    n = sc$n,
    p_x = sc$p_x,
    p_y = sc$p_y,
    n_boot = n_boot,
    fast_time = fast_time,
    slow_time = slow_time,
    speedup = slow_time / fast_time,
    singular_value_max_abs = max(rep_df$singular_value_max_abs),
    slow_trailing_sv_max_abs = max(rep_df$slow_trailing_sv_max_abs),
    center_x_max_abs = max(rep_df$center_x_max_abs),
    center_y_max_abs = max(rep_df$center_y_max_abs),
    cross_recon_max_abs = max(rep_df$cross_recon_max_abs),
    aligned_loading_x_max_abs = max(rep_df$aligned_loading_x_max_abs),
    aligned_loading_y_max_abs = max(rep_df$aligned_loading_y_max_abs),
    aligned_score_x_max_abs = max(rep_df$aligned_score_x_max_abs),
    aligned_score_y_max_abs = max(rep_df$aligned_score_y_max_abs),
    variable_estimate_max_abs = stability_df$variable_estimate_max_abs,
    variable_stability_max_abs = stability_df$variable_stability_max_abs,
    score_estimate_max_abs = stability_df$score_estimate_max_abs,
    score_lower_max_abs = stability_df$score_lower_max_abs,
    score_upper_max_abs = stability_df$score_upper_max_abs,
    subspace_mean_max_abs = stability_df$subspace_mean_max_abs,
    subspace_max_max_abs = stability_df$subspace_max_max_abs,
    stringsAsFactors = FALSE
  )

  list(summary = summary_df, reps = transform(rep_df, scenario = sc$name))
}

run_study <- function(n_boot, outdir) {
  scenarios <- scenario_grid()
  summary_rows <- vector("list", length(scenarios))
  rep_rows <- vector("list", length(scenarios))

  for (i in seq_along(scenarios)) {
    res <- run_scenario(scenarios[[i]], n_boot = n_boot, seed = 5000L + i)
    summary_rows[[i]] <- res$summary
    rep_rows[[i]] <- res$reps
  }

  summary_df <- do.call(rbind, summary_rows)
  rep_df <- do.call(rbind, rep_rows)

  utils::write.csv(summary_df, file.path(outdir, "cross_core_bootstrap_summary.csv"), row.names = FALSE)
  utils::write.csv(rep_df, file.path(outdir, "cross_core_bootstrap_per_rep.csv"), row.names = FALSE)

  cat("\nCross core bootstrap summary:\n")
  print(summary_df, row.names = FALSE)

  invisible(list(summary = summary_df, reps = rep_df))
}

run_study(n_boot = n_boot, outdir = outdir)
