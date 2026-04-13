#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg  <- grep("^--file=", args_full, value = TRUE)
script    <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else ""
root      <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)

source(file.path(root, "R", "bench_generators.R"), local = TRUE)

center_cols <- function(X) {
  sweep(X, 2L, colMeans(X), "-")
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

permute_columns <- function(X) {
  apply(X, 2L, sample)
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

engine_shortcut <- function(permuted_residual) {
  tail_ratio_from_s2(svd(permuted_residual)$d^2, 1L)
}

observed_original <- function(X, step) {
  tail_ratio_from_s2(svd(center_cols(X))$d^2, step)
}

observed_collapsed <- function(X, step) {
  tail_ratio_from_s2(svd(deflate_rank1_steps(X, step - 1L))$d^2, 1L)
}

time_eval <- function(permuted_residuals, step, fn) {
  values <- numeric(length(permuted_residuals))
  start <- proc.time()[["elapsed"]]
  for (i in seq_along(permuted_residuals)) {
    values[[i]] <- fn(permuted_residuals[[i]], step)
  }
  elapsed <- proc.time()[["elapsed"]] - start
  list(values = values, elapsed = elapsed)
}

scenario_matrix <- function(kind, seed) {
  switch(
    kind,
    null_small = bench_oneblock_null(
      n = 80, p = 15, noise = "gaussian", seed = seed
    )$X,
    shadowing_medium = bench_oneblock_shadowing(
      n = 160, p = 40, root_profile = c(10, 6, 3, 1.6), noise = "gaussian", seed = seed
    )$X,
    shadowing_hetero = bench_oneblock_shadowing(
      n = 180, p = 45, root_profile = c(8, 5, 2.5, 1.5), noise = "heteroscedastic", seed = seed
    )$X,
    wide_medium = bench_oneblock_shadowing(
      n = 90, p = 60, root_profile = c(12, 8, 4), noise = "gaussian", seed = seed
    )$X,
    stop(sprintf("Unknown scenario: %s", kind), call. = FALSE)
  )
}

run_scenario <- function(kind, B = 48L, n_rep = 3L, max_step = 4L) {
  rows <- vector("list", n_rep * max_step)
  idx  <- 1L

  for (rep_id in seq_len(n_rep)) {
    X <- scenario_matrix(kind, seed = 1000L + rep_id)
    max_valid_step <- min(max_step, min(dim(X)) - 1L)

    for (step in seq_len(max_valid_step)) {
      residual <- deflate_rank1_steps(X, step - 1L)
      permuted <- replicate(B, permute_columns(residual), simplify = FALSE)

      observed_gap <- abs(observed_original(X, step) - observed_collapsed(X, step))

      original  <- time_eval(permuted, step, vitale_p3_original)
      collapsed <- time_eval(permuted, step, vitale_p3_collapsed)
      shortcut  <- time_eval(permuted, step, function(mat, ignored_step) {
        engine_shortcut(mat)
      })

      rows[[idx]] <- data.frame(
        scenario = kind,
        rep      = rep_id,
        n        = nrow(X),
        p        = ncol(X),
        step     = step,
        B        = B,
        observed_gap         = observed_gap,
        exact_max_abs_error   = max(abs(original$values - collapsed$values)),
        exact_mean_abs_error  = mean(abs(original$values - collapsed$values)),
        shortcut_max_abs_error  = max(abs(original$values - shortcut$values)),
        shortcut_mean_abs_error = mean(abs(original$values - shortcut$values)),
        elapsed_original      = original$elapsed,
        elapsed_collapsed     = collapsed$elapsed,
        elapsed_shortcut      = shortcut$elapsed,
        speedup_collapsed     = original$elapsed / max(collapsed$elapsed, .Machine$double.eps),
        speedup_shortcut      = original$elapsed / max(shortcut$elapsed, .Machine$double.eps),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows[seq_len(idx - 1L)])
}

summarise_runs <- function(results) {
  key <- interaction(results$scenario, results$step, drop = TRUE)
  parts <- split(results, key)
  out <- lapply(parts, function(df) {
    data.frame(
      scenario = df$scenario[[1L]],
      step     = df$step[[1L]],
      n        = df$n[[1L]],
      p        = df$p[[1L]],
      B        = df$B[[1L]],
      observed_gap          = max(df$observed_gap),
      exact_max_abs_error   = max(df$exact_max_abs_error),
      exact_mean_abs_error  = mean(df$exact_mean_abs_error),
      shortcut_max_abs_error  = mean(df$shortcut_max_abs_error),
      shortcut_mean_abs_error = mean(df$shortcut_mean_abs_error),
      median_speedup_collapsed = stats::median(df$speedup_collapsed),
      median_speedup_shortcut  = stats::median(df$speedup_shortcut),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

main <- function() {
  scenarios <- c("null_small", "shadowing_medium", "shadowing_hetero", "wide_medium")
  results <- do.call(
    rbind,
    lapply(scenarios, function(kind) run_scenario(kind = kind, B = 48L, n_rep = 3L, max_step = 4L))
  )
  summary <- summarise_runs(results)

  cat("Vitale P3 parity benchmark\n")
  cat("==========================\n\n")
  cat("Columns:\n")
  cat("- observed_gap: original observed tail-ratio vs. deflated-residual formulation\n")
  cat("- exact_*: original projected-null P3 vs. collapsed exact tail-ratio on the same permuted residuals\n")
  cat("- shortcut_*: original projected-null P3 vs. current engine null shortcut\n")
  cat("- speedup_*: original projected-null runtime divided by comparison runtime\n\n")

  print(summary, row.names = FALSE, digits = 6)

  invisible(summary)
}

if (identical(environment(), globalenv())) {
  main()
}
