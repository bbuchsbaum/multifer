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

run_bench <- identical(arg_value("run-bench", "true"), "true")
outdir <- arg_value("outdir", file.path(root, "tools", "results", "mature_guardrails"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  pkgload::load_all(root, quiet = TRUE)
})

run_speed_pair <- function(size, seed = 1L) {
  dat <- bench_speed_agreement(size = size, seed = seed)

  t_fast <- system.time({
    fast <- infer_plsc(
      dat$X, dat$Y,
      B = 49L, R = 8L, alpha = 0.05, seed = seed,
      fast_path = "auto"
    )
  })[["elapsed"]]

  t_slow <- system.time({
    slow <- infer_plsc(
      dat$X, dat$Y,
      B = 49L, R = 8L, alpha = 0.05, seed = seed,
      fast_path = "off"
    )
  })[["elapsed"]]

  data.frame(
    suite = paste0("speed_", size),
    fast_time = t_fast,
    slow_time = t_slow,
    speedup = t_slow / t_fast,
    selected_agree = identical(fast$units$selected, slow$units$selected),
    p_value_max_abs = max(abs(fast$component_tests$p_value - slow$component_tests$p_value), na.rm = TRUE),
    statistic_max_abs = max(abs(fast$component_tests$statistic - slow$component_tests$statistic), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

cat("Running mature-family guardrail tests...\n")
testthat::test_local(
  root,
  filter = "exact-path-regression|calibration|core-update-cross|vitale-p3-parity|infer-wrappers|auto-subspace",
  reporter = "summary"
)

if (run_bench) {
  bench_df <- do.call(rbind, list(
    run_speed_pair("medium", seed = 41L),
    run_speed_pair("large", seed = 43L)
  ))
  utils::write.csv(bench_df, file.path(outdir, "mature_guardrails_speed.csv"), row.names = FALSE)
  cat("\nMature-family speed guardrails:\n")
  print(bench_df, row.names = FALSE)
}
