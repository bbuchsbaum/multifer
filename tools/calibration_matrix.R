#!/usr/bin/env Rscript
# Offline calibration-matrix generator for multifer v1 (multifer-hor.1).
#
# Writes tools/results/calibration_matrix/summary.csv with one row per
# (engine, alpha, n, B, n_sim) regime, reporting the empirical
# rejection rate under a true null. The package's calibration test
# (tests/testthat/test-calibration.R) reads the CSV and asserts each
# row falls inside a k-SE binomial band; it does NOT regenerate the
# matrix so CI time stays bounded.
#
# Usage:
#   Rscript tools/calibration_matrix.R
#   Rscript tools/calibration_matrix.R --n-sim=200 --outdir=tools/results/calibration_matrix
#
# The defaults below are sized for an author-laptop run that is
# informative enough for the binomial band to mean something
# without being punishing to rerun.
#
# Covered families (one row per regime):
#   - oneblock variance          (svd_oneblock)
#   - cross covariance           (cross_svd covariance)
#   - cross correlation paired   (cross_svd correlation, paired_rows)

suppressMessages({
  # When run from the package root, prefer devtools::load_all so the
  # script tracks the working tree without an install step.  Fall
  # back to library(multifer) for CRAN-style environments.
  if (file.exists("DESCRIPTION") &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(multifer)
  }
})

# ---- argument parsing --------------------------------------------------------

parse_args <- function(args) {
  defaults <- list(
    n_sim   = 80L,
    outdir  = "tools/results/calibration_matrix",
    alphas  = "0.01,0.05",
    ns      = "60,150,400",
    Bs      = "199,499",
    verbose = TRUE
  )
  for (a in args) {
    if (!grepl("^--", a)) next
    kv <- sub("^--", "", a)
    key <- sub("=.*$", "", kv)
    val <- sub("^[^=]*=?", "", kv)
    if (!nzchar(val)) val <- "TRUE"
    defaults[[gsub("-", "_", key, fixed = TRUE)]] <- val
  }
  defaults$n_sim  <- as.integer(defaults$n_sim)
  defaults$alphas <- as.numeric(strsplit(defaults$alphas, ",")[[1L]])
  defaults$ns     <- as.integer(strsplit(defaults$ns, ",")[[1L]])
  defaults$Bs     <- as.integer(strsplit(defaults$Bs, ",")[[1L]])
  defaults
}

# ---- per-regime harness ------------------------------------------------------

run_regime <- function(engine, alpha, n, B, n_sim, seed_offset) {
  sim_fn <- switch(
    engine,
    "oneblock_variance" = function(i) {
      bench_oneblock_null(n = n, p = max(8L, n %/% 6L),
                          seed = seed_offset + i)$X
    },
    "cross_covariance" = function(i) {
      d <- bench_cross_null(
        n = n, p_x = max(5L, n %/% 10L), p_y = max(4L, n %/% 12L),
        within_rank_x = 2L, within_rank_y = 2L,
        seed = seed_offset + i
      )
      list(X = d$X, Y = d$Y)
    },
    "cross_correlation_paired" = function(i) {
      d <- bench_cross_null(
        n = n, p_x = max(5L, n %/% 10L), p_y = max(4L, n %/% 12L),
        within_rank_x = 2L, within_rank_y = 2L,
        seed = seed_offset + i
      )
      list(X = d$X, Y = d$Y)
    },
    stop("unknown engine: ", engine, call. = FALSE)
  )

  infer_fn <- switch(
    engine,
    "oneblock_variance" = function(dat, alpha, seed) {
      infer(
        adapter = "svd_oneblock", data = dat,
        geometry = "oneblock", relation = "variance",
        B = B, R = 0L, alpha = alpha, seed = seed
      )
    },
    "cross_covariance" = function(dat, alpha, seed) {
      infer(
        adapter = "cross_svd", data = dat,
        geometry = "cross", relation = "covariance",
        B = B, R = 0L, alpha = alpha, seed = seed
      )
    },
    "cross_correlation_paired" = function(dat, alpha, seed) {
      infer(
        adapter = "cross_svd", data = dat,
        geometry = "cross", relation = "correlation",
        B = B, R = 0L, alpha = alpha, seed = seed
      )
    }
  )

  rejections <- 0L
  for (i in seq_len(n_sim)) {
    dat <- sim_fn(i)
    res <- infer_fn(dat, alpha = alpha, seed = seed_offset + 500L + i)
    if (sum(res$units$selected) >= 1L) rejections <- rejections + 1L
  }
  list(rejections = rejections,
       rate       = rejections / n_sim)
}

# ---- main --------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  outdir <- args$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  engines <- c("oneblock_variance",
               "cross_covariance",
               "cross_correlation_paired")

  rows <- list()
  row_idx <- 0L
  total <- length(engines) * length(args$alphas) *
    length(args$ns) * length(args$Bs)

  t0 <- Sys.time()
  seed_base <- 100000L

  for (engine in engines) {
    for (alpha in args$alphas) {
      for (n in args$ns) {
        for (B in args$Bs) {
          row_idx <- row_idx + 1L
          seed_offset <- seed_base +
            1000L * match(engine, engines) +
            100L * match(alpha, args$alphas) +
            10L * match(n, args$ns) +
            match(B, args$Bs)
          r <- run_regime(engine, alpha, n, B,
                          n_sim = args$n_sim,
                          seed_offset = seed_offset)
          rows[[row_idx]] <- data.frame(
            engine = engine,
            alpha = alpha,
            n = n,
            B = B,
            n_sim = args$n_sim,
            rejections = r$rejections,
            rate = r$rate,
            seed_offset = seed_offset,
            stringsAsFactors = FALSE
          )
          if (isTRUE(as.logical(args$verbose))) {
            cat(sprintf(
              "[%3d/%3d] %s  alpha=%.3f n=%4d B=%4d  rate=%.4f  (%d/%d)\n",
              row_idx, total, engine, alpha, n, B, r$rate,
              r$rejections, args$n_sim
            ))
          }
        }
      }
    }
  }

  summary_df <- do.call(rbind, rows)
  out_path <- file.path(outdir, "summary.csv")
  write.csv(summary_df, out_path, row.names = FALSE)

  cat(sprintf("\nWrote %s (%d rows) in %.1fs\n",
              out_path, nrow(summary_df),
              as.numeric(Sys.time() - t0, units = "secs")))
}

if (sys.nframe() == 0L) {
  main()
}
