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

n_perm <- as.integer(arg_value("n-perm", "120"))
n_boot <- as.integer(arg_value("n-boot", "120"))
outdir <- arg_value("outdir", file.path(root, "tools", "results", "procrustes_bias"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  pkgload::load_all(root, quiet = TRUE)
})

clamp01 <- function(x) pmin(pmax(x, 0), 1)

orthonormal_matrix <- function(n, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  qr.Q(qr(matrix(stats::rnorm(n * k), nrow = n, ncol = k)))
}

simulate_cross_latent <- function(n, p_x, p_y, signal = c(4, 3), noise_x = 1, noise_y = 1,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- length(signal)
  T_scores <- scale(matrix(stats::rnorm(n * k), nrow = n, ncol = k), center = TRUE, scale = FALSE)
  T_scores <- qr.Q(qr(T_scores))
  Wx_true <- orthonormal_matrix(p_x, k)
  Wy_true <- orthonormal_matrix(p_y, k)

  X_signal <- T_scores %*% diag(signal, nrow = k) %*% t(Wx_true)
  Y_signal <- T_scores %*% diag(signal, nrow = k) %*% t(Wy_true)

  X <- X_signal + matrix(stats::rnorm(n * p_x, sd = noise_x), nrow = n, ncol = p_x)
  Y <- Y_signal + matrix(stats::rnorm(n * p_y, sd = noise_y), nrow = n, ncol = p_y)

  list(
    X = X,
    Y = Y,
    Wx_true = Wx_true,
    Wy_true = Wy_true,
    T = T_scores
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

make_cross_recipe <- function() {
  adapter <- adapter_cross_svd()
  infer_recipe(geometry = "cross", relation = "covariance", adapter = adapter)
}

procrustes_q <- function(Vref, Vb) {
  perm <- match_components(Vref, Vb)
  Vb_perm <- Vb[, perm, drop = FALSE]
  M <- crossprod(Vref, Vb_perm)
  sv <- svd(M)
  list(
    perm = perm,
    Q = sv$v %*% t(sv$u),
    Vb_perm = Vb_perm
  )
}

mcintosh_rotated_d <- function(d, perm, Q) {
  d_perm <- d[perm]
  sqrt(colSums((diag(d_perm, nrow = length(d_perm), ncol = length(d_perm)) %*% Q)^2))
}

axis_metrics <- function(V_aligned, V_true, V_ref) {
  p_truth <- match_components(V_true, V_aligned)
  V_truth_perm <- V_aligned[, p_truth, drop = FALSE]
  V_truth_sign <- align_sign(V_true, V_truth_perm)
  truth_overlap <- clamp01(abs(colSums(V_true * V_truth_sign)))
  truth_angles <- acos(truth_overlap)

  p_ref <- match_components(V_ref, V_aligned)
  V_ref_perm <- V_aligned[, p_ref, drop = FALSE]
  V_ref_sign <- align_sign(V_ref, V_ref_perm)
  ref_overlap <- clamp01(abs(colSums(V_ref * V_ref_sign)))
  ref_angles <- acos(ref_overlap)

  G <- abs(crossprod(V_true, V_aligned))^2
  offdiag_share <- if (sum(G) > 0) {
    truth_diag <- abs(crossprod(V_true, V_truth_sign))^2
    (sum(truth_diag) - sum(diag(truth_diag))) / sum(truth_diag)
  } else {
    NA_real_
  }

  list(
    truth_angles = truth_angles,
    ref_angles = ref_angles,
    subspace_angle_max = max(principal_angles(V_true, V_aligned)),
    offdiag_share = offdiag_share
  )
}

null_distortion_summary <- function(dat, scenario_name, align_domain = c("X", "Y"),
                                    n_perm = 120L, seed = 1L) {
  align_domain <- match.arg(align_domain)
  recipe <- make_cross_recipe()
  adapter <- recipe$adapter
  fit0 <- adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))
  Vref <- adapter$loadings(fit0, align_domain)[, 1:2, drop = FALSE]

  shifts <- matrix(NA_real_, nrow = n_perm, ncol = 3L)
  colnames(shifts) <- c("lv1_shift", "lv2_shift", "tail_ratio_shift")
  mixing <- numeric(n_perm)

  set.seed(seed)
  for (b in seq_len(n_perm)) {
    perm <- sample.int(nrow(dat$Y))
    fitb <- adapter$refit(fit0, list(X = dat$X, Y = dat$Y[perm, , drop = FALSE], relation = "covariance"))
    Vb <- adapter$loadings(fitb, align_domain)[, 1:2, drop = FALSE]
    d_raw <- fitb$d[1:2]

    pq <- procrustes_q(Vref, Vb)
    d_rot <- mcintosh_rotated_d(d_raw, pq$perm, pq$Q)

    shifts[b, 1L] <- d_rot[1L] - d_raw[pq$perm][1L]
    shifts[b, 2L] <- d_rot[2L] - d_raw[pq$perm][2L]
    shifts[b, 3L] <- (d_rot[1L]^2 / sum(d_rot^2)) -
      (d_raw[pq$perm][1L]^2 / sum(d_raw[pq$perm]^2))
    mixing[b] <- 1 - mean(diag(pq$Q^2))
  }

  data.frame(
    scenario = scenario_name,
    align_domain = align_domain,
    n_perm = n_perm,
    lv1_mean_shift = mean(shifts[, 1L]),
    lv2_mean_shift = mean(shifts[, 2L]),
    tail_ratio_mean_shift = mean(shifts[, 3L]),
    lv1_abs_shift_mean = mean(abs(shifts[, 1L])),
    lv2_abs_shift_mean = mean(abs(shifts[, 2L])),
    mixing_mean = mean(mixing),
    stringsAsFactors = FALSE
  )
}

bootstrap_alignment_summary <- function(dat, scenario_name, method = c("sign", "procrustes"),
                                        n_boot = 120L, seed = 1L) {
  method <- match.arg(method)
  recipe <- make_cross_recipe()
  adapter <- recipe$adapter
  fit0 <- adapter$refit(NULL, list(X = dat$X, Y = dat$Y, relation = "covariance"))
  units <- form_units(adapter$roots(fit0))

  artifact <- bootstrap_fits(
    recipe = recipe,
    adapter = adapter,
    data = list(X = dat$X, Y = dat$Y),
    original_fit = fit0,
    units = units,
    R = n_boot,
    method_align = method,
    seed = seed,
    store_aligned_scores = FALSE
  )

  domains <- c("X", "Y")
  truths <- list(X = dat$Wx_true[, 1:2, drop = FALSE], Y = dat$Wy_true[, 1:2, drop = FALSE])
  refs <- list(X = adapter$loadings(fit0, "X")[, 1:2, drop = FALSE],
               Y = adapter$loadings(fit0, "Y")[, 1:2, drop = FALSE])

  rows <- vector("list", length(domains))
  for (i in seq_along(domains)) {
    d <- domains[[i]]
    per_rep <- vector("list", n_boot)
    arr <- array(NA_real_, dim = c(nrow(truths[[d]]), 2L, n_boot))
    for (b in seq_len(n_boot)) {
      Vb <- artifact$reps[[b]]$aligned_loadings[[d]][, 1:2, drop = FALSE]
      per_rep[[b]] <- axis_metrics(Vb, truths[[d]], refs[[d]])
      arr[, , b] <- Vb
    }

    truth_angle_mat <- do.call(rbind, lapply(per_rep, function(x) x$truth_angles))
    ref_angle_mat <- do.call(rbind, lapply(per_rep, function(x) x$ref_angles))
    offdiag <- vapply(per_rep, `[[`, numeric(1L), "offdiag_share")
    subspace <- vapply(per_rep, `[[`, numeric(1L), "subspace_angle_max")
    loading_sd_mean <- mean(apply(arr, c(1L, 2L), stats::sd))

    rows[[i]] <- data.frame(
      scenario = scenario_name,
      method = method,
      domain = d,
      n_boot = n_boot,
      truth_angle1_mean = mean(truth_angle_mat[, 1L]),
      truth_angle2_mean = mean(truth_angle_mat[, 2L]),
      reference_angle1_mean = mean(ref_angle_mat[, 1L]),
      reference_angle2_mean = mean(ref_angle_mat[, 2L]),
      offdiag_share_mean = mean(offdiag),
      subspace_angle_max_mean = mean(subspace),
      loading_sd_mean = loading_sd_mean,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

run_study <- function(n_perm = 120L, n_boot = 120L, outdir) {
  scenarios <- scenario_grid()

  null_rows <- vector("list", length(scenarios) * 2L)
  boot_rows <- vector("list", length(scenarios) * 2L)
  pos_null <- 1L
  pos_boot <- 1L

  for (i in seq_along(scenarios)) {
    sc <- scenarios[[i]]
    dat <- simulate_cross_latent(
      n = sc$n,
      p_x = sc$p_x,
      p_y = sc$p_y,
      signal = sc$signal,
      noise_x = sc$noise_x,
      noise_y = sc$noise_y,
      seed = 1000L + i
    )

    for (align_domain in c("X", "Y")) {
      null_rows[[pos_null]] <- null_distortion_summary(
        dat,
        scenario_name = sc$name,
        align_domain = align_domain,
        n_perm = n_perm,
        seed = 2000L + i * 10L + if (align_domain == "X") 1L else 2L
      )
      pos_null <- pos_null + 1L
    }

    for (method in c("sign", "procrustes")) {
      boot_rows[[pos_boot]] <- bootstrap_alignment_summary(
        dat,
        scenario_name = sc$name,
        method = method,
        n_boot = n_boot,
        seed = 3000L + i * 10L + if (method == "sign") 1L else 2L
      )
      pos_boot <- pos_boot + 1L
    }
  }

  null_df <- do.call(rbind, null_rows[seq_len(pos_null - 1L)])
  boot_df <- do.call(rbind, boot_rows[seq_len(pos_boot - 1L)])
  sign_df <- boot_df[boot_df$method == "sign", , drop = FALSE]
  proc_df <- boot_df[boot_df$method == "procrustes", , drop = FALSE]
  rownames(sign_df) <- paste(sign_df$scenario, sign_df$domain, sep = "::")
  rownames(proc_df) <- paste(proc_df$scenario, proc_df$domain, sep = "::")
  shared <- intersect(rownames(sign_df), rownames(proc_df))
  compare_df <- data.frame(
    scenario = sign_df[shared, "scenario"],
    domain = sign_df[shared, "domain"],
    delta_truth_angle1 = proc_df[shared, "truth_angle1_mean"] - sign_df[shared, "truth_angle1_mean"],
    delta_truth_angle2 = proc_df[shared, "truth_angle2_mean"] - sign_df[shared, "truth_angle2_mean"],
    delta_reference_angle1 = proc_df[shared, "reference_angle1_mean"] - sign_df[shared, "reference_angle1_mean"],
    delta_reference_angle2 = proc_df[shared, "reference_angle2_mean"] - sign_df[shared, "reference_angle2_mean"],
    delta_loading_sd = proc_df[shared, "loading_sd_mean"] - sign_df[shared, "loading_sd_mean"],
    delta_subspace_angle = proc_df[shared, "subspace_angle_max_mean"] - sign_df[shared, "subspace_angle_max_mean"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  utils::write.csv(null_df, file.path(outdir, "procrustes_null_distortion.csv"), row.names = FALSE)
  utils::write.csv(boot_df, file.path(outdir, "procrustes_bootstrap_alignment.csv"), row.names = FALSE)
  utils::write.csv(compare_df, file.path(outdir, "procrustes_method_deltas.csv"), row.names = FALSE)

  cat("\nNull distortion summary:\n")
  print(null_df, row.names = FALSE)
  cat("\nBootstrap alignment summary:\n")
  print(boot_df, row.names = FALSE)
  cat("\nMethod deltas (procrustes - sign):\n")
  print(compare_df, row.names = FALSE)

  invisible(list(null = null_df, bootstrap = boot_df, deltas = compare_df))
}

run_study(n_perm = n_perm, n_boot = n_boot, outdir = outdir)
