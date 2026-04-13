.interp_quantile_sorted <- function(sorted_x, p) {
  n <- length(sorted_x)
  if (n == 0L) {
    return(NA_real_)
  }
  if (n == 1L) {
    return(sorted_x[1L])
  }

  h <- (n - 1) * p + 1
  lo <- floor(h)
  hi <- ceiling(h)

  if (lo == hi) {
    return(sorted_x[lo])
  }

  sorted_x[lo] + (h - lo) * (sorted_x[hi] - sorted_x[lo])
}

.row_quantile_pair <- function(mat, probs) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }

  n <- nrow(mat)
  if (n == 0L) {
    return(list(lower = numeric(0L), upper = numeric(0L)))
  }

  # Fast vectorized path: matrixStats::rowQuantiles runs a single C loop
  # over rows and matches base::quantile(type = 7) interpolation, which is
  # identical to .interp_quantile_sorted below. Only usable when every row
  # is finite; a single row with NA/NaN/Inf falls back to the per-row loop.
  if (!anyNA(mat) && all(is.finite(mat))) {
    q <- matrixStats::rowQuantiles(mat, probs = as.numeric(probs),
                                   na.rm = FALSE)
    # rowQuantiles returns a matrix with one column per prob; guard against
    # the single-column degeneration that matrixStats sometimes drops.
    if (is.matrix(q)) {
      return(list(lower = q[, 1L], upper = q[, 2L]))
    }
    # single-row case: q is a length-2 vector.
    return(list(lower = q[1L], upper = q[2L]))
  }

  # Slow path: per-row loop with NA/NaN/Inf tolerance.
  lower <- rep(NA_real_, n)
  upper <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    v <- mat[i, ]
    v <- v[is.finite(v)]
    if (length(v) == 0L) {
      next
    }
    v <- sort(v)
    lower[i] <- .interp_quantile_sorted(v, probs[1L])
    upper[i] <- .interp_quantile_sorted(v, probs[2L])
  }
  list(lower = lower, upper = upper)
}

.row_mean_sd <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }

  if (nrow(mat) == 0L) {
    return(list(mean = numeric(0L), sd = numeric(0L)))
  }

  if (!anyNA(mat) && all(is.finite(mat))) {
    mean <- rowMeans(mat)
    if (ncol(mat) < 2L) {
      sd <- rep(0, nrow(mat))
    } else {
      centered <- mat - mean
      sd <- sqrt(pmax(0, rowSums(centered * centered) / (ncol(mat) - 1L)))
    }
    return(list(mean = mean, sd = sd))
  }

  counts <- rowSums(is.finite(mat))
  mean <- rowMeans(mat, na.rm = TRUE)
  mean[counts == 0L] <- NA_real_

  sd <- rep(0, nrow(mat))
  has_var <- counts >= 2L
  if (any(has_var)) {
    centered <- mat
    centered[!is.finite(centered)] <- NA_real_
    centered <- centered - mean
    sd[has_var] <- sqrt(rowSums(centered[has_var, , drop = FALSE]^2, na.rm = TRUE) /
                          (counts[has_var] - 1L))
  }
  sd[counts == 0L] <- NA_real_

  list(mean = mean, sd = sd)
}
