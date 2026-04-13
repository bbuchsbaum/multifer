#' Top-k singular values
#'
#' Returns the top `k` singular values of `X`. Uses `RSpectra::svds()`
#' when the package is installed and `k` is small enough to make the
#' Lanczos iteration faster than a full `base::svd()`; otherwise falls
#' back to `base::svd()$d[seq_len(k)]`.
#'
#' The ladder's Monte Carlo draws only need the top few singular values
#' of a permuted residual -- the collapsed Vitale P3 null statistic is
#' `s[step]^2 / sum(s[step:end]^2)`, which on a permuted residual
#' equals `s[step]^2 / (||M_b||_F^2 - sum(s[1:(step-1)]^2))`. Only the
#' top `step` singular values are needed per draw, so computing a full
#' SVD is wasteful for small `step`. Switching to a partial SVD gives
#' a ~10x per-draw speedup at step=1 on n=500, p=300 inputs.
#'
#' @param X Numeric matrix.
#' @param k Positive integer. Number of singular values to return.
#' @return Numeric vector of length `min(k, min(dim(X)))`, sorted
#'   descending.
#' @keywords internal
#' @noRd
top_singular_values <- function(X, k) {
  k <- as.integer(k)
  mn <- min(dim(X))
  if (k >= mn) {
    return(base::svd(X)$d)
  }
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    # RSpectra::svds requires k < min(dim(X)) strictly; guard anyway.
    k_use <- min(k, mn - 1L)
    if (k_use < 1L) return(base::svd(X)$d)
    res <- tryCatch(
      RSpectra::svds(X, k = k_use, nu = 0L, nv = 0L),
      error = function(e) NULL
    )
    if (!is.null(res) && length(res$d) == k_use) {
      return(res$d)
    }
  }
  base::svd(X)$d[seq_len(k)]
}
