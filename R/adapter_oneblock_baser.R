#' Adapter: base::svd for (oneblock, variance)
#'
#' Wraps a raw `base::svd()` of a centered matrix. The fit object is a list
#' with fields `u` (left singular vectors), `d` (singular values), `v` (right
#' singular vectors), and `center` (column means used during centering). Every
#' hook reads from and returns objects of this same shape.
#'
#' @param adapter_id Character scalar. Registry key.
#' @param adapter_version Character scalar version string.
#'
#' @return A `multifer_adapter` object.
#' @export
adapter_svd <- function(adapter_id = "svd_oneblock", adapter_version = "0.0.1") {
  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),
    roots = function(x, ...) x$d^2,

    scores = function(x, domain = NULL, ...) {
      x$u %*% diag(x$d, nrow = length(x$d))
    },

    loadings = function(x, domain = NULL, ...) x$v,

    truncate = function(x, k, ...) {
      k <- min(k, length(x$d))
      list(
        u      = x$u[, seq_len(k), drop = FALSE],
        d      = x$d[seq_len(k)],
        v      = x$v[, seq_len(k), drop = FALSE],
        center = x$center
      )
    },

    residualize = function(x, k, data, ...) {
      k   <- min(k, length(x$d))
      idx <- seq_len(k)
      # rank-k reconstruction: U_k D_k V_k' + column means
      recon <- x$u[, idx, drop = FALSE] %*%
               diag(x$d[idx], nrow = k) %*%
               t(x$v[, idx, drop = FALSE])
      recon <- sweep(recon, 2L, x$center, "+")
      data - recon
    },

    refit = function(x, new_data, ...) {
      # Accept either a raw matrix or a list with $X.
      if (is.list(new_data) && !is.null(new_data$X)) {
        new_data <- new_data$X
      }
      ctr    <- colMeans(new_data)
      Xc     <- sweep(new_data, 2L, ctr, "-")
      fit    <- svd(Xc)
      fit$center <- ctr
      fit
    },

    null_action = function(x, data, ...) {
      # Column-wise permutation of the residual matrix.
      apply(data, 2L, sample)
    },

    component_stat = function(x, data, k, ...) {
      # Relative tail-ratio lambda_k / sum_{q >= k} lambda_q on the SVD of
      # the centered `data`. Mirrors the Part 1 section 1 algebraic
      # simplification of the Vitale P3 statistic.
      Xc <- sweep(data, 2L, colMeans(data), "-")
      s  <- svd(Xc)$d
      s2 <- s^2
      if (k > length(s2) || sum(s2[k:length(s2)]) == 0) return(NA_real_)
      s2[k] / sum(s2[k:length(s2)])
    },

    validity_level       = "conditional",
    declared_assumptions = c("rows_exchangeable", "centered_data"),
    checked_assumptions  = list()
  )
}


#' Adapter: stats::prcomp for (oneblock, variance)
#'
#' Wraps `stats::prcomp()` with centering enabled. The fit object is the list
#' returned by `stats::prcomp()`, which carries `sdev`, `rotation`, `x`, and
#' `center`.
#'
#' @param adapter_id Character scalar. Registry key.
#' @param adapter_version Character scalar version string.
#'
#' @return A `multifer_adapter` object.
#' @export
adapter_prcomp <- function(adapter_id = "prcomp_oneblock",
                           adapter_version = "0.0.1") {
  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),
    roots = function(x, ...) x$sdev^2,

    scores = function(x, domain = NULL, ...) x$x,

    loadings = function(x, domain = NULL, ...) x$rotation,

    truncate = function(x, k, ...) {
      k <- min(k, length(x$sdev))
      out          <- x
      out$sdev     <- x$sdev[seq_len(k)]
      out$rotation <- x$rotation[, seq_len(k), drop = FALSE]
      if (!is.null(x$x)) {
        out$x <- x$x[, seq_len(k), drop = FALSE]
      }
      out
    },

    residualize = function(x, k, data, ...) {
      k   <- min(k, length(x$sdev))
      idx <- seq_len(k)
      scr <- x$x[, idx, drop = FALSE]
      rot <- x$rotation[, idx, drop = FALSE]
      # data - scores[, 1:k] %*% t(rotation[, 1:k])
      # Note: prcomp centers internally; scores are in centered space.
      # Reconstruct in original space by adding back center.
      recon <- scr %*% t(rot)
      ctr   <- if (is.null(x$center)) 0 else x$center
      data - sweep(recon, 2L, ctr, "+")
    },

    refit = function(x, new_data, ...) {
      if (is.list(new_data) && !is.null(new_data$X)) {
        new_data <- new_data$X
      }
      stats::prcomp(new_data, center = TRUE, scale. = FALSE)
    },

    null_action = function(x, data, ...) apply(data, 2L, sample),

    component_stat = function(x, data, k, ...) {
      s2 <- stats::prcomp(data, center = TRUE, scale. = FALSE)$sdev^2
      if (k > length(s2) || sum(s2[k:length(s2)]) == 0) return(NA_real_)
      s2[k] / sum(s2[k:length(s2)])
    },

    validity_level       = "conditional",
    declared_assumptions = c("rows_exchangeable"),
    checked_assumptions  = list()
  )
}


#' Register the base-R oneblock reference adapters
#'
#' Called from the package `.onLoad` hook. Safe to call repeatedly.
#'
#' @keywords internal
register_oneblock_baser_adapters <- function() {
  register_infer_adapter(
    adapter_id = "svd_oneblock",
    adapter    = adapter_svd(),
    overwrite  = TRUE
  )
  register_infer_adapter(
    adapter_id = "prcomp_oneblock",
    adapter    = adapter_prcomp(),
    overwrite  = TRUE
  )
}
