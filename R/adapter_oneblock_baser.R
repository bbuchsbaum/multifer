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

    # Phase 1.5 fast path (§17 / §30 rank 2). Builds a core representation
    # from the original fit and the centered original data matrix. The
    # core stores (U, d, V, center) so that row-bootstrap resamples can
    # be solved by an n x k inner SVD instead of an n x p refit.
    core = function(x, data, ...) {
      if (is.list(data) && !is.null(data$X)) {
        data <- data$X
      }
      list(
        U      = x$u,
        d      = x$d,
        V      = x$v,
        center = x$center,
        n      = nrow(data),
        k      = length(x$d)
      )
    },

    # Core-space update. `indices` is an integer vector of length n_orig
    # (the bootstrap resample indices). Returns a fit-shaped object
    # (list with u, d, v, center) matching `refit()`'s return contract,
    # so downstream `loadings()` and `scores()` hooks are unchanged.
    update_core = function(core_obj, indices = NULL, ...) {
      if (is.null(indices)) {
        stop("`update_core` for adapter_svd requires integer `indices`.",
             call. = FALSE)
      }
      U  <- core_obj$U
      d  <- core_obj$d
      V  <- core_obj$V
      ctr_orig <- core_obj$center

      U_idx  <- U[indices, , drop = FALSE]
      u_bar  <- base::colMeans(U_idx)
      U_tilde <- U_idx - base::rep(u_bar, each = nrow(U_idx))
      # M = U_tilde %*% diag(d) is an n x k inner matrix.
      M  <- base::sweep(U_tilde, 2L, d, `*`)
      sv <- base::svd(M)
      # SVD of centered resampled = (sv$u, sv$d, V %*% sv$v).
      v_new <- V %*% sv$v
      # Compose the new column-center in original variable space:
      # (original center) + V %*% d %*% u_bar
      new_ctr <- ctr_orig + as.numeric(V %*% (d * u_bar))
      list(
        u      = sv$u,
        d      = sv$d,
        v      = v_new,
        center = new_ctr
      )
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

    # Phase 1.5 fast path. prcomp stores (sdev, rotation, x, center); we
    # reconstruct the thin SVD (u, d, v) from those fields so the shared
    # Fisher-style trick applies. Return a prcomp-shaped list so the
    # `loadings()` and `scores()` hooks are unchanged.
    core = function(x, data, ...) {
      if (is.list(data) && !is.null(data$X)) data <- data$X
      n <- nrow(data)
      # prcomp(..., scale.=FALSE): sdev = sqrt(sum(centered^2) / (n-1)) per col.
      # Recover the SVD: d = sdev * sqrt(n-1), u = x / d (where x are scores),
      # v = rotation.
      d <- x$sdev * sqrt(max(1L, n - 1L))
      # Guard against zero singular values when reconstructing u.
      d_inv <- ifelse(d > 0, 1 / d, 0)
      u <- x$x %*% diag(d_inv, nrow = length(d_inv))
      list(
        U      = u,
        d      = d,
        V      = x$rotation,
        center = if (is.null(x$center)) rep(0, ncol(data)) else x$center,
        n      = n,
        k      = length(d)
      )
    },

    update_core = function(core_obj, indices = NULL, ...) {
      if (is.null(indices)) {
        stop("`update_core` for adapter_prcomp requires integer `indices`.",
             call. = FALSE)
      }
      U <- core_obj$U; d <- core_obj$d; V <- core_obj$V
      ctr_orig <- core_obj$center

      U_idx   <- U[indices, , drop = FALSE]
      u_bar   <- base::colMeans(U_idx)
      U_tilde <- U_idx - base::rep(u_bar, each = nrow(U_idx))
      M       <- base::sweep(U_tilde, 2L, d, `*`)
      sv      <- base::svd(M)
      v_new   <- V %*% sv$v
      new_ctr <- ctr_orig + as.numeric(V %*% (d * u_bar))
      n_new   <- length(indices)

      # Return a prcomp-shaped object so adapter$loadings/scores work
      # unchanged. sdev uses the same n-1 convention as stats::prcomp.
      structure(
        list(
          sdev     = sv$d / sqrt(max(1L, n_new - 1L)),
          rotation = v_new,
          center   = new_ctr,
          scale    = FALSE,
          x        = base::sweep(sv$u, 2L, sv$d, `*`)
        ),
        class = "prcomp"
      )
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
