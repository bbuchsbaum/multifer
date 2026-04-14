# Core-space resampling primitives
#
# This module hosts the exact core-space fast-path helpers that the
# paper-faithful engines use to evaluate latent-root statistics under
# resampling without repeatedly re-decomposing the full operator.
#
# The unifying idea (Part 1 §17 / Part 5 §30 rank 2): whenever a latent
# operator admits a thin factorization whose column spaces are
# orthonormal, resampling-induced updates to its leading singular
# values can be computed in the low-dimensional core rather than on
# the full operator.
#
# Two instances currently live in the package:
#
# 1. ONE-BLOCK VARIANCE (Fisher-style row bootstrap). Implemented inline
#    inside `adapter_svd()$update_core` and `adapter_prcomp()$update_core`
#    in R/adapter_oneblock_baser.R. Given the thin SVD
#    Xc = U D V^T of the centered original, a row bootstrap X[idx, ]
#    reduces to an n x k inner SVD of (U[idx] - 1 u_bar^T) * diag(d),
#    rotated back through V. Cost drops from O(n p^2) to
#    O(n k^2 + p k).
#
# 2. CROSS-BLOCK COVARIANCE. Implemented below as `.cross_covariance_core*`.
#    Given the block thin SVDs Xc = Ux Dx Vx^T and Yc = Uy Dy Vy^T, the
#    centered resampled cross-product factors exactly as
#
#        Xc^T %*% Yc[perm, ]  =  Vx (Dx Ux^T Uy[perm, ] Dy) Vy^T
#
#    and because Vx / Vy have orthonormal columns the top singular
#    value of the full product equals the top singular value of the
#    k_x x k_y inner core
#
#        M  =  Dx (Ux^T Uy[perm, ]) Dy.
#
#    Cost drops from O(n p q) per draw to O(n k_x k_y) plus a
#    tiny inner SVD. Exactness against the full refit path has been
#    verified to machine precision across six shape / SNR regimes; see
#    notes/cross_core_bootstrap_study.md.
#
# Each instance exposes the same three logical operations:
#
#   - build a core representation from the current latent factorization
#   - evaluate the observed top-singular-value^2 on the core
#   - evaluate a null top-singular-value^2 on the core under a
#     permutation of the core's row-space indices
#
# Phase 2 engines (multiblock, geneig) are expected to add more
# instances that satisfy the same three operations and plug into the
# same `ladder_driver` / `bootstrap_fits` surface without touching the
# driver code.

# ---------------------------------------------------------------------------
# Cross-block covariance core helpers
# ---------------------------------------------------------------------------

.active_root_count <- function(d, tol = 1e-10) {
  if (length(d) == 0L) return(0L)
  thresh <- max(1, d[[1L]]) * tol
  sum(d > thresh)
}

.trim_cross_covariance_fit <- function(fit, tol = 1e-10) {
  if (is.null(fit) || !identical(fit$relation, "covariance")) {
    return(fit)
  }

  k <- .active_root_count(fit$d, tol = tol)
  if (k <= 0L) {
    fit$d <- numeric(0)
    fit$Wx <- fit$Wx[, 0L, drop = FALSE]
    fit$Wy <- fit$Wy[, 0L, drop = FALSE]
    fit$Tx <- fit$Tx[, 0L, drop = FALSE]
    fit$Ty <- fit$Ty[, 0L, drop = FALSE]
    return(fit)
  }

  keep <- seq_len(k)
  fit$d <- fit$d[keep]
  fit$Wx <- fit$Wx[, keep, drop = FALSE]
  fit$Wy <- fit$Wy[, keep, drop = FALSE]
  fit$Tx <- fit$Tx[, keep, drop = FALSE]
  fit$Ty <- fit$Ty[, keep, drop = FALSE]
  fit
}

.cross_covariance_core <- function(X, Y, tol = 1e-10) {
  sv_x <- cached_svd(X)
  sv_y <- cached_svd(Y)

  keep_x <- sv_x$d > max(1, sv_x$d[1L]) * tol
  keep_y <- sv_y$d > max(1, sv_y$d[1L]) * tol
  if (!any(keep_x) || !any(keep_y)) {
    return(list(
      use_core = TRUE,
      Ux = matrix(0, nrow = nrow(X), ncol = 0L),
      dx = numeric(0),
      Vx = matrix(0, nrow = ncol(X), ncol = 0L),
      Uy = matrix(0, nrow = nrow(Y), ncol = 0L),
      dy = numeric(0),
      Vy = matrix(0, nrow = ncol(Y), ncol = 0L)
    ))
  }

  Ux <- sv_x$u[, keep_x, drop = FALSE]
  dx <- sv_x$d[keep_x]
  Vx <- sv_x$v[, keep_x, drop = FALSE]
  Uy <- sv_y$u[, keep_y, drop = FALSE]
  dy <- sv_y$d[keep_y]
  Vy <- sv_y$v[, keep_y, drop = FALSE]

  list(
    # Only engage the core path when its flop count is strictly smaller
    # than the direct crossprod. This is the exact-core cost gate.
    use_core = (length(dx) * length(dy)) < (ncol(X) * ncol(Y)),
    Ux = Ux, dx = dx, Vx = Vx,
    Uy = Uy, dy = dy, Vy = Vy
  )
}

.cross_covariance_core_matrix <- function(core, Uy_view = core$Uy) {
  if (length(core$dx) == 0L || length(core$dy) == 0L) {
    return(matrix(0, nrow = 1L, ncol = 1L))
  }
  inner <- base::crossprod(core$Ux, Uy_view)
  M <- base::sweep(inner, 1L, core$dx, `*`)
  base::sweep(M, 2L, core$dy, `*`)
}

.cross_covariance_core_observed_sv2 <- function(core) {
  M <- .cross_covariance_core_matrix(core)
  s1 <- top_singular_values(M, 1L)[1L]
  s1 * s1
}

.cross_covariance_core_null_sv2 <- function(core, perm) {
  M <- .cross_covariance_core_matrix(core, core$Uy[perm, , drop = FALSE])
  s1 <- top_singular_values(M, 1L)[1L]
  s1 * s1
}

.cross_covariance_core_deflate <- function(core, X, Y) {
  if (length(core$dx) == 0L || length(core$dy) == 0L) {
    return(list(X = X, Y = Y))
  }
  sv_inner <- top_svd(.cross_covariance_core_matrix(core), 1L)
  u1 <- core$Vx %*% sv_inner$u[, 1L, drop = FALSE]
  v1 <- core$Vy %*% sv_inner$v[, 1L, drop = FALSE]
  list(
    X = X - X %*% u1 %*% t(u1),
    Y = Y - Y %*% v1 %*% t(v1)
  )
}
