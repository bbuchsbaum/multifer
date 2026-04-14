#' Adapter: thin-SVD cross-block adapter for (cross, covariance) and (cross, correlation)
#'
#' `adapter_cross_svd()` wraps a thin SVD of the cross-product matrix and
#' declares support for BOTH `(cross, covariance)` AND `(cross, correlation)`.
#' This dual declaration is the canonical example of the strict-dispatch
#' ambiguity story (Part 5 section 35): calling `infer_recipe()` with
#' `geometry = "cross"` and no `relation` on this adapter will error under
#' `strict = TRUE`.
#'
#' @section Fit object convention:
#' The fit object produced and consumed by every hook is a plain list with
#' the following fields:
#'
#' - `relation` -- character scalar, either `"covariance"` or
#'   `"correlation"`. Set by `refit()` from `new_data$relation`.
#' - `d` -- numeric vector of latent roots (singular values of the
#'   relation-appropriate cross matrix, length r).
#' - `Wx` -- numeric matrix (p x r). X weight vectors (right singular
#'   vectors of the cross matrix, or right singular vectors of the
#'   orthonormal cross matrix for correlation mode).
#' - `Wy` -- numeric matrix (q x r). Y weight vectors.
#' - `Tx` -- numeric matrix (n x r). X scores: `X %*% Wx`.
#' - `Ty` -- numeric matrix (n x r). Y scores: `Y %*% Wy`.
#' - `center_x` -- numeric vector (p). Column means of X used for
#'   centering.
#' - `center_y` -- numeric vector (q). Column means of Y used for
#'   centering.
#'
#' The `relation` field governs which math path is taken: covariance mode
#' is the PLSC path (`svd(t(X) %*% Y)`) and correlation mode is the CCA
#' path (`svd(t(Qx) %*% Qy)` where Qx, Qy are thin QR orthonormal
#' factors). Under row permutation P, the permuted covariance root is
#' `s_a(t(X) %*% P %*% Y)` and the permuted correlation root is
#' `s_a(t(Qx) %*% P %*% Qy)`. Phase 1.5 can exploit the small-problem
#' structure this induces; in Phase 1 (refit-first) we use the full
#' recomputation path.
#'
#' @param adapter_id Character scalar. Registry key.
#' @param adapter_version Character scalar version string.
#'
#' @return A `multifer_adapter` object.
#' @export
adapter_cross_svd <- function(adapter_id = "cross_svd",
                              adapter_version = "0.0.1") {
  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(geometry = "cross", relation = "covariance",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability")),
      list(geometry = "cross", relation = "correlation",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),

    roots = function(x, ...) x$d,

    scores = function(x, domain = c("X", "Y"), ...) {
      domain <- match.arg(domain)
      if (domain == "X") x$Tx else x$Ty
    },

    loadings = function(x, domain = c("X", "Y"), ...) {
      domain <- match.arg(domain)
      if (domain == "X") x$Wx else x$Wy
    },

    truncate = function(x, k, ...) {
      k <- min(k, length(x$d))
      idx <- seq_len(k)
      list(
        relation = x$relation,
        d        = x$d[idx],
        Wx       = x$Wx[, idx, drop = FALSE],
        Wy       = x$Wy[, idx, drop = FALSE],
        Tx       = x$Tx[, idx, drop = FALSE],
        Ty       = x$Ty[, idx, drop = FALSE],
        center_x = x$center_x,
        center_y = x$center_y
      )
    },

    residualize = function(x, k, data, ...) {
      # data is list(X = ., Y = .).
      # This reference adapter keeps a simple adapter-local residualize()
      # hook so the capability gate is satisfied for both declared
      # relations. The correlation engine used by infer_cca()/infer()
      # does NOT rely on this approximation for supported CCA designs:
      # run_cross_ladder() performs the QR-whitened stepwise deflation
      # directly for relation = "correlation". This hook therefore
      # remains a generic reference implementation and fallback, not the
      # canonical CCA ladder path.
      k   <- min(k, length(x$d))
      idx <- seq_len(k)
      Wx_k <- x$Wx[, idx, drop = FALSE]
      Wy_k <- x$Wy[, idx, drop = FALSE]
      Xd <- data$X - data$X %*% Wx_k %*% t(Wx_k)
      Yd <- data$Y - data$Y %*% Wy_k %*% t(Wy_k)
      list(X = Xd, Y = Yd)
    },

    refit = function(x, new_data, ...) {
      # new_data is list(X = ., Y = ., relation = "covariance"|"correlation").
      # If new_data$relation is absent, fall back to x$relation if available.
      rel <- new_data$relation
      if (is.null(rel) && !is.null(x) && !is.null(x$relation)) {
        rel <- x$relation
      }
      if (is.null(rel)) {
        stop("adapter_cross_svd refit: `new_data$relation` must be set to \"covariance\" or \"correlation\".",
             call. = FALSE)
      }
      if (!(rel %in% c("covariance", "correlation"))) {
        stop(sprintf(
          "adapter_cross_svd refit: unknown relation '%s'. Must be \"covariance\" or \"correlation\".",
          rel
        ), call. = FALSE)
      }

      X <- new_data$X
      Y <- new_data$Y

      # Center both blocks.
      cx <- colMeans(X)
      cy <- colMeans(Y)
      Xc <- sweep(X, 2L, cx, "-")
      Yc <- sweep(Y, 2L, cy, "-")

      if (rel == "covariance") {
        # PLSC path: SVD of the cross-product matrix t(X) %*% Y.
        # crossprod(Xc, Yc) is (p x q).
        # sv$u is (p x r): X weight vectors.
        # sv$v is (q x r): Y weight vectors.
        # d[a] = s_a(Xc^t Yc).
        sv <- svd(crossprod(Xc, Yc))
        Wx <- sv$u   # p x r
        Wy <- sv$v   # q x r
        Tx <- Xc %*% Wx
        Ty <- Yc %*% Wy
      } else {
        # CCA path: SVD of orthonormal cross-product t(Qx) %*% Qy.
        # From the thin QR: Xc = Qx %*% Rx, so Qx is (n x p), Rx is (p x p).
        # crossprod(Qx, Qy) is (p x q).
        # sv$u is (p x r) in Qx column space (X canonical directions in QR
        # basis); sv$v is (q x r) in Qy column space.
        # Scores in original space: Tx = Qx %*% sv$u (n x r).
        # Weights in original variable space: Wx = Rx^{-1} %*% sv$u (p x r)
        # so that Xc %*% Wx = Qx %*% sv$u = Tx.
        # d[a] = s_a(Qx^t Qy) = a-th canonical correlation.
        qx  <- qr(Xc)
        qy  <- qr(Yc)
        Qx  <- qr.Q(qx)
        Qy  <- qr.Q(qy)
        sv  <- svd(crossprod(Qx, Qy))
        Tx  <- Qx %*% sv$u   # n x r -- scores in original space
        Ty  <- Qy %*% sv$v   # n x r
        # Weights: solve Rx %*% Wx = sv$u  =>  Wx = backsolve(Rx, sv$u)
        Rx  <- qr.R(qx)
        Ry  <- qr.R(qy)
        Wx  <- backsolve(Rx, sv$u)   # p x r
        Wy  <- backsolve(Ry, sv$v)   # q x r
      }

      fit <- list(
        relation = rel,
        d        = sv$d,
        Wx       = Wx,
        Wy       = Wy,
        Tx       = Tx,
        Ty       = Ty,
        center_x = cx,
        center_y = cy
      )
      .trim_cross_covariance_fit(fit)
    },

    # Phase 1.5 fast path (§17 / §30 rank 2).
    #
    # Covariance mode: the existing D-weighted identity lands the
    # collapsed Vitale P3 null in a k_x x k_y inner SVD.
    #
    # Correlation mode (multifer-9u9.1.3): the per-replicate column
    # rescaling DOES admit a clean k_x x k_y update via the whitening
    # identity
    #
    #     Qx_boot^T Qy_boot = B_x^{-1/2} (Ux_tilde^T Uy_tilde) B_y^{-1/2}
    #
    # where B_x = Ux_tilde^T Ux_tilde and B_y = Uy_tilde^T Uy_tilde are
    # the k_x x k_x and k_y x k_y Gram matrices of the centered-and-
    # resampled block singular vectors. Canonical correlations are the
    # singular values of the whitened small matrix. Scores and weights
    # lift back to the original variable space via Vx / Vy. This
    # identity is exact to machine precision and is the path taken when
    # core_obj$relation == "correlation".
    #
    # Exactness of the correlation fast path requires the core to carry
    # the FULL column rank of the centered blocks -- any rank truncation
    # changes the whitening and produces a biased canonical correlation.
    # .truncate_core() in R/bootstrap.R therefore leaves correlation-mode
    # cores alone; callers that need the truncation speed-up must stay on
    # the covariance path.
    core = function(x, data, ...) {
      X <- data$X; Y <- data$Y
      cx <- base::colMeans(X); cy <- base::colMeans(Y)
      Xc <- base::sweep(X, 2L, cx, "-")
      Yc <- base::sweep(Y, 2L, cy, "-")
      sv_x <- base::svd(Xc)
      sv_y <- base::svd(Yc)
      rel <- if (!is.null(x) && !is.null(x$relation)) x$relation else "covariance"
      list(
        relation = rel,
        Ux       = sv_x$u, dx = sv_x$d, Vx = sv_x$v,
        Uy       = sv_y$u, dy = sv_y$d, Vy = sv_y$v,
        center_x = cx,
        center_y = cy,
        n        = nrow(X)
      )
    },

    update_core = function(core_obj, indices = NULL, ...) {
      if (is.null(indices)) {
        stop("`update_core` for adapter_cross_svd requires integer `indices`.",
             call. = FALSE)
      }
      Ux <- core_obj$Ux; dx <- core_obj$dx; Vx <- core_obj$Vx
      Uy <- core_obj$Uy; dy <- core_obj$dy; Vy <- core_obj$Vy
      rel <- if (!is.null(core_obj$relation)) core_obj$relation else "covariance"

      Ux_idx <- Ux[indices, , drop = FALSE]
      Uy_idx <- Uy[indices, , drop = FALSE]
      n_new  <- nrow(Ux_idx)

      ux_bar <- base::colMeans(Ux_idx)
      uy_bar <- base::colMeans(Uy_idx)
      Ux_tilde <- Ux_idx - base::rep(ux_bar, each = n_new)
      Uy_tilde <- Uy_idx - base::rep(uy_bar, each = n_new)

      if (rel == "covariance") {
        # DUx[i, j] = Ux_tilde[i, j] * dx[j]; same for DUy.
        DUx <- base::sweep(Ux_tilde, 2L, dx, `*`)
        DUy <- base::sweep(Uy_tilde, 2L, dy, `*`)

        # M = D_x (Ux_tilde^T Uy_tilde) D_y  ==  crossprod(DUx, DUy), k_x x k_y.
        M  <- base::crossprod(DUx, DUy)
        sv <- base::svd(M)

        # The centered resampled cross-product factors as Vx M Vy^T, so its
        # SVD is (Vx %*% sv$u) * diag(sv$d) * (Vy %*% sv$v)^T.
        Wx_new <- Vx %*% sv$u                 # p_x x r
        Wy_new <- Vy %*% sv$v                 # p_y x r
        d_new  <- sv$d

        # Scores: X_boot_c %*% Wx_new = Ux_tilde Dx Vx^T (Vx A) = Ux_tilde Dx A
        #                              = DUx %*% sv$u   (since Vx has orthonormal cols).
        Tx_new <- DUx %*% sv$u
        Ty_new <- DUy %*% sv$v
      } else {
        # correlation path: whitened k_x x k_y update.
        Bx  <- base::crossprod(Ux_tilde)
        By  <- base::crossprod(Uy_tilde)
        Cxy <- base::crossprod(Ux_tilde, Uy_tilde)

        # Symmetric matrix square-root inverse via eigendecomposition.
        # Guard: near-zero or slightly negative eigenvalues arise from
        # two sources:
        #   (a) bootstrap resamples that do not span the centered block
        #       column space (too few unique rows relative to k_x or k_y);
        #   (b) floating-point loss of PSD in the symmetrisation step.
        # Either way the B^{-1/2} whitening identity breaks down. We signal
        # this with a classed condition so the caller can fall back to a
        # full refit for this replicate. The fallback lives in
        # bootstrap.R (see the tryCatch around adapter$update_core) and is
        # tested in tests/testthat/test-correlation-fallback.R.
        # Do NOT catch errors other than multifer_core_rank_deficient here;
        # unrelated failures must propagate.
        ex <- base::eigen((Bx + base::t(Bx)) / 2, symmetric = TRUE)
        ey <- base::eigen((By + base::t(By)) / 2, symmetric = TRUE)
        tol_core <- max(1, ex$values[1], ey$values[1]) * sqrt(.Machine$double.eps)
        if (any(ex$values <= tol_core) || any(ey$values <= tol_core)) {
          cond <- structure(
            class = c("multifer_core_rank_deficient", "error", "condition"),
            list(
              message = "adapter_cross_svd update_core: correlation-mode bootstrap sample is rank-deficient in the core representation; caller should fall back to refit for this replicate.",
              call    = NULL
            )
          )
          stop(cond)
        }
        Bx_inv_sqrt <- ex$vectors %*%
          base::sweep(base::t(ex$vectors), 1L, 1 / sqrt(ex$values), `*`)
        By_inv_sqrt <- ey$vectors %*%
          base::sweep(base::t(ey$vectors), 1L, 1 / sqrt(ey$values), `*`)

        M  <- Bx_inv_sqrt %*% Cxy %*% By_inv_sqrt        # k_x x k_y
        sv <- base::svd(M)

        # Weights in original variable space:
        #   Wx = Vx %*% diag(1/dx) %*% Bx_inv_sqrt %*% sv$u
        Wx_new <- Vx %*% base::sweep(Bx_inv_sqrt %*% sv$u, 1L, dx, `/`)
        Wy_new <- Vy %*% base::sweep(By_inv_sqrt %*% sv$v, 1L, dy, `/`)
        d_new  <- sv$d

        # Scores: Tx = Qx_boot %*% sv$u = Ux_tilde Bx_inv_sqrt sv$u.
        Tx_new <- Ux_tilde %*% (Bx_inv_sqrt %*% sv$u)
        Ty_new <- Uy_tilde %*% (By_inv_sqrt %*% sv$v)
      }

      new_cx <- core_obj$center_x + as.numeric(Vx %*% (dx * ux_bar))
      new_cy <- core_obj$center_y + as.numeric(Vy %*% (dy * uy_bar))

      fit <- list(
        relation = rel,
        d        = d_new,
        Wx       = Wx_new,
        Wy       = Wy_new,
        Tx       = Tx_new,
        Ty       = Ty_new,
        center_x = new_cx,
        center_y = new_cy
      )
      .trim_cross_covariance_fit(fit)
    },

    null_action = function(x, data, ...) {
      # Break pairing by permuting rows of Y.
      perm <- sample(seq_len(nrow(data$Y)))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },

    component_stat = function(x, data, k, ...) {
      # Relative tail-ratio analogue for cross-block shapes:
      # s_k^2 / sum_{q >= k} s_q^2 where s are singular values of the
      # relation-appropriate cross matrix. `data` is the DEFLATED
      # cross-block data at step k, list(X = ., Y = .).
      # Dispatch on x$relation (set by refit() from the recipe's resolved
      # relation).
      if (is.null(x$relation)) {
        stop("adapter_cross_svd: `relation` must be set in the fit object before component_stat is called",
             call. = FALSE)
      }
      s <- switch(
        x$relation,
        covariance  = svd(crossprod(data$X, data$Y))$d,
        correlation = svd(crossprod(qr.Q(qr(data$X)), qr.Q(qr(data$Y))))$d,
        stop("unknown relation: ", x$relation, call. = FALSE)
      )
      if (k > length(s) || sum(s[k:length(s)]^2) == 0) return(NA_real_)
      s[k]^2 / sum(s[k:length(s)]^2)
    },

    validity_level       = "conditional",
    declared_assumptions = c("paired_rows", "centered_blocks"),
    checked_assumptions  = .cross_baser_checks()
  )
}


#' Adapter: stats::cancor for (cross, correlation)
#'
#' Wraps `stats::cancor()`. Declares ONLY `(cross, correlation)`, so it has
#' no ambiguity under strict dispatch. This is the canonical default for
#' [infer_cca()]. Use `adapter_cross_svd()` when you need both covariance
#' and correlation modes from a single reference adapter.
#'
#' @section Fit object convention:
#' The fit object produced and consumed by every hook is a plain list with
#' the following fields:
#'
#' - `cor` -- numeric vector of canonical correlations.
#' - `xcoef` -- numeric matrix (p x r). X canonical coefficients.
#' - `ycoef` -- numeric matrix (q x r). Y canonical coefficients.
#' - `xcoef_scores` -- numeric matrix (n x r). `X %*% xcoef`.
#' - `ycoef_scores` -- numeric matrix (n x r). `Y %*% ycoef`.
#' - `xcenter` -- numeric vector (p). X column means.
#' - `ycenter` -- numeric vector (q). Y column means.
#'
#' @param adapter_id Character scalar. Registry key.
#' @param adapter_version Character scalar version string.
#'
#' @return A `multifer_adapter` object.
#' @export
adapter_cancor <- function(adapter_id = "cancor_cross",
                           adapter_version = "0.0.1") {
  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(geometry = "cross", relation = "correlation",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),

    roots = function(x, ...) x$cor,

    scores = function(x, domain = c("X", "Y"), ...) {
      domain <- match.arg(domain)
      if (domain == "X") x$xcoef_scores else x$ycoef_scores
    },

    loadings = function(x, domain = c("X", "Y"), ...) {
      domain <- match.arg(domain)
      if (domain == "X") x$xcoef else x$ycoef
    },

    truncate = function(x, k, ...) {
      k <- min(k, length(x$cor))
      idx <- seq_len(k)
      list(
        cor          = x$cor[idx],
        xcoef        = x$xcoef[, idx, drop = FALSE],
        ycoef        = x$ycoef[, idx, drop = FALSE],
        xcoef_scores = x$xcoef_scores[, idx, drop = FALSE],
        ycoef_scores = x$ycoef_scores[, idx, drop = FALSE],
        xcenter      = x$xcenter,
        ycenter      = x$ycenter
      )
    },

    residualize = function(x, k, data, ...) {
      # Deflate by projecting X and Y onto the orthogonal complement of
      # their first-k canonical directions. This removes the rank-k
      # contribution captured by the leading canonical pairs.
      k    <- min(k, length(x$cor))
      idx  <- seq_len(k)
      Wx_k <- x$xcoef[, idx, drop = FALSE]
      Wy_k <- x$ycoef[, idx, drop = FALSE]
      Xd <- data$X - data$X %*% Wx_k %*% t(Wx_k)
      Yd <- data$Y - data$Y %*% Wy_k %*% t(Wy_k)
      list(X = Xd, Y = Yd)
    },

    refit = function(x, new_data, ...) {
      cc <- stats::cancor(new_data$X, new_data$Y)
      list(
        cor          = cc$cor,
        xcoef        = cc$xcoef,
        ycoef        = cc$ycoef,
        xcoef_scores = new_data$X %*% cc$xcoef,
        ycoef_scores = new_data$Y %*% cc$ycoef,
        xcenter      = cc$xcenter,
        ycenter      = cc$ycenter
      )
    },

    null_action = function(x, data, ...) {
      # Break pairing by permuting rows of Y.
      perm <- sample(seq_len(nrow(data$Y)))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },

    component_stat = function(x, data, k, ...) {
      # Relative tail-ratio: rho_k^2 / sum_{q >= k} rho_q^2.
      cc  <- stats::cancor(data$X, data$Y)
      rho <- cc$cor
      if (k > length(rho) || sum(rho[k:length(rho)]^2) == 0) return(NA_real_)
      rho[k]^2 / sum(rho[k:length(rho)]^2)
    },

    validity_level       = "conditional",
    declared_assumptions = c("paired_rows"),
    checked_assumptions  = .cross_baser_checks()
  )
}


# -----------------------------------------------------------------------------
# Executable validity checks shared by adapter_cross_svd() and adapter_cancor().
#
# Both cross adapters assume the data argument is a list with X and Y,
# that nrow(X) == nrow(Y) (paired rows), and that both blocks are
# finite numeric matrices of sane shape. Strict-mode infer() fails
# fast with a named reason if any check fires.
# -----------------------------------------------------------------------------

.is_cross_correlation_recipe <- function(recipe) {
  !is.null(recipe) &&
    inherits(recipe, "multifer_infer_recipe") &&
    identical(recipe$shape$geometry$kind, "cross") &&
    identical(recipe$shape$relation$kind, "correlation")
}

.cross_recipe_design <- function(recipe) {
  if (is.null(recipe) || is.null(recipe$shape) || is.null(recipe$shape$design)) {
    return(NULL)
  }
  recipe$shape$design
}

#' @noRd
.cross_baser_checks <- function() {
  list(
    list(
      name   = "cross_data_is_list_with_xy",
      detail = "cross `data` must be a list containing numeric matrices `X` and `Y`",
      check  = function(data, ...) {
        is.list(data) &&
          !is.null(data$X) && !is.null(data$Y) &&
          is.matrix(data$X) && is.matrix(data$Y) &&
          is.numeric(data$X) && is.numeric(data$Y)
      }
    ),
    list(
      name   = "cross_paired_rows",
      detail = "cross `data$X` and `data$Y` must have the same number of rows (paired-row design)",
      check  = function(data, ...) {
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)
        nrow(data$X) == nrow(data$Y)
      }
    ),
    list(
      name   = "cross_blocks_are_finite",
      detail = "cross `data$X` and `data$Y` must not contain NA, NaN, or Inf",
      check  = function(data, ...) {
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)
        all(is.finite(data$X)) && all(is.finite(data$Y))
      }
    ),
    list(
      name   = "cross_min_dimensions",
      detail = "cross blocks must each have at least 2 rows and at least 1 column",
      check  = function(data, ...) {
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)
        nrow(data$X) >= 2L &&
          ncol(data$X) >= 1L && ncol(data$Y) >= 1L
      }
    ),
    list(
      name   = "cross_columns_have_variance",
      detail = "every column of cross blocks `X` and `Y` must have non-zero sample variance",
      check  = function(data, ...) {
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)
        if (!all(is.finite(data$X)) || !all(is.finite(data$Y))) return(TRUE)
        if (nrow(data$X) < 2L || nrow(data$Y) < 2L) return(TRUE)
        vx <- apply(data$X, 2L, stats::var)
        vy <- apply(data$Y, 2L, stats::var)
        all(vx > 0) && all(vy > 0)
      }
    ),
    list(
      name   = "cross_blocks_full_column_rank",
      detail = paste0(
        "centered cross blocks `X` and `Y` must have full numerical column rank ",
        "(collinear or duplicated columns break the canonical-correlation whitening)"
      ),
      check  = function(data, ...) {
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)
        if (!all(is.finite(data$X)) || !all(is.finite(data$Y))) return(TRUE)
        if (nrow(data$X) < 2L || nrow(data$Y) < 2L) return(TRUE)
        # Center then check rank via QR. Use a generous tolerance relative
        # to the block magnitudes so mild near-collinearity (common in
        # standardized biomedical data) does not fire the check.
        Xc <- sweep(data$X, 2L, colMeans(data$X), "-")
        Yc <- sweep(data$Y, 2L, colMeans(data$Y), "-")
        tol_x <- max(1, max(abs(Xc))) * sqrt(.Machine$double.eps)
        tol_y <- max(1, max(abs(Yc))) * sqrt(.Machine$double.eps)
        rx <- qr(Xc, tol = tol_x)$rank
        ry <- qr(Yc, tol = tol_y)$rank
        rx == ncol(Xc) && ry == ncol(Yc)
      }
    ),
    list(
      name   = "cross_correlation_sample_size",
      detail = paste0(
        "correlation-mode whitening requires enough effective rows: ",
        "for paired/blocked designs n must exceed p + q; for nuisance-adjusted ",
        "designs n - rank(Z) must exceed p + q"
      ),
      check  = function(data, recipe = NULL, ...) {
        if (!.is_cross_correlation_recipe(recipe)) return(TRUE)
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)

        n <- nrow(data$X)
        p <- ncol(data$X)
        q <- ncol(data$Y)
        design <- .cross_recipe_design(recipe)
        design_kind <- design$kind %||% "paired_rows"

        if (identical(design_kind, "nuisance_adjusted")) {
          Z <- design$Z
          if (is.null(Z) || !is.matrix(Z) || !is.numeric(Z)) return(TRUE)
          r <- qr(Z, tol = 1e-10)$rank
          n_eff <- n - r
          return(list(
            passed = n_eff > (p + q),
            detail = sprintf(
              "nuisance-adjusted correlation requires n - rank(Z) = %d to exceed p + q = %d.",
              n_eff, p + q
            )
          ))
        }

        list(
          passed = n > (p + q),
          detail = sprintf(
            "correlation-mode whitening requires n = %d to exceed p + q = %d.",
            n, p + q
          )
        )
      }
    ),
    list(
      name   = "cross_nuisance_design_rank",
      detail = paste0(
        "nuisance-adjusted correlation requires a full-column-rank Z and ",
        "enough residual degrees of freedom for whitening"
      ),
      check  = function(data, recipe = NULL, ...) {
        if (!.is_cross_correlation_recipe(recipe)) return(TRUE)
        design <- .cross_recipe_design(recipe)
        if (is.null(design) || !identical(design$kind, "nuisance_adjusted")) {
          return(TRUE)
        }
        Z <- design$Z
        if (is.null(Z) || !is.matrix(Z) || !is.numeric(Z)) return(TRUE)
        if (!is.list(data) || is.null(data$X) || is.null(data$Y)) return(TRUE)
        if (!is.matrix(data$X) || !is.matrix(data$Y)) return(TRUE)

        r <- qr(Z, tol = 1e-10)$rank
        z_full_rank <- identical(r, ncol(Z))
        n_eff <- nrow(data$X) - r
        p_plus_q <- ncol(data$X) + ncol(data$Y)

        if (!z_full_rank) {
          return(list(
            passed = FALSE,
            detail = sprintf(
              "nuisance matrix Z has numerical rank %d but %d columns; full column rank is required.",
              r, ncol(Z)
            )
          ))
        }

        list(
          passed = n_eff > p_plus_q,
          detail = sprintf(
            "nuisance-adjusted correlation leaves n - rank(Z) = %d residual rows, which must exceed p + q = %d.",
            n_eff, p_plus_q
          )
        )
      }
    ),
    list(
      name   = "cross_grouped_design_consistency",
      detail = paste0(
        "grouped correlation designs require at least one exchangeable block ",
        "with size > 1, and nuisance-adjusted grouped designs must retain such a block ",
        "after residual-basis row reduction"
      ),
      check  = function(data, recipe = NULL, ...) {
        if (!.is_cross_correlation_recipe(recipe)) return(TRUE)
        design <- .cross_recipe_design(recipe)
        if (is.null(design)) return(TRUE)

        groups <- NULL
        if (identical(design$kind, "blocked_rows")) {
          groups <- design$groups
        } else if (identical(design$kind, "nuisance_adjusted") && !is.null(design$groups)) {
          groups <- design$groups
        } else {
          return(TRUE)
        }

        groups <- as.factor(groups)
        group_sizes <- table(groups)
        if (!any(group_sizes > 1L)) {
          return(list(
            passed = FALSE,
            detail = "grouped correlation design is singleton-only; restricted permutation would be the identity."
          ))
        }

        if (!identical(design$kind, "nuisance_adjusted")) {
          return(TRUE)
        }

        Z <- design$Z
        if (is.null(Z) || !is.matrix(Z) || !is.numeric(Z)) return(TRUE)
        r <- qr(Z, tol = 1e-10)$rank
        n_keep <- nrow(Z) - r
        if (n_keep <= 0L) {
          return(list(
            passed = FALSE,
            detail = "nuisance-adjusted grouped correlation leaves no residual rows after projecting out Z."
          ))
        }

        keep_idx <- tryCatch(
          .theil_keep_indices(groups, n_keep),
          error = function(e) e
        )
        if (inherits(keep_idx, "error")) {
          return(list(
            passed = FALSE,
            detail = paste0("group structure is incompatible with the residual-basis reduction: ",
                            conditionMessage(keep_idx))
          ))
        }

        reduced_sizes <- table(droplevels(groups[keep_idx]))
        list(
          passed = any(reduced_sizes > 1L),
          detail = "after residual-basis reduction, grouped correlation must retain at least one exchangeable block with size > 1."
        )
      }
    )
  )
}


#' Register the base-R cross reference adapters
#'
#' Called from the package `.onLoad` hook. Safe to call repeatedly.
#'
#' @keywords internal
register_cross_baser_adapters <- function() {
  register_infer_adapter(
    adapter_id = "cross_svd",
    adapter    = adapter_cross_svd(),
    overwrite  = TRUE
  )
  register_infer_adapter(
    adapter_id = "cancor_cross",
    adapter    = adapter_cancor(),
    overwrite  = TRUE
  )
}
