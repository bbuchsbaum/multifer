#' Run the cross sequential deflation engine
#'
#' Implements the test ladder for (cross, covariance) and
#' (cross, correlation) recipes via the shared
#' \code{\link{ladder_driver}}. Both relations share the same scaffold
#' but use different per-step statistics:
#'
#' \itemize{
#'   \item covariance: squared leading singular value of the centered
#'     cross-product \code{X^t Y} on the deflated residual.
#'   \item correlation: squared leading singular value of the
#'     orthonormal cross-product \code{Q_x^t Q_y} on the deflated
#'     residual, where \code{Q_x} and \code{Q_y} are the thin QR
#'     factors of the centered residual blocks.
#' }
#'
#' For covariance-mode cross problems, the top singular value of the
#' deflated cross-block matrix at step \code{a} is the \code{a}-th
#' latent root of the original problem, so root-strength testing on the
#' deflated problem gives the intended stepwise ladder.
#'
#' For correlation-mode cross problems, valid multi-root testing requires
#' stepwise residualization in the whitened CCA space. `multifer`
#' implements that ladder for the plain paired-row design, and for
#' nuisance-adjusted designs it first projects both blocks into an
#' exchangeability-preserving residual basis before running the same
#' stepwise ladder. More complex designs still fall back to the
#' conservative first-root cap.
#'
#' Phase 1 is refit-first; the Part 2 section 9 core-space update is
#' a Phase 1.5 optimization.
#'
#' @param recipe A compiled \code{multifer_infer_recipe} with geometry
#'   \code{"cross"} and relation \code{"covariance"} or
#'   \code{"correlation"}.
#' @param X Numeric matrix, n x p.
#' @param Y Numeric matrix, n x q. Must have the same number of rows
#'   as \code{X}.
#' @param B Per-rung cap on Monte Carlo draws.
#' @param B_total Optional integer global Monte Carlo budget shared
#'   across ladder rungs. Defaults to `B * max_steps`.
#' @param batch_size Positive integer, Besag-Clifford batch size within
#'   each rung. Default `32L`.
#' @param alpha Significance threshold.
#' @param max_steps Maximum ladder rungs. Default
#'   \code{min(nrow, ncol(X), ncol(Y)) - 1}, capped at 50.
#' @param seed Integer or NULL.
#'
#' @return A list with \code{units}, \code{component_tests},
#'   \code{roots_observed}, and \code{ladder_result}, mirroring
#'   \code{\link{run_oneblock_ladder}}.
#'
#' @export
run_cross_ladder <- function(recipe,
                             X,
                             Y,
                             B              = 1000L,
                             B_total        = NULL,
                             batch_size     = 32L,
                             alpha          = 0.05,
                             max_steps      = NULL,
                             seed           = NULL,
                             cross_rank_cap = 50L) {

  orth_basis <- function(M, tol = 1e-10) {
    q <- qr(M, tol = tol)
    r <- q$rank
    if (r <= 0L) {
      return(matrix(0, nrow = nrow(M), ncol = 0L))
    }
    qr.Q(q, complete = FALSE)[, seq_len(r), drop = FALSE]
  }

  residual_basis <- function(Z, n, tol = 1e-10) {
    if (!is.matrix(Z) || !is.numeric(Z)) {
      stop("For nuisance-adjusted correlation designs, `design$Z` must be a numeric matrix.",
           call. = FALSE)
    }
    if (nrow(Z) != n) {
      stop("For nuisance-adjusted correlation designs, `design$Z` must have one row per observation.",
           call. = FALSE)
    }
    qz <- qr(Z, tol = tol)
    r <- qz$rank
    if (r >= n) {
      return(matrix(0, nrow = n, ncol = 0L))
    }
    q_full <- qr.Q(qz, complete = TRUE)
    q_full[, seq.int(r + 1L, n), drop = FALSE]
  }

  ## --- validate recipe --------------------------------------------------------

  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "cross") {
    stop(paste0("`recipe` geometry must be \"cross\"; got \"",
                recipe$shape$geometry$kind, "\"."), call. = FALSE)
  }
  rel_kind <- recipe$shape$relation$kind
  if (!(rel_kind %in% c("covariance", "correlation"))) {
    stop(paste0("`recipe` relation must be \"covariance\" or \"correlation\"; got \"",
                rel_kind, "\"."), call. = FALSE)
  }

  ## --- validate X, Y ----------------------------------------------------------

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("`X` and `Y` must not contain NA, NaN, or Inf.", call. = FALSE)
  }
  if (nrow(X) < 2L || ncol(X) < 1L || ncol(Y) < 1L) {
    stop("`X` must have at least 2 rows; `X` and `Y` need at least 1 column each.",
         call. = FALSE)
  }

  ## --- center both blocks -----------------------------------------------------

  Xc <- sweep(X, 2L, colMeans(X), "-")
  Yc <- sweep(Y, 2L, colMeans(Y), "-")
  design_kind <- recipe$shape$design$kind
  allow_multiroot_correlation <- rel_kind == "correlation" &&
    design_kind %in% c("paired_rows", "nuisance_adjusted")

  ## --- relation-specific cross statistic --------------------------------------
  # cross_stat(X_residual, Y_residual) returns the FULL singular-value
  # vector of the relation-appropriate cross matrix. observed_stat_fn and
  # null_stat_fn both call this and apply the same tail-ratio.

  cross_stat <- if (rel_kind == "covariance") {
    function(Xr, Yr) {
      cached_svd(crossprod(Xr, Yr))$d
    }
  } else {
    function(Qx, Qy) {
      if (ncol(Qx) == 0L || ncol(Qy) == 0L) {
        return(0)
      }
      cached_svd(crossprod(Qx, Qy))$d
    }
  }

  ## --- full observed root vector (for form_units) -----------------------------

  qx_full <- NULL
  qy_full <- NULL
  if (rel_kind == "covariance") {
    s_full <- cross_stat(Xc, Yc)
  } else {
    if (design_kind == "nuisance_adjusted") {
      Qz <- residual_basis(recipe$shape$design$Z, nrow(Xc))
      X_corr <- crossprod(Qz, Xc)
      Y_corr <- crossprod(Qz, Yc)
    } else {
      X_corr <- Xc
      Y_corr <- Yc
    }
    qx_full <- orth_basis(X_corr)
    qy_full <- orth_basis(Y_corr)
    s_full <- cross_stat(qx_full, qy_full)
  }
  roots_observed <- s_full^2
  zero_tol <- max(1, sum(roots_observed)) * .Machine$double.eps

  ## --- max_steps default ------------------------------------------------------

  if (is.null(max_steps)) {
    if (rel_kind == "correlation" && allow_multiroot_correlation) {
      max_steps <- min(length(roots_observed), 50L)
    } else {
      max_steps <- min(min(nrow(Xc), ncol(Xc), ncol(Yc)) - 1L, 50L)
    }
  }
  max_steps <- as.integer(max(1L, max_steps))
  if (rel_kind == "correlation" && !allow_multiroot_correlation) {
    max_steps <- min(max_steps, 1L)
  }

  ## --- callbacks --------------------------------------------------------------

  # Build the cross matrix directly and compute (top-1 SV)^2.
  # Only the leading singular value is needed for the test statistic --
  # we hand the matrix to top_singular_values() which routes through
  # RSpectra::svds(, k=1) when available, giving a large speedup over
  # a full SVD on the cross-product.
  cross_matrix <- if (rel_kind == "covariance") {
    function(Xr, Yr) base::crossprod(Xr, Yr)
  } else {
    function(Qx, Qy) {
      if (ncol(Qx) == 0L || ncol(Qy) == 0L) {
        return(matrix(0, nrow = 1L, ncol = 1L))
      }
      base::crossprod(Qx, Qy)
    }
  }

  # Phase 1.5 follow-up: rank-capped core approximation for rung 1.
  #
  # For (cross, covariance) the top singular value of Xc^T · Yc[perm,] is
  # the quantity we need per null draw. If we thin-SVD both blocks once
  # before the ladder, Xc = Ux Dx Vx^T and Yc = Uy Dy Vy^T, then
  #   Xc^T · Yc[perm,]  =  Vx Dx (Ux^T · P Uy) Dy Vy^T
  # where P permutes rows. Since Vx / Vy have orthonormal columns, the
  # top singular value of that product equals the top singular value of
  #   M  =  Dx · (Ux^T · Uy[perm,]) · Dy
  # which is a k_x x k_y inner matrix. Truncating the thin SVDs to
  # cross_rank_cap makes the per-draw work O(k_x * n * k_y) instead of
  # O(n * p * q); with the default k = 50 on n=500, p=250, q=200 this is
  # a ~20x drop in flops on the per-draw hot spot.
  #
  # The approximation only fires on rung 1 because deeper rungs see a
  # deflated (Xd, Yd) that we do not have a cheap cached SVD for. Most
  # ladder runs stop at rung 1 under null-ish data, so this is where
  # the savings live.
  use_fast_rung1 <- rel_kind == "covariance" &&
                    is.numeric(cross_rank_cap) &&
                    length(cross_rank_cap) == 1L &&
                    cross_rank_cap >= 1L

  core_Ux <- NULL; core_Dx <- NULL
  core_Uy <- NULL; core_Dy <- NULL
  if (use_fast_rung1) {
    k_max_x <- min(as.integer(cross_rank_cap), nrow(Xc) - 1L, ncol(Xc))
    k_max_y <- min(as.integer(cross_rank_cap), nrow(Yc) - 1L, ncol(Yc))
    if (k_max_x >= 1L && k_max_y >= 1L) {
      sv_x <- top_svd(Xc, k_max_x)
      sv_y <- top_svd(Yc, k_max_y)
      core_Ux <- sv_x$u
      core_Dx <- sv_x$d
      core_Uy <- sv_y$u
      core_Dy <- sv_y$d
    } else {
      use_fast_rung1 <- FALSE
    }
  }

  rung1_top_sv2 <- function(Uy_view) {
    # M = diag(Dx) · (Ux^T · Uy_view) · diag(Dy), a k_x x k_y inner matrix.
    # Return the top singular value squared of M.
    inner <- base::crossprod(core_Ux, Uy_view)          # k_x x k_y
    M     <- base::sweep(inner, 1L, core_Dx, `*`)       # scale rows
    M     <- base::sweep(M,     2L, core_Dy, `*`)       # scale columns
    s1    <- top_singular_values(M, 1L)[1L]
    s1 * s1
  }

  observed_stat_fn <- function(step, data) {
    if (use_fast_rung1 && step == 1L) {
      return(rung1_top_sv2(core_Uy))
    }
    if (rel_kind == "covariance") {
      M <- cross_matrix(data$X, data$Y)
    } else {
      M <- cross_matrix(data$Qx, data$Qy)
    }
    s1 <- top_singular_values(M, 1L)[1L]
    s1 * s1
  }

  null_stat_fn <- function(step, data) {
    if (use_fast_rung1 && step == 1L) {
      perm <- base::sample.int(nrow(core_Uy))
      return(rung1_top_sv2(core_Uy[perm, , drop = FALSE]))
    }
    if (rel_kind == "covariance") {
      perm <- base::sample.int(nrow(data$Y))
      Yp   <- data$Y[perm, , drop = FALSE]
      M    <- cross_matrix(data$X, Yp)
    } else {
      perm <- base::sample.int(nrow(data$Qy))
      Qyp  <- data$Qy[perm, , drop = FALSE]
      M    <- cross_matrix(data$Qx, Qyp)
    }
    s1   <- top_singular_values(M, 1L)[1L]
    s1 * s1
  }

  deflate_fn <- function(step, data) {
    # Remove the rank-1 contribution of the top cross-pair from BOTH
    # blocks. For covariance: take the top left/right singular vectors
    # of crossprod(X, Y) and subtract their projection from X and Y.
    # For correlation: same idea but using the whitened SVD. Only the
    # top factor is needed, so route through the partial-SVD helper.
    if (rel_kind == "covariance") {
      sv <- top_svd(crossprod(data$X, data$Y), 1L)
      u1 <- sv$u[, 1L, drop = FALSE]
      v1 <- sv$v[, 1L, drop = FALSE]
      Xn <- data$X - data$X %*% u1 %*% t(u1)
      Yn <- data$Y - data$Y %*% v1 %*% t(v1)
    } else {
      if (ncol(data$Qx) == 0L || ncol(data$Qy) == 0L) {
        return(list(Qx = data$Qx, Qy = data$Qy))
      }
      sv <- top_svd(crossprod(data$Qx, data$Qy), 1L)
      tx <- data$Qx %*% sv$u[, 1L, drop = FALSE]
      ty <- data$Qy %*% sv$v[, 1L, drop = FALSE]
      Qx_next <- orth_basis(data$Qx - tx %*% crossprod(tx, data$Qx))
      Qy_next <- orth_basis(data$Qy - ty %*% crossprod(ty, data$Qy))
      return(list(Qx = Qx_next, Qy = Qy_next))
    }
    if (sum(Xn^2) + sum(Yn^2) <= zero_tol) {
      return(list(
        X = matrix(0, nrow = nrow(data$X), ncol = ncol(data$X)),
        Y = matrix(0, nrow = nrow(data$Y), ncol = ncol(data$Y))
      ))
    }
    list(X = Xn, Y = Yn)
  }

  ## --- run the ladder ---------------------------------------------------------

  ladder_result <- ladder_driver(
    observed_stat_fn = observed_stat_fn,
    null_stat_fn     = null_stat_fn,
    deflate_fn       = deflate_fn,
    initial_data     = if (rel_kind == "covariance") {
      list(X = Xc, Y = Yc)
    } else {
      list(Qx = qx_full, Qy = qy_full)
    },
    max_steps        = max_steps,
    B                = B,
    B_total          = B_total,
    batch_size       = batch_size,
    alpha            = alpha,
    seed             = seed
  )

  ## --- selected vector + units ------------------------------------------------

  rejected_through <- ladder_result$rejected_through
  n_roots          <- length(roots_observed)
  selected         <- logical(n_roots)
  if (rejected_through >= 1L) {
    selected[seq_len(rejected_through)] <- TRUE
  }

  units <- form_units(roots_observed, selected = selected)

  ## --- component_tests data.frame ---------------------------------------------

  sr <- ladder_result$step_results
  component_tests <- data.frame(
    step          = vapply(sr, function(x) x$step,          integer(1L)),
    observed_stat = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value       = vapply(sr, function(x) x$p_value,       double(1L)),
    mc_se         = vapply(sr, function(x) x$mc_se,         double(1L)),
    r             = vapply(sr, function(x) x$r,             integer(1L)),
    B             = vapply(sr, function(x) x$B,             integer(1L)),
    selected      = vapply(sr, function(x) x$selected,      logical(1L)),
    stringsAsFactors = FALSE
  )

  list(
    units           = units,
    component_tests = component_tests,
    roots_observed  = roots_observed,
    ladder_result   = ladder_result
  )
}
