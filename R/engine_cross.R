#' Run the cross sequential deflation engine
#'
#' Implements the test ladder for (cross, covariance) and
#' (cross, correlation) recipes via the shared
#' \code{\link{ladder_driver}}. Both relations share the same scaffold
#' but use different per-step statistics:
#'
#' \itemize{
#'   \item covariance: top singular value of the centered cross-product
#'     \code{X^t Y} on the deflated residual, normalized by the squared
#'     Frobenius mass of the remaining singular values.
#'   \item correlation: top singular value of the orthonormal
#'     cross-product \code{Q_x^t Q_y} on the deflated residual, where
#'     \code{Q_x} and \code{Q_y} are the thin QR factors of the
#'     centered residual blocks.
#' }
#'
#' Both use the Part 1 section 1 simplification: the top singular value
#' of the deflated cross-block matrix at step \code{a} is the a-th
#' latent root of the original problem.
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
                             B          = 1000L,
                             B_total    = NULL,
                             batch_size = 32L,
                             alpha      = 0.05,
                             max_steps  = NULL,
                             seed       = NULL) {

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

  ## --- relation-specific cross statistic --------------------------------------
  # cross_stat(X_residual, Y_residual) returns the FULL singular-value
  # vector of the relation-appropriate cross matrix. observed_stat_fn and
  # null_stat_fn both call this and apply the same tail-ratio.

  cross_stat <- if (rel_kind == "covariance") {
    function(Xr, Yr) {
      cached_svd(crossprod(Xr, Yr))$d
    }
  } else {
    function(Xr, Yr) {
      # Thin QR for whitening; fall back to centered if rank-deficient.
      qx <- qr(Xr)
      qy <- qr(Yr)
      Qx <- qr.Q(qx)
      Qy <- qr.Q(qy)
      cached_svd(crossprod(Qx, Qy))$d
    }
  }

  ## --- full observed root vector (for form_units) -----------------------------

  s_full <- cross_stat(Xc, Yc)
  roots_observed <- s_full^2
  zero_tol <- max(1, sum(roots_observed)) * .Machine$double.eps

  ## --- max_steps default ------------------------------------------------------

  if (is.null(max_steps)) {
    max_steps <- min(min(nrow(Xc), ncol(Xc), ncol(Yc)) - 1L, 50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  ## --- callbacks --------------------------------------------------------------

  observed_stat_fn <- function(step, data) {
    s  <- cross_stat(data$X, data$Y)
    s2 <- s^2
    total <- sum(s2)
    if (total <= zero_tol || length(s2) == 0L) return(0)
    s2[1L] / total
  }

  null_stat_fn <- function(step, data) {
    perm <- base::sample.int(nrow(data$Y))
    Yp   <- data$Y[perm, , drop = FALSE]
    s    <- cross_stat(data$X, Yp)
    s2   <- s^2
    total <- sum(s2)
    if (total <= zero_tol || length(s2) == 0L) return(0)
    s2[1L] / total
  }

  deflate_fn <- function(step, data) {
    # Remove the rank-1 contribution of the top cross-pair from BOTH
    # blocks. For covariance: take the top left/right singular vectors
    # of crossprod(X, Y) and subtract their projection from X and Y.
    # For correlation: same idea but using the whitened SVD.
    if (rel_kind == "covariance") {
      sv <- svd(crossprod(data$X, data$Y))
      u1 <- sv$u[, 1L, drop = FALSE]
      v1 <- sv$v[, 1L, drop = FALSE]
      Xn <- data$X - data$X %*% u1 %*% t(u1)
      Yn <- data$Y - data$Y %*% v1 %*% t(v1)
    } else {
      qx <- qr(data$X)
      qy <- qr(data$Y)
      Qx <- qr.Q(qx)
      Qy <- qr.Q(qy)
      sv <- svd(crossprod(Qx, Qy))
      # Map whitened directions back to original variable space.
      Wx <- backsolve(qr.R(qx), sv$u[, 1L, drop = FALSE])
      Wy <- backsolve(qr.R(qy), sv$v[, 1L, drop = FALSE])
      # Normalize columns so the projector is well-defined.
      norm_or_one <- function(w) {
        nm <- sqrt(sum(w^2))
        if (nm < 1e-12) w else w / nm
      }
      Wx <- apply(Wx, 2L, norm_or_one)
      Wy <- apply(Wy, 2L, norm_or_one)
      if (!is.matrix(Wx)) Wx <- matrix(Wx, ncol = 1L)
      if (!is.matrix(Wy)) Wy <- matrix(Wy, ncol = 1L)
      Xn <- data$X - data$X %*% Wx %*% t(Wx)
      Yn <- data$Y - data$Y %*% Wy %*% t(Wy)
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
    initial_data     = list(X = Xc, Y = Yc),
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
