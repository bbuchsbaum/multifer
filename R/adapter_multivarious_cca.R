#' multifer adapter for multivarious::cca
#'
#' A Tier-2 adapter wrapping \code{multivarious::cca} for the
#' (cross, correlation) shape. Requires the multivarious package at
#' invocation time only -- constructing and registering the adapter
#' does not itself require multivarious to be loaded.
#'
#' This is a SINGLE-relation adapter. It declares only correlation
#' because canonical correlation is fundamentally a correlation-mode
#' method. Users wanting covariance-mode inference on cross data should
#' use \code{adapter_cross_svd()} or \code{adapter_multivarious_plsc()}.
#'
#' Declares \code{component_significance}, \code{variable_stability},
#' \code{score_stability}, and \code{subspace_stability} for
#' (cross, correlation). Does NOT claim \code{variable_significance}
#' (deferred to Phase 3 per Part 5 section 38).
#'
#' @param adapter_id Character, default \code{"multivarious_cca"}.
#' @param adapter_version Character, default \code{"0.0.1"}.
#' @param ncomp Integer, number of components to keep when refitting.
#'   Default \code{NULL} lets the constructor pick a sensible value.
#' @param preproc_x Preprocessor for the X block. Default \code{NULL}
#'   delegates to \code{multivarious::cca()}.
#' @param preproc_y Preprocessor for the Y block. Default \code{NULL}
#'   delegates to \code{multivarious::cca()}.
#' @param lambda Shared ridge shrinkage passed through to
#'   \code{multivarious::cca()} when block-specific lambdas are not
#'   supplied.
#' @param lambda_x Optional X-block ridge shrinkage.
#' @param lambda_y Optional Y-block ridge shrinkage.
#'
#' @return A \code{multifer_adapter}.
#' @export
adapter_multivarious_cca <- function(adapter_id = "multivarious_cca",
                                     adapter_version = "0.0.1",
                                     ncomp = NULL,
                                     preproc_x = NULL,
                                     preproc_y = NULL,
                                     lambda = 1e-4,
                                     lambda_x = NULL,
                                     lambda_y = NULL) {

  require_multivarious <- function() {
    if (!requireNamespace("multivarious", quietly = TRUE)) {
      stop("multivarious is required by adapter_multivarious_cca but is not installed.",
           call. = FALSE)
    }
  }

  fit_cca <- function(new_data) {
    require_multivarious()
    if (is.null(new_data$X) || is.null(new_data$Y)) {
      stop("adapter_multivarious_cca: new_data must be a list with X and Y.",
           call. = FALSE)
    }

    args <- list(
      X = new_data$X,
      Y = new_data$Y,
      ncomp = ncomp,
      lambda = lambda
    )
    if (!is.null(preproc_x)) args$preproc_x <- preproc_x
    if (!is.null(preproc_y)) args$preproc_y <- preproc_y
    if (!is.null(lambda_x)) args$lambda_x <- lambda_x
    if (!is.null(lambda_y)) args$lambda_y <- lambda_y

    do.call(multivarious::cca, args)
  }

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(geometry = "cross", relation = "correlation",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),

    roots = function(x, ...) {
      require_multivarious()
      x$cor
    },

    scores = function(x, domain = c("X", "Y"), ...) {
      require_multivarious()
      domain <- match.arg(domain)
      multivarious::scores(x, block = domain)
    },

    loadings = function(x, domain = c("X", "Y"), ...) {
      require_multivarious()
      domain <- match.arg(domain)
      stats::coef(x, source = domain)
    },

    truncate = function(x, k, ...) {
      require_multivarious()
      multivarious::truncate(x, k)
    },

    residualize = function(x, k, data, ...) {
      require_multivarious()
      k <- min(k, length(x$cor))
      idx <- seq_len(k)
      Wx_k <- stats::coef(x, source = "X")[, idx, drop = FALSE]
      Wy_k <- stats::coef(x, source = "Y")[, idx, drop = FALSE]
      list(
        X = data$X - data$X %*% Wx_k %*% t(Wx_k),
        Y = data$Y - data$Y %*% Wy_k %*% t(Wy_k)
      )
    },

    refit = function(x, new_data, ...) {
      fit_cca(new_data)
    },

    null_action = function(x, data, ...) {
      perm <- sample(seq_len(nrow(data$Y)))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },

    component_stat = function(x, data, k, ...) {
      fit <- fit_cca(data)
      rho <- fit$cor
      if (k > length(rho) || sum(rho[k:length(rho)]^2) == 0) return(NA_real_)
      rho[k]^2 / sum(rho[k:length(rho)]^2)
    },

    validity_level       = "conditional",
    declared_assumptions = c("paired_rows"),
    checked_assumptions  = .cross_baser_checks()
  )
}


#' Register the multivarious::cca adapter when multivarious is installed
#'
#' Called from the package \code{.onLoad} hook. Skips registration silently
#' if \code{multivarious} is not installed.
#'
#' @keywords internal
register_multivarious_cca_adapter <- function() {
  if (!requireNamespace("multivarious", quietly = TRUE)) {
    return(invisible(NULL))
  }
  register_infer_adapter(
    adapter_id = "multivarious_cca",
    adapter    = adapter_multivarious_cca(),
    overwrite  = TRUE
  )
  invisible(NULL)
}
