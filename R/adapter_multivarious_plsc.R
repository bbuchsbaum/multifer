#' multifer adapter for multivarious::plsc
#'
#' A Tier-2 adapter wrapping \code{multivarious::plsc} for the
#' (cross, covariance) shape. Requires the multivarious package at
#' invocation time only -- constructing and registering the adapter
#' does not itself require multivarious to be loaded.
#'
#' This is a SINGLE-relation adapter. It declares only covariance
#' because plsc is fundamentally a covariance method. Users wanting
#' correlation-based inference on cross data should use
#' \code{adapter_cross_svd()} or \code{adapter_cancor()} instead.
#' Because only one relation is declared, the strict-dispatch
#' ambiguity error never fires when \code{infer()} is called on this
#' adapter without an explicit relation -- the single declared
#' relation is auto-used.
#'
#' Declares \code{component_significance}, \code{variable_stability},
#' \code{score_stability}, and \code{subspace_stability} for
#' (cross, covariance). Does NOT claim \code{variable_significance}
#' (deferred to Phase 3 per Part 5 section 38).
#'
#' @param adapter_id Character, default \code{"multivarious_plsc"}.
#' @param adapter_version Character, default \code{"0.0.1"}.
#' @param ncomp Integer, number of components to keep when refitting.
#'   Default \code{NULL} lets the constructor pick a sensible value
#'   based on data shape.
#'
#' @return A \code{multifer_adapter}.
#' @export
adapter_multivarious_plsc <- function(adapter_id = "multivarious_plsc",
                                      adapter_version = "0.0.1",
                                      ncomp = NULL) {

  require_multivarious <- function() {
    if (!requireNamespace("multivarious", quietly = TRUE)) {
      stop("multivarious is required by adapter_multivarious_plsc but is not installed.",
           call. = FALSE)
    }
  }

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "cross",
    capabilities    = capability_matrix(
      list(geometry = "cross", relation = "covariance",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),

    roots = function(x, ...) {
      # plsc latent roots are the singular values of the cross product
      # X^t Y, stored as `singvals`. We expose them squared so the
      # convention matches the variance shape (eigenvalue scale).
      require_multivarious()
      x$singvals^2
    },

    scores = function(x, domain = c("X", "Y"), ...) {
      require_multivarious()
      domain <- match.arg(domain)
      multivarious::scores(x, source = domain)
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
      # data is a list with X and Y. Deflate both blocks by their
      # rank-k projection onto the first k plsc loadings.
      require_multivarious()
      Wx <- stats::coef(x, source = "X")[, seq_len(k), drop = FALSE]
      Wy <- stats::coef(x, source = "Y")[, seq_len(k), drop = FALSE]
      list(
        X = data$X - data$X %*% Wx %*% t(Wx),
        Y = data$Y - data$Y %*% Wy %*% t(Wy)
      )
    },

    refit = function(x, new_data, ...) {
      require_multivarious()
      if (is.null(new_data$X) || is.null(new_data$Y)) {
        stop("adapter_multivarious_plsc: new_data must be a list with X and Y.",
             call. = FALSE)
      }
      k <- if (is.null(ncomp)) {
        max(1L, min(nrow(new_data$X),
                    ncol(new_data$X),
                    ncol(new_data$Y)) - 1L)
      } else {
        ncomp
      }
      multivarious::plsc(new_data$X, new_data$Y, ncomp = k)
    },

    null_action = function(x, data, ...) {
      # Row permutation of Y only: breaks the cross-block pairing.
      perm <- sample(seq_len(nrow(data$Y)))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },

    component_stat = function(x, data, k, ...) {
      # Covariance-mode tail ratio: singular values of X^t Y.
      s <- svd(crossprod(data$X, data$Y))$d
      if (k > length(s) || sum(s[k:length(s)]^2) == 0) return(NA_real_)
      s[k]^2 / sum(s[k:length(s)]^2)
    },

    validity_level       = "conditional",
    declared_assumptions = c("paired_rows", "centered_blocks"),
    checked_assumptions  = list()
  )
}


#' Register the multivarious::plsc adapter when multivarious is installed
#'
#' Called from the package \code{.onLoad} hook. Skips registration silently
#' if \code{multivarious} is not installed.
#'
#' @keywords internal
register_multivarious_plsc_adapter <- function() {
  if (!requireNamespace("multivarious", quietly = TRUE)) {
    return(invisible(NULL))
  }
  register_infer_adapter(
    adapter_id = "multivarious_plsc",
    adapter    = adapter_multivarious_plsc(),
    overwrite  = TRUE
  )
  invisible(NULL)
}
