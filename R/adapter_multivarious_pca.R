#' multifer adapter for multivarious::pca
#'
#' A Tier-2 adapter wrapping \code{multivarious::pca} for the
#' (oneblock, variance) shape. Requires the multivarious package
#' at invocation time only -- constructing and registering the
#' adapter does not itself require multivarious to be loaded.
#'
#' Declares \code{component_significance}, \code{variable_stability},
#' \code{score_stability}, and \code{subspace_stability} for
#' (oneblock, variance). Does NOT claim \code{variable_significance}
#' (deferred to Phase 3 per Part 5 section 38).
#'
#' @param adapter_id Character, default \code{"multivarious_pca"}.
#' @param adapter_version Character, default \code{"0.0.1"}.
#' @param ncomp Integer, number of components to keep in refitted
#'   fits. Default \code{NULL} lets \code{multivarious::pca} choose.
#'
#' @return A \code{multifer_adapter}.
#' @export
adapter_multivarious_pca <- function(adapter_id = "multivarious_pca",
                                     adapter_version = "0.0.1",
                                     ncomp = NULL) {

  # Hooks dispatch to multivarious generics at invocation time.
  require_multivarious <- function() {
    if (!requireNamespace("multivarious", quietly = TRUE)) {
      stop("multivarious is required by adapter_multivarious_pca but is not installed.",
           call. = FALSE)
    }
  }

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = adapter_version,
    shape_kinds     = "oneblock",
    capabilities    = capability_matrix(
      list(geometry = "oneblock", relation = "variance",
           targets  = c("component_significance", "variable_stability",
                        "score_stability", "subspace_stability"))
    ),

    roots = function(x, ...) {
      require_multivarious()
      multivarious::sdev(x)^2
    },

    scores = function(x, domain = NULL, ...) {
      require_multivarious()
      multivarious::scores(x)
    },

    loadings = function(x, domain = NULL, ...) {
      require_multivarious()
      # coef() on a bi_projector returns the V (right singular vectors / loadings)
      stats::coef(x)
    },

    truncate = function(x, k, ...) {
      require_multivarious()
      multivarious::truncate(x, k)
    },

    residualize = function(x, k, data, ...) {
      # Subtract rank-k reconstruction using multivarious generics.
      require_multivarious()
      V <- stats::coef(x)[, seq_len(k), drop = FALSE]
      S <- multivarious::scores(x)[, seq_len(k), drop = FALSE]
      data - tcrossprod(S, V)
    },

    refit = function(x, new_data, ...) {
      require_multivarious()
      if (is.list(new_data) && !is.null(new_data$X)) {
        new_data <- new_data$X
      }
      k <- if (is.null(ncomp)) min(dim(new_data)) - 1L else ncomp
      multivarious::pca(new_data, ncomp = k)
    },

    null_action = function(x, data, ...) {
      # Column-wise permutation; no multivarious needed.
      apply(data, 2L, sample)
    },

    component_stat = function(x, data, k, ...) {
      # Relative tail-ratio analog for the variance shape:
      # lambda_k / sum_{q >= k} lambda_q on SVD of data.
      s  <- svd(data)$d
      s2 <- s^2
      if (k > length(s2) || sum(s2[k:length(s2)]) == 0) return(NA_real_)
      s2[k] / sum(s2[k:length(s2)])
    },

    validity_level       = "conditional",
    declared_assumptions = c("rows_exchangeable"),
    checked_assumptions  = list()
  )
}


#' Register the multivarious::pca adapter when multivarious is installed
#'
#' Called from the package \code{.onLoad} hook. Skips registration silently
#' if \code{multivarious} is not installed, so users on a minimal R install
#' never see a failure.
#'
#' @keywords internal
register_multivarious_pca_adapter <- function() {
  if (!requireNamespace("multivarious", quietly = TRUE)) {
    return(invisible(NULL))
  }
  register_infer_adapter(
    adapter_id = "multivarious_pca",
    adapter    = adapter_multivarious_pca(),
    overwrite  = TRUE
  )
  invisible(NULL)
}
