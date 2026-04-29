#' multifer adapter for multivarious::pca
#'
#' A Tier-2 adapter wrapping \code{multivarious::pca} for the
#' (oneblock, variance) shape. \pkg{multivarious} is a hard dependency
#' of \pkg{multifer}; this adapter is the default engine behind
#' \code{\link[=infer_pca]{infer_pca()}}.
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

    feature_evidence_action = function(x,
                                       data,
                                       units,
                                       design,
                                       statistic = c("loading",
                                                     "squared_loading",
                                                     "subspace_norm"),
                                       orientation = c("auto", "signed",
                                                       "unsigned",
                                                       "subspace_norm"),
                                       R = NULL,
                                       seed = NULL,
                                       ...) {
      require_multivarious()
      if (!exists("bootstrap_pca", envir = asNamespace("multivarious"),
                  inherits = FALSE)) {
        stop("multivarious::bootstrap_pca() is required for PCA native feature evidence.",
             call. = FALSE)
      }
      if (!inherits(units, "multifer_units")) {
        stop("`units` must be a multifer_units table.", call. = FALSE)
      }
      statistic <- match.arg(statistic)
      orientation <- match.arg(orientation)
      members <- attr(units, "members")
      V_all <- stats::coef(x)
      available <- ncol(V_all)
      keep <- vapply(members, function(m) all(m <= available), logical(1L))
      if (!all(keep)) {
        units <- units[keep, , drop = FALSE]
        members <- members[keep]
        attr(units, "members") <- members
      }
      max_member <- max(unlist(members, use.names = FALSE), 0L)
      if (max_member < 1L) {
        return(infer_feature_evidence())
      }
      nboot <- if (is.null(R)) 500L else as.integer(R)
      boot <- multivarious::bootstrap_pca(
        x,
        nboot = nboot,
        k = max_member,
        seed = seed,
        ...
      )
      V <- V_all[, seq_len(max_member), drop = FALSE]
      reps <- lapply(seq_len(boot$nboot), function(b) {
        list(aligned_loadings = list(X = V %*% boot$Ab_array[, , b]))
      })
      artifact <- structure(
        list(
          reps = reps,
          R = boot$nboot,
          method_align = "subspace",
          domains = "X",
          seed = seed
        ),
        class = "multifer_bootstrap_artifact"
      )
      adapter_stub <- structure(
        list(loadings = function(fit, domain = NULL, ...) stats::coef(fit)),
        class = "multifer_adapter"
      )
      out <- feature_evidence_from_bootstrap(
        artifact = artifact,
        units = units,
        original_fit = x,
        adapter = adapter_stub,
        statistic = statistic,
        scope = "unit",
        k = "all",
        orientation = orientation,
        normalize = "none"
      )
      out$method <- "conditional_subspace_bootstrap"
      out$validity_level <- "conditional"
      out$calibration <- "bootstrap"
      out$p_value <- NA_real_
      out$p_adjusted <- NA_real_
      out
    },

    validity_level       = "conditional",
    declared_assumptions = c("rows_exchangeable"),
    checked_assumptions  = .oneblock_baser_checks()
  )
}


#' Register the multivarious::pca adapter
#'
#' Called from the package \code{.onLoad} hook.
#'
#' @keywords internal
register_multivarious_pca_adapter <- function() {
  register_infer_adapter(
    adapter_id = "multivarious_pca",
    adapter    = adapter_multivarious_pca(),
    overwrite  = TRUE
  )
  invisible(NULL)
}
