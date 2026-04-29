#' Infer PCA-family significance and stability
#'
#' @description
#' **Maturity: mature.** Paper-backed exact path for one-block variance
#' inference, defended by executable calibration and parity evidence.
#'
#' Thin convenience wrapper over [infer()] for one-block variance models.
#' The default adapter is `"multivarious_pca"`. Base-R reference adapters
#' such as `"prcomp_oneblock"` and `"svd_oneblock"` remain available, but
#' callers must request them explicitly with `adapter = `.
#'
#' If you have already fit a model via `multivarious::pca()` (or any
#' other compatible fitter), pass it through `model = ` and `infer_pca()`
#' will reuse your fit directly instead of refitting the original data.
#'
#' @param X Numeric matrix with observations in rows and variables in columns.
#' @param adapter Adapter id or object. Defaults to `"multivarious_pca"`.
#' @param ... Additional arguments forwarded to [infer()]. In particular
#'   `model = ` accepts a pre-fit projector (e.g. `multivarious::pca(X)`)
#'   so you can fit once and re-infer.
#'
#' @return An [infer_result].
#' @export
infer_pca <- function(X,
                      adapter = "multivarious_pca",
                      ...) {
  infer(
    adapter = adapter,
    data = X,
    geometry = "oneblock",
    relation = "variance",
    ...
  )
}

#' Infer PLSC-family significance and stability
#'
#' @description
#' **Maturity: mature.** Paper-backed exact path for two-block
#' covariance inference, including an exact core-space bootstrap
#' identity.
#'
#' Thin convenience wrapper over [infer()] for two-block covariance models.
#' The default adapter is `"multivarious_plsc"`. The `"cross_svd"`
#' reference adapter remains available, but callers must request it
#' explicitly with `adapter = "cross_svd"`.
#'
#' If you have already fit a model via `multivarious::plsc()`, pass it
#' through `model = ` and `infer_plsc()` will reuse your fit directly.
#'
#' @param X Numeric matrix for the first block.
#' @param Y Numeric matrix for the second block. Must have the same
#'   number of rows as `X`.
#' @param adapter Adapter id or object. Defaults to `"multivarious_plsc"`.
#' @param ... Additional arguments forwarded to [infer()]. In particular
#'   `model = ` accepts a pre-fit projector (e.g. `multivarious::plsc(X, Y)`).
#'
#' @return An [infer_result].
#' @export
infer_plsc <- function(X,
                       Y,
                       adapter = "multivarious_plsc",
                       ...) {
  infer(
    adapter = adapter,
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "covariance",
    ...
  )
}

#' Infer PLSR-family significance and stability
#'
#' @description
#' **Maturity: narrow.** Shipped today via `plsr_refit` only. The
#' inferential target is cross-fitted held-out predictive gain, not a
#' cross-covariance root. Broader predictive-cross families (reduced
#' rank regression, canonical ridge, kernel variants) are not part of
#' the v1 public wrapper surface.
#'
#' Thin convenience wrapper over [infer()] for two-block predictive
#' models. The inferential target is **cross-fitted held-out predictive
#' gain**, not an `X^T Y` cross-covariance root. In other words,
#' `infer_plsr()` asks whether the next latent predictor improves
#' out-of-sample prediction of `Y` beyond the previous ones.
#'
#' The default adapter is `"plsr_refit"`, a Tier-1 wrapper around
#' `pls::plsr()`. For the predictive-relation doctrine and the
#' registration-time cross-fit gate, see `notes/plsr_maturity_contract.md`
#' and
#' `vignette("writing-adapters")`, especially the section
#' "The predictive relation is stricter on purpose".
#'
#' The public predictive surface is intentionally narrow in v1:
#' `infer_plsr()` is the shipped path, while broader predictive-cross
#' families remain future work.
#'
#' As elsewhere in v1, `variable_significance` remains out of scope.
#' Component p-values are significance outputs; variable / score /
#' subspace summaries remain bootstrap stability outputs.
#'
#' @param X Numeric matrix for the predictor block.
#' @param Y Numeric matrix for the response block. Must have the same
#'   number of rows as `X`.
#' @param adapter Adapter id or object. Defaults to `"plsr_refit"`.
#' @param ... Additional arguments forwarded to [infer()].
#'
#' @return An [infer_result].
#' @export
infer_plsr <- function(X,
                       Y,
                       adapter = "plsr_refit",
                       ...) {
  dots <- list(...)

  if (!is.null(dots$relation) && !identical(dots$relation, "predictive")) {
    stop(
      paste0(
        "`infer_plsr()` uses relation = \"predictive\". ",
        "Do not request covariance-style inference through this wrapper."
      ),
      call. = FALSE
    )
  }
  if (!is.null(dots$geometry) && !identical(dots$geometry, "cross")) {
    stop(
      "`infer_plsr()` uses geometry = \"cross\" and does not accept another geometry.",
      call. = FALSE
    )
  }

  dots$relation <- NULL
  dots$geometry <- NULL

  do.call(
    infer,
    c(
      list(
        adapter = adapter,
        data = list(X = X, Y = Y),
        geometry = "cross",
        relation = "predictive"
      ),
      dots
    )
  )
}

#' Infer CCA-family significance and stability
#'
#' @description
#' **Maturity: mature** on the shipped support matrix (paired-row,
#' common-`Z` nuisance-adjusted, nuisance-adjusted within-block, and
#' `blocked_rows(groups)` designs); first-root-only outside that
#' matrix.
#'
#' Thin convenience wrapper over [infer()] for two-block correlation
#' models. The shipped multi-root support boundary covers the paired-row
#' design and the supported nuisance-adjusted variants described in the
#' package README and vignettes. Outside that support matrix, `multifer`
#' conservatively caps the ladder to first-root inference rather than
#' making unsupported stepwise claims.
#'
#' The default adapter is `"cancor_cross"`, the dedicated
#' correlation-only wrapper around [stats::cancor()]. This keeps the
#' default CCA path single-relation and avoids the strict-dispatch
#' ambiguity that is intentional in the dual-relation `"cross_svd"`
#' reference adapter.
#'
#' `infer_cca()` still dispatches into the same correlation engine:
#' QR-whitened stepwise deflation for paired rows and supported
#' nuisance-adjusted designs, with conservative first-root capping for
#' richer structured designs. Use `adapter = "cross_svd"` when you
#' explicitly want the dual covariance/correlation reference adapter.
#'
#' @param X Numeric matrix for the first block.
#' @param Y Numeric matrix for the second block. Must have the same
#'   number of rows as `X`.
#' @param adapter Adapter id or object. Defaults to `"cancor_cross"`.
#' @param ... Additional arguments forwarded to [infer()].
#'
#' @return An [infer_result].
#' @export
infer_cca <- function(X,
                      Y,
                      adapter = "cancor_cross",
                      ...) {
  infer(
    adapter = adapter,
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "correlation",
    ...
  )
}

#' Infer discriminant-root significance for LDA-family models
#'
#' @description
#' **Maturity: narrow.** Shipped today for a deliberately limited
#' public surface: only discriminant-root significance for LDA via
#' label permutation and B-metric deflation. Broader metric-weighted
#' or contrastive generalized-eigen families are not part of the v1
#' public wrapper surface.
#'
#' Thin convenience wrapper over [infer()] for `(geneig,
#' generalized_eigen)` models driven by class labels. The default adapter
#' is `"lda_refit"`, which wraps `MASS::lda()` and tests the ordered
#' discriminant roots, not specific variables and not pairwise class
#' contrasts.
#'
#' The estimand is the number of discriminant roots that rise above
#' noise. The ladder is therefore naturally short: with `K` classes,
#' there are at most `K - 1` non-zero discriminant roots to test.
#' The MASS-backed v1 adapter also requires full within-class column
#' rank before fitting; p > n - K and collinear within-class designs are
#' rejected as validity failures rather than stabilized silently.
#'
#' For the current engine doctrine and the B-metric deflation rule, see
#' `notes/engine_geneig_spec.md`.
#'
#' The public geneig surface is intentionally narrow in v1:
#' `infer_lda()` exposes discriminant-root significance, while broader
#' generalized-eigen bootstrap/stability workflows remain outside the
#' wrapper-level contract for now.
#'
#' @param X Numeric matrix with observations in rows and variables in columns.
#' @param labels Factor of class labels with one entry per row of `X`.
#' @param targets Requested inferential targets. Defaults to
#'   `"component_significance"` because the geneig bootstrap/stability
#'   path is not yet part of the public wrapper surface.
#' @param adapter Adapter id or object. Defaults to `"lda_refit"`.
#' @param ... Additional arguments forwarded to [infer()].
#'
#' @return An [infer_result].
#' @export
infer_lda <- function(X,
                      labels,
                      targets = "component_significance",
                      adapter = "lda_refit",
                      ...) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.factor(labels)) {
    stop("`labels` must be a factor.", call. = FALSE)
  }
  if (length(labels) != nrow(X)) {
    stop("`labels` must have one value per row of `X`.", call. = FALSE)
  }

  infer(
    adapter = adapter,
    data = list(X = X, y = labels),
    geometry = "geneig",
    relation = "generalized_eigen",
    targets = targets,
    ...
  )
}
