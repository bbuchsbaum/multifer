#' Infer PCA-family significance and stability
#'
#' Thin convenience wrapper over [infer()] for one-block variance models.
#' The default adapter is `"multivarious_pca"`, which wraps
#' `multivarious::pca()` and is the recommended fitting path for
#' PCA-family inference in `multifer`. The base-R adapters
#' `"svd_oneblock"` and `"prcomp_oneblock"` remain available as reference
#' examples and as fall-back fitting routines when `multivarious` is not
#' convenient.
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
#' Thin convenience wrapper over [infer()] for two-block covariance models.
#' The default adapter is `"multivarious_plsc"`, which wraps
#' `multivarious::plsc()` and is the recommended fitting path for
#' PLSC-family inference. The base-R adapter `"cross_svd"` remains
#' available as a reference example and fall-back fitting routine.
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

#' Infer CCA-family significance and stability
#'
#' Thin convenience wrapper over [infer()] for two-block correlation
#' models. Current validity is strongest for the paired-row design and
#' supported nuisance-adjusted variants described in the package README
#' and vignettes.
#'
#' CCA inference currently uses the `"cross_svd"` reference adapter in
#' correlation mode: it ships the QR-whitened stepwise deflation and
#' the design-aware null actions (`paired_rows`, `nuisance_adjusted`,
#' `blocked_rows`). A `multivarious`-backed CCA adapter is a natural
#' extension point but has not been written yet; see
#' `notes/package_vision.md` for the roadmap.
#'
#' @param X Numeric matrix for the first block.
#' @param Y Numeric matrix for the second block. Must have the same
#'   number of rows as `X`.
#' @param adapter Adapter id or object. Defaults to `"cross_svd"`.
#' @param ... Additional arguments forwarded to [infer()].
#'
#' @return An [infer_result].
#' @export
infer_cca <- function(X,
                      Y,
                      adapter = "cross_svd",
                      ...) {
  infer(
    adapter = adapter,
    data = list(X = X, Y = Y),
    geometry = "cross",
    relation = "correlation",
    ...
  )
}
