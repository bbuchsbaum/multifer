#' Infer PCA-family significance and stability
#'
#' Thin convenience wrapper over [infer()] for one-block variance models.
#'
#' @param X Numeric matrix with observations in rows and variables in columns.
#' @param adapter Adapter id or object. Defaults to `"prcomp_oneblock"`.
#' @param ... Additional arguments forwarded to [infer()].
#'
#' @return An [infer_result].
#' @export
infer_pca <- function(X,
                      adapter = "prcomp_oneblock",
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
#'
#' @param X Numeric matrix for the first block.
#' @param Y Numeric matrix for the second block. Must have the same number of rows as `X`.
#' @param adapter Adapter id or object. Defaults to `"cross_svd"`.
#' @param ... Additional arguments forwarded to [infer()].
#'
#' @return An [infer_result].
#' @export
infer_plsc <- function(X,
                       Y,
                       adapter = "cross_svd",
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
#' Thin convenience wrapper over [infer()] for two-block correlation models.
#' Current validity is strongest for the paired-row design and supported
#' nuisance-adjusted variants described in the package README and vignettes.
#'
#' @param X Numeric matrix for the first block.
#' @param Y Numeric matrix for the second block. Must have the same number of rows as `X`.
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
