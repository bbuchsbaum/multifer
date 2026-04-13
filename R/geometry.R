#' Geometry of a latent operator
#'
#' The first component of a typed shape. Declares the geometric family
#' of the latent operator being inferred on, independent of the relation
#' (variance / covariance / correlation / generalized eigenvalue) and
#' independent of the experimental design.
#'
#' @param kind One of `"oneblock"`, `"cross"`, `"multiblock"`, `"geneig"`.
#'
#' @return An object of class `multifer_geometry` (and a `kind`-specific
#'   subclass).
#'
#' @details
#' - `oneblock` — a single data block, e.g. PCA, kernel PCA.
#' - `cross` — two paired blocks, e.g. PLS, CCA, reduced-rank regression.
#' - `multiblock` — three or more blocks with shared structure,
#'   e.g. JIVE, AJIVE, multiblock projectors.
#' - `geneig` — generalized-eigen / contrastive / metric-weighted fits,
#'   e.g. cPCA++, discriminant analysis.
#'
#' @export
geometry <- function(kind) {
  valid <- c("oneblock", "cross", "multiblock", "geneig")
  if (!is.character(kind) || length(kind) != 1L || is.na(kind)) {
    stop("`kind` must be a single non-NA character string.", call. = FALSE)
  }
  if (!(kind %in% valid)) {
    stop(sprintf(
      "Unknown geometry `%s`. Must be one of: %s.",
      kind, paste(valid, collapse = ", ")
    ), call. = FALSE)
  }
  structure(
    list(kind = kind),
    class = c(paste0("multifer_geometry_", kind), "multifer_geometry")
  )
}

#' @export
print.multifer_geometry <- function(x, ...) {
  cat("<multifer_geometry: ", x$kind, ">\n", sep = "")
  invisible(x)
}

#' Test whether an object is a multifer geometry
#' @param x Object to test.
#' @export
is_geometry <- function(x) inherits(x, "multifer_geometry")
