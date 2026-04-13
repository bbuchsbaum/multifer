#' Relation of a latent operator
#'
#' The second component of a typed shape. Declares which quantity the
#' latent roots represent. This is the piece that disambiguates a `cross`
#' geometry into PLS-like (`covariance`) vs. CCA-like (`correlation`)
#' inference, which is the most common silent-failure mode that strict
#' dispatch is designed to prevent.
#'
#' @param kind One of `"variance"`, `"covariance"`, `"correlation"`,
#'   `"generalized_eigen"`.
#'
#' @return An object of class `multifer_relation` (and a `kind`-specific
#'   subclass).
#'
#' @export
relation <- function(kind) {
  valid <- c("variance", "covariance", "correlation", "generalized_eigen")
  if (!is.character(kind) || length(kind) != 1L || is.na(kind)) {
    stop("`kind` must be a single non-NA character string.", call. = FALSE)
  }
  if (!(kind %in% valid)) {
    stop(sprintf(
      "Unknown relation `%s`. Must be one of: %s.",
      kind, paste(valid, collapse = ", ")
    ), call. = FALSE)
  }
  structure(
    list(kind = kind),
    class = c(paste0("multifer_relation_", kind), "multifer_relation")
  )
}

#' @export
print.multifer_relation <- function(x, ...) {
  cat("<multifer_relation: ", x$kind, ">\n", sep = "")
  invisible(x)
}

#' Test whether an object is a multifer relation
#' @param x Object to test.
#' @export
is_relation <- function(x) inherits(x, "multifer_relation")
