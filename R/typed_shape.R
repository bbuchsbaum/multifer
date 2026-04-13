#' Typed shape — the central inferential object
#'
#' A typed shape is the triple `(geometry, relation, design)`. It is the
#' minimum information needed to decide whether an inference procedure
#' is valid on a given model and dataset, and to compile a recipe that
#' selects the appropriate null action, test statistic, and adapter
#' hooks.
#'
#' The purpose of making this a typed object (rather than a loose list
#' of options) is to prevent silent ambiguity: a `cross` geometry with
#' no declared relation could equally mean PLS-like or CCA-like
#' inference, and those have different nulls. Typed shapes force the
#' caller to be explicit, and `strict = TRUE` recipe compilation errors
#' on unspecified relations.
#'
#' @param geometry A `multifer_geometry` object (see [geometry()]).
#' @param relation A `multifer_relation` object (see [relation()]).
#' @param design   A `multifer_design` object (see [exchangeable_rows()],
#'   [paired_rows()], [blocked_rows()], [nuisance_adjusted()]).
#'
#' @return An object of class `multifer_typed_shape`.
#'
#' @details
#' This constructor performs **structural** validation only: it checks
#' that each slot is of the right S3 class. Compatibility between
#' geometry and relation (e.g. `cross` + `variance` being unusual) is
#' **not** checked here — that belongs in recipe compilation where the
#' adapter capability matrix is available.
#'
#' @export
typed_shape <- function(geometry, relation, design) {
  if (!is_geometry(geometry)) {
    stop("`geometry` must be a multifer_geometry (see geometry()).",
         call. = FALSE)
  }
  if (!is_relation(relation)) {
    stop("`relation` must be a multifer_relation (see relation()).",
         call. = FALSE)
  }
  if (!is_design(design)) {
    stop("`design` must be a multifer_design (see exchangeable_rows(), etc.).",
         call. = FALSE)
  }
  structure(
    list(geometry = geometry, relation = relation, design = design),
    class = "multifer_typed_shape"
  )
}

#' @export
print.multifer_typed_shape <- function(x, ...) {
  cat("<multifer_typed_shape>\n",
      "  geometry: ", x$geometry$kind, "\n",
      "  relation: ", x$relation$kind, "\n",
      "  design:   ", x$design$kind, "\n",
      sep = "")
  invisible(x)
}

#' Test whether an object is a multifer typed shape
#' @param x Object to test.
#' @export
is_typed_shape <- function(x) inherits(x, "multifer_typed_shape")
