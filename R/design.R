#' Designs: the third component of a typed shape
#'
#' A design object encodes the exchangeability / nuisance structure of
#' the data. It determines which null action is valid (permutation,
#' sign-flip, block permutation, nuisance-aware residualization) and
#' therefore feeds directly into recipe compilation and validity
#' checking.
#'
#' The v1 vocabulary is deliberately small. Richer designs (longitudinal,
#' hierarchical, structural) belong in later phases and should extend
#' these constructors rather than replace them.
#'
#' @name design
NULL

new_design <- function(kind, ...) {
  structure(
    c(list(kind = kind), list(...)),
    class = c(paste0("multifer_design_", kind), "multifer_design")
  )
}

#' Exchangeable rows
#'
#' Rows are iid-exchangeable. The default design for one-block
#' decompositions without grouping structure.
#'
#' @return A `multifer_design` object.
#' @export
exchangeable_rows <- function() {
  new_design("exchangeable_rows")
}

#' Paired rows
#'
#' Rows of two blocks are paired (each row in block X corresponds to the
#' same observation in block Y). The default design for cross-block
#' inference. The valid null action breaks the pairing in exactly one
#' block.
#'
#' @return A `multifer_design` object.
#' @export
paired_rows <- function() {
  new_design("paired_rows")
}

#' Blocked rows
#'
#' Rows come in exchangeable blocks. Permutation is valid only *within*
#' each block, and sign-flips may be applied block-wise.
#'
#' @param groups A factor, integer, or character vector of length `n`
#'   giving the block label of each row.
#'
#' @return A `multifer_design` object.
#' @export
blocked_rows <- function(groups) {
  if (is.null(groups)) {
    stop("`groups` must not be NULL.", call. = FALSE)
  }
  if (!(is.factor(groups) || is.integer(groups) ||
        is.character(groups) || is.numeric(groups))) {
    stop("`groups` must be a factor, integer, numeric, or character vector.",
         call. = FALSE)
  }
  if (length(groups) == 0L) {
    stop("`groups` must have length >= 1.", call. = FALSE)
  }
  if (anyNA(groups)) {
    stop("`groups` must not contain NA.", call. = FALSE)
  }
  new_design("blocked_rows", groups = as.factor(groups))
}

#' Nuisance-adjusted design
#'
#' Marks that inference is conditional on a nuisance covariate matrix
#' `Z`. The recipe compiler is expected to use an exchangeability-
#' preserving basis transform (Winkler-style) before applying any null
#' action. The design object itself does not perform the transform — it
#' only records that one is required.
#'
#' @param Z A numeric matrix of nuisance covariates, with one row per
#'   observation.
#'
#' @return A `multifer_design` object.
#' @export
nuisance_adjusted <- function(Z) {
  if (is.null(Z)) {
    stop("`Z` must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix.", call. = FALSE)
  }
  if (anyNA(Z)) {
    stop("`Z` must not contain NA.", call. = FALSE)
  }
  new_design("nuisance_adjusted", Z = Z)
}

#' @export
print.multifer_design <- function(x, ...) {
  cat("<multifer_design: ", x$kind, sep = "")
  if (x$kind == "blocked_rows") {
    cat(" | ", length(levels(x$groups)), " blocks, n = ",
        length(x$groups), sep = "")
  } else if (x$kind == "nuisance_adjusted") {
    cat(" | Z: ", nrow(x$Z), "x", ncol(x$Z), sep = "")
  }
  cat(">\n")
  invisible(x)
}

#' Test whether an object is a multifer design
#' @param x Object to test.
#' @export
is_design <- function(x) inherits(x, "multifer_design")
