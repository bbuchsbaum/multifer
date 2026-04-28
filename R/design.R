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

#' Clustered rows
#'
#' Rows belong to sampling units such as subjects. Bootstrap resampling
#' should sample whole clusters with replacement and carry all rows from
#' each sampled cluster. Optional strata restrict cluster sampling within
#' fixed strata.
#'
#' @param clusters A factor, integer, numeric, or character vector of
#'   length `n` giving the cluster label of each row.
#' @param strata Optional factor, integer, numeric, or character vector
#'   of length `n` giving the stratum label of each row. Each cluster
#'   must belong to exactly one stratum.
#'
#' @return A `multifer_design` object.
#' @export
clustered_rows <- function(clusters, strata = NULL) {
  if (is.null(clusters)) {
    stop("`clusters` must not be NULL.", call. = FALSE)
  }
  if (!(is.factor(clusters) || is.integer(clusters) ||
        is.character(clusters) || is.numeric(clusters))) {
    stop("`clusters` must be a factor, integer, numeric, or character vector.",
         call. = FALSE)
  }
  if (length(clusters) == 0L) {
    stop("`clusters` must have length >= 1.", call. = FALSE)
  }
  if (anyNA(clusters)) {
    stop("`clusters` must not contain NA.", call. = FALSE)
  }

  clusters <- as.factor(clusters)

  if (!is.null(strata)) {
    if (!(is.factor(strata) || is.integer(strata) ||
          is.character(strata) || is.numeric(strata))) {
      stop("`strata` must be a factor, integer, numeric, or character vector.",
           call. = FALSE)
    }
    if (length(strata) != length(clusters)) {
      stop("`strata` must have the same length as `clusters`.", call. = FALSE)
    }
    if (anyNA(strata)) {
      stop("`strata` must not contain NA.", call. = FALSE)
    }
    strata <- as.factor(strata)

    by_cluster <- split(strata, clusters)
    spans <- vapply(by_cluster, function(x) length(unique(x)) > 1L, logical(1L))
    if (any(spans)) {
      stop("Each cluster must belong to exactly one stratum.", call. = FALSE)
    }
  }

  new_design("clustered_rows", clusters = clusters, strata = strata)
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
#' @param groups Optional factor, integer, numeric, or character vector
#'   of length `nrow(Z)` giving exchangeability blocks for structured
#'   nuisance-adjusted designs. When supplied, row permutations are
#'   restricted within blocks after the residual-basis transform.
#'
#' @return A `multifer_design` object.
#' @export
nuisance_adjusted <- function(Z, groups = NULL) {
  if (is.null(Z)) {
    stop("`Z` must not be NULL.", call. = FALSE)
  }
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("`Z` must be a numeric matrix.", call. = FALSE)
  }
  if (anyNA(Z)) {
    stop("`Z` must not contain NA.", call. = FALSE)
  }
  if (!is.null(groups)) {
    if (!(is.factor(groups) || is.integer(groups) ||
          is.character(groups) || is.numeric(groups))) {
      stop("`groups` must be a factor, integer, numeric, or character vector.",
           call. = FALSE)
    }
    if (length(groups) != nrow(Z)) {
      stop("`groups` must have length nrow(Z).", call. = FALSE)
    }
    if (anyNA(groups)) {
      stop("`groups` must not contain NA.", call. = FALSE)
    }
    groups <- as.factor(groups)
  }
  new_design("nuisance_adjusted", Z = Z, groups = groups)
}

#' @export
print.multifer_design <- function(x, ...) {
  cat("<multifer_design: ", x$kind, sep = "")
  if (x$kind == "blocked_rows") {
    cat(" | ", length(levels(x$groups)), " blocks, n = ",
        length(x$groups), sep = "")
  } else if (x$kind == "clustered_rows") {
    cat(" | ", length(levels(x$clusters)), " clusters, n = ",
        length(x$clusters), sep = "")
    if (!is.null(x$strata)) {
      cat(" | ", length(levels(x$strata)), " strata", sep = "")
    }
  } else if (x$kind == "nuisance_adjusted") {
    cat(" | Z: ", nrow(x$Z), "x", ncol(x$Z), sep = "")
    if (!is.null(x$groups)) {
      cat(" | ", length(levels(x$groups)), " blocks", sep = "")
    }
  }
  cat(">\n")
  invisible(x)
}

#' Test whether an object is a multifer design
#' @param x Object to test.
#' @export
is_design <- function(x) inherits(x, "multifer_design")
