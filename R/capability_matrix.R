#' Valid capability targets
#'
#' The set of inferential targets a multifer adapter may claim to
#' support. Matches Part 5 section 38 and the frozen `infer_result`
#' sub-block names.
#'
#' @return Character vector of valid target names.
#' @export
valid_capability_targets <- function() {
  c("component_significance",
    "variable_stability",
    "score_stability",
    "subspace_stability",
    "variable_significance")
}

#' Build a capability matrix
#'
#' A machine-readable grid of `(geometry, relation, target) ->
#' supported?` that every adapter declares. The two places it is
#' consulted are:
#'
#' 1. **Adapter registration.** An adapter cannot *claim* a
#'    `(geometry, relation, target)` triple unless it provides the
#'    hooks needed to compute it (Part 5 section 36, point 1). This
#'    check is enforced by the adapter constructor in a later phase;
#'    here we only define the data structure.
#' 2. **Recipe compilation.** `infer_recipe()` refuses to compile a
#'    recipe whose requested target is not declared supported
#'    (Part 5 section 36, point 2).
#'
#' @param ... Zero or more entry lists, each with fields `geometry`
#'   (character scalar), `relation` (character scalar), and `targets`
#'   (character vector of valid target names, see
#'   [valid_capability_targets()]).
#'
#' @return An object of class `multifer_capability_matrix`. Internally
#'   a data.frame with columns `geometry`, `relation`, `target`.
#'
#' @details
#' The constructor normalizes the input into a long-format data.frame
#' where each row is one supported `(geometry, relation, target)`
#' triple. Geometry and relation values are validated against the
#' typed shape vocabulary: unknown values are rejected.
#'
#' @examples
#' cm <- capability_matrix(
#'   list(geometry = "oneblock", relation = "variance",
#'        targets = c("component_significance", "variable_stability")),
#'   list(geometry = "cross", relation = "covariance",
#'        targets = c("component_significance", "subspace_stability")),
#'   list(geometry = "cross", relation = "correlation",
#'        targets = "component_significance")
#' )
#' supports(cm, "cross", "correlation", "component_significance")
#' capabilities_for(cm, "cross", "covariance")
#'
#' @export
capability_matrix <- function(...) {
  entries <- list(...)

  valid_geom <- c("oneblock", "cross", "multiblock", "geneig", "adapter")
  valid_rel  <- c("variance", "covariance", "correlation",
                  "generalized_eigen", "predictive")
  valid_tgt  <- valid_capability_targets()

  rows_geom <- character(0)
  rows_rel  <- character(0)
  rows_tgt  <- character(0)

  for (i in seq_along(entries)) {
    e <- entries[[i]]
    if (!is.list(e) || is.null(e$geometry) || is.null(e$relation) ||
        is.null(e$targets)) {
      stop(sprintf(
        "entry %d: must be a list with fields `geometry`, `relation`, `targets`.",
        i
      ), call. = FALSE)
    }
    if (!(is.character(e$geometry) && length(e$geometry) == 1L &&
          e$geometry %in% valid_geom)) {
      stop(sprintf(
        "entry %d: `geometry` must be one of %s.",
        i, paste(valid_geom, collapse = ", ")
      ), call. = FALSE)
    }
    if (!(is.character(e$relation) && length(e$relation) == 1L &&
          e$relation %in% valid_rel)) {
      stop(sprintf(
        "entry %d: `relation` must be one of %s.",
        i, paste(valid_rel, collapse = ", ")
      ), call. = FALSE)
    }
    if (identical(e$relation, "predictive") &&
        !(e$geometry %in% c("cross", "adapter"))) {
      stop(sprintf(
        paste0("entry %d: relation 'predictive' requires geometry ",
               "'cross' or 'adapter'; ",
               "got geometry '%s'."),
        i, e$geometry
      ), call. = FALSE)
    }
    if (!is.character(e$targets) || length(e$targets) == 0L) {
      stop(sprintf(
        "entry %d: `targets` must be a non-empty character vector.",
        i
      ), call. = FALSE)
    }
    bad <- setdiff(e$targets, valid_tgt)
    if (length(bad) > 0L) {
      stop(sprintf(
        "entry %d: unknown target(s): %s. Valid: %s.",
        i, paste(bad, collapse = ", "), paste(valid_tgt, collapse = ", ")
      ), call. = FALSE)
    }

    t <- unique(e$targets)
    rows_geom <- c(rows_geom, rep(e$geometry, length(t)))
    rows_rel  <- c(rows_rel,  rep(e$relation, length(t)))
    rows_tgt  <- c(rows_tgt,  t)
  }

  df <- data.frame(
    geometry = rows_geom,
    relation = rows_rel,
    target   = rows_tgt,
    stringsAsFactors = FALSE
  )
  df <- unique(df)

  structure(df, class = c("multifer_capability_matrix", "data.frame"))
}

#' Test whether an object is a capability matrix
#' @param x Object to test.
#' @export
is_capability_matrix <- function(x) {
  inherits(x, "multifer_capability_matrix")
}

#' Query whether a capability matrix supports a triple
#'
#' @param cm A `multifer_capability_matrix`.
#' @param geometry Character scalar.
#' @param relation Character scalar.
#' @param target Character scalar in [valid_capability_targets()].
#'
#' @return `TRUE` if the triple is declared supported, `FALSE` otherwise.
#'
#' @export
supports <- function(cm, geometry, relation, target) {
  if (!is_capability_matrix(cm)) {
    stop("`cm` must be a multifer_capability_matrix.", call. = FALSE)
  }
  if (!(is.character(geometry) && length(geometry) == 1L &&
        is.character(relation) && length(relation) == 1L &&
        is.character(target) && length(target) == 1L)) {
    stop("`geometry`, `relation`, `target` must each be a single string.",
         call. = FALSE)
  }
  any(cm$geometry == geometry &
      cm$relation == relation &
      cm$target   == target)
}

#' List all supported targets for a (geometry, relation) pair
#'
#' @param cm A `multifer_capability_matrix`.
#' @param geometry Character scalar.
#' @param relation Character scalar.
#'
#' @return Character vector of target names (possibly empty).
#' @export
capabilities_for <- function(cm, geometry, relation) {
  if (!is_capability_matrix(cm)) {
    stop("`cm` must be a multifer_capability_matrix.", call. = FALSE)
  }
  mask <- cm$geometry == geometry & cm$relation == relation
  unique(cm$target[mask])
}

#' @export
print.multifer_capability_matrix <- function(x, ...) {
  cat("<multifer_capability_matrix: ", nrow(x), " entries>\n", sep = "")
  if (nrow(x) == 0L) {
    cat("  (empty)\n")
    return(invisible(x))
  }
  pairs <- unique(x[, c("geometry", "relation")])
  for (i in seq_len(nrow(pairs))) {
    g <- pairs$geometry[i]
    r <- pairs$relation[i]
    tgts <- capabilities_for(x, g, r)
    cat("  ", g, " + ", r, ":\n", sep = "")
    for (t in tgts) cat("    - ", t, "\n", sep = "")
  }
  invisible(x)
}
