#' Compile a typed shape against an adapter's capability matrix
#'
#' `infer_recipe()` is the compile-time validity gate (Part 5 section 36,
#' point 2). It takes a typed shape -- the `(geometry, relation, design)`
#' triple -- and a registered adapter, checks that every requested target
#' is declared supported for that `(geometry, relation)` pair, and returns
#' a compiled recipe object that downstream inference can execute without
#' repeating these checks.
#'
#' @section Strict dispatch (Part 5 section 35):
#' With `strict = TRUE` (default) the constructor errors on silent
#' ambiguity: if `relation` is `NULL` and the adapter supports multiple
#' relations for the requested `geometry`, the caller must disambiguate.
#' With `strict = FALSE` the ambiguity is resolved automatically (first
#' declared relation wins), a `warning()` is emitted, and the recipe's
#' `validity_level` is downgraded one step.
#'
#' @param shape A `multifer_typed_shape` (see [typed_shape()]).  When
#'   supplied, the individual `geometry`, `relation`, and `design`
#'   arguments are ignored.
#' @param geometry Character scalar -- one of `"oneblock"`, `"cross"`,
#'   `"multiblock"`, `"geneig"`, `"adapter"`.  Used when `shape` is `NULL`.
#' @param relation Character scalar -- one of `"variance"`,
#'   `"covariance"`, `"correlation"`, `"generalized_eigen"`,
#'   `"predictive"`.  Used when `shape` is `NULL`.  May be `NULL` if the
#'   adapter has exactly one relation for the given geometry (strict mode
#'   errors if it has more).
#' @param design A `multifer_design` object (see [exchangeable_rows()],
#'   [paired_rows()], [blocked_rows()], [nuisance_adjusted()]).  Used
#'   when `shape` is `NULL`.  If `NULL`, a sensible default is chosen
#'   based on `geometry`: `"oneblock"`, `"multiblock"`, and `"geneig"`
#'   get [exchangeable_rows()]; `"cross"` gets [paired_rows()].
#' @param targets Character vector of requested inferential targets, or
#'   the special sentinel `"default"`.  The default set is every target
#'   the adapter declares for this `(geometry, relation)` pair, excluding
#'   `"variable_significance"` (deferred to Phase 3, Part 5 section 38).
#' @param adapter Either a `multifer_adapter` object or a character
#'   scalar naming a registered adapter (looked up via
#'   [get_infer_adapter()]).
#' @param strict Logical scalar.  `TRUE` (default) enforces strict
#'   dispatch (Part 5 section 35); `FALSE` relaxes the ambiguous-relation
#'   rule only -- unsupported targets and `variable_significance` still
#'   error unconditionally.
#'
#' @return An object of class `multifer_infer_recipe` with fields:
#'   `problem`, `shape`, `targets`, `adapter_id`, `adapter`,
#'   `validity_level`, `strict`, `downgrade_reason`,
#'   `declared_assumptions`, `checked_assumptions`, `call`. The
#'   `problem` field is the first-class compiled inferential object and
#'   is the preferred place to inspect the inferred shape, target
#'   family, engine kind, and validity metadata in one place.
#'
#' @seealso [infer_problem()]
#'
#' @export
infer_recipe <- function(shape    = NULL,
                         geometry = NULL,
                         relation = NULL,
                         design   = NULL,
                         targets  = "default",
                         adapter  = NULL,
                         strict   = TRUE) {

  the_call <- match.call()

  # -- resolve adapter ----------------------------------------------------------
  if (is.character(adapter)) {
    if (length(adapter) != 1L || is.na(adapter) || nchar(adapter) == 0L) {
      stop("`adapter` must be a non-empty single string or a multifer_adapter.",
           call. = FALSE)
    }
    adapter <- get_infer_adapter(adapter)
  }
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter object or a registered adapter id string.",
         call. = FALSE)
  }

  # -- resolve shape vs. individual components ----------------------------------
  if (!is.null(shape)) {
    if (!is_typed_shape(shape)) {
      stop("`shape` must be a multifer_typed_shape (see typed_shape()).",
           call. = FALSE)
    }
    # When shape is provided, individual components are ignored.
    geom_kind <- shape$geometry$kind
    rel_kind  <- shape$relation$kind
  } else {
    # geometry is required when shape is NULL
    if (is.null(geometry)) {
      stop("Either `shape` or `geometry` must be supplied.", call. = FALSE)
    }
    if (!is.character(geometry) || length(geometry) != 1L || is.na(geometry)) {
      stop("`geometry` must be a single non-NA string.", call. = FALSE)
    }
    valid_geom <- c("oneblock", "cross", "multiblock", "geneig", "adapter")
    if (!(geometry %in% valid_geom)) {
      stop(sprintf("`geometry` must be one of: %s.",
                   paste(valid_geom, collapse = ", ")),
           call. = FALSE)
    }
    geom_kind <- geometry
  }

  if (!is.null(shape) &&
      identical(rel_kind, "predictive") &&
      !(geom_kind %in% c("cross", "adapter"))) {
    stop(
      sprintf(
        "Relation 'predictive' requires geometry 'cross' or 'adapter' during recipe compilation; got '%s'.",
        geom_kind
      ),
      call. = FALSE
    )
  }
  if (!is.null(shape) &&
      identical(geom_kind, "geneig") &&
      !identical(rel_kind, "generalized_eigen")) {
    stop(
      "Geometry 'geneig' requires relation 'generalized_eigen' during recipe compilation.",
      call. = FALSE
    )
  }
  if (is.null(shape) &&
      !is.null(relation) &&
      identical(relation, "predictive") &&
      !(geom_kind %in% c("cross", "adapter"))) {
    stop(
      sprintf(
        "Relation 'predictive' requires geometry 'cross' or 'adapter' during recipe compilation; got '%s'.",
        geom_kind
      ),
      call. = FALSE
    )
  }
  if (is.null(shape) &&
      !is.null(relation) &&
      identical(geom_kind, "geneig") &&
      !identical(relation, "generalized_eigen")) {
    stop(
      "Geometry 'geneig' requires relation 'generalized_eigen' during recipe compilation.",
      call. = FALSE
    )
  }

  # -- strict dispatch rule 4: geometry must be in adapter$shape_kinds ---------
  if (!(geom_kind %in% adapter$shape_kinds)) {
    stop(sprintf(
      "Geometry '%s' is not in this adapter's shape_kinds (%s).",
      geom_kind, paste(adapter$shape_kinds, collapse = ", ")
    ), call. = FALSE)
  }

  # -- resolve relation (when shape is NULL) ------------------------------------
  downgrade_reason <- NA_character_
  validity_level   <- adapter$validity_level
  vl_order         <- c("exact", "conditional", "asymptotic", "heuristic")

  if (is.null(shape)) {
    cm <- adapter_capabilities(adapter)

    if (is.null(relation)) {
      # Find all relations the adapter declares for this geometry
      cm_geom   <- cm[cm$geometry == geom_kind, , drop = FALSE]
      avail_rel <- unique(cm_geom$relation)

      if (length(avail_rel) == 0L) {
        stop(sprintf(
          "Adapter '%s' declares no capabilities for geometry '%s'.",
          adapter$adapter_id, geom_kind
        ), call. = FALSE)
      }

      if (length(avail_rel) > 1L) {
        # Strict dispatch rule 1: ambiguous relation
        if (isTRUE(strict)) {
          stop(sprintf(
            paste0("Ambiguous relation: adapter '%s' supports multiple relations ",
                   "for geometry '%s': %s. ",
                   "Specify `relation` explicitly to avoid silent ambiguity ",
                   "(Part 5 section 35)."),
            adapter$adapter_id, geom_kind,
            paste(sort(avail_rel), collapse = ", ")
          ), call. = FALSE)
        } else {
          # Non-strict: auto-pick first, warn, downgrade
          picked_rel       <- avail_rel[1L]
          downgrade_reason <- sprintf(
            "Auto-picked relation '%s' from multiple candidates (%s) because strict = FALSE.",
            picked_rel, paste(sort(avail_rel), collapse = ", ")
          )
          warning(sprintf(
            paste0("Ambiguous relation auto-resolved to '%s' (strict = FALSE). ",
                   "Adapter '%s' supports: %s. ",
                   "validity_level downgraded. ",
                   "Set `relation` explicitly to suppress this warning."),
            picked_rel, adapter$adapter_id,
            paste(sort(avail_rel), collapse = ", ")
          ), call. = FALSE)
          rel_kind <- picked_rel
          # Downgrade validity_level one step
          pos <- match(validity_level, vl_order)
          if (!is.na(pos) && pos < length(vl_order)) {
            validity_level <- vl_order[pos + 1L]
          }
        }
      } else {
        rel_kind <- avail_rel[1L]
      }
    } else {
      # relation was supplied
      if (!is.character(relation) || length(relation) != 1L || is.na(relation)) {
        stop("`relation` must be a single non-NA string.", call. = FALSE)
      }
      valid_rel <- c("variance", "covariance", "correlation",
                     "generalized_eigen", "predictive")
      if (!(relation %in% valid_rel)) {
        stop(sprintf("`relation` must be one of: %s.",
                     paste(valid_rel, collapse = ", ")),
             call. = FALSE)
      }
      rel_kind <- relation
    }

    # default design based on geometry
    if (is.null(design)) {
      design <- if (geom_kind == "cross") paired_rows() else exchangeable_rows()
    }
    if (!is_design(design)) {
      stop("`design` must be a multifer_design (see exchangeable_rows(), etc.).",
           call. = FALSE)
    }

    shape <- typed_shape(
      geometry = geometry(geom_kind),
      relation = relation(rel_kind),
      design   = design
    )
  } else {
    rel_kind <- shape$relation$kind
  }

  if (identical(geom_kind, "geneig") && !identical(rel_kind, "generalized_eigen")) {
    stop(
      "Geometry 'geneig' requires relation 'generalized_eigen' during recipe compilation.",
      call. = FALSE
    )
  }
  if (identical(rel_kind, "predictive") &&
      !(geom_kind %in% c("cross", "adapter"))) {
    stop(
      sprintf(
        "Relation 'predictive' requires geometry 'cross' or 'adapter' during recipe compilation; got '%s'.",
        geom_kind
      ),
      call. = FALSE
    )
  }

  # -- resolve targets ----------------------------------------------------------
  cm <- adapter_capabilities(adapter)

  # strict dispatch rule 3: variable_significance always errors
  if (is.character(targets) && length(targets) > 0L &&
      !identical(targets, "default")) {
    if ("variable_significance" %in% targets) {
      stop(paste0(
        "Target 'variable_significance' is not permitted in v1. ",
        "Variable significance is deferred to Phase 3 (Part 5 section 38). ",
        "Remove 'variable_significance' from `targets`."
      ), call. = FALSE)
    }
  }

  if (identical(targets, "default") || (length(targets) == 1L && targets == "default")) {
    # Expand to all supported targets for this (geometry, relation), minus variable_significance
    all_supported <- capabilities_for(cm, geom_kind, rel_kind)
    resolved_targets <- setdiff(all_supported, "variable_significance")
    if (length(resolved_targets) == 0L) {
      stop(sprintf(
        "Adapter '%s' declares no supported targets for (%s, %s).",
        adapter$adapter_id, geom_kind, rel_kind
      ), call. = FALSE)
    }
  } else {
    if (!is.character(targets) || length(targets) == 0L) {
      stop("`targets` must be a non-empty character vector or \"default\".",
           call. = FALSE)
    }
    resolved_targets <- unique(targets)
  }

  # strict dispatch rule 2: target not supported errors regardless of strict
  supported_for_pair <- capabilities_for(cm, geom_kind, rel_kind)
  unsupported <- setdiff(resolved_targets, supported_for_pair)
  if (length(unsupported) > 0L) {
    stop(sprintf(
      paste0("Target(s) not supported by adapter '%s' for (%s, %s).\n",
             "  Requested:  %s\n",
             "  Supported:  %s"),
      adapter$adapter_id, geom_kind, rel_kind,
      paste(resolved_targets, collapse = ", "),
      if (length(supported_for_pair) > 0L)
        paste(supported_for_pair, collapse = ", ")
      else
        "(none)"
    ), call. = FALSE)
  }

  # -- assemble -----------------------------------------------------------------
  problem <- infer_problem(
    shape = shape,
    targets = resolved_targets,
    adapter_id = adapter$adapter_id,
    adapter_version = adapter$adapter_version,
    validity_level = validity_level,
    strict = isTRUE(strict),
    downgrade_reason = downgrade_reason,
    call = the_call
  )

  structure(
    list(
      problem              = problem,
      shape                = problem$shape,
      targets              = problem$targets,
      adapter_id           = problem$adapter_id,
      adapter              = adapter,
      validity_level       = problem$validity_level,
      strict               = problem$strict,
      downgrade_reason     = problem$downgrade_reason,
      declared_assumptions = adapter$declared_assumptions,
      checked_assumptions  = adapter$checked_assumptions,
      call                 = the_call
    ),
    class = "multifer_infer_recipe"
  )
}

#' Test whether an object is a compiled infer recipe
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_infer_recipe <- function(x) inherits(x, "multifer_infer_recipe")

#' @export
print.multifer_infer_recipe <- function(x, ...) {
  cat("<multifer_infer_recipe>\n")
  problem <- .recipe_problem(x)
  cat("  geometry:      ", problem$shape$geometry$kind, "\n", sep = "")
  cat("  relation:      ", problem$shape$relation$kind, "\n", sep = "")
  cat("  design:        ", problem$shape$design$kind,   "\n", sep = "")
  cat("  target_family: ", problem$target_family, "\n", sep = "")
  cat("  engine:        ", problem$engine_kind, "\n", sep = "")
  cat("  adapter:       ", problem$adapter_id, "\n", sep = "")
  cat("  targets:       ", paste(problem$targets, collapse = ", "), "\n", sep = "")
  cat("  validity:      ", problem$validity_level, "\n", sep = "")
  cat("  strict:        ", problem$strict, "\n", sep = "")
  if (!is.na(problem$downgrade_reason)) {
    cat("  downgrade:     ", problem$downgrade_reason, "\n", sep = "")
  }
  invisible(x)
}
