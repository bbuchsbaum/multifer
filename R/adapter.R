#' Adapter for plugging an external model into the multifer inference engine
#'
#' An `infer_adapter` is the contract object that lets external packages (or
#' advanced users) register a model class with `multifer` without inheriting
#' from any `multivarious` class.  It carries three things:
#'
#' 1. **Capability declarations** -- a `multifer_capability_matrix` listing
#'    every `(geometry, relation, target)` triple this adapter claims to
#'    support.
#' 2. **Hook functions** -- R functions the engine calls at runtime.
#' 3. **Validity metadata** -- `validity_level`, `declared_assumptions`, and
#'    `checked_assumptions` (Part 5 section 36).
#'
#' @section Registration-time gate (Part 5 section 36, point 1):
#' The constructor enforces that every claimed capability is backed by the
#' required hooks.  Claiming a capability without the hooks is a registration
#' error that fires immediately, before any data is touched.
#'
#' @param adapter_id Non-empty character scalar.  Unique identifier used as
#'   the registry key.
#' @param adapter_version Character scalar, e.g. `"0.0.1"`.
#' @param shape_kinds Character vector.  Subset of
#'   `c("oneblock", "cross", "multiblock", "geneig")`.  Declares which
#'   geometry families this adapter handles.
#' @param capabilities A `multifer_capability_matrix` (see
#'   [capability_matrix()]).  Every `(geometry, relation, target)` triple
#'   this adapter claims to support.
#' @param ... Hook functions, passed by name.  See Details for the full list.
#' @param validity_level One of `c("exact", "conditional", "asymptotic",
#'   "heuristic")`.
#' @param declared_assumptions Character vector.  Trust-the-user assumptions
#'   (noise model, design correctness, etc.).
#' @param checked_assumptions Named list of runtime-testable assumption checks.
#'   Each element should be a list with fields `name` (character),
#'   `check` (function(x, data) returning logical), and `detail` (character).
#'
#' @details
#' **Hook functions** (all passed via `...`, all optional except as required by
#' the capability gate):
#'
#' - `roots(x)` -- numeric vector of latent roots.
#' - `scores(x, domain)` -- matrix of scores.
#' - `loadings(x, domain)` -- matrix of loadings.
#' - `truncate(x, k)` -- truncated fit.
#' - `residualize(x, k, data)` -- deflated residual after removing k components.
#' - `refit(x, new_data)` -- refit on perturbed data (slow-path fallback).
#' - `core(x, data)` -- optional fast-path core representation.
#' - `update_core(core_obj, ...)` -- optional fast-path core update.
#' - `align(xb, xref)` -- alignment (sign / legacy Procrustes / subspace).
#' - `null_action(x, data)` -- generate one null resample.
#' - `component_stat(x, data, k)` -- per-component test statistic.
#' - `variable_stat(x, data, domain, k)` -- variable-level statistic.
#' - `score_stat(x, data, domain, k)` -- score-level statistic.
#'
#' **Capability gate rules** (Part 5 section 36, point 1):
#'
#' - `component_significance` requires `null_action`, `component_stat`, AND
#'   `residualize`.
#' - `variable_stability` requires at least one of `{refit, core+update_core}`
#'   AND (`variable_stat` OR `loadings`).
#' - `score_stability` requires at least one of `{refit, core+update_core}`
#'   AND `scores`.
#' - `subspace_stability` requires at least one of `{refit, core+update_core}`
#'   AND `loadings`.
#' - `variable_significance` is excluded from v1 (Part 5 section 38) and
#'   always errors.
#'
#' @return An object of class `multifer_adapter`.
#'
#' @export
infer_adapter <- function(adapter_id,
                          adapter_version = "0.0.1",
                          shape_kinds,
                          capabilities,
                          ...,
                          validity_level,
                          declared_assumptions = character(0),
                          checked_assumptions = list()) {

  # -- scalar string checks ---------------------------------------------------
  if (!is.character(adapter_id) || length(adapter_id) != 1L ||
      is.na(adapter_id) || nchar(adapter_id) == 0L) {
    stop("`adapter_id` must be a non-empty single string.", call. = FALSE)
  }
  if (!is.character(adapter_version) || length(adapter_version) != 1L ||
      is.na(adapter_version)) {
    stop("`adapter_version` must be a single string.", call. = FALSE)
  }

  # -- shape_kinds ------------------------------------------------------------
  valid_geom <- c("oneblock", "cross", "multiblock", "geneig")
  if (!is.character(shape_kinds) || length(shape_kinds) == 0L) {
    stop("`shape_kinds` must be a non-empty character vector.", call. = FALSE)
  }
  bad_geom <- setdiff(shape_kinds, valid_geom)
  if (length(bad_geom) > 0L) {
    stop(sprintf(
      "`shape_kinds` contains unknown geometry/geometries: %s. Must be subset of: %s.",
      paste(bad_geom, collapse = ", "),
      paste(valid_geom, collapse = ", ")
    ), call. = FALSE)
  }

  # -- capabilities -----------------------------------------------------------
  if (!is_capability_matrix(capabilities)) {
    stop("`capabilities` must be a multifer_capability_matrix (see capability_matrix()).",
         call. = FALSE)
  }

  # Every geometry in capabilities must appear in shape_kinds
  cap_geoms <- unique(capabilities$geometry)
  missing_geom <- setdiff(cap_geoms, shape_kinds)
  if (length(missing_geom) > 0L) {
    stop(sprintf(
      "capabilities reference geometry/geometries not listed in `shape_kinds`: %s.",
      paste(missing_geom, collapse = ", ")
    ), call. = FALSE)
  }

  # -- validity_level ---------------------------------------------------------
  valid_vl <- c("exact", "conditional", "asymptotic", "heuristic")
  if (!is.character(validity_level) || length(validity_level) != 1L ||
      is.na(validity_level) || !(validity_level %in% valid_vl)) {
    stop(sprintf(
      "`validity_level` must be one of: %s.",
      paste(valid_vl, collapse = ", ")
    ), call. = FALSE)
  }

  # -- collect hooks ----------------------------------------------------------
  hooks_raw <- list(...)
  valid_hooks <- c("roots", "scores", "loadings", "truncate", "residualize",
                   "refit", "core", "update_core", "align", "null_action",
                   "component_stat", "variable_stat", "score_stat")
  bad_hooks <- setdiff(names(hooks_raw), valid_hooks)
  if (length(bad_hooks) > 0L) {
    stop(sprintf(
      "Unknown hook(s) in `...`: %s. Valid hooks: %s.",
      paste(bad_hooks, collapse = ", "),
      paste(valid_hooks, collapse = ", ")
    ), call. = FALSE)
  }
  # Normalise: named list with NULL for absent hooks
  hooks <- stats::setNames(
    lapply(valid_hooks, function(h) hooks_raw[[h]]),
    valid_hooks
  )

  # -- variable_significance excluded from v1 (Part 5 section 38) ------------
  if (any(capabilities$target == "variable_significance")) {
    stop(paste0(
      "Claiming `variable_significance` is not permitted in v1. ",
      "Variable significance is deferred to Phase 3 (Part 5 section 38). ",
      "Remove `variable_significance` from `capabilities` to register this adapter."
    ), call. = FALSE)
  }

  # -- registration-time capability gate (Part 5 section 36, point 1) --------
  claimed_targets <- unique(capabilities$target)

  # helper to check if fast path is available: core + update_core both present
  has_fast_path <- !is.null(hooks[["core"]]) && !is.null(hooks[["update_core"]])
  has_slow_path <- !is.null(hooks[["refit"]])
  has_perturbation_path <- has_fast_path || has_slow_path

  if ("component_significance" %in% claimed_targets) {
    missing_cs <- character(0)
    if (is.null(hooks[["null_action"]]))    missing_cs <- c(missing_cs, "null_action")
    if (is.null(hooks[["component_stat"]])) missing_cs <- c(missing_cs, "component_stat")
    if (is.null(hooks[["residualize"]]))    missing_cs <- c(missing_cs, "residualize")
    if (length(missing_cs) > 0L) {
      stop(sprintf(
        "Claiming `component_significance` requires hooks: %s. Missing: %s.",
        "null_action, component_stat, residualize",
        paste(missing_cs, collapse = ", ")
      ), call. = FALSE)
    }
  }

  if ("variable_stability" %in% claimed_targets) {
    # Requires a perturbation path (refit or core/update_core) AND
    # variable_stat or loadings.
    # Rationale (Part 5 section 36): variable stability is computed by
    # perturbing the data and summarising loading variation; both a
    # perturbation mechanism and a way to extract loadings are required.
    if (!has_perturbation_path) {
      stop(paste0(
        "Claiming `variable_stability` requires at least one perturbation hook: ",
        "`refit` or both `core` and `update_core`. None found."
      ), call. = FALSE)
    }
    if (is.null(hooks[["variable_stat"]]) && is.null(hooks[["loadings"]])) {
      stop(paste0(
        "Claiming `variable_stability` requires `variable_stat` or `loadings`. ",
        "Neither is provided."
      ), call. = FALSE)
    }
  }

  if ("score_stability" %in% claimed_targets) {
    if (!has_perturbation_path) {
      stop(paste0(
        "Claiming `score_stability` requires at least one perturbation hook: ",
        "`refit` or both `core` and `update_core`. None found."
      ), call. = FALSE)
    }
    if (is.null(hooks[["scores"]])) {
      stop(paste0(
        "Claiming `score_stability` requires the `scores` hook. It is not provided."
      ), call. = FALSE)
    }
  }

  if ("subspace_stability" %in% claimed_targets) {
    if (!has_perturbation_path) {
      stop(paste0(
        "Claiming `subspace_stability` requires at least one perturbation hook: ",
        "`refit` or both `core` and `update_core`. None found."
      ), call. = FALSE)
    }
    if (is.null(hooks[["loadings"]])) {
      stop(paste0(
        "Claiming `subspace_stability` requires the `loadings` hook. ",
        "It is not provided."
      ), call. = FALSE)
    }
  }

  # -- declared_assumptions ---------------------------------------------------
  if (!is.character(declared_assumptions)) {
    stop("`declared_assumptions` must be a character vector.", call. = FALSE)
  }

  # -- checked_assumptions ----------------------------------------------------
  if (!is.list(checked_assumptions)) {
    stop("`checked_assumptions` must be a list.", call. = FALSE)
  }

  # -- assemble ---------------------------------------------------------------
  structure(
    c(
      list(
        adapter_id           = adapter_id,
        adapter_version      = adapter_version,
        shape_kinds          = shape_kinds,
        capabilities         = capabilities,
        validity_level       = validity_level,
        declared_assumptions = declared_assumptions,
        checked_assumptions  = checked_assumptions
      ),
      hooks
    ),
    class = "multifer_adapter"
  )
}

#' Test whether an object is a multifer adapter
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_infer_adapter <- function(x) inherits(x, "multifer_adapter")

#' Return the capability matrix of an adapter
#'
#' @param adapter A `multifer_adapter` object.
#' @return A `multifer_capability_matrix`.
#' @export
adapter_capabilities <- function(adapter) {
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  adapter$capabilities
}

#' Query whether an adapter supports a specific (geometry, relation, target) triple
#'
#' Thin wrapper over [supports()] using the adapter's capability matrix.
#'
#' @param adapter A `multifer_adapter` object.
#' @param geometry Character scalar.
#' @param relation Character scalar.
#' @param target Character scalar.
#' @return Logical scalar.
#' @export
adapter_supports <- function(adapter, geometry, relation, target) {
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  supports(adapter$capabilities, geometry, relation, target)
}

#' @export
print.multifer_adapter <- function(x, ...) {
  n_caps <- nrow(x$capabilities)
  cat("<multifer_adapter>\n")
  cat("  id:           ", x$adapter_id, "\n", sep = "")
  cat("  version:      ", x$adapter_version, "\n", sep = "")
  cat("  shape_kinds:  ", paste(x$shape_kinds, collapse = ", "), "\n", sep = "")
  cat("  validity:     ", x$validity_level, "\n", sep = "")
  cat("  capabilities: ", n_caps, " declared triple(s)\n", sep = "")
  invisible(x)
}

# ==============================================================================
# Registry
# ==============================================================================

# Package-level adapter registry.  Non-exported environment.
.adapter_registry <- new.env(parent = emptyenv())

#' Register an infer adapter
#'
#' Stores a `multifer_adapter` in the package-level registry keyed by its
#' `adapter_id`.  Errors if the id is already registered unless
#' `overwrite = TRUE`.
#'
#' @param adapter_id Character scalar.  Must match `adapter$adapter_id`.
#' @param adapter A `multifer_adapter` object (see [infer_adapter()]).
#' @param overwrite Logical.  If `FALSE` (default), errors on duplicate ids.
#'   Set to `TRUE` to replace an existing registration.
#'
#' @return Invisibly returns `adapter_id`.
#' @export
register_infer_adapter <- function(adapter_id, adapter, overwrite = FALSE) {
  if (!is.character(adapter_id) || length(adapter_id) != 1L ||
      is.na(adapter_id) || nchar(adapter_id) == 0L) {
    stop("`adapter_id` must be a non-empty single string.", call. = FALSE)
  }
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!isTRUE(overwrite) && exists(adapter_id, envir = .adapter_registry,
                                   inherits = FALSE)) {
    stop(sprintf(
      "Adapter '%s' is already registered. Use `overwrite = TRUE` to replace it.",
      adapter_id
    ), call. = FALSE)
  }
  assign(adapter_id, adapter, envir = .adapter_registry)
  invisible(adapter_id)
}

#' Look up a registered adapter by id
#'
#' @param adapter_id Character scalar.
#' @return The `multifer_adapter` registered under that id.
#' @export
get_infer_adapter <- function(adapter_id) {
  if (!is.character(adapter_id) || length(adapter_id) != 1L) {
    stop("`adapter_id` must be a single string.", call. = FALSE)
  }
  if (!exists(adapter_id, envir = .adapter_registry, inherits = FALSE)) {
    stop(sprintf("No adapter registered with id '%s'.", adapter_id),
         call. = FALSE)
  }
  get(adapter_id, envir = .adapter_registry, inherits = FALSE)
}

#' List all registered adapter ids
#'
#' @return Character vector of registered adapter ids (possibly empty).
#' @export
list_infer_adapters <- function() {
  ls(.adapter_registry)
}

#' Remove a registered adapter
#'
#' @param adapter_id Character scalar.  Errors if not registered.
#' @return Invisibly returns `adapter_id`.
#' @export
unregister_infer_adapter <- function(adapter_id) {
  if (!is.character(adapter_id) || length(adapter_id) != 1L) {
    stop("`adapter_id` must be a single string.", call. = FALSE)
  }
  if (!exists(adapter_id, envir = .adapter_registry, inherits = FALSE)) {
    stop(sprintf("No adapter registered with id '%s'.", adapter_id),
         call. = FALSE)
  }
  rm(list = adapter_id, envir = .adapter_registry)
  invisible(adapter_id)
}

#' Clear all registered adapters (test helper)
#'
#' Removes every entry from the package adapter registry.  Intended for use
#' in `setup.R` or at the top of test files that exercise the registry, so
#' that tests do not interfere with each other.
#'
#' This function is exported so test files can call it without `:::`.
#'
#' @return Invisibly returns `NULL`.
#' @export
clear_adapter_registry <- function() {
  rm(list = ls(.adapter_registry), envir = .adapter_registry)
  invisible(NULL)
}
