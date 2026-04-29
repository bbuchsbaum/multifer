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
#'   `c("oneblock", "cross", "multiblock", "geneig", "adapter")`.
#'   Declares which geometry families this adapter handles.
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
#' - `domains(x, data)` -- optional character vector of loading/score domains.
#' - `project_scores(x, data, domain)` -- optional adapter-native projection of
#'   original data through a bootstrap fit for score stability.
#' - `truncate(x, k)` -- truncated fit.
#' - `residualize(x, k, data)` -- deflated residual after removing k components.
#' - `refit(x, new_data)` -- refit on perturbed data (slow-path fallback).
#' - `component_engine` -- optional character scalar. Set to `"adapter"` to
#'   route component-significance through adapter-owned `component_stat`,
#'   `null_action`, and `residualize` hooks for geometries that otherwise have
#'   a built-in reference engine.
#' - `bootstrap_action(x, data, design, replicate)` -- optional adapter-owned
#'   perturbation hook returning a replicate fit or replicate data.
#' - `core(x, data)` -- optional fast-path core representation.
#' - `update_core(core_obj, ...)` -- optional fast-path core update.
#' - `align(...)` -- optional adapter-owned bootstrap alignment for loadings
#'   and scores.
#' - `null_action(x, data)` -- generate one null resample.
#' - `component_stat(x, data, k, split = NULL)` -- per-component test statistic.
#' - `predict_response(x, new_data, k = NULL)` -- fitted responses for
#'   predictive cross-family adapters.
#' - `variable_stat(x, data, domain, k)` -- variable-level statistic.
#' - `score_stat(x, data, domain, k)` -- score-level statistic.
#'
#' **Capability gate rules** (Part 5 section 36, point 1):
#'
#' - `component_significance` requires `null_action`, `component_stat`, AND
#'   `residualize`. For `(geneig, generalized_eigen)`, the `residualize`
#'   hook must explicitly declare B-metric deflation via
#'   `attr(residualize, "b_metric") <- TRUE` or
#'   `attr(residualize, "delegates_geneig_deflation") <- TRUE`. For
#'   `(cross, predictive)`, the adapter must also provide `refit`,
#'   `predict_response`, and a split-aware `component_stat(..., split = NULL)`.
#'   For `(adapter, predictive)`, the adapter must provide `refit` and
#'   a split-aware `component_stat(..., split = NULL)`; prediction lives
#'   inside the adapter-owned statistic for opaque payloads.
#' - `variable_stability` requires at least one of
#'   `{refit, bootstrap_action, core+update_core}` AND
#'   (`variable_stat` OR `loadings`).
#' - `score_stability` requires at least one of
#'   `{refit, bootstrap_action, core+update_core}`, `loadings`, AND
#'   (`scores` OR `project_scores`). For `geometry = "adapter"`,
#'   `project_scores` is required because multifer cannot infer scores from
#'   opaque adapter-owned payloads.
#' - `subspace_stability` requires at least one of
#'   `{refit, bootstrap_action, core+update_core}` AND `loadings`.
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

  .has_formal_arg <- function(fn, arg) {
    is.function(fn) && arg %in% names(formals(fn))
  }

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
  valid_geom <- c("oneblock", "cross", "multiblock", "geneig", "adapter")
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
  valid_hooks <- c("roots", "scores", "loadings", "domains", "project_scores",
                   "truncate", "residualize", "refit", "component_engine",
                   "bootstrap_action", "core", "update_core", "align",
                   "null_action", "component_stat", "predict_response",
                   "variable_stat", "score_stat")
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

  if (!is.null(hooks[["component_engine"]])) {
    if (!is.character(hooks[["component_engine"]]) ||
        length(hooks[["component_engine"]]) != 1L ||
        is.na(hooks[["component_engine"]]) ||
        !(hooks[["component_engine"]] %in% c("default", "adapter"))) {
      stop("`component_engine` must be either \"default\" or \"adapter\".",
           call. = FALSE)
    }
  }

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
  has_adapter_bootstrap <- !is.null(hooks[["bootstrap_action"]])
  has_perturbation_path <- has_fast_path || has_slow_path || has_adapter_bootstrap

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

  geneig_component_claim <- any(
    capabilities$geometry == "geneig" &
      capabilities$relation == "generalized_eigen" &
      capabilities$target == "component_significance"
  )
  if (geneig_component_claim) {
    residualize_hook <- hooks[["residualize"]]
    has_b_metric_marker <- isTRUE(attr(residualize_hook, "b_metric")) ||
      isTRUE(attr(residualize_hook, "delegates_geneig_deflation"))
    if (!has_b_metric_marker) {
      stop(
        paste0(
          "Claiming (geneig, generalized_eigen, component_significance) requires ",
          "a residualize hook marked with `b_metric = TRUE` or ",
          "`delegates_geneig_deflation = TRUE`. ",
          "B-metric deflation required - see notes/engine_geneig_spec.md."
        ),
        call. = FALSE
      )
    }
  }

  predictive_component_claim <- any(
    capabilities$relation == "predictive" &
      capabilities$target == "component_significance"
  )
  predictive_cross_component_claim <- any(
    capabilities$geometry == "cross" &
      capabilities$relation == "predictive" &
      capabilities$target == "component_significance"
  )
  predictive_adapter_component_claim <- any(
    capabilities$geometry == "adapter" &
      capabilities$relation == "predictive" &
      capabilities$target == "component_significance"
  )
  if (predictive_component_claim) {
    missing_predictive <- character(0)
    if (is.null(hooks[["refit"]])) {
      missing_predictive <- c(missing_predictive, "refit")
    }
    if (predictive_cross_component_claim && is.null(hooks[["predict_response"]])) {
      missing_predictive <- c(missing_predictive, "predict_response")
    }
    if (length(missing_predictive) > 0L || !.has_formal_arg(hooks[["component_stat"]], "split")) {
      split_msg <- if (.has_formal_arg(hooks[["component_stat"]], "split")) {
        NULL
      } else {
        "`component_stat` with a formal `split` argument"
      }
      required <- c(
        "refit",
        if (predictive_cross_component_claim) "predict_response",
        split_msg
      )
      required <- required[!vapply(required, is.null, logical(1))]
      detail <- c(
        if (length(missing_predictive) > 0L) {
          paste0("missing: ", paste(missing_predictive, collapse = ", "))
        },
        if (!.has_formal_arg(hooks[["component_stat"]], "split")) {
          "component_stat is not split-aware"
        }
      )
      stop(
        paste0(
          "Predictive admissibility rule: claiming ",
          if (predictive_adapter_component_claim && !predictive_cross_component_claim) {
            "(adapter, predictive, component_significance)"
          } else if (predictive_cross_component_claim && !predictive_adapter_component_claim) {
            "(cross, predictive, component_significance)"
          } else {
            "predictive component_significance"
          },
          " requires ",
          paste(required, collapse = ", "),
          ". ",
          paste(detail, collapse = "; "),
          ". See notes/engine_predictive_spec.md for the gate contract."
        ),
        call. = FALSE
      )
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
        "`refit`, `bootstrap_action`, or both `core` and `update_core`. None found."
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
    adapter_score_claim <- any(
      capabilities$geometry == "adapter" &
        capabilities$target == "score_stability"
    )
    if (!has_perturbation_path) {
      stop(paste0(
        "Claiming `score_stability` requires at least one perturbation hook: ",
        "`refit`, `bootstrap_action`, or both `core` and `update_core`. None found."
      ), call. = FALSE)
    }
    if (is.null(hooks[["loadings"]])) {
      stop(paste0(
        "Claiming `score_stability` requires the `loadings` hook for bootstrap alignment. ",
        "It is not provided."
      ), call. = FALSE)
    }
    if (is.null(hooks[["scores"]]) && is.null(hooks[["project_scores"]])) {
      stop(paste0(
        "Claiming `score_stability` requires `scores` or `project_scores`. ",
        "Neither is provided."
      ), call. = FALSE)
    }
    if (adapter_score_claim && is.null(hooks[["project_scores"]])) {
      stop(paste0(
        "Claiming `score_stability` for `geometry = \"adapter\"` requires ",
        "`project_scores`. Opaque adapter geometry cannot use the scores-only ",
        "fallback because multifer does not interpret adapter-owned payloads."
      ), call. = FALSE)
    }
  }

  if ("subspace_stability" %in% claimed_targets) {
    if (!has_perturbation_path) {
      stop(paste0(
        "Claiming `subspace_stability` requires at least one perturbation hook: ",
        "`refit`, `bootstrap_action`, or both `core` and `update_core`. None found."
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
#' By default returns a plain character vector so call sites like
#' `"my_adapter" %in% list_infer_adapters()` keep working. Pass
#' `details = TRUE` to get a data.frame keyed by `adapter_id` with the
#' canonical `maturity` label (from [multifer_adapter_maturity()]) plus
#' the distinct `(geometry, relation)` rows the adapter claims. That
#' extended form is the one the drift-guard test and `print` paths
#' consume.
#'
#' @param details Logical. `FALSE` (default) returns a character vector
#'   of adapter ids, preserving the old contract. `TRUE` returns a
#'   data.frame with columns `adapter_id`, `maturity`, `geometry`,
#'   `relation`.
#'
#' @return Character vector or data.frame; see `details`.
#' @seealso [multifer_adapter_maturity()]
#' @export
list_infer_adapters <- function(details = FALSE) {
  ids <- ls(.adapter_registry)
  if (!isTRUE(details)) {
    return(ids)
  }
  if (length(ids) == 0L) {
    return(data.frame(
      adapter_id = character(0),
      maturity = character(0),
      geometry = character(0),
      relation = character(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(ids, function(id) {
    a <- get(id, envir = .adapter_registry, inherits = FALSE)
    caps <- a$capabilities
    pairs <- unique(caps[, c("geometry", "relation")])
    data.frame(
      adapter_id = id,
      maturity = multifer_adapter_maturity(id),
      geometry = pairs$geometry,
      relation = pairs$relation,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
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
