#' Compiled inferential problem object
#'
#' `infer_problem()` is the first-class representation of the inferential
#' question that `multifer` has compiled from a typed shape and an adapter's
#' declared capabilities. It keeps the package's ontology in one inspectable
#' place rather than spreading it across recipe fields and dispatch locals.
#'
#' The object records:
#'
#' - the typed shape being inferred on,
#' - the requested inferential targets,
#' - the target family (`"latent_root"` or `"predictive_gain"`),
#' - the engine kind selected for execution,
#' - the adapter identity and compiled validity level,
#' - whether strict dispatch was enforced and whether any downgrade occurred.
#'
#' @param shape A `multifer_typed_shape`.
#' @param targets Character vector of resolved inferential targets.
#' @param adapter_id Character scalar naming the adapter.
#' @param adapter_version Optional adapter version string.
#' @param validity_level One of `c("exact", "conditional", "asymptotic",
#'   "heuristic")`.
#' @param strict Logical scalar recording whether strict dispatch was enforced.
#' @param downgrade_reason Optional character scalar. Use `NA_character_` when
#'   no downgrade occurred.
#' @param target_family Optional character scalar. When `NULL`, derived from
#'   the shape relation.
#' @param engine_kind Optional character scalar. When `NULL`, derived from the
#'   shape and relation.
#' @param call Optional matched call or other provenance object.
#'
#' @return An object of class `multifer_infer_problem`.
#' @export
infer_problem <- function(shape,
                          targets,
                          adapter_id,
                          adapter_version = NA_character_,
                          validity_level,
                          strict = TRUE,
                          downgrade_reason = NA_character_,
                          target_family = NULL,
                          engine_kind = NULL,
                          call = NULL) {
  if (!is_typed_shape(shape)) {
    stop("`shape` must be a multifer_typed_shape.", call. = FALSE)
  }
  if (!is.character(targets) || length(targets) == 0L) {
    stop("`targets` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.character(adapter_id) || length(adapter_id) != 1L ||
      is.na(adapter_id) || nchar(adapter_id) == 0L) {
    stop("`adapter_id` must be a non-empty single string.", call. = FALSE)
  }
  valid_vl <- c("exact", "conditional", "asymptotic", "heuristic")
  if (!is.character(validity_level) || length(validity_level) != 1L ||
      is.na(validity_level) || !(validity_level %in% valid_vl)) {
    stop(sprintf(
      "`validity_level` must be one of: %s.",
      paste(valid_vl, collapse = ", ")
    ), call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be a single non-NA logical.", call. = FALSE)
  }

  if (is.null(target_family)) {
    target_family <- .infer_target_family(shape$relation$kind)
  }
  if (is.null(engine_kind)) {
    engine_kind <- .infer_engine_kind(shape)
  }

  structure(
    list(
      shape = shape,
      targets = unique(as.character(targets)),
      target_family = as.character(target_family),
      engine_kind = as.character(engine_kind),
      adapter_id = adapter_id,
      adapter_version = as.character(adapter_version),
      validity_level = validity_level,
      strict = as.logical(strict),
      downgrade_reason = downgrade_reason,
      call = call
    ),
    class = "multifer_infer_problem"
  )
}

.infer_target_family <- function(relation_kind) {
  if (identical(relation_kind, "predictive")) {
    "predictive_gain"
  } else {
    "latent_root"
  }
}

.infer_engine_kind <- function(shape) {
  geom_kind <- shape$geometry$kind
  rel_kind  <- shape$relation$kind

  if (identical(geom_kind, "cross") && identical(rel_kind, "predictive")) {
    return("predictive")
  }
  if (identical(geom_kind, "adapter") && identical(rel_kind, "predictive")) {
    return("adapter_predictive")
  }
  if (identical(geom_kind, "oneblock")) {
    return("oneblock")
  }
  if (identical(geom_kind, "cross")) {
    return("cross")
  }
  if (identical(geom_kind, "geneig")) {
    return("geneig")
  }

  geom_kind
}

.recipe_problem <- function(recipe) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a multifer_infer_recipe.", call. = FALSE)
  }
  if (!is.null(recipe$problem) && is_infer_problem(recipe$problem)) {
    return(recipe$problem)
  }

  infer_problem(
    shape = recipe$shape,
    targets = recipe$targets,
    adapter_id = recipe$adapter_id,
    adapter_version = recipe$adapter$adapter_version %||% NA_character_,
    validity_level = recipe$validity_level,
    strict = recipe$strict,
    downgrade_reason = recipe$downgrade_reason,
    call = recipe$call
  )
}

#' Test whether an object is a compiled inferential problem
#'
#' @param x Object to test.
#' @return Logical scalar.
#' @export
is_infer_problem <- function(x) inherits(x, "multifer_infer_problem")

#' @export
print.multifer_infer_problem <- function(x, ...) {
  cat("<multifer_infer_problem>\n")
  cat("  geometry:      ", x$shape$geometry$kind, "\n", sep = "")
  cat("  relation:      ", x$shape$relation$kind, "\n", sep = "")
  cat("  design:        ", x$shape$design$kind, "\n", sep = "")
  cat("  target_family: ", x$target_family, "\n", sep = "")
  cat("  engine:        ", x$engine_kind, "\n", sep = "")
  cat("  adapter:       ", x$adapter_id, "\n", sep = "")
  cat("  targets:       ", paste(x$targets, collapse = ", "), "\n", sep = "")
  cat("  validity:      ", x$validity_level, "\n", sep = "")
  cat("  strict:        ", x$strict, "\n", sep = "")
  if (!is.na(x$downgrade_reason)) {
    cat("  downgrade:     ", x$downgrade_reason, "\n", sep = "")
  }
  invisible(x)
}
