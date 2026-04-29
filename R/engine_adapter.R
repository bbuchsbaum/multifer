#' Run the generic adapter-owned sequential deflation engine
#'
#' Executes the latent-root ladder through adapter hooks. This engine is for
#' model families whose component statistic, null action, or residualization is
#' not the package's built-in Euclidean oneblock/cross reference calculation.
#' The adapter owns `component_stat()`, `null_action()`, and `residualize()`;
#' the engine only supplies the ladder control flow.
#'
#' At each rung, hook data are the current residual data object. The engine
#' calls `component_stat(fit, data, 1L)` and `residualize(fit, 1L, data)`, so
#' adapters should interpret `k = 1L` as "the leading component of the current
#' residual", not "remove the first global k components again".
#'
#' @param recipe A compiled `multifer_infer_recipe`.
#' @param adapter A `multifer_adapter`.
#' @param data Data object accepted by the adapter hooks.
#' @param validate_data Optional function used to validate the data object at
#'   the original, null, and residual rungs.
#' @param labels Optional named character list with `statistic`, `null`, and
#'   `estimand` labels for result provenance.
#' @param original_fit Optional pre-fit object for `data`.
#' @inheritParams run_oneblock_ladder
#'
#' @return A list with the same fields as `run_oneblock_ladder()`.
run_adapter_ladder <- function(recipe,
                               adapter,
                               data,
                               B             = 1000L,
                               B_total       = NULL,
                               batch_size    = 32L,
                               alpha         = 0.05,
                               max_steps     = NULL,
                               seed          = NULL,
                               auto_subspace = TRUE,
                               tie_threshold = 0.01,
                               original_fit  = NULL,
                               validate_data = NULL,
                               labels        = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!is.null(validate_data) && !is.function(validate_data)) {
    stop("`validate_data` must be NULL or a function.", call. = FALSE)
  }

  if (!is.null(validate_data)) {
    validate_data(data)
  }

  if (is.null(original_fit)) {
    original_fit <- adapter$refit(NULL, data)
  }
  roots_observed <- adapter$roots(original_fit)
  if (!is.numeric(roots_observed) || length(roots_observed) == 0L ||
      any(!is.finite(roots_observed))) {
    stop("`adapter$roots()` must return a non-empty finite numeric vector.",
         call. = FALSE)
  }

  if (is.null(max_steps)) {
    max_steps <- min(length(roots_observed), 50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  rel_kind <- recipe$shape$relation$kind
  current_step <- NA_integer_
  current_fit <- NULL
  get_fit <- function(step, rung_data) {
    if (!identical(current_step, as.integer(step)) || is.null(current_fit)) {
      current_fit <<- if (identical(step, 1L)) {
        original_fit
      } else {
        adapter$refit(original_fit, rung_data)
      }
      current_step <<- as.integer(step)
    }
    current_fit
  }
  reset_fit <- function() {
    current_step <<- NA_integer_
    current_fit <<- NULL
  }

  observed_stat_fn <- function(step, rung_data) {
    fit <- get_fit(step, rung_data)
    stat <- if (identical(rel_kind, "predictive")) {
      adapter$component_stat(fit, rung_data, 1L, split = NULL)
    } else {
      adapter$component_stat(fit, rung_data, 1L)
    }
    .multifer_scalar_stat(stat)
  }

  null_stat_fn <- function(step, rung_data) {
    fit <- get_fit(step, rung_data)
    null_data <- adapter$null_action(fit, rung_data)
    if (!is.null(validate_data)) {
      validate_data(null_data)
    }
    stat <- if (identical(rel_kind, "predictive")) {
      adapter$component_stat(fit, null_data, 1L, split = NULL)
    } else {
      adapter$component_stat(fit, null_data, 1L)
    }
    .multifer_scalar_stat(stat)
  }

  deflate_fn <- function(step, rung_data) {
    fit <- get_fit(step, rung_data)
    out <- adapter$residualize(fit, 1L, rung_data)
    if (!is.null(validate_data)) {
      validate_data(out)
    }
    reset_fit()
    out
  }

  ladder_result <- ladder_driver(
    observed_stat_fn = observed_stat_fn,
    null_stat_fn     = null_stat_fn,
    deflate_fn       = deflate_fn,
    initial_data     = data,
    max_steps        = max_steps,
    B                = B,
    B_total          = B_total,
    batch_size       = batch_size,
    alpha            = alpha,
    seed             = seed
  )

  rejected_through <- ladder_result$rejected_through
  selected <- logical(length(roots_observed))
  if (rejected_through >= 1L) {
    selected[seq_len(min(rejected_through, length(selected)))] <- TRUE
  }

  units <- form_units(
    roots_observed,
    selected        = selected,
    group_near_ties = isTRUE(auto_subspace),
    tie_threshold   = tie_threshold
  )

  sr <- ladder_result$step_results
  component_tests <- data.frame(
    step           = vapply(sr, function(x) x$step,          integer(1L)),
    observed_stat  = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value        = vapply(sr, function(x) x$p_value,       double(1L)),
    mc_se          = vapply(sr, function(x) x$mc_se,         double(1L)),
    r              = vapply(sr, function(x) x$r,             integer(1L)),
    B              = vapply(sr, function(x) x$B,             integer(1L)),
    selected       = vapply(sr, function(x) x$selected,      logical(1L)),
    stringsAsFactors = FALSE
  )

  geom_kind <- recipe$shape$geometry$kind
  default_statistic <- if (identical(rel_kind, "predictive")) {
    sprintf("adapter predictive component_stat on %s data", geom_kind)
  } else {
    sprintf("adapter component_stat on %s data", geom_kind)
  }
  default_estimand <- if (identical(rel_kind, "predictive")) {
    "adapter predictive gain"
  } else {
    "adapter latent roots"
  }
  labels <- labels %||% list()
  list(
    units           = units,
    component_tests = component_tests,
    roots_observed  = roots_observed,
    ladder_result   = ladder_result,
    labels          = list(
      statistic = labels$statistic %||%
        default_statistic,
      null      = labels$null %||%
        sprintf("adapter null_action on %s data", geom_kind),
      estimand  = labels$estimand %||% default_estimand
    ),
    original_fit = original_fit
  )
}

.validate_oneblock_data <- function(data) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("For oneblock geometry, `data` must be a numeric matrix.",
         call. = FALSE)
  }
  if (any(!is.finite(data))) {
    stop("Oneblock data must not contain NA, NaN, or Inf.", call. = FALSE)
  }
  if (nrow(data) < 2L || ncol(data) < 2L) {
    stop("Oneblock data must have at least 2 rows and 2 columns.",
         call. = FALSE)
  }
  invisible(data)
}
