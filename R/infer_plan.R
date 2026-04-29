# Internal executable plan compiler.
#
# The public adapter contract is declarative capabilities plus hooks. Runtime
# engines should consume only resolved callbacks. This first plan compiler is
# intentionally narrow: it covers the adapter-driven ladder path and exists so
# relation-specific hook details are absorbed before the ladder runs.

compile_infer_plan <- function(recipe,
                               adapter = NULL,
                               validate_data = NULL,
                               labels = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  adapter <- adapter %||% recipe$adapter
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!is.null(validate_data) && !is.function(validate_data)) {
    stop("`validate_data` must be NULL or a function.", call. = FALSE)
  }

  rel_kind <- recipe$shape$relation$kind
  component_stat <- if (identical(rel_kind, "predictive")) {
    function(fit, data, k) adapter$component_stat(fit, data, k, split = NULL)
  } else {
    function(fit, data, k) adapter$component_stat(fit, data, k)
  }

  validate_cb <- if (is.null(validate_data)) {
    function(data) invisible(data)
  } else {
    validate_data
  }

  structure(
    list(
      recipe = recipe,
      adapter = adapter,
      fit = function(previous_fit, data) adapter$refit(previous_fit, data),
      roots = function(fit) adapter$roots(fit),
      component_stat = component_stat,
      null_action = function(fit, data) adapter$null_action(fit, data),
      residualize = function(fit, k, data) adapter$residualize(fit, k, data),
      validate_data = validate_cb,
      labels = .infer_plan_labels(recipe, labels)
    ),
    class = "multifer_infer_plan"
  )
}

is_infer_plan <- function(x) inherits(x, "multifer_infer_plan")

.infer_plan_labels <- function(recipe, labels = NULL) {
  labels <- labels %||% list()
  geom_kind <- recipe$shape$geometry$kind
  rel_kind <- recipe$shape$relation$kind

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

  list(
    statistic = labels$statistic %||% default_statistic,
    null = labels$null %||% sprintf("adapter null_action on %s data", geom_kind),
    estimand = labels$estimand %||% default_estimand
  )
}
