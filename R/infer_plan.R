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

compile_oneblock_ladder_plan <- function(recipe, X, max_steps = NULL) {
  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "oneblock") {
    stop(
      paste0("`recipe` geometry must be \"oneblock\"; got \"",
             recipe$shape$geometry$kind, "\"."),
      call. = FALSE
    )
  }
  if (recipe$shape$relation$kind != "variance") {
    stop(
      paste0("`recipe` relation must be \"variance\"; got \"",
             recipe$shape$relation$kind, "\"."),
      call. = FALSE
    )
  }

  .validate_oneblock_data(X)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  roots_observed <- cached_svd(Xc)$d^2
  zero_tol <- max(1, sum(Xc^2)) * .Machine$double.eps

  if (is.null(max_steps)) {
    max_steps <- min(min(dim(Xc)) - 1L, 50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  observed_stat_fn <- function(step, data) {
    s_top <- top_singular_values(data, 1L)[1L]
    total <- sum(data * data)
    if (total <= zero_tol) return(0)
    (s_top * s_top) / total
  }

  null_stat_fn <- function(step, data) {
    perm <- base::apply(data, 2L, base::sample)
    s <- top_singular_values(perm, step)
    if (length(s) < step) return(0)
    s2 <- s * s
    total_F2 <- sum(data * data)
    tail_total <- total_F2 - sum(s2[seq_len(step - 1L)])
    if (tail_total <= zero_tol) return(0)
    s2[step] / tail_total
  }

  deflate_fn <- function(step, data) {
    sv <- top_svd(data, 1L)
    resid <- data - sv$d[1L] *
      (sv$u[, 1L, drop = FALSE] %*% t(sv$v[, 1L, drop = FALSE]))
    if (sum(resid^2) <= zero_tol) {
      return(matrix(0, nrow = nrow(data), ncol = ncol(data)))
    }
    resid
  }

  structure(
    list(
      initial_data = Xc,
      max_steps = max_steps,
      roots_observed = roots_observed,
      observed_stat_fn = observed_stat_fn,
      null_stat_fn = null_stat_fn,
      deflate_fn = deflate_fn,
      labels = list(
        statistic = "Vitale P3 tail-ratio on latent variance roots",
        null = "column permutation of residual matrix",
        estimand = "latent variance roots"
      )
    ),
    class = "multifer_ladder_plan"
  )
}

is_ladder_plan <- function(x) inherits(x, "multifer_ladder_plan")

.ladder_plan_result <- function(plan,
                                ladder_result,
                                auto_subspace = TRUE,
                                tie_threshold = 0.01) {
  if (!is_ladder_plan(plan)) {
    stop("`plan` must be a multifer_ladder_plan.", call. = FALSE)
  }
  rejected_through <- ladder_result$rejected_through
  roots_observed <- plan$roots_observed
  selected <- logical(length(roots_observed))
  if (rejected_through >= 1L) {
    selected[seq_len(rejected_through)] <- TRUE
  }

  units <- form_units(
    roots_observed,
    selected = selected,
    group_near_ties = isTRUE(auto_subspace),
    tie_threshold = tie_threshold
  )

  sr <- ladder_result$step_results
  component_tests <- data.frame(
    step = vapply(sr, function(x) x$step, integer(1L)),
    observed_stat = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value = vapply(sr, function(x) x$p_value, double(1L)),
    mc_se = vapply(sr, function(x) x$mc_se, double(1L)),
    r = vapply(sr, function(x) x$r, integer(1L)),
    B = vapply(sr, function(x) x$B, integer(1L)),
    selected = vapply(sr, function(x) x$selected, logical(1L)),
    stringsAsFactors = FALSE
  )

  list(
    units = units,
    component_tests = component_tests,
    roots_observed = roots_observed,
    ladder_result = ladder_result,
    labels = plan$labels
  )
}
