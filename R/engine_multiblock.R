#' Run the generic multiblock sequential deflation engine
#'
#' Executes the latent-root ladder for adapters that declare
#' `geometry = "multiblock"`. Unlike the oneblock and cross reference
#' engines, this path is deliberately adapter-driven: the adapter owns
#' refitting, null generation, component statistics, and residualization.
#' At each rung the hook data are the current residual block list, so
#' `component_stat()` and `residualize()` are called with `k = 1L` to act
#' on the leading component of that current residual. The ladder step remains
#' the global rung label in the returned results.
#'
#' @param recipe A compiled `multifer_infer_recipe` with geometry
#'   `"multiblock"`.
#' @param adapter A `multifer_adapter` whose hooks support the recipe.
#' @param data A named or unnamed list of aligned numeric matrix blocks.
#' @param original_fit Optional pre-fit object for `data`. When supplied,
#'   the observed roots and first ladder rung use this fit instead of
#'   refitting `data`.
#' @param max_steps Maximum number of ladder rungs. Defaults to the number
#'   of roots returned by the adapter, capped at 50.
#' @inheritParams run_oneblock_ladder
#'
#' @return A list with the same fields as `run_oneblock_ladder()`.
#' @export
run_multiblock_ladder <- function(recipe,
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
                                  original_fit  = NULL) {

  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "multiblock") {
    stop(
      paste0("`recipe` geometry must be \"multiblock\"; got \"",
             recipe$shape$geometry$kind, "\"."),
      call. = FALSE
    )
  }
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }

  .validate_multiblock_data(data)

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

  current_step <- NA_integer_
  current_fit <- NULL
  get_fit <- function(step, block_data) {
    if (!identical(current_step, as.integer(step)) || is.null(current_fit)) {
      current_fit <<- if (identical(step, 1L)) {
        original_fit
      } else {
        adapter$refit(original_fit, block_data)
      }
      current_step <<- as.integer(step)
    }
    current_fit
  }
  reset_fit <- function() {
    current_step <<- NA_integer_
    current_fit <<- NULL
  }

  observed_stat_fn <- function(step, block_data) {
    fit <- get_fit(step, block_data)
    .multifer_scalar_stat(adapter$component_stat(fit, block_data, 1L))
  }

  null_stat_fn <- function(step, block_data) {
    fit <- get_fit(step, block_data)
    null_data <- adapter$null_action(fit, block_data)
    .validate_multiblock_data(null_data)
    .multifer_scalar_stat(adapter$component_stat(fit, null_data, 1L))
  }

  deflate_fn <- function(step, block_data) {
    fit <- get_fit(step, block_data)
    out <- adapter$residualize(fit, 1L, block_data)
    .validate_multiblock_data(out)
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

  list(
    units           = units,
    component_tests = component_tests,
    roots_observed  = roots_observed,
    ladder_result   = ladder_result,
    labels          = list(
      statistic = "adapter component_stat on multiblock data",
      null      = "adapter null_action on aligned block-list data",
      estimand  = "adapter latent roots"
    )
  )
}

.multifer_scalar_stat <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`adapter$component_stat()` must return one finite numeric scalar.",
         call. = FALSE)
  }
  as.numeric(x)
}

.validate_multiblock_data <- function(data, min_blocks = 2L) {
  if (!is.list(data) || length(data) < min_blocks) {
    stop("For multiblock geometry, `data` must be a list of matrix blocks.",
         call. = FALSE)
  }
  bad <- vapply(data, function(x) !is.matrix(x) || !is.numeric(x), logical(1L))
  if (any(bad)) {
    stop("Every multiblock data block must be a numeric matrix.",
         call. = FALSE)
  }
  if (any(vapply(data, function(x) any(!is.finite(x)), logical(1L)))) {
    stop("Multiblock data blocks must not contain NA, NaN, or Inf.",
         call. = FALSE)
  }
  n <- vapply(data, nrow, integer(1L))
  if (any(n < 2L) || length(unique(n)) != 1L) {
    stop("All multiblock data blocks must have the same number of rows >= 2.",
         call. = FALSE)
  }
  if (any(vapply(data, ncol, integer(1L)) < 1L)) {
    stop("Every multiblock data block must have at least one column.",
         call. = FALSE)
  }
  invisible(data)
}

.multiblock_domain_names <- function(data) {
  .validate_multiblock_data(data)
  nm <- names(data)
  if (is.null(nm)) {
    nm <- character(length(data))
  }
  missing <- is.na(nm) | !nzchar(nm)
  nm[missing] <- paste0("block", which(missing))
  if (any(duplicated(nm))) {
    stop("Multiblock block names must be unique when used as domains.",
         call. = FALSE)
  }
  nm
}

.adapter_domains <- function(adapter, fit = NULL, data = NULL, geom_kind = NULL) {
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }

  out <- NULL
  if (!is.null(adapter$domains)) {
    out <- adapter$domains(fit, data = data)
  } else if (identical(geom_kind, "oneblock")) {
    out <- "X"
  } else if (identical(geom_kind, "cross")) {
    out <- c("X", "Y")
  } else if (identical(geom_kind, "multiblock")) {
    out <- .multiblock_domain_names(data)
  }

  if (!is.character(out) || length(out) == 0L ||
      anyNA(out) || any(!nzchar(out)) || any(duplicated(out))) {
    stop("Adapter domains must be a non-empty character vector of unique names.",
         call. = FALSE)
  }
  out
}

.resample_multiblock_data <- function(data, indices) {
  .validate_multiblock_data(data)
  out <- lapply(data, function(block) block[indices, , drop = FALSE])
  names(out) <- names(data)
  out
}
