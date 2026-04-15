#' Run the predictive sequential deflation engine
#'
#' Implements the test ladder for `(cross, predictive)` recipes via the
#' shared [ladder_driver()]. The observed statistic at each rung is the
#' cross-fitted incremental predictive gain of the leading latent
#' predictor on the current deflated data. Gains are computed out of
#' sample over user-supplied folds (or a deterministic balanced K-fold
#' partition when `folds = NULL`).
#'
#' Predictive relations are ordered by rung rather than by descending
#' latent-root magnitude. The returned `roots_observed` vector therefore
#' records the held-out gain sequence actually evaluated by the ladder,
#' while unit formation uses a strictly decreasing proxy so the shared
#' [form_units()] machinery can preserve rung order without implying a
#' latent-root interpretation.
#'
#' @param recipe A compiled `multifer_infer_recipe`. Must have geometry
#'   `"cross"` and relation `"predictive"`.
#' @param X Numeric matrix, `n x p`.
#' @param Y Numeric matrix, `n x q`. Must have the same number of rows
#'   as `X`.
#' @param folds Optional fold assignment vector of length `n`. When
#'   `NULL`, a balanced fold assignment is created internally.
#' @param n_folds Positive integer. Used only when `folds = NULL`.
#'   Default `5L`.
#' @param B Per-rung cap on Monte Carlo draws.
#' @param B_total Optional integer global Monte Carlo budget shared
#'   across ladder rungs. Defaults to `B * max_steps`.
#' @param batch_size Positive integer, Besag-Clifford batch size within
#'   each rung. Default `32L`.
#' @param alpha Significance threshold.
#' @param max_steps Maximum ladder rungs. Default
#'   `min(nrow(X) - 1, ncol(X), ncol(Y), 50L)`.
#' @param seed Integer or NULL.
#' @param auto_subspace Logical. Passed through to [form_units()].
#' @param tie_threshold Positive numeric. Passed through to [form_units()].
#'
#' @return A list with `units`, `component_tests`, `roots_observed`, and
#'   `ladder_result`, mirroring [run_oneblock_ladder()] and
#'   [run_cross_ladder()].
#'
#' @export
run_predictive_ladder <- function(recipe,
                                  X,
                                  Y,
                                  folds         = NULL,
                                  n_folds       = 5L,
                                  B             = 1000L,
                                  B_total       = NULL,
                                  batch_size    = 32L,
                                  alpha         = 0.05,
                                  max_steps     = NULL,
                                  seed          = NULL,
                                  auto_subspace = TRUE,
                                  tie_threshold = 0.01) {

  if (!is_infer_recipe(recipe)) {
    stop("`recipe` must be a compiled multifer_infer_recipe.", call. = FALSE)
  }
  if (recipe$shape$geometry$kind != "cross") {
    stop(
      paste0("`recipe` geometry must be \"cross\"; got \"",
             recipe$shape$geometry$kind, "\"."),
      call. = FALSE
    )
  }
  if (recipe$shape$relation$kind != "predictive") {
    stop(
      paste0("`recipe` relation must be \"predictive\"; got \"",
             recipe$shape$relation$kind, "\"."),
      call. = FALSE
    )
  }

  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("`Y` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y)) {
    stop("`X` and `Y` must have the same number of rows.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("`X` and `Y` must not contain NA, NaN, or Inf.", call. = FALSE)
  }
  if (nrow(X) < 4L || ncol(X) < 1L || ncol(Y) < 1L) {
    stop(
      "`X` and `Y` must have at least 4 rows and at least 1 column each.",
      call. = FALSE
    )
  }

  adapter <- recipe$adapter
  if (is.null(adapter$refit) || is.null(adapter$predict_response) ||
      is.null(adapter$residualize) || is.null(adapter$null_action)) {
    stop(
      paste0(
        "Predictive engine requires adapter hooks `refit`, `predict_response`, ",
        "`residualize`, and `null_action`."
      ),
      call. = FALSE
    )
  }

  current_data <- list(
    X = sweep(X, 2L, colMeans(X), "-"),
    Y = sweep(Y, 2L, colMeans(Y), "-")
  )
  zero_tol <- max(1, sum(current_data$Y^2)) * .Machine$double.eps

  fold_ids <- .fold_resolve_ids(
    n = nrow(current_data$X),
    folds = folds,
    n_folds = n_folds,
    seed = seed
  )

  if (is.null(max_steps)) {
    max_steps <- min(nrow(current_data$X) - 1L,
                     ncol(current_data$X),
                     ncol(current_data$Y),
                     50L)
  }
  max_steps <- as.integer(max(1L, max_steps))

  fit_cache <- new.env(parent = emptyenv())
  get_step_fit <- function(step, data) {
    key <- as.character(step)
    fit <- fit_cache[[key]]
    if (!is.null(fit)) {
      return(fit)
    }
    fit <- adapter$refit(NULL, data)
    fit_cache[[key]] <- fit
    fit
  }

  observed_stat_fn <- function(step, data) {
    .predictive_incremental_gain(adapter, data, fold_ids, k = 1L,
                                 zero_tol = zero_tol)
  }

  null_stat_fn <- function(step, data) {
    fit <- get_step_fit(step, data)
    null_data <- adapter$null_action(fit, data)
    .predictive_validate_payload(null_data, template = data)
    .predictive_incremental_gain(adapter, null_data, fold_ids, k = 1L,
                                 zero_tol = zero_tol)
  }

  deflate_fn <- function(step, data) {
    fit <- get_step_fit(step, data)
    next_data <- adapter$residualize(fit, 1L, data)
    .predictive_validate_payload(next_data, template = data)

    if (sum(next_data$X^2) + sum(next_data$Y^2) <= zero_tol) {
      return(list(
        X = matrix(0, nrow = nrow(data$X), ncol = ncol(data$X)),
        Y = matrix(0, nrow = nrow(data$Y), ncol = ncol(data$Y))
      ))
    }

    next_data
  }

  ladder_result <- ladder_driver(
    observed_stat_fn = observed_stat_fn,
    null_stat_fn     = null_stat_fn,
    deflate_fn       = deflate_fn,
    initial_data     = current_data,
    max_steps        = max_steps,
    B                = B,
    B_total          = B_total,
    batch_size       = batch_size,
    alpha            = alpha,
    seed             = seed
  )

  roots_observed <- vapply(
    ladder_result$step_results,
    function(x) as.numeric(x$observed_stat),
    double(1L)
  )
  n_steps <- length(roots_observed)
  proxy_roots <- rev(seq_len(max(1L, n_steps)))

  selected <- logical(n_steps)
  if (ladder_result$rejected_through >= 1L) {
    selected[seq_len(ladder_result$rejected_through)] <- TRUE
  }

  units <- form_units(
    roots = proxy_roots,
    selected = selected,
    group_near_ties = isTRUE(auto_subspace),
    tie_threshold = tie_threshold
  )

  sr <- ladder_result$step_results
  component_tests <- data.frame(
    step          = vapply(sr, function(x) x$step,          integer(1L)),
    observed_stat = vapply(sr, function(x) x$observed_stat, double(1L)),
    p_value       = vapply(sr, function(x) x$p_value,       double(1L)),
    mc_se         = vapply(sr, function(x) x$mc_se,         double(1L)),
    r             = vapply(sr, function(x) x$r,             integer(1L)),
    B             = vapply(sr, function(x) x$B,             integer(1L)),
    selected      = vapply(sr, function(x) x$selected,      logical(1L)),
    stringsAsFactors = FALSE
  )

  list(
    units           = units,
    component_tests = component_tests,
    roots_observed  = roots_observed,
    ladder_result   = ladder_result,
    labels          = list(
      statistic = "cross-fit incremental held-out predictive gain",
      null      = "row permutation of Y (adapter-supplied null_action)",
      estimand  = "held-out predictive gain"
    )
  )
}

.predictive_resolve_folds <- function(n, folds = NULL, n_folds = 5L, seed = NULL) {
  .fold_resolve_ids(n = n, folds = folds, n_folds = n_folds, seed = seed)
}

.predictive_subset_data <- function(data, idx) {
  .fold_subset_row_aligned(data, idx)
}

.predictive_validate_payload <- function(data, template = NULL) {
  if (!is.list(data) || is.null(data$X) || is.null(data$Y)) {
    stop("Predictive engine payloads must be lists with X and Y.", call. = FALSE)
  }
  if (!is.matrix(data$X) || !is.numeric(data$X) ||
      !is.matrix(data$Y) || !is.numeric(data$Y)) {
    stop("Predictive engine payloads must contain numeric matrices X and Y.",
         call. = FALSE)
  }
  if (nrow(data$X) != nrow(data$Y)) {
    stop("Predictive engine payloads must keep X and Y row counts aligned.",
         call. = FALSE)
  }
  if (any(!is.finite(data$X)) || any(!is.finite(data$Y))) {
    stop("Predictive engine payloads must not contain non-finite values.",
         call. = FALSE)
  }
  if (!is.null(template)) {
    if (nrow(data$X) != nrow(template$X) ||
        ncol(data$X) != ncol(template$X) ||
        nrow(data$Y) != nrow(template$Y) ||
        ncol(data$Y) != ncol(template$Y)) {
      stop(
        "Predictive engine payloads must preserve the current X/Y dimensions.",
        call. = FALSE
      )
    }
  }
  invisible(data)
}

.predictive_total_gain <- function(adapter, data, folds, k, zero_tol) {
  if (k <= 0L) {
    return(0)
  }

  y_centered <- sweep(data$Y, 2L, colMeans(data$Y), "-")
  tss <- sum(y_centered^2)
  if (tss <= zero_tol) {
    return(0)
  }

  fold_press <- .fold_map(
    data = data,
    fold_ids = folds,
    subset_data = .predictive_subset_data,
    fold_fun = function(train_data, test_data, split, adapter, k, y_cols) {
      fit <- adapter$refit(NULL, train_data)
      pred <- adapter$predict_response(fit, test_data, k = k)
      if (is.vector(pred)) {
        pred <- matrix(pred, ncol = y_cols)
      }
      if (!is.matrix(pred) || !is.numeric(pred) ||
          nrow(pred) != length(split$test) || ncol(pred) != y_cols) {
        stop(
          paste0(
            "`predict_response` must return a numeric matrix with dimensions ",
            "length(test_idx) x ncol(Y)."
          ),
          call. = FALSE
        )
      }
      sum((test_data$Y - pred)^2)
    },
    adapter = adapter,
    k = k,
    y_cols = ncol(data$Y)
  )

  press <- sum(unlist(fold_press, use.names = FALSE))

  1 - (press / tss)
}

.predictive_incremental_gain <- function(adapter, data, folds, k = 1L, zero_tol) {
  gain_k <- .predictive_total_gain(adapter, data, folds, k = k, zero_tol = zero_tol)
  gain_prev <- .predictive_total_gain(adapter, data, folds, k = k - 1L, zero_tol = zero_tol)
  max(0, gain_k - gain_prev)
}
