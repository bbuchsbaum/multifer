#' Feature importance from a bootstrap artifact
#'
#' Summarize feature-level importance from an existing
#' \code{\link{bootstrap_fits}()} artifact. This is a descriptive
#' bootstrap consumer: it does not refit models, does not run another
#' bootstrap, and does not report p-values.
#'
#' @param artifact A \code{multifer_bootstrap_artifact}.
#' @param units A \code{multifer_units} table.
#' @param k Unit selection. One of \code{"selected"}, \code{"all"}, a
#'   character vector of unit ids, or an integer vector of unit rows.
#'   \code{"selected"} means \code{units$selected == TRUE}.
#' @param weights One of \code{"equal"}, \code{"root"}, \code{"root2"},
#'   or a numeric vector matching either all units or the selected units.
#' @param roots Optional numeric vector of component roots used for
#'   \code{"root"} and \code{"root2"} weights.
#' @param normalize One of \code{"block"}, \code{"global"}, or
#'   \code{"none"}.
#' @param top_m Optional positive integer. When supplied,
#'   \code{top_m_frequency} is the fraction of bootstrap replicates in
#'   which the feature rank is at most \code{top_m}.
#' @param probs Numeric length-2 interval probabilities.
#' @param scope One of \code{"aggregate"}, \code{"unit"}, or \code{"both"}.
#' @param original_fit Optional original adapter fit. When supplied with
#'   \code{adapter}, the returned \code{estimate} is computed from the
#'   original fit; otherwise it is the bootstrap mean.
#' @param adapter Optional \code{multifer_adapter} used with
#'   \code{original_fit}.
#'
#' @return A \code{multifer_feature_importance} data frame.
#' @export
feature_importance_from_bootstrap <- function(artifact,
                                              units,
                                              k = c("selected", "all"),
                                              weights = c("equal", "root", "root2"),
                                              roots = NULL,
                                              normalize = c("block", "global", "none"),
                                              top_m = NULL,
                                              probs = c(0.025, 0.975),
                                              scope = c("aggregate", "unit", "both"),
                                              original_fit = NULL,
                                              adapter = NULL) {
  if (!inherits(artifact, "multifer_bootstrap_artifact")) {
    stop("`artifact` must be a multifer_bootstrap_artifact.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }
  if (!is.null(original_fit) && !inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be supplied when `original_fit` is supplied.",
         call. = FALSE)
  }
  if (!is.numeric(probs) || length(probs) != 2L || anyNA(probs) ||
      any(probs < 0) || any(probs > 1) || probs[1L] >= probs[2L]) {
    stop("`probs` must be an increasing numeric vector of length 2 in [0, 1].",
         call. = FALSE)
  }
  if (!is.null(top_m)) {
    if (!is.numeric(top_m) || length(top_m) != 1L || is.na(top_m) ||
        top_m != as.integer(top_m) || top_m < 1L) {
      stop("`top_m` must be a positive integer or NULL.", call. = FALSE)
    }
    top_m <- as.integer(top_m)
  }

  normalize <- match.arg(normalize)
  scope <- match.arg(scope)
  unit_idx <- .fi_resolve_units(units, k)
  if (length(unit_idx) == 0L || artifact$R < 1L) {
    return(.fi_empty())
  }

  if (is.null(roots) && !is.null(original_fit) && !is.null(adapter$roots)) {
    roots <- adapter$roots(original_fit)
  }
  alpha <- .fi_unit_weights(weights, units, unit_idx, roots)

  first_loadings <- artifact$reps[[1L]]$aligned_loadings
  template <- .fi_template(first_loadings, units, unit_idx, scope)
  if (nrow(template) == 0L) {
    return(.fi_empty())
  }

  vals <- matrix(NA_real_, nrow = nrow(template), ncol = artifact$R)
  ranks <- matrix(NA_real_, nrow = nrow(template), ncol = artifact$R)
  for (b in seq_len(artifact$R)) {
    v <- .fi_values(
      artifact$reps[[b]]$aligned_loadings,
      template = template,
      units = units,
      unit_idx = unit_idx,
      alpha = alpha,
      normalize = normalize
    )
    vals[, b] <- v
    ranks[, b] <- .fi_ranks(v, template)
  }

  observed_vals <- NULL
  if (!is.null(original_fit)) {
    original_loadings <- stats::setNames(
      lapply(artifact$domains, function(d) adapter$loadings(original_fit, d)),
      artifact$domains
    )
    observed_vals <- .fi_values(
      original_loadings,
      template = template,
      units = units,
      unit_idx = unit_idx,
      alpha = alpha,
      normalize = normalize
    )
  }

  val_stats <- .row_mean_sd(vals)
  rank_stats <- .row_mean_sd(ranks)
  val_ci <- .row_quantile_pair(vals, probs)
  rank_ci <- .row_quantile_pair(ranks, probs)
  estimate <- if (is.null(observed_vals)) val_stats$mean else observed_vals
  eps <- 1e-12
  stability <- 1 / (1 + val_stats$sd / (abs(estimate) + eps))

  top_freq <- rep(NA_real_, nrow(template))
  if (!is.null(top_m)) {
    top_freq <- rowMeans(ranks <= top_m, na.rm = TRUE)
  }

  out <- data.frame(
    scope = template$scope,
    domain = template$domain,
    variable = template$variable,
    unit_id = template$unit_id,
    estimate = as.numeric(estimate),
    lower = as.numeric(val_ci$lower),
    upper = as.numeric(val_ci$upper),
    std_error = as.numeric(val_stats$sd),
    stability = as.numeric(stability),
    rank_mean = as.numeric(rank_stats$mean),
    rank_lower = as.numeric(rank_ci$lower),
    rank_upper = as.numeric(rank_ci$upper),
    top_m_frequency = as.numeric(top_freq),
    weighting = .fi_weight_label(weights),
    normalization = normalize,
    method = "bootstrap_descriptive",
    stringsAsFactors = FALSE
  )
  structure(out, class = c("multifer_feature_importance", "data.frame"))
}

#' Feature importance from an adapter fit
#'
#' Compute observed feature-level importance from a single fitted adapter.
#' The default statistic is the same squared-loading statistic used by
#' \code{\link{feature_importance_from_bootstrap}()}. Adapters can instead
#' provide a family-specific \code{variable_stat()} hook and request
#' \code{statistic = "adapter"}.
#'
#' @param fit Fitted adapter object.
#' @param adapter A \code{multifer_adapter}.
#' @param units A \code{multifer_units} table.
#' @param data Original data passed to \code{adapter$variable_stat()} when
#'   \code{statistic = "adapter"}. The loading-based statistic does not use
#'   this argument.
#' @param k Unit selection. One of \code{"selected"}, \code{"all"}, a
#'   character vector of unit ids, or an integer vector of unit rows.
#' @param weights One of \code{"equal"}, \code{"root"}, \code{"root2"},
#'   or a numeric vector matching either all units or the selected units.
#' @param roots Optional numeric vector of component roots used for
#'   \code{"root"} and \code{"root2"} weights. If omitted and
#'   \code{adapter$roots()} is available, roots are read from \code{fit}.
#' @param normalize One of \code{"block"}, \code{"global"}, or
#'   \code{"none"}.
#' @param scope One of \code{"aggregate"}, \code{"unit"}, or \code{"both"}.
#' @param statistic One of \code{"loadings"} or \code{"adapter"}.
#'
#' @return A \code{multifer_feature_importance} data frame.
#' @export
feature_importance_from_fit <- function(fit,
                                        adapter,
                                        units,
                                        data = NULL,
                                        k = c("selected", "all"),
                                        weights = c("equal", "root", "root2"),
                                        roots = NULL,
                                        normalize = c("block", "global", "none"),
                                        scope = c("aggregate", "unit", "both"),
                                        statistic = c("loadings", "adapter")) {
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }

  normalize <- match.arg(normalize)
  scope <- match.arg(scope)
  statistic <- match.arg(statistic)
  unit_idx <- .fi_resolve_units(units, k)
  if (length(unit_idx) == 0L) {
    return(.fi_empty())
  }

  if (is.null(roots) && is.function(adapter$roots)) {
    roots <- adapter$roots(fit)
  }
  alpha <- .fi_unit_weights(weights, units, unit_idx, roots)
  domains <- .fi_fit_domains(adapter, fit, data)

  if (identical(statistic, "loadings")) {
    if (!is.function(adapter$loadings)) {
      stop("`statistic = \"loadings\"` requires an adapter `loadings` hook.",
           call. = FALSE)
    }
    loadings <- stats::setNames(
      lapply(domains, function(d) adapter$loadings(fit, d)),
      domains
    )
    template <- .fi_template(loadings, units, unit_idx, scope)
    if (nrow(template) == 0L) {
      return(.fi_empty())
    }
    values <- .fi_values(
      loadings,
      template = template,
      units = units,
      unit_idx = unit_idx,
      alpha = alpha,
      normalize = normalize
    )
    return(.fi_observed_output(
      template = template,
      values = values,
      weights = weights,
      normalize = normalize,
      method = "fit_loadings"
    ))
  }

  if (!is.function(adapter$variable_stat)) {
    stop("`statistic = \"adapter\"` requires an adapter `variable_stat` hook.",
         call. = FALSE)
  }
  statistics <- .fi_adapter_stat_matrices(adapter, fit, data, domains, units, unit_idx)
  template <- .fi_template(statistics, units, unit_idx, scope)
  if (nrow(template) == 0L) {
    return(.fi_empty())
  }
  values <- .fi_stat_values(
    statistics = statistics,
    template = template,
    unit_idx = unit_idx,
    alpha = alpha,
    normalize = normalize
  )
  .fi_observed_output(
    template = template,
    values = values,
    weights = weights,
    normalize = normalize,
    method = "fit_adapter_variable_stat"
  )
}

#' Null-calibrated p-values for feature importance
#'
#' Compute Monte Carlo p-values for unsigned aggregate feature importance
#' using an adapter's design-valid \code{null_action}. This is not a
#' bootstrap procedure.
#'
#' @param recipe A \code{multifer_infer_recipe}.
#' @param adapter A \code{multifer_adapter}.
#' @param data Original row-aligned data.
#' @param units A \code{multifer_units} table.
#' @param original_fit Optional original fit. If \code{NULL}, the adapter
#'   is refit on \code{data}.
#' @param B Positive integer number of null draws.
#' @param k,weights,roots,normalize Passed to the same importance
#'   statistic used by \code{\link{feature_importance_from_bootstrap}()}.
#' @param adjust One of \code{"maxT"}, \code{"BH"}, or \code{"none"}.
#' @param seed Optional integer seed.
#'
#' @details This sidecar API uses the adapter's \code{domains()},
#'   \code{loadings()}, \code{null_action()}, and \code{refit()} hooks. When
#'   a geometry's null payload is not already shaped like the adapter's
#'   \code{refit()} input, the adapter can provide a \code{refit_data()} hook
#'   to translate it.
#'
#' @return A \code{multifer_feature_importance_pvalues} data frame.
#' @export
feature_importance_pvalues <- function(recipe,
                                       adapter,
                                       data,
                                       units,
                                       original_fit = NULL,
                                       B = 1000L,
                                       k = c("selected", "all"),
                                       weights = c("equal", "root", "root2"),
                                       roots = NULL,
                                       normalize = c("block", "global", "none"),
                                       adjust = c("maxT", "BH", "none"),
                                       seed = NULL) {
  if (!inherits(recipe, "multifer_infer_recipe")) {
    stop("`recipe` must be a multifer_infer_recipe.", call. = FALSE)
  }
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) ||
      B != as.integer(B) || B < 1L) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  normalize <- match.arg(normalize)
  adjust <- match.arg(adjust)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  plan <- compile_variable_importance_plan(recipe, adapter, data, model = original_fit)
  if (is.null(original_fit)) {
    original_fit <- plan$original_fit()
  }
  if (is.null(roots) && !is.null(adapter$roots)) {
    roots <- adapter$roots(original_fit)
  }

  unit_idx <- .fi_resolve_units(units, k)
  if (length(unit_idx) == 0L) {
    return(.fi_empty_pvalues())
  }
  domains <- plan$domains(original_fit)
  alpha <- .fi_unit_weights(weights, units, unit_idx, roots)
  original_loadings <- stats::setNames(
    lapply(domains, function(d) adapter$loadings(original_fit, d)),
    domains
  )
  template <- .fi_template(original_loadings, units, unit_idx, "aggregate")
  observed <- .fi_values(
    original_loadings,
    template = template,
    units = units,
    unit_idx = unit_idx,
    alpha = alpha,
    normalize = normalize
  )

  null_mat <- matrix(NA_real_, nrow = length(observed), ncol = B)
  for (b in seq_len(B)) {
    payload <- plan$null_payload(original_fit)
    null_fit <- plan$refit_null(original_fit, payload)
    null_loadings <- stats::setNames(
      lapply(domains, function(d) adapter$loadings(null_fit, d)),
      domains
    )
    null_mat[, b] <- .fi_values(
      null_loadings,
      template = template,
      units = units,
      unit_idx = unit_idx,
      alpha = alpha,
      normalize = normalize
    )
  }

  exceed <- rowSums(null_mat >= observed, na.rm = TRUE)
  p <- (1 + exceed) / (B + 1)
  max_null <- apply(null_mat, 2L, max, na.rm = TRUE)
  p_adj <- switch(
    adjust,
    maxT = (1 + vapply(observed, function(x) sum(max_null >= x), integer(1L))) / (B + 1),
    BH = stats::p.adjust(p, method = "BH"),
    none = p
  )
  mc_se <- sqrt(p * (1 - p) / B)
  null_label <- plan$null_label

  out <- data.frame(
    domain = template$domain,
    variable = template$variable,
    statistic = as.numeric(observed),
    p_value = as.numeric(p),
    p_adjusted = as.numeric(p_adj),
    adjustment = adjust,
    B = B,
    mc_uncertainty = as.numeric(mc_se),
    null_label = null_label,
    validity_level = recipe$validity_level,
    stringsAsFactors = FALSE
  )
  structure(out, class = c("multifer_feature_importance_pvalues", "data.frame"))
}

#' @export
print.multifer_feature_importance <- function(x, ...) {
  cat("<multifer_feature_importance>\n")
  cat("  rows:   ", nrow(x), "\n", sep = "")
  if (nrow(x) > 0L) {
    cat("  method: ", unique(x$method)[1L], "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.multifer_feature_importance_pvalues <- function(x, ...) {
  cat("<multifer_feature_importance_pvalues>\n")
  cat("  rows: ", nrow(x), "\n", sep = "")
  if (nrow(x) > 0L) {
    cat("  B:    ", unique(x$B)[1L], "\n", sep = "")
  }
  invisible(x)
}

.fi_empty <- function() {
  structure(
    data.frame(
      scope = character(0), domain = character(0), variable = character(0),
      unit_id = character(0), estimate = numeric(0), lower = numeric(0),
      upper = numeric(0), std_error = numeric(0), stability = numeric(0),
      rank_mean = numeric(0), rank_lower = numeric(0), rank_upper = numeric(0),
      top_m_frequency = numeric(0), weighting = character(0),
      normalization = character(0), method = character(0),
      stringsAsFactors = FALSE
    ),
    class = c("multifer_feature_importance", "data.frame")
  )
}

.fi_empty_pvalues <- function() {
  structure(
    data.frame(
      domain = character(0), variable = character(0), statistic = numeric(0),
      p_value = numeric(0), p_adjusted = numeric(0), adjustment = character(0),
      B = integer(0), mc_uncertainty = numeric(0), null_label = character(0),
      validity_level = character(0), stringsAsFactors = FALSE
    ),
    class = c("multifer_feature_importance_pvalues", "data.frame")
  )
}

.fi_validate_pvalue_hooks <- function(adapter) {
  required <- c("loadings", "null_action", "refit")
  missing <- required[!vapply(required, function(h) is.function(adapter[[h]]), logical(1L))]
  if (length(missing) > 0L) {
    stop(sprintf(
      "feature_importance_pvalues() requires adapter hook(s): %s.",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(adapter)
}

.fi_fit_domains <- function(adapter, fit, data) {
  geom_kind <- "adapter"
  shape_kinds <- adapter$shape_kinds
  if (is.character(shape_kinds) && length(shape_kinds) > 0L) {
    geom_kind <- shape_kinds[[1L]]
  }
  .adapter_domains(adapter, fit = fit, data = data, geom_kind = geom_kind)
}

.fi_observed_output <- function(template, values, weights, normalize, method) {
  ranks <- .fi_ranks(values, template)
  out <- data.frame(
    scope = template$scope,
    domain = template$domain,
    variable = template$variable,
    unit_id = template$unit_id,
    estimate = as.numeric(values),
    lower = NA_real_,
    upper = NA_real_,
    std_error = NA_real_,
    stability = NA_real_,
    rank_mean = as.numeric(ranks),
    rank_lower = as.numeric(ranks),
    rank_upper = as.numeric(ranks),
    top_m_frequency = NA_real_,
    weighting = .fi_weight_label(weights),
    normalization = normalize,
    method = method,
    stringsAsFactors = FALSE
  )
  structure(out, class = c("multifer_feature_importance", "data.frame"))
}

.fi_adapter_stat_matrices <- function(adapter, fit, data, domains, units, unit_idx) {
  members <- attr(units, "members")
  stats::setNames(lapply(domains, function(d) {
    q <- NULL
    vars <- NULL
    for (u in seq_along(unit_idx)) {
      unit_row <- unit_idx[[u]]
      k <- as.integer(members[[unit_row]])
      if (length(k) == 0L) {
        stop("Selected units must contain at least one component member.",
             call. = FALSE)
      }
      stat <- adapter$variable_stat(fit, data, domain = d, k = k)
      if (!is.numeric(stat) || is.matrix(stat) || length(stat) == 0L ||
          any(!is.finite(stat)) || any(stat < 0)) {
        stop(
          "Adapter `variable_stat` must return a non-empty finite non-negative numeric vector.",
          call. = FALSE
        )
      }
      stat_vars <- names(stat)
      if (is.null(stat_vars)) {
        stat_vars <- as.character(seq_along(stat))
      }
      if (is.null(q)) {
        vars <- stat_vars
        q <- matrix(NA_real_, nrow = length(stat), ncol = length(unit_idx))
        rownames(q) <- vars
      } else if (length(stat) != nrow(q) || !identical(stat_vars, vars)) {
        stop(
          "Adapter `variable_stat` must return the same named variables for every selected unit within a domain.",
          call. = FALSE
        )
      }
      q[, u] <- as.numeric(stat)
    }
    q
  }), domains)
}

.fi_stat_values <- function(statistics, template, unit_idx, alpha, normalize) {
  out <- numeric(nrow(template))
  for (d in names(statistics)) {
    q <- statistics[[d]]
    vars <- .fi_vars(q)
    agg <- as.numeric(q %*% alpha / sum(alpha))
    idx <- template$domain == d & template$scope == "aggregate"
    out[idx] <- agg
    idx_unit <- which(template$domain == d & template$scope == "unit")
    if (length(idx_unit) > 0L) {
      for (row in idx_unit) {
        u <- match(template$unit_index[row], unit_idx)
        out[row] <- q[template$variable[row] == vars, u]
      }
    }
  }
  .fi_normalize(out, template, normalize)
}

compile_variable_importance_plan <- function(recipe, adapter, data, model = NULL) {
  .fi_validate_pvalue_hooks(adapter)
  geom_kind <- recipe$shape$geometry$kind
  rel_kind  <- recipe$shape$relation$kind
  refit_data <- .compile_fi_refit_data(adapter, geom_kind, rel_kind, data)
  list(
    geom_kind    = geom_kind,
    rel_kind     = rel_kind,
    domains      = function(fit) {
      .adapter_domains(adapter, fit = fit, data = data, geom_kind = geom_kind)
    },
    original_fit = function() {
      .infer_original_fit(adapter, data, geom_kind, rel_kind, model)
    },
    null_payload = function(fit) adapter$null_action(fit, data),
    refit_null   = function(fit, null_payload) {
      adapter$refit(fit, refit_data(fit, null_payload))
    },
    null_label = .fi_null_label(geom_kind)
  )
}

.compile_fi_refit_data <- function(adapter, geom_kind, rel_kind, data) {
  if (is.function(adapter$refit_data)) {
    return(function(fit, null_payload) {
      adapter$refit_data(fit, null_payload, data)
    })
  }
  .default_fi_refit_data_translator(geom_kind, rel_kind)
}

.default_fi_refit_data_translator <- function(geom_kind, rel_kind) {
  if (identical(geom_kind, "oneblock")) {
    return(function(fit, null_payload) {
      if (is.list(null_payload) && !is.null(null_payload$X)) {
        null_payload$X
      } else {
        null_payload
      }
    })
  }

  if (identical(geom_kind, "cross")) {
    force(rel_kind)
    return(function(fit, null_payload) {
      if (!is.list(null_payload) || is.null(null_payload$X) ||
          is.null(null_payload$Y)) {
        stop(
          "feature_importance_pvalues: cross null payloads must be lists with `X` and `Y`, or the adapter must provide `refit_data`.",
          call. = FALSE
        )
      }
      list(X = null_payload$X, Y = null_payload$Y, relation = rel_kind)
    })
  }

  if (identical(geom_kind, "multiblock")) {
    return(function(fit, null_payload) {
      .validate_multiblock_data(null_payload)
      null_payload
    })
  }

  stop(sprintf(
    "feature_importance_pvalues: geometry '%s' requires the adapter to provide a `refit_data` hook to translate null payloads back to refit input.",
    geom_kind
  ), call. = FALSE)
}

.fi_null_label <- function(geom_kind) {
  switch(
    geom_kind,
    oneblock = "adapter_null_action_oneblock",
    cross = "adapter_null_action_cross",
    multiblock = "adapter_null_action_multiblock",
    geneig = "adapter_null_action_geneig",
    adapter = "adapter_null_action",
    "adapter_null_action"
  )
}

.fi_resolve_units <- function(units, k) {
  n <- nrow(units)
  if (n == 0L) return(integer(0))
  if (is.character(k) && length(k) > 1L &&
      any(k %in% c("selected", "all"))) {
    k <- match.arg(k, c("selected", "all"))
  }
  if (is.character(k) && length(k) == 1L && k == "selected") {
    return(which(units$selected))
  }
  if (is.character(k) && length(k) == 1L && k == "all") {
    return(seq_len(n))
  }
  if (is.character(k)) {
    idx <- match(k, units$unit_id)
    if (anyNA(idx)) {
      stop("Unknown unit id in `k`.", call. = FALSE)
    }
    return(as.integer(idx))
  }
  if (is.numeric(k) && all(k == as.integer(k), na.rm = TRUE)) {
    idx <- as.integer(k)
    if (anyNA(idx) || any(idx < 1L) || any(idx > n)) {
      stop("Numeric `k` values must be valid unit row indices.", call. = FALSE)
    }
    return(idx)
  }
  stop("`k` must be \"selected\", \"all\", unit ids, or unit row indices.",
       call. = FALSE)
}

.fi_weight_label <- function(weights) {
  if (is.numeric(weights)) "custom" else as.character(weights)[1L]
}

.fi_unit_weights <- function(weights, units, unit_idx, roots) {
  if (is.numeric(weights)) {
    if (length(weights) == nrow(units)) {
      out <- weights[unit_idx]
    } else if (length(weights) == length(unit_idx)) {
      out <- weights
    } else {
      stop("Numeric `weights` must match all units or the selected units.",
           call. = FALSE)
    }
  } else {
    weights <- match.arg(weights, c("equal", "root", "root2"))
    if (weights == "equal") {
      out <- rep(1, length(unit_idx))
    } else {
      if (is.null(roots) || !is.numeric(roots)) {
        stop("`roots` must be supplied for root-based weights.", call. = FALSE)
      }
      members <- attr(units, "members")
      out <- vapply(unit_idx, function(i) {
        m <- as.integer(members[[i]])
        if (any(m > length(roots))) {
          stop("`roots` is shorter than the selected unit members.",
               call. = FALSE)
        }
        if (weights == "root") sum(roots[m]) else sum(roots[m]^2)
      }, numeric(1L))
    }
  }
  if (any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
    stop("Resolved `weights` must be finite, non-negative, and have positive sum.",
         call. = FALSE)
  }
  as.numeric(out)
}

.fi_template <- function(loadings, units, unit_idx, scope) {
  domains <- names(loadings)
  include_aggregate <- scope %in% c("aggregate", "both")
  include_unit <- scope %in% c("unit", "both")
  rows <- vector("list", 0L)
  pos <- 1L
  for (d in domains) {
    L <- loadings[[d]]
    vars <- rownames(L)
    if (is.null(vars)) vars <- as.character(seq_len(nrow(L)))
    if (include_aggregate) {
      rows[[pos]] <- data.frame(
        scope = "aggregate", domain = d, variable = vars,
        unit_id = NA_character_, unit_index = NA_integer_,
        stringsAsFactors = FALSE
      )
      pos <- pos + 1L
    }
    if (include_unit) {
      for (i in unit_idx) {
        rows[[pos]] <- data.frame(
          scope = "unit", domain = d, variable = vars,
          unit_id = units$unit_id[i], unit_index = i,
          stringsAsFactors = FALSE
        )
        pos <- pos + 1L
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(scope = character(0), domain = character(0),
                      variable = character(0), unit_id = character(0),
                      unit_index = integer(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

.fi_values <- function(loadings, template, units, unit_idx, alpha, normalize) {
  members <- attr(units, "members")
  q_by_domain <- lapply(names(loadings), function(d) {
    L <- loadings[[d]]
    q <- matrix(NA_real_, nrow = nrow(L), ncol = length(unit_idx))
    for (u in seq_along(unit_idx)) {
      m <- as.integer(members[[unit_idx[u]]])
      if (length(m) == 0L) {
        stop("Selected units must contain at least one component member.",
             call. = FALSE)
      }
      if (any(m < 1L) || any(m > ncol(L))) {
        stop(
          sprintf(
            "Selected unit '%s' has component members outside the available loading columns for domain '%s'.",
            units$unit_id[unit_idx[u]], d
          ),
          call. = FALSE
        )
      }
      q[, u] <- rowSums(L[, m, drop = FALSE]^2)
    }
    q
  })
  names(q_by_domain) <- names(loadings)

  out <- numeric(nrow(template))
  for (d in names(q_by_domain)) {
    q <- q_by_domain[[d]]
    agg <- as.numeric(q %*% alpha / sum(alpha))
    idx <- template$domain == d & template$scope == "aggregate"
    out[idx] <- agg
    idx_unit <- which(template$domain == d & template$scope == "unit")
    if (length(idx_unit) > 0L) {
      for (row in idx_unit) {
        u <- match(template$unit_index[row], unit_idx)
        out[row] <- q[template$variable[row] == .fi_vars(loadings[[d]]), u]
      }
    }
  }
  .fi_normalize(out, template, normalize)
}

.fi_vars <- function(L) {
  vars <- rownames(L)
  if (is.null(vars)) vars <- as.character(seq_len(nrow(L)))
  vars
}

.fi_normalize <- function(x, template, normalize) {
  if (normalize == "none" || length(x) == 0L) return(x)
  out <- x
  if (normalize == "global") {
    total <- sum(out, na.rm = TRUE)
    if (is.finite(total) && total > 0) out <- out / total
    return(out)
  }
  groups <- paste(template$scope, template$domain, template$unit_id, sep = "\r")
  for (g in unique(groups)) {
    idx <- groups == g
    total <- sum(out[idx], na.rm = TRUE)
    if (is.finite(total) && total > 0) {
      out[idx] <- out[idx] / total
    }
  }
  out
}

.fi_ranks <- function(values, template) {
  out <- rep(NA_real_, length(values))
  groups <- paste(template$scope, template$domain, template$unit_id, sep = "\r")
  for (g in unique(groups)) {
    idx <- which(groups == g)
    out[idx] <- rank(-values[idx], ties.method = "average", na.last = "keep")
  }
  out
}
