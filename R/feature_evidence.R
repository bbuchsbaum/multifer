#' Build the `feature_evidence` table
#'
#' `feature_evidence` is a sidecar table for feature-to-unit evidence and
#' aggregate feature importance. It is deliberately separate from
#' `variable_stability`: bootstrap ratios and null-calibrated p-values are
#' evidence summaries, not stability scores.
#'
#' @param scope Character vector. Typical values are `"unit"` and
#'   `"aggregate"`.
#' @param domain Character vector of domains or block labels.
#' @param feature Character vector of feature identifiers.
#' @param unit_id Character vector of unit ids, or `NA` for aggregate rows.
#' @param unit_set Character vector describing the unit set.
#' @param statistic Character vector naming the statistic.
#' @param estimate,bootstrap_mean,bias,std_error,lower,upper,z Numeric columns.
#' @param z_label Character vector describing `z`.
#' @param p_value,p_adjusted Numeric p-value columns. Use `NA` when no null
#'   calibration was run.
#' @param adjustment Character vector naming the multiplicity adjustment.
#' @param calibration Character vector, one of `"bootstrap"`, `"null_refit"`,
#'   `"analytic"`, or `"none"`.
#' @param interval_method Character vector, one of `"percentile"`, `"normal"`,
#'   or `"none"`.
#' @param null_label Character vector describing the null.
#' @param method Character vector describing the evidence method.
#' @param orientation Character vector, one of `"signed"`, `"unsigned"`,
#'   `"subspace_norm"`, or `"adapter_invariant"`.
#' @param identifiable Logical vector.
#' @param validity_level Character vector, one of `"exact"`, `"conditional"`,
#'   `"asymptotic"`, or `"heuristic"`.
#' @param warnings Character vector with row-level caveats.
#'
#' @return A data frame with class `multifer_feature_evidence`.
#' @export
infer_feature_evidence <- function(scope = character(0),
                                   domain = character(0),
                                   feature = character(0),
                                   unit_id = character(0),
                                   unit_set = character(0),
                                   statistic = character(0),
                                   estimate = numeric(0),
                                   bootstrap_mean = numeric(0),
                                   bias = numeric(0),
                                   std_error = numeric(0),
                                   lower = numeric(0),
                                   upper = numeric(0),
                                   z = numeric(0),
                                   z_label = character(0),
                                   p_value = numeric(0),
                                   p_adjusted = numeric(0),
                                   adjustment = character(0),
                                   calibration = character(0),
                                   interval_method = character(0),
                                   null_label = character(0),
                                   method = character(0),
                                   orientation = character(0),
                                   identifiable = logical(0),
                                   validity_level = character(0),
                                   warnings = character(0)) {
  n <- length(scope)
  lens <- c(
    length(domain), length(feature), length(unit_id), length(unit_set),
    length(statistic), length(estimate), length(bootstrap_mean), length(bias),
    length(std_error), length(lower), length(upper), length(z),
    length(z_label), length(p_value), length(p_adjusted), length(adjustment),
    length(calibration), length(interval_method), length(null_label),
    length(method), length(orientation), length(identifiable),
    length(validity_level), length(warnings)
  )
  if (!all(lens == n)) {
    stop("all `infer_feature_evidence()` arguments must have the same length.",
         call. = FALSE)
  }
  if (n > 0L) {
    .fe_validate_values(scope, c("unit", "aggregate", "subspace", "score", "block"),
                        "scope")
    .fe_validate_values(orientation,
                        c("signed", "unsigned", "subspace_norm",
                          "adapter_invariant"),
                        "orientation")
    .fe_validate_values(calibration,
                        c("bootstrap", "null_refit", "analytic", "none"),
                        "calibration")
    .fe_validate_values(interval_method, c("percentile", "normal", "none"),
                        "interval_method")
    .fe_validate_values(validity_level,
                        c("exact", "conditional", "asymptotic", "heuristic"),
                        "validity_level")
  }
  structure(
    data.frame(
      scope = as.character(scope),
      domain = as.character(domain),
      feature = as.character(feature),
      unit_id = as.character(unit_id),
      unit_set = as.character(unit_set),
      statistic = as.character(statistic),
      estimate = as.numeric(estimate),
      bootstrap_mean = as.numeric(bootstrap_mean),
      bias = as.numeric(bias),
      std_error = as.numeric(std_error),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      z = as.numeric(z),
      z_label = as.character(z_label),
      p_value = as.numeric(p_value),
      p_adjusted = as.numeric(p_adjusted),
      adjustment = as.character(adjustment),
      calibration = as.character(calibration),
      interval_method = as.character(interval_method),
      null_label = as.character(null_label),
      method = as.character(method),
      orientation = as.character(orientation),
      identifiable = as.logical(identifiable),
      validity_level = as.character(validity_level),
      warnings = as.character(warnings),
      stringsAsFactors = FALSE
    ),
    class = c("multifer_feature_evidence", "data.frame")
  )
}

#' Feature evidence from a bootstrap artifact
#'
#' Summarize feature-to-unit evidence from an existing
#' `multifer_bootstrap_artifact`. Component-wise signed rows report bootstrap
#' ratios, not p-values. Null-calibrated p-values require
#' [feature_evidence_pvalues()].
#'
#' @param artifact A `multifer_bootstrap_artifact`.
#' @param units A `multifer_units` table.
#' @param original_fit Optional original adapter fit. When supplied with
#'   `adapter`, estimates are computed from the original fit; otherwise
#'   bootstrap means are used as estimates.
#' @param adapter Optional `multifer_adapter`.
#' @param statistic One of `"loading"`, `"squared_loading"`, or
#'   `"subspace_norm"`.
#' @param scope One of `"unit"`, `"aggregate"`, or `"both"`.
#' @param k Unit selection. One of `"selected"`, `"all"`, unit ids, or unit
#'   row indices.
#' @param orientation One of `"auto"`, `"signed"`, `"unsigned"`, or
#'   `"subspace_norm"`.
#' @param weights,roots Weighting controls for aggregate evidence.
#' @param normalize One of `"block"`, `"global"`, or `"none"`.
#' @param probs Numeric length-2 interval probabilities.
#'
#' @return A `multifer_feature_evidence` data frame.
#' @export
feature_evidence_from_bootstrap <- function(artifact,
                                            units,
                                            original_fit = NULL,
                                            adapter = NULL,
                                            statistic = c("loading",
                                                          "squared_loading",
                                                          "subspace_norm"),
                                            scope = c("unit", "aggregate", "both"),
                                            k = c("selected", "all"),
                                            orientation = c("auto", "signed",
                                                            "unsigned",
                                                            "subspace_norm"),
                                            weights = c("equal", "root", "root2"),
                                            roots = NULL,
                                            normalize = c("block", "global",
                                                          "none"),
                                            probs = c(0.025, 0.975)) {
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

  statistic <- match.arg(statistic)
  scope <- match.arg(scope)
  orientation <- match.arg(orientation)
  normalize <- match.arg(normalize)

  unit_idx <- .fi_resolve_units(units, k)
  if (length(unit_idx) == 0L || artifact$R < 1L || length(artifact$reps) == 0L) {
    return(.fe_empty())
  }

  domains <- artifact$domains
  if (is.null(domains)) {
    domains <- names(artifact$reps[[1L]]$aligned_loadings)
  }
  first_loadings <- artifact$reps[[1L]]$aligned_loadings[domains]
  template <- .fe_template(first_loadings, units, unit_idx, scope)
  if (nrow(template) == 0L) {
    return(.fe_empty())
  }

  if (is.null(roots) && !is.null(original_fit) && is.function(adapter$roots)) {
    roots <- adapter$roots(original_fit)
  }
  alpha <- .fi_unit_weights(weights, units, unit_idx, roots)
  spec <- .fe_resolve_template_spec(
    template = template,
    units = units,
    unit_idx = unit_idx,
    statistic = statistic,
    orientation = orientation
  )

  vals <- matrix(NA_real_, nrow = nrow(template), ncol = artifact$R)
  for (b in seq_len(artifact$R)) {
    vals[, b] <- .fe_values(
      loadings = artifact$reps[[b]]$aligned_loadings,
      template = template,
      units = units,
      unit_idx = unit_idx,
      alpha = alpha,
      statistic = spec$statistic,
      normalize = normalize
    )
  }

  observed <- NULL
  if (!is.null(original_fit)) {
    original_loadings <- stats::setNames(
      lapply(domains, function(d) adapter$loadings(original_fit, d)),
      domains
    )
    observed <- .fe_values(
      loadings = original_loadings,
      template = template,
      units = units,
      unit_idx = unit_idx,
      alpha = alpha,
      statistic = spec$statistic,
      normalize = normalize
    )
  }

  stats <- .row_mean_sd(vals)
  ci <- .row_quantile_pair(vals, probs)
  estimate <- if (is.null(observed)) stats$mean else observed
  se <- stats$sd
  z <- estimate / se
  z[!is.finite(z)] <- NA_real_
  bias <- stats$mean - estimate

  infer_feature_evidence(
    scope = template$scope,
    domain = template$domain,
    feature = template$feature,
    unit_id = template$unit_id,
    unit_set = spec$unit_set,
    statistic = spec$statistic,
    estimate = estimate,
    bootstrap_mean = stats$mean,
    bias = bias,
    std_error = se,
    lower = ci$lower,
    upper = ci$upper,
    z = z,
    z_label = rep("estimate / bootstrap_se", nrow(template)),
    p_value = rep(NA_real_, nrow(template)),
    p_adjusted = rep(NA_real_, nrow(template)),
    adjustment = rep("none", nrow(template)),
    calibration = rep("bootstrap", nrow(template)),
    interval_method = rep("percentile", nrow(template)),
    null_label = rep(NA_character_, nrow(template)),
    method = .fe_method_label(spec$statistic),
    orientation = spec$orientation,
    identifiable = spec$identifiable,
    validity_level = rep("conditional", nrow(template)),
    warnings = spec$warnings
  )
}

#' Null-calibrated p-values for feature evidence
#'
#' Compute feature evidence with p-values from an adapter null action. Ordinary
#' bootstrap ratios are intentionally not converted into p-values.
#'
#' @param recipe A `multifer_infer_recipe`.
#' @param adapter A `multifer_adapter`.
#' @param data Original data.
#' @param units A `multifer_units` table.
#' @param original_fit Optional original fit.
#' @param statistic One of `"squared_loading"`, `"subspace_norm"`, or
#'   `"adapter"` for an adapter-owned nonnegative statistic.
#' @param scope One of `"unit"` or `"aggregate"`.
#' @param B Positive integer number of null draws.
#' @param k,weights,roots,normalize Passed to the statistic resolver.
#' @param adjust One of `"BH"`, `"holm"`, `"none"`, or `"maxT"`.
#' @param seed Optional integer seed.
#'
#' @return A `multifer_feature_evidence` data frame.
#' @export
feature_evidence_pvalues <- function(recipe,
                                     adapter,
                                     data,
                                     units,
                                     original_fit = NULL,
                                     statistic = c("squared_loading",
                                                   "subspace_norm",
                                                   "adapter"),
                                     scope = c("unit", "aggregate"),
                                     B = 1000L,
                                     k = c("selected", "all"),
                                     weights = c("equal", "root", "root2"),
                                     roots = NULL,
                                     normalize = c("block", "global", "none"),
                                     adjust = c("BH", "holm", "none", "maxT"),
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
  statistic <- match.arg(statistic)
  scope <- match.arg(scope)
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

  plan <- compile_variable_importance_plan(recipe, adapter, data,
                                           model = original_fit)
  if (is.null(original_fit)) {
    original_fit <- plan$original_fit()
  }
  if (is.null(roots) && is.function(adapter$roots)) {
    roots <- adapter$roots(original_fit)
  }

  unit_idx <- .fi_resolve_units(units, k)
  if (length(unit_idx) == 0L) {
    return(.fe_empty())
  }
  domains <- plan$domains(original_fit)
  alpha <- .fi_unit_weights(weights, units, unit_idx, roots)
  if (identical(statistic, "adapter")) {
    observed_statistics <- .fe_adapter_stat_matrices(
      adapter = adapter,
      fit = original_fit,
      data = data,
      domains = domains,
      units = units,
      unit_idx = unit_idx,
      statistic = statistic,
      orientation = "adapter_invariant"
    )
    template <- .fe_template(observed_statistics, units, unit_idx, scope)
    spec <- .fe_resolve_template_spec(
      template = template,
      units = units,
      unit_idx = unit_idx,
      statistic = statistic,
      orientation = "adapter_invariant"
    )
    observed <- .fe_stat_values(
      statistics = observed_statistics,
      template = template,
      unit_idx = unit_idx,
      alpha = alpha,
      normalize = normalize
    )
    stat_for_null <- function(fit, stat_data) {
      statistics <- .fe_adapter_stat_matrices(
        adapter = adapter,
        fit = fit,
        data = stat_data,
        domains = domains,
        units = units,
        unit_idx = unit_idx,
        statistic = statistic,
        orientation = "adapter_invariant"
      )
      .fe_stat_values(statistics, template, unit_idx, alpha, normalize)
    }
  } else {
    original_loadings <- stats::setNames(
      lapply(domains, function(d) adapter$loadings(original_fit, d)),
      domains
    )
    template <- .fe_template(original_loadings, units, unit_idx, scope)
    spec <- .fe_resolve_template_spec(
      template = template,
      units = units,
      unit_idx = unit_idx,
      statistic = statistic,
      orientation = "auto"
    )
    observed <- .fe_values(
      original_loadings, template, units, unit_idx, alpha,
      statistic = spec$statistic, normalize = normalize
    )
    stat_for_null <- function(fit, stat_data) {
      null_loadings <- stats::setNames(
        lapply(domains, function(d) adapter$loadings(fit, d)),
        domains
      )
      .fe_values(null_loadings, template, units, unit_idx, alpha,
                 statistic = spec$statistic, normalize = normalize)
    }
  }

  null_mat <- matrix(NA_real_, nrow = length(observed), ncol = B)
  for (b in seq_len(B)) {
    payload <- plan$null_payload(original_fit)
    null_fit <- plan$refit_null(original_fit, payload)
    null_mat[, b] <- stat_for_null(null_fit, data)
  }

  exceed <- rowSums(null_mat >= observed, na.rm = TRUE)
  p <- (1 + exceed) / (B + 1)
  p_adj <- switch(
    adjust,
    maxT = {
      max_null <- apply(null_mat, 2L, max, na.rm = TRUE)
      (1 + vapply(observed, function(x) sum(max_null >= x), integer(1L))) /
        (B + 1)
    },
    BH = stats::p.adjust(p, method = "BH"),
    holm = stats::p.adjust(p, method = "holm"),
    none = p
  )

  infer_feature_evidence(
    scope = template$scope,
    domain = template$domain,
    feature = template$feature,
    unit_id = template$unit_id,
    unit_set = spec$unit_set,
    statistic = spec$statistic,
    estimate = observed,
    bootstrap_mean = rep(NA_real_, length(observed)),
    bias = rep(NA_real_, length(observed)),
    std_error = rep(NA_real_, length(observed)),
    lower = rep(NA_real_, length(observed)),
    upper = rep(NA_real_, length(observed)),
    z = rep(NA_real_, length(observed)),
    z_label = rep("none", length(observed)),
    p_value = p,
    p_adjusted = p_adj,
    adjustment = rep(adjust, length(observed)),
    calibration = rep("null_refit", length(observed)),
    interval_method = rep("none", length(observed)),
    null_label = rep(plan$null_label, length(observed)),
    method = .fe_method_label(spec$statistic),
    orientation = spec$orientation,
    identifiable = spec$identifiable,
    validity_level = rep(recipe$validity_level, length(observed)),
    warnings = spec$warnings
  )
}

#' Feature evidence from an adapter-native action
#'
#' Calls an adapter's optional `feature_evidence_action` hook and normalizes
#' the result to the `multifer_feature_evidence` class. This is the extension
#' point for method-native procedures such as conditional PCA subspace
#' bootstrap, PLSC bootstrap ratios, VIP shortcuts, or other family-specific
#' evidence calculations.
#'
#' @param adapter A `multifer_adapter`.
#' @param fit Original adapter fit.
#' @param data Original data.
#' @param units A `multifer_units` table.
#' @param design Optional design object passed to the hook.
#' @param statistic,orientation Evidence request labels passed to the hook.
#' @param R Optional number of bootstrap/evidence draws.
#' @param seed Optional integer seed.
#' @param ... Additional hook arguments.
#'
#' @return A `multifer_feature_evidence` data frame.
#' @export
feature_evidence_from_adapter <- function(adapter,
                                          fit,
                                          data,
                                          units,
                                          design = NULL,
                                          statistic = "adapter",
                                          orientation = "adapter_invariant",
                                          R = NULL,
                                          seed = NULL,
                                          ...) {
  if (!inherits(adapter, "multifer_adapter")) {
    stop("`adapter` must be a multifer_adapter.", call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }
  if (!is.function(adapter$feature_evidence_action)) {
    stop("`adapter` must provide a `feature_evidence_action` hook.",
         call. = FALSE)
  }

  out <- adapter$feature_evidence_action(
    fit,
    data = data,
    units = units,
    design = design,
    statistic = statistic,
    orientation = orientation,
    R = R,
    seed = seed,
    ...
  )
  .normalize_feature_evidence_action(out)
}

#' @export
print.multifer_feature_evidence <- function(x, ...) {
  cat("<multifer_feature_evidence>\n")
  cat("  rows: ", nrow(x), "\n", sep = "")
  if (nrow(x) > 0L) {
    cat("  calibration: ", unique(x$calibration)[1L], "\n", sep = "")
  }
  invisible(x)
}

.normalize_feature_evidence_action <- function(out) {
  if (inherits(out, "multifer_feature_evidence")) {
    return(out)
  }
  if (!is.data.frame(out)) {
    stop("`feature_evidence_action()` must return a data frame or multifer_feature_evidence.",
         call. = FALSE)
  }
  required <- names(infer_feature_evidence())
  missing <- setdiff(required, names(out))
  if (length(missing) > 0L) {
    stop(sprintf(
      "`feature_evidence_action()` output is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  do.call(infer_feature_evidence, out[required])
}

.fe_empty <- function() {
  infer_feature_evidence()
}

.fe_validate_values <- function(x, allowed, what) {
  bad <- setdiff(stats::na.omit(unique(x)), allowed)
  if (length(bad) > 0L) {
    stop(sprintf("unknown %s value(s): %s", what, paste(bad, collapse = ", ")),
         call. = FALSE)
  }
}

.fe_template <- function(loadings, units, unit_idx, scope) {
  include_aggregate <- scope %in% c("aggregate", "both")
  include_unit <- scope %in% c("unit", "both")
  rows <- vector("list", 0L)
  pos <- 1L
  for (d in names(loadings)) {
    L <- loadings[[d]]
    vars <- .fi_vars(L)
    if (include_aggregate) {
      rows[[pos]] <- data.frame(
        scope = "aggregate",
        domain = d,
        feature = vars,
        unit_id = NA_character_,
        unit_index = NA_integer_,
        stringsAsFactors = FALSE
      )
      pos <- pos + 1L
    }
    if (include_unit) {
      for (i in unit_idx) {
        rows[[pos]] <- data.frame(
          scope = "unit",
          domain = d,
          feature = vars,
          unit_id = units$unit_id[i],
          unit_index = i,
          stringsAsFactors = FALSE
        )
        pos <- pos + 1L
      }
    }
  }
  if (length(rows) == 0L) {
    return(data.frame(scope = character(0), domain = character(0),
                      feature = character(0), unit_id = character(0),
                      unit_index = integer(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

.fe_resolve_template_spec <- function(template,
                                      units,
                                      unit_idx,
                                      statistic,
                                      orientation) {
  members <- attr(units, "members")
  n <- nrow(template)
  stat <- rep(statistic, n)
  orient <- rep(NA_character_, n)
  identifiable <- rep(FALSE, n)
  unit_set <- rep(.fe_unit_set_label(units, unit_idx), n)
  warnings <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    if (template$scope[i] == "aggregate") {
      if (identical(statistic, "loading")) {
        stop("`statistic = \"loading\"` is not defined for aggregate scope.",
             call. = FALSE)
      }
      orient[i] <- if (identical(statistic, "adapter")) {
        "adapter_invariant"
      } else if (identical(statistic, "subspace_norm")) {
        "subspace_norm"
      } else {
        "unsigned"
      }
      identifiable[i] <- FALSE
      next
    }

    u <- template$unit_index[i]
    m <- as.integer(members[[u]])
    is_singleton <- length(m) == 1L && isTRUE(units$identifiable[u])
    identifiable[i] <- is_singleton

    if (identical(statistic, "loading")) {
      if (orientation == "signed" && !is_singleton) {
        stop(sprintf(
          "Signed feature evidence is not defined for non-identifiable unit '%s'. Use orientation = \"subspace_norm\" or \"auto\".",
          units$unit_id[u]
        ), call. = FALSE)
      }
      if (orientation == "auto" && !is_singleton) {
        stat[i] <- "subspace_norm"
        orient[i] <- "subspace_norm"
        warnings[i] <- "signed loading not identifiable; used subspace_norm"
      } else {
        orient[i] <- "signed"
      }
    } else if (identical(statistic, "squared_loading")) {
      orient[i] <- if (length(m) > 1L || !isTRUE(units$identifiable[u])) {
        "subspace_norm"
      } else {
        "unsigned"
      }
    } else if (identical(statistic, "subspace_norm")) {
      orient[i] <- "subspace_norm"
    } else if (identical(statistic, "adapter")) {
      orient[i] <- "adapter_invariant"
    }
  }

  list(
    statistic = stat,
    orientation = orient,
    identifiable = identifiable,
    unit_set = unit_set,
    warnings = warnings
  )
}

.fe_values <- function(loadings,
                       template,
                       units,
                       unit_idx,
                       alpha,
                       statistic,
                       normalize) {
  members <- attr(units, "members")
  out <- rep(NA_real_, nrow(template))

  for (d in names(loadings)) {
    L <- loadings[[d]]
    vars <- .fi_vars(L)
    q <- matrix(NA_real_, nrow = nrow(L), ncol = length(unit_idx))
    for (u in seq_along(unit_idx)) {
      m <- as.integer(members[[unit_idx[u]]])
      .fe_check_members(m, L, units$unit_id[unit_idx[u]], d)
      q[, u] <- rowSums(L[, m, drop = FALSE]^2)
    }

    idx <- which(template$domain == d & template$scope == "aggregate")
    if (length(idx) > 0L) {
      agg <- as.numeric(q %*% alpha / sum(alpha))
      out[idx] <- agg[match(template$feature[idx], vars)]
    }

    idx <- which(template$domain == d & template$scope == "unit")
    for (row in idx) {
      u_pos <- match(template$unit_index[row], unit_idx)
      m <- as.integer(members[[template$unit_index[row]]])
      .fe_check_members(m, L, template$unit_id[row], d)
      feature_row <- match(template$feature[row], vars)
      st <- statistic[row]
      out[row] <- switch(
        st,
        loading = L[feature_row, m[1L]],
        squared_loading = sum(L[feature_row, m]^2),
        subspace_norm = sqrt(sum(L[feature_row, m]^2)),
        stop(sprintf("unknown feature statistic: %s", st), call. = FALSE)
      )
      if (identical(st, "squared_loading") && length(m) > 1L) {
        out[row] <- q[feature_row, u_pos]
      }
    }
  }

  if (any(template$scope == "aggregate") || any(statistic != "loading")) {
    norm_template <- template
    names(norm_template)[names(norm_template) == "feature"] <- "variable"
    out <- .fi_normalize(out, norm_template, normalize)
  }
  out
}

.fe_adapter_stat_matrices <- function(adapter,
                                      fit,
                                      data,
                                      domains,
                                      units,
                                      unit_idx,
                                      statistic,
                                      orientation) {
  if (!is.function(adapter$feature_stat) && !is.function(adapter$variable_stat)) {
    stop("`statistic = \"adapter\"` requires adapter `feature_stat` or `variable_stat`.",
         call. = FALSE)
  }
  members <- attr(units, "members")
  stats::setNames(lapply(domains, function(d) {
    q <- NULL
    vars <- NULL
    for (u in seq_along(unit_idx)) {
      unit_row <- unit_idx[[u]]
      m <- as.integer(members[[unit_row]])
      if (length(m) == 0L) {
        stop("Selected units must contain at least one component member.",
             call. = FALSE)
      }
      stat <- if (is.function(adapter$feature_stat)) {
        adapter$feature_stat(
          fit,
          data = data,
          domain = d,
          unit_id = units$unit_id[unit_row],
          members = m,
          statistic = statistic,
          orientation = orientation,
          scope = "unit"
        )
      } else {
        adapter$variable_stat(fit, data, domain = d, k = m)
      }
      if (!is.numeric(stat) || is.matrix(stat) || length(stat) == 0L ||
          any(!is.finite(stat)) || any(stat < 0)) {
        stop(
          "Adapter feature statistic must return a non-empty finite non-negative numeric vector.",
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
          "Adapter feature statistic must return the same named features for every selected unit within a domain.",
          call. = FALSE
        )
      }
      q[, u] <- as.numeric(stat)
    }
    q
  }), domains)
}

.fe_stat_values <- function(statistics, template, unit_idx, alpha, normalize) {
  stat_template <- template
  names(stat_template)[names(stat_template) == "feature"] <- "variable"
  .fi_stat_values(
    statistics = statistics,
    template = stat_template,
    unit_idx = unit_idx,
    alpha = alpha,
    normalize = normalize
  )
}

.fe_check_members <- function(m, L, unit_id, domain) {
  if (length(m) == 0L) {
    stop("Selected units must contain at least one component member.",
         call. = FALSE)
  }
  if (any(m < 1L) || any(m > ncol(L))) {
    stop(
      sprintf(
        "Selected unit '%s' has component members outside the available loading columns for domain '%s'.",
        unit_id, domain
      ),
      call. = FALSE
    )
  }
}

.fe_unit_set_label <- function(units, unit_idx) {
  if (length(unit_idx) == 0L) {
    return("")
  }
  if (identical(unit_idx, which(units$selected))) {
    return("selected")
  }
  if (identical(unit_idx, seq_len(nrow(units)))) {
    return("all")
  }
  paste(units$unit_id[unit_idx], collapse = "+")
}

.fe_method_label <- function(statistic) {
  ifelse(
    statistic == "loading",
    "aligned_bootstrap_loading",
    ifelse(
      statistic == "squared_loading",
      "aligned_bootstrap_squared_loading",
      ifelse(
        statistic == "adapter",
        "adapter_feature_stat",
        "aligned_bootstrap_subspace_norm"
      )
    )
  )
}
