# helper-stub-adapters.R
#
# Testthat helper file -- automatically sourced before any test runs.
# Provides stub hook functions and factory functions that build minimal
# adapters and canned infer_result objects for contract-proof tests.
#
# Design notes:
#   - Hook bodies are trivially minimal (return empty structures) so they
#     satisfy the registration-time existence check in R/adapter.R.
#   - Factory functions return adapter objects; they do NOT register or
#     clear the registry -- that is the caller's responsibility.
#   - make_stub_cross_adapter() accepts a `relations` argument so 9qu.9
#     can build an adapter that declares two relations, triggering the
#     strict-dispatch ambiguity error when relation is omitted.

# ---------------------------------------------------------------------------
# Stub hook functions
# ---------------------------------------------------------------------------

# Returns an empty numeric vector as a placeholder for latent roots.
stub_hook_roots <- function(x, ...) numeric(0)

# Returns an empty matrix as a placeholder for scores.
stub_hook_scores <- function(x, domain = NULL, ...) matrix(0, 0, 0)

# Returns an empty matrix as a placeholder for loadings.
stub_hook_loadings <- function(x, domain = NULL, ...) matrix(0, 0, 0)

# Returns x unchanged (truncation no-op).
stub_hook_truncate <- function(x, k, ...) x

# Returns data unchanged (residualization no-op).
stub_hook_residualize <- function(x, k, data, ...) data

# Returns x unchanged (refit no-op).
stub_hook_refit <- function(x, new_data, ...) x

# Returns data unchanged (null action no-op).
stub_hook_null_action <- function(x, data, ...) data

# Returns scalar 0 as a placeholder component test statistic.
stub_hook_component_stat <- function(x, data, k, ...) 0

# Returns scalar 0 as a placeholder variable-level statistic.
stub_hook_variable_stat <- function(x, data, domain, k, ...) 0

# Returns scalar 0 as a placeholder score-level statistic.
stub_hook_score_stat <- function(x, data, domain, k, ...) 0

# Returns replicate loadings unchanged (alignment no-op).
stub_hook_align <- function(reference_loadings, replicate_loadings, ...) {
  replicate_loadings
}

# ---------------------------------------------------------------------------
# Factory: oneblock stub adapter (PCA-like)
# ---------------------------------------------------------------------------

# Builds a stub adapter that declares the four v1 targets for
# (oneblock, variance): component_significance, variable_stability,
# score_stability, subspace_stability.
# variable_significance is intentionally excluded (blocked in v1).
#
# Returns a multifer_adapter object ready to pass to register_infer_adapter().
# Does NOT register or touch the registry.
make_stub_oneblock_adapter <- function(adapter_id = "stub_pca",
                                       validity_level = "conditional") {
  caps <- capability_matrix(
    list(
      geometry = "oneblock",
      relation = "variance",
      targets  = c("component_significance",
                   "variable_stability",
                   "score_stability",
                   "subspace_stability")
    )
  )

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = "0.0.1",
    shape_kinds     = "oneblock",
    capabilities    = caps,
    # hooks required by the four v1 targets
    roots          = stub_hook_roots,
    scores         = stub_hook_scores,
    loadings       = stub_hook_loadings,
    truncate       = stub_hook_truncate,
    residualize    = stub_hook_residualize,
    refit          = stub_hook_refit,
    null_action    = stub_hook_null_action,
    component_stat = stub_hook_component_stat,
    variable_stat  = stub_hook_variable_stat,
    score_stat     = stub_hook_score_stat,
    align          = stub_hook_align,
    validity_level = validity_level,
    declared_assumptions = c("gaussian_noise", "correct_rank")
  )
}

# ---------------------------------------------------------------------------
# Factory: cross-block stub adapter
# ---------------------------------------------------------------------------

# Builds a stub adapter for cross-block inference.
# `relations` controls which relations are declared in the capability matrix.
# When length(relations) == 1, a recipe compiled with no explicit relation
# argument will resolve unambiguously (no strict-dispatch error).
# When length(relations) >= 2 (e.g. c("covariance", "correlation")), a recipe
# compiled with no relation argument under strict = TRUE will error -- this is
# the ambiguity case that issue 9qu.9 exercises.
#
# Returns a multifer_adapter object ready to pass to register_infer_adapter().
# Does NOT register or touch the registry.
make_stub_cross_adapter <- function(adapter_id = "stub_cross",
                                    relations = c("covariance", "correlation"),
                                    validity_level = "conditional") {
  # Build one capability-matrix entry per relation, each declaring all four
  # v1 targets.
  entries <- lapply(relations, function(rel) {
    list(
      geometry = "cross",
      relation = rel,
      targets  = c("component_significance",
                   "variable_stability",
                   "score_stability",
                   "subspace_stability")
    )
  })

  caps <- do.call(capability_matrix, entries)

  infer_adapter(
    adapter_id      = adapter_id,
    adapter_version = "0.0.1",
    shape_kinds     = "cross",
    capabilities    = caps,
    roots          = stub_hook_roots,
    scores         = stub_hook_scores,
    loadings       = stub_hook_loadings,
    truncate       = stub_hook_truncate,
    residualize    = stub_hook_residualize,
    refit          = stub_hook_refit,
    null_action    = stub_hook_null_action,
    component_stat = stub_hook_component_stat,
    variable_stat  = stub_hook_variable_stat,
    score_stat     = stub_hook_score_stat,
    align          = stub_hook_align,
    validity_level = validity_level,
    declared_assumptions = c("gaussian_noise", "correct_rank")
  )
}

# ---------------------------------------------------------------------------
# Factory: canned infer_result for cross-block contract-proof tests
# ---------------------------------------------------------------------------

# Builds a fully valid infer_result for a cross-block case.
# variable_stability and score_stability rows span BOTH domains "X" and "Y",
# exercising the multi-domain path in the schema.
#
# unit_ids: character vector of unit ids, each treated as "component" type
# adapter_id: character scalar used in the provenance block
make_stub_cross_infer_result <- function(unit_ids = c("u1", "u2"),
                                         adapter_id = "stub_cross") {
  n_units <- length(unit_ids)

  # units table: all components, identifiable, selected
  units <- infer_units(
    unit_id      = unit_ids,
    unit_type    = rep("component", n_units),
    members      = lapply(seq_len(n_units), function(i) as.integer(i)),
    identifiable = rep(TRUE, n_units),
    selected     = rep(TRUE, n_units)
  )

  # component_tests: one row per unit, canned statistic and p-value
  component_tests <- infer_component_tests(
    unit_id        = unit_ids,
    statistic      = seq(5.0, by = -1.0, length.out = n_units),
    p_value        = rep(0.01, n_units),
    mc_uncertainty = rep(0.002, n_units),
    stopped_at     = rep(500L, n_units),
    null_label     = rep("sign_flip", n_units),
    validity_level = rep("conditional", n_units)
  )

  # variable_stability: 2 variables per unit in BOTH domains "X" and "Y"
  n_vars   <- 2L
  domains  <- c("X", "Y")
  vs_unit_ids <- rep(rep(unit_ids, each = n_vars), times = length(domains))
  vs_vars     <- rep(rep(paste0("v", seq_len(n_vars)), times = n_units),
                     times = length(domains))
  vs_domains  <- rep(domains, each = n_units * n_vars)
  n_vs        <- length(vs_unit_ids)
  variable_stability <- infer_variable_stability(
    unit_id            = vs_unit_ids,
    domain             = vs_domains,
    variable           = vs_vars,
    estimate           = rep(c(0.7, 0.4), times = n_units * length(domains)),
    stability          = rep(c(0.8, 0.6), times = n_units * length(domains)),
    selection_freq     = rep(NA_real_, n_vs),
    weight_sensitivity = rep(NA_real_, n_vs)
  )

  # score_stability: 2 observations per unit in BOTH domains "X" and "Y"
  n_obs       <- 2L
  ss_unit_ids <- rep(rep(unit_ids, each = n_obs), times = length(domains))
  ss_domains  <- rep(domains, each = n_units * n_obs)
  n_ss        <- length(ss_unit_ids)
  score_stability <- infer_score_stability(
    unit_id     = ss_unit_ids,
    domain      = ss_domains,
    observation = rep(seq_len(n_obs), times = n_units * length(domains)),
    estimate    = rep(c(0.3, -0.1), times = n_units * length(domains)),
    lower       = rep(c(0.1, -0.3), times = n_units * length(domains)),
    upper       = rep(c(0.5,  0.1), times = n_units * length(domains)),
    leverage    = rep(NA_real_, n_ss)
  )

  # subspace_stability: one row per unit
  subspace_stability <- infer_subspace_stability(
    unit_id              = unit_ids,
    principal_angle_mean = rep(0.04, n_units),
    principal_angle_max  = rep(0.09, n_units),
    alignment_method     = rep("procrustes", n_units),
    stability_label      = rep("stable", n_units)
  )

  assumptions <- infer_assumptions(
    declared = c("gaussian_noise", "correct_rank"),
    checked  = list()
  )

  mc <- infer_mc(
    rng_seed          = 99L,
    rng_kind          = "Mersenne-Twister",
    stopping_boundary = "none",
    batch_schedule    = 500L,
    stop_iteration    = NA_integer_,
    total_draws_used  = 500L,
    exceedance_counts = stats::setNames(rep(1L, n_units), unit_ids)
  )

  cost <- infer_cost(
    full_data_ops    = 1L,
    core_updates     = 1L,
    mc_budget_spent  = stats::setNames(rep(500L, n_units), unit_ids),
    cache_hits       = 0.0,
    wall_time_phases = c(fit = 0.02, mc = 0.06)
  )

  provenance <- infer_provenance(
    adapter_id      = adapter_id,
    adapter_version = "0.0.1",
    capabilities    = c("component_significance",
                        "variable_stability",
                        "score_stability",
                        "subspace_stability"),
    call            = NULL
  )

  infer_result(
    units              = units,
    component_tests    = component_tests,
    variable_stability = variable_stability,
    score_stability    = score_stability,
    subspace_stability = subspace_stability,
    assumptions        = assumptions,
    mc                 = mc,
    cost               = cost,
    provenance         = provenance
  )
}

# ---------------------------------------------------------------------------
# Factory: canned infer_result for oneblock contract-proof tests
# ---------------------------------------------------------------------------

# Builds a fully valid infer_result with canned (non-real) values.
# All downstream tables key off unit_ids drawn from the `units` table.
# Suitable for testing schema shape, foreign-key validation, and print
# output without running any real inference engine.
#
# unit_ids: character vector of unit ids, each treated as "component" type
# adapter_id: character scalar used in the provenance block
make_stub_oneblock_infer_result <- function(unit_ids = c("u1", "u2"),
                                            adapter_id = "stub_pca") {
  n_units <- length(unit_ids)

  # units table: all components, identifiable, selected
  units <- infer_units(
    unit_id      = unit_ids,
    unit_type    = rep("component", n_units),
    members      = lapply(seq_len(n_units), function(i) as.integer(i)),
    identifiable = rep(TRUE, n_units),
    selected     = rep(TRUE, n_units)
  )

  # component_tests: one row per unit, canned statistic and p-value
  component_tests <- infer_component_tests(
    unit_id        = unit_ids,
    statistic      = seq(5.0, by = -1.0, length.out = n_units),
    p_value        = rep(0.01, n_units),
    mc_uncertainty = rep(0.002, n_units),
    stopped_at     = rep(500L, n_units),
    null_label     = rep("sign_flip", n_units),
    validity_level = rep("conditional", n_units)
  )

  # variable_stability: 3 variables per unit in domain "X"
  n_vars <- 3L
  vs_unit_ids <- rep(unit_ids, each = n_vars)
  vs_vars     <- rep(paste0("v", seq_len(n_vars)), times = n_units)
  variable_stability <- infer_variable_stability(
    unit_id            = vs_unit_ids,
    domain             = rep("X", length(vs_unit_ids)),
    variable           = vs_vars,
    estimate           = rep(c(0.8, 0.5, 0.3), times = n_units),
    stability          = rep(c(0.9, 0.7, 0.6), times = n_units),
    selection_freq     = rep(NA_real_, length(vs_unit_ids)),
    weight_sensitivity = rep(NA_real_, length(vs_unit_ids))
  )

  # score_stability: 2 observations per unit in domain "X"
  n_obs <- 2L
  ss_unit_ids <- rep(unit_ids, each = n_obs)
  score_stability <- infer_score_stability(
    unit_id     = ss_unit_ids,
    domain      = rep("X", length(ss_unit_ids)),
    observation = rep(seq_len(n_obs), times = n_units),
    estimate    = rep(c(0.2, -0.1), times = n_units),
    lower       = rep(c(0.0, -0.3), times = n_units),
    upper       = rep(c(0.4,  0.1), times = n_units),
    leverage    = rep(NA_real_, length(ss_unit_ids))
  )

  # subspace_stability: one row per unit
  subspace_stability <- infer_subspace_stability(
    unit_id              = unit_ids,
    principal_angle_mean = rep(0.05, n_units),
    principal_angle_max  = rep(0.10, n_units),
    alignment_method     = rep("sign", n_units),
    stability_label      = rep("stable", n_units)
  )

  assumptions <- infer_assumptions(
    declared = c("gaussian_noise", "correct_rank"),
    checked  = list()
  )

  mc <- infer_mc(
    rng_seed          = 42L,
    rng_kind          = "Mersenne-Twister",
    stopping_boundary = "none",
    batch_schedule    = 500L,
    stop_iteration    = NA_integer_,
    total_draws_used  = 500L,
    exceedance_counts = stats::setNames(rep(1L, n_units), unit_ids)
  )

  cost <- infer_cost(
    full_data_ops    = 1L,
    core_updates     = 0L,
    mc_budget_spent  = stats::setNames(rep(500L, n_units), unit_ids),
    cache_hits       = 0.0,
    wall_time_phases = c(fit = 0.01, mc = 0.05)
  )

  provenance <- infer_provenance(
    adapter_id      = adapter_id,
    adapter_version = "0.0.1",
    capabilities    = c("component_significance",
                        "variable_stability",
                        "score_stability",
                        "subspace_stability"),
    call            = NULL
  )

  infer_result(
    units              = units,
    component_tests    = component_tests,
    variable_stability = variable_stability,
    score_stability    = score_stability,
    subspace_stability = subspace_stability,
    assumptions        = assumptions,
    mc                 = mc,
    cost               = cost,
    provenance         = provenance
  )
}
