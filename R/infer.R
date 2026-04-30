#' Top-level inference dispatcher
#'
#' The user-facing entry point for the shipped inference engine.
#' Given an adapter and data, runs the full pipeline:
#' \enumerate{
#'   \item resolve / build a \code{\link{typed_shape}};
#'   \item compile a recipe via \code{\link{infer_recipe}()} with
#'     strict dispatch by default;
#'   \item compile an executable plan and run its engine callback;
#'   \item fit (or use the supplied) original fit, then run
#'     \code{\link{bootstrap_fits}()} to produce a perturbation
#'     artifact;
#'   \item invoke \code{\link{variable_stability_from_bootstrap}()},
#'     \code{\link{score_stability_from_bootstrap}()}, and
#'     \code{\link{subspace_stability_from_bootstrap}()};
#'   \item assemble the frozen \code{\link{infer_result}()} schema.
#' }
#'
#' The default path runs the exact refit-first engine for component testing and
#' bootstrap summaries. When adapters expose \code{core()} and
#' \code{update_core()} hooks, \code{fast_path = "auto"} uses the Phase 1.5
#' core-update path for bootstrap refits; \code{fast_path = "off"} keeps the
#' refit path for agreement checks. Sequential Monte Carlo uses a global budget
#' allocator when \code{B_total} is supplied, and \code{parallel} can fan out
#' bootstrap replicates through the \code{mirai} backend.
#'
#' @param adapter Either a \code{multifer_adapter} object or a string
#'   id registered with \code{\link{register_infer_adapter}()}.
#' @param data For oneblock, a numeric matrix. For cross, a list with
#'   elements \code{X} and \code{Y}. For multiblock, a list of aligned
#'   numeric matrix blocks. For geneig, a data payload accepted by the
#'   adapter's operator constructor (currently the LDA wrapper uses a
#'   list with \code{X} and \code{y}). For adapter-owned geometry, any
#'   payload accepted by the adapter hooks.
#' @param ... Reserved for future use.
#' @param recipe Optional pre-compiled
#'   \code{multifer_infer_recipe}; if supplied, \code{geometry},
#'   \code{relation}, \code{design}, \code{targets}, and \code{strict}
#'   are ignored.
#' @param geometry Character scalar for the geometry kind. Passed to
#'   \code{\link{infer_recipe}()}.
#' @param relation Character scalar for the relation kind.
#' @param design A \code{multifer_design} object.
#' @param targets Character vector of requested targets, or
#'   \code{"default"} to expand to all supported targets except
#'   \code{variable_significance}. Geneig wrappers currently default to
#'   \code{"component_significance"} until the public stability path is
#'   added.
#' @param strict Logical, passed to \code{\link{infer_recipe}()}.
#'   Default \code{TRUE}.
#' @param B Integer, per-rung cap on Monte Carlo draws.
#' @param B_total Optional integer, global Monte Carlo budget shared
#'   across ladder rungs. Defaults to `B * max_steps`.
#' @param mc_batch_size Positive integer, Besag-Clifford batch size
#'   within each rung. Default `32L`.
#' @param R Integer, number of bootstrap replicates.
#' @param alpha Numeric, significance threshold for the test ladder.
#' @param model Optional pre-fit object. When NULL, \code{adapter$refit(NULL, data)}
#'   is called to produce the original fit.
#' @param seed Optional integer seed.
#' @param parallel One of `"sequential"`, `"mirai"`, `"auto"`. Controls
#'   the bootstrap fan-out backend. `"sequential"` keeps the single-stream
#'   Phase 1 loop. `"mirai"` / `"auto"` fan out per-replicate tasks to a
#'   mirai daemon pool with deterministic per-task seeds.
#' @param fast_path One of `"auto"`, `"off"`. `"auto"` uses the Phase 1.5
#'   core()/update_core() fast path when the adapter exposes it.
#'   `"off"` forces the refit path, useful for correctness validation.
#' @param auto_subspace Logical. When `TRUE` (default), near-tied
#'   roots are automatically bundled into a subspace unit via
#'   [form_units()] `group_near_ties = TRUE`. Subspace units cause
#'   `print.infer_result()` to lead with principal-angle stability.
#' @param tie_threshold Positive numeric. Relative-gap threshold used
#'   when `auto_subspace = TRUE`. Default `0.01`.
#' @param return_bundle Logical. When `TRUE`, return an `infer_bundle`
#'   containing the frozen `infer_result` plus reusable artifacts such as the
#'   bootstrap artifact. Default `FALSE`.
#'
#' @return An \code{\link{infer_result}}.
#' @export
infer <- function(adapter,
                  data,
                  ...,
                  recipe   = NULL,
                  geometry = NULL,
                  relation = NULL,
                  design   = NULL,
                  targets  = "default",
                  strict   = TRUE,
                  B             = 1000L,
                  B_total       = NULL,
                  mc_batch_size = 32L,
                  R             = 500L,
                  alpha         = 0.05,
                  model         = NULL,
                  seed          = NULL,
                  parallel      = c("sequential", "mirai", "auto"),
                  fast_path     = c("auto", "off"),
                  auto_subspace = TRUE,
                  tie_threshold = 0.01,
                  return_bundle = FALSE) {

  call <- match.call()
  parallel  <- match.arg(parallel)
  fast_path <- match.arg(fast_path)
  t_start_total <- proc.time()[["elapsed"]]

  ## --- reset thin-SVD cache for this call ------------------------------------
  # Cache hit rate is reported in the $cost block; see Phase 1.5 §30 item 1.
  svd_cache_reset()

  ## --- resolve adapter --------------------------------------------------------

  if (is.character(adapter) && length(adapter) == 1L) {
    adapter_obj <- get_infer_adapter(adapter)
  } else if (inherits(adapter, "multifer_adapter")) {
    adapter_obj <- adapter
  } else {
    stop("`adapter` must be a multifer_adapter or a registered adapter id.",
         call. = FALSE)
  }

  ## --- compile recipe ---------------------------------------------------------

  if (is.null(recipe)) {
    recipe <- infer_recipe(
      geometry = geometry,
      relation = relation,
      design   = design,
      targets  = targets,
      adapter  = adapter_obj,
      strict   = strict
    )
  } else if (!inherits(recipe, "multifer_infer_recipe")) {
    stop("`recipe` must be a multifer_infer_recipe.", call. = FALSE)
  }

  execution_plan <- compile_infer_execution_plan(
    recipe = recipe,
    adapter = adapter_obj,
    data = data,
    model = model
  )
  problem <- execution_plan$problem
  resolved_targets <- execution_plan$targets

  ## --- executable validity checks --------------------------------------------
  # Every mature adapter ships concrete checked_assumptions; strict mode
  # (the default) errors on any violation before the engine starts.
  check_results <- run_adapter_checks(
    adapter_obj,
    data,
    recipe = recipe,
    strict = strict
  )

  ## --- engine: sequential deflation ladder ------------------------------------

  t_engine_start <- proc.time()[["elapsed"]]

  engine_out <- execution_plan$run_engine(
    B = B,
    B_total = B_total,
    batch_size = mc_batch_size,
    alpha = alpha,
    seed = seed,
    auto_subspace = auto_subspace,
    tie_threshold = tie_threshold
  )

  units <- engine_out$units

  ## --- component_tests for the frozen schema ----------------------------------

  step_results <- engine_out$ladder_result$step_results

  # Build per-unit component_tests rows, one per unit. For tested rungs,
  # use the ladder p-value; for untested rungs, p_value = NA_real_.
  rejected <- engine_out$ladder_result$rejected_through
  n_units  <- nrow(units)
  unit_ids <- units$unit_id
  steps_tested <- length(step_results)

  ct_unit_id  <- character(0)
  ct_stat     <- numeric(0)
  ct_p        <- numeric(0)
  ct_se       <- numeric(0)
  ct_stopped  <- integer(0)
  ct_null_lab <- character(0)
  ct_validity <- character(0)

  # The per-rung null label is the engine's own null label. Dispatch no
  # longer inspects geometry or design to guess what null the engine
  # ran -- that was the architectural smell multifer-5ow was filed
  # against. Engines populate labels$null via infer_result.R's
  # engine-provenance contract; here we just read it.
  null_label <- engine_out$labels$null %||% NA_character_
  validity   <- adapter_obj$validity_level
  if (!is.null(recipe$downgrade_reason) && !is.na(recipe$downgrade_reason)) {
    validity <- recipe$validity_level
  }

  for (i in seq_len(n_units)) {
    if (i <= steps_tested) {
      sr <- step_results[[i]]
      ct_unit_id  <- c(ct_unit_id, unit_ids[i])
      ct_stat     <- c(ct_stat,    as.numeric(sr$observed_stat))
      ct_p        <- c(ct_p,       as.numeric(sr$p_value))
      ct_se       <- c(ct_se,      as.numeric(sr$mc_se))
      ct_stopped  <- c(ct_stopped, as.integer(sr$B))
      ct_null_lab <- c(ct_null_lab, null_label)
      ct_validity <- c(ct_validity, validity)
    }
  }

  component_tests <- infer_component_tests(
    unit_id        = ct_unit_id,
    statistic      = ct_stat,
    p_value        = ct_p,
    mc_uncertainty = ct_se,
    stopped_at     = ct_stopped,
    null_label     = ct_null_lab,
    validity_level = ct_validity
  )

  t_engine_end <- proc.time()[["elapsed"]]

  ## --- bootstrap stream + stability consumers ---------------------------------

  needs_bootstrap <- any(c("variable_stability",
                           "score_stability",
                           "subspace_stability") %in% resolved_targets)

  variable_stab <- infer_variable_stability()
  score_stab    <- infer_score_stability()
  subspace_stab <- infer_subspace_stability()
  artifact <- NULL
  original_fit_for_bundle <- NULL

  t_boot_start  <- proc.time()[["elapsed"]]
  t_boot_end    <- t_boot_start

  if (needs_bootstrap && R > 0L) {
    if (!isTRUE(execution_plan$bootstrap_supported)) {
      stop(execution_plan$bootstrap_error, call. = FALSE)
    }

    original_fit <- execution_plan$original_fit(engine_out)
    original_fit_for_bundle <- original_fit

    artifact <- bootstrap_fits(
      recipe       = recipe,
      adapter      = adapter_obj,
      data         = data,
      original_fit = original_fit,
      units        = units,
      R            = R,
      method_align = "sign",
      seed         = seed,
      parallel     = parallel,
      fast_path    = fast_path,
      store_aligned_scores = "score_stability" %in% resolved_targets &&
        !is.null(adapter_obj$project_scores)
    )

    if ("variable_stability" %in% resolved_targets) {
      variable_stab <- variable_stability_from_bootstrap(artifact, units)
    }
    if ("score_stability" %in% resolved_targets) {
      score_stab <- score_stability_from_bootstrap(artifact, data, units)
    }
    if ("subspace_stability" %in% resolved_targets) {
      subspace_stab <- subspace_stability_from_bootstrap(
        artifact, original_fit, adapter_obj, units
      )
    }

    t_boot_end <- proc.time()[["elapsed"]]
  }

  ## --- assumptions + mc + cost + provenance -----------------------------------

  assumptions <- infer_assumptions(
    declared = adapter_obj$declared_assumptions,
    checked  = check_results
  )

  ladder_res_meta    <- engine_out$ladder_result
  total_draws_used   <- as.integer(ladder_res_meta$total_draws_used %||%
                                   sum(vapply(step_results,
                                              function(sr) as.integer(sr$drawn %||% sr$B %||% 0L),
                                              integer(1L))))
  batch_schedule_agg <- as.integer(ladder_res_meta$batch_schedule %||% integer(0L))
  stop_boundary_lab  <- ladder_res_meta$stopping_boundary %||% "fixed_B"

  engine_labels <- engine_out$labels %||% list()
  label_or_na <- function(key) {
    val <- engine_labels[[key]]
    if (is.null(val) || length(val) == 0L) {
      return(NA_character_)
    }
    as.character(val[[1L]])
  }

  mc_block <- infer_mc(
    rng_seed          = if (is.null(seed)) NA_integer_ else as.integer(seed),
    rng_kind          = RNGkind()[1L],
    stopping_boundary = stop_boundary_lab,
    batch_schedule    = batch_schedule_agg,
    stop_iteration    = total_draws_used,
    total_draws_used  = total_draws_used,
    exceedance_counts = vapply(
      step_results,
      function(sr) as.integer(sr$r),
      integer(1L)
    ),
    statistic_label = label_or_na("statistic"),
    null_label      = label_or_na("null"),
    estimand_label  = label_or_na("estimand")
  )

  t_total_end <- proc.time()[["elapsed"]]
  cache_rate  <- svd_cache_rate()
  if (is.na(cache_rate)) cache_rate <- 0
  core_updates_used <- if (exists("artifact", inherits = FALSE) &&
                           isTRUE(artifact$used_fast_path)) {
    as.integer(R)
  } else {
    0L
  }
  cost_block <- infer_cost(
    full_data_ops    = 1L,
    core_updates     = core_updates_used,
    mc_budget_spent  = total_draws_used,
    cache_hits       = cache_rate,
    wall_time_phases = c(
      engine = t_engine_end - t_engine_start,
      mc     = t_engine_end - t_engine_start,
      boot   = t_boot_end - t_boot_start,
      total  = t_total_end - t_start_total
    )
  )

  provenance <- infer_provenance(
    adapter_id      = adapter_obj$adapter_id,
    adapter_version = adapter_obj$adapter_version,
    capabilities    = paste0(
      problem$shape$geometry$kind, "/",
      problem$shape$relation$kind, ":",
      paste(resolved_targets, collapse = "+")
    ),
    call = call
  )

  result <- infer_result(
    units              = units,
    component_tests    = component_tests,
    variable_stability = variable_stab,
    score_stability    = score_stab,
    subspace_stability = subspace_stab,
    assumptions        = assumptions,
    mc                 = mc_block,
    cost               = cost_block,
    provenance         = provenance
  )

  if (isTRUE(return_bundle)) {
    return(infer_bundle_result(
      result = result,
      artifacts = list(
        bootstrap = artifact,
        original_fit = original_fit_for_bundle,
        engine = engine_out,
        recipe = recipe,
        adapter = adapter_obj
      )
    ))
  }

  result
}

#' Return inference results with reusable artifacts
#'
#' `infer_bundle()` is a convenience wrapper around [infer()] that keeps the
#' frozen `infer_result` unchanged while exposing artifacts useful to sidecar
#' consumers such as [feature_evidence_from_bootstrap()].
#'
#' @param ... Passed to [infer()].
#'
#' @return A list with `result` and `artifacts`, class `multifer_infer_bundle`.
#' @export
infer_bundle <- function(...) {
  infer(..., return_bundle = TRUE)
}

infer_bundle_result <- function(result, artifacts) {
  if (!inherits(result, "infer_result")) {
    stop("`result` must be an infer_result.", call. = FALSE)
  }
  structure(
    list(result = result, artifacts = artifacts),
    class = c("multifer_infer_bundle", "list")
  )
}

.infer_original_fit <- function(adapter_obj, data, geom_kind, rel_kind, model = NULL) {
  if (!is.null(model)) {
    return(model)
  }
  if (identical(geom_kind, "cross")) {
    adapter_obj$refit(NULL, list(X = data$X, Y = data$Y, relation = rel_kind))
  } else {
    adapter_obj$refit(NULL, data)
  }
}
