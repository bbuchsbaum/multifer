#' The frozen `infer_result` schema
#'
#' `infer_result()` is the output contract of every `infer()` call.
#' It is deliberately **unit-centered**: every downstream table keys off
#' a `unit_id` drawn from the `units` table, not a raw component index.
#' A unit is either a singleton component (identifiable) or a subspace
#' (grouped tied / near-tied roots). This is what makes tied-root
#' handling honest by default (Part 5 section 37, section 44).
#'
#' This schema is **frozen** at the end of Phase 0. Adding or removing
#' fields requires an explicit design review.
#'
#' @section Fields:
#'
#' - `units` -- data.frame with columns `unit_id`, `unit_type`
#'   (`"component"` or `"subspace"`), `members` (list-column of integer
#'   vectors), `identifiable` (logical), `selected` (logical).
#' - `component_tests` -- per-unit significance row: statistic, p-value,
#'   Monte Carlo uncertainty, null label, validity level.
#' - `variable_stability` -- keyed by `unit_id` and `domain`.
#' - `score_stability` -- keyed by `unit_id` and `domain`.
#' - `subspace_stability` -- keyed by `unit_id`.
#' - `assumptions` -- list with `declared` and `checked`.
#' - `mc` -- reproducibility log for the Monte Carlo layer (Part 5 section 41).
#' - `cost` -- cost accounting block (Part 4 section 32).
#' - `provenance` -- adapter id, version, capability flags, call.
#'
#' @name infer_result
NULL

# ---- units ------------------------------------------------------------------

#' Build the `units` table
#'
#' @param unit_id Character vector of unique unit identifiers.
#' @param unit_type Character vector, each element in
#'   `c("component", "subspace")`.
#' @param members A list of integer vectors; element `i` gives the raw
#'   component indices that belong to unit `i`. For `"component"` units
#'   this is a length-1 integer vector; for `"subspace"` units it is
#'   length >= 2.
#' @param identifiable Logical vector. `TRUE` means sign or legacy
#'   Procrustes alignment resolves the unit. Subspace units are typically
#'   `FALSE`.
#' @param selected Logical vector. `TRUE` means the sequential test
#'   ladder accepted this unit as significant.
#'
#' @return A data.frame with an `multifer_units` S3 class.
#' @export
infer_units <- function(unit_id = character(0),
                        unit_type = character(0),
                        members = list(),
                        identifiable = logical(0),
                        selected = logical(0)) {
  n <- length(unit_id)
  if (length(unit_type) != n || length(members) != n ||
      length(identifiable) != n || length(selected) != n) {
    stop("all `infer_units()` arguments must have the same length.",
         call. = FALSE)
  }
  if (n > 0L) {
    if (anyDuplicated(unit_id)) {
      stop("`unit_id` must be unique.", call. = FALSE)
    }
    bad_type <- setdiff(unit_type, c("component", "subspace"))
    if (length(bad_type) > 0L) {
      stop(sprintf("unknown unit_type(s): %s",
                   paste(bad_type, collapse = ", ")), call. = FALSE)
    }
    if (!all(vapply(members, function(m) is.integer(m) || is.numeric(m),
                    logical(1L)))) {
      stop("`members` must be a list of integer vectors.", call. = FALSE)
    }
  }
  structure(
    data.frame(
      unit_id = as.character(unit_id),
      unit_type = as.character(unit_type),
      identifiable = as.logical(identifiable),
      selected = as.logical(selected),
      stringsAsFactors = FALSE
    ),
    members = members,
    class = c("multifer_units", "data.frame")
  )
}

#' @export
`$.multifer_units` <- function(x, name) {
  if (identical(name, "members")) attr(x, "members") else NextMethod()
}

# ---- component_tests --------------------------------------------------------

#' Build the `component_tests` table
#'
#' One row per unit in `units`, carrying the sequential test result.
#'
#' @param unit_id Character vector referencing `units$unit_id`.
#' @param statistic Observed test statistic (numeric).
#' @param p_value Monte Carlo p-value `(1 + r) / (B + 1)` (numeric).
#' @param mc_uncertainty Monte Carlo standard error of `p_value`.
#' @param stopped_at Integer, the MC iteration at which sequential
#'   stopping terminated (or total draws if run to completion).
#' @param null_label Character, name of the null action used.
#' @param validity_level Character, one of
#'   `c("exact", "conditional", "asymptotic", "heuristic")`.
#'
#' @return A data.frame with an `multifer_component_tests` S3 class.
#' @export
infer_component_tests <- function(unit_id = character(0),
                                  statistic = numeric(0),
                                  p_value = numeric(0),
                                  mc_uncertainty = numeric(0),
                                  stopped_at = integer(0),
                                  null_label = character(0),
                                  validity_level = character(0)) {
  n <- length(unit_id)
  lens <- c(length(statistic), length(p_value), length(mc_uncertainty),
            length(stopped_at), length(null_label), length(validity_level))
  if (!all(lens == n)) {
    stop("all `infer_component_tests()` arguments must have the same length.",
         call. = FALSE)
  }
  if (n > 0L) {
    bad <- setdiff(validity_level,
                   c("exact", "conditional", "asymptotic", "heuristic"))
    if (length(bad) > 0L) {
      stop(sprintf("unknown validity_level(s): %s",
                   paste(bad, collapse = ", ")), call. = FALSE)
    }
  }
  structure(
    data.frame(
      unit_id = as.character(unit_id),
      statistic = as.numeric(statistic),
      p_value = as.numeric(p_value),
      mc_uncertainty = as.numeric(mc_uncertainty),
      stopped_at = as.integer(stopped_at),
      null_label = as.character(null_label),
      validity_level = as.character(validity_level),
      stringsAsFactors = FALSE
    ),
    class = c("multifer_component_tests", "data.frame")
  )
}

# ---- variable_stability -----------------------------------------------------

#' Build the `variable_stability` table
#'
#' Keyed by `(unit_id, domain, variable)`. Carries stability measures,
#' not p-values. Variable significance is **not** part of v1
#' (Part 5 section 38).
#'
#' @param unit_id Character vector referencing `units$unit_id`.
#' @param domain Character vector, the block label (e.g. `"X"`, `"Y"`,
#'   `"block_1"`).
#' @param variable Character vector, variable name or id within the
#'   domain.
#' @param estimate Point estimate of the loading (numeric).
#' @param stability Stability measure (numeric, typically in `[0, 1]`).
#' @param selection_freq Optional numeric; `NA_real_` if not used.
#' @param weight_sensitivity Optional numeric; `NA_real_` if not used.
#'
#' @return A data.frame with an `multifer_variable_stability` S3 class.
#' @export
infer_variable_stability <- function(unit_id = character(0),
                                     domain = character(0),
                                     variable = character(0),
                                     estimate = numeric(0),
                                     stability = numeric(0),
                                     selection_freq = numeric(0),
                                     weight_sensitivity = numeric(0)) {
  n <- length(unit_id)
  lens <- c(length(domain), length(variable), length(estimate),
            length(stability), length(selection_freq),
            length(weight_sensitivity))
  if (!all(lens == n)) {
    stop("all `infer_variable_stability()` arguments must have the same length.",
         call. = FALSE)
  }
  structure(
    data.frame(
      unit_id = as.character(unit_id),
      domain = as.character(domain),
      variable = as.character(variable),
      estimate = as.numeric(estimate),
      stability = as.numeric(stability),
      selection_freq = as.numeric(selection_freq),
      weight_sensitivity = as.numeric(weight_sensitivity),
      stringsAsFactors = FALSE
    ),
    class = c("multifer_variable_stability", "data.frame")
  )
}

# ---- score_stability --------------------------------------------------------

#' Build the `score_stability` table
#'
#' Keyed by `(unit_id, domain, observation)`.
#'
#' @param unit_id Character vector referencing `units$unit_id`.
#' @param domain Character vector of block labels.
#' @param observation Integer vector, row index within the domain.
#' @param estimate Numeric score point estimate.
#' @param lower Numeric lower bound of the score interval.
#' @param upper Numeric upper bound of the score interval.
#' @param leverage Optional numeric leverage / influence measure.
#'
#' @return A data.frame with an `multifer_score_stability` S3 class.
#' @export
infer_score_stability <- function(unit_id = character(0),
                                  domain = character(0),
                                  observation = integer(0),
                                  estimate = numeric(0),
                                  lower = numeric(0),
                                  upper = numeric(0),
                                  leverage = numeric(0)) {
  n <- length(unit_id)
  lens <- c(length(domain), length(observation), length(estimate),
            length(lower), length(upper), length(leverage))
  if (!all(lens == n)) {
    stop("all `infer_score_stability()` arguments must have the same length.",
         call. = FALSE)
  }
  structure(
    data.frame(
      unit_id = as.character(unit_id),
      domain = as.character(domain),
      observation = as.integer(observation),
      estimate = as.numeric(estimate),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      leverage = as.numeric(leverage),
      stringsAsFactors = FALSE
    ),
    class = c("multifer_score_stability", "data.frame")
  )
}

# ---- subspace_stability -----------------------------------------------------

#' Build the `subspace_stability` table
#'
#' One row per unit, summarizing principal-angle stability.
#'
#' @param unit_id Character vector referencing `units$unit_id`.
#' @param principal_angle_mean Mean principal angle across bootstrap /
#'   perturbation replicates (numeric, radians).
#' @param principal_angle_max Max principal angle (numeric, radians).
#' @param alignment_method Character, one of
#'   `c("sign", "procrustes", "subspace")`. `procrustes` is retained as
#'   a legacy label for backward compatibility.
#' @param stability_label Character free-text label (`"stable"`,
#'   `"tied"`, `"unstable"`, etc.).
#'
#' @return A data.frame with an `multifer_subspace_stability` S3 class.
#' @export
infer_subspace_stability <- function(unit_id = character(0),
                                     principal_angle_mean = numeric(0),
                                     principal_angle_max = numeric(0),
                                     alignment_method = character(0),
                                     stability_label = character(0)) {
  n <- length(unit_id)
  lens <- c(length(principal_angle_mean), length(principal_angle_max),
            length(alignment_method), length(stability_label))
  if (!all(lens == n)) {
    stop("all `infer_subspace_stability()` arguments must have the same length.",
         call. = FALSE)
  }
  if (n > 0L) {
    bad <- setdiff(alignment_method, c("sign", "procrustes", "subspace"))
    if (length(bad) > 0L) {
      stop(sprintf("unknown alignment_method(s): %s",
                   paste(bad, collapse = ", ")), call. = FALSE)
    }
  }
  structure(
    data.frame(
      unit_id = as.character(unit_id),
      principal_angle_mean = as.numeric(principal_angle_mean),
      principal_angle_max = as.numeric(principal_angle_max),
      alignment_method = as.character(alignment_method),
      stability_label = as.character(stability_label),
      stringsAsFactors = FALSE
    ),
    class = c("multifer_subspace_stability", "data.frame")
  )
}

# ---- sidecar blocks: assumptions / mc / cost / provenance -------------------

#' Build an `assumptions` block
#'
#' @param declared Character vector of assumptions the adapter declares
#'   (trust-the-user: noise model, design correctness).
#' @param checked Named list of runtime checks; each element should be a
#'   list with `passed` (logical) and `detail` (character).
#'
#' @return A list with class `multifer_assumptions`.
#' @export
infer_assumptions <- function(declared = character(0),
                              checked = list()) {
  if (!is.character(declared)) {
    stop("`declared` must be a character vector.", call. = FALSE)
  }
  if (!is.list(checked)) {
    stop("`checked` must be a list.", call. = FALSE)
  }
  structure(
    list(declared = declared, checked = checked),
    class = "multifer_assumptions"
  )
}

#' Build the Monte Carlo reproducibility log
#'
#' Implements Part 5 section 41.
#'
#' The three `*_label` slots are **engine-provenance labels**. Each
#' engine is responsible for populating them so that `infer_result`
#' never needs top-level dispatch logic to describe what statistic,
#' null action, or estimand a given result reflects. The contract is:
#'
#' - `statistic_label` names the observed rung statistic (for example
#'   `"Vitale P3 tail-ratio on latent variance roots"` for the oneblock
#'   engine, or `"cross-fit incremental R^2"` for the predictive
#'   engine).
#' - `null_label` names the null action the engine applied to form
#'   null draws. It is a single engine-level label and is distinct
#'   from `component_tests$null_label`, which may vary per rung.
#' - `estimand_label` names the ordered object the ladder tests (for
#'   example `"latent variance roots"`, `"canonical correlations"`,
#'   `"held-out predictive gain"`). `print.infer_result` shows this
#'   label in its header so predictive and geneig results never read
#'   like a generic latent-root printout.
#'
#' @param rng_seed Integer or character.
#' @param rng_kind Character (e.g. `"Mersenne-Twister"`).
#' @param stopping_boundary Numeric or character describing the
#'   anytime-valid stopping boundary.
#' @param batch_schedule Integer vector, draws per batch.
#' @param stop_iteration Integer.
#' @param total_draws_used Integer.
#' @param exceedance_counts Named integer vector, one per unit.
#' @param statistic_label Character, the engine's label for the
#'   observed rung statistic. Defaults to `NA_character_`.
#' @param null_label Character, the engine's label for the null
#'   action. Defaults to `NA_character_`.
#' @param estimand_label Character, the engine's label for the
#'   ordered estimand the ladder tests. Defaults to `NA_character_`.
#'
#' @return A list with class `multifer_mc`.
#' @export
infer_mc <- function(rng_seed = NA_integer_,
                     rng_kind = NA_character_,
                     stopping_boundary = NA_character_,
                     batch_schedule = integer(0),
                     stop_iteration = NA_integer_,
                     total_draws_used = NA_integer_,
                     exceedance_counts = integer(0),
                     statistic_label = NA_character_,
                     null_label = NA_character_,
                     estimand_label = NA_character_) {
  structure(
    list(
      rng_seed = rng_seed,
      rng_kind = rng_kind,
      stopping_boundary = stopping_boundary,
      batch_schedule = as.integer(batch_schedule),
      stop_iteration = as.integer(stop_iteration),
      total_draws_used = as.integer(total_draws_used),
      exceedance_counts = exceedance_counts,
      statistic_label = as.character(statistic_label),
      null_label = as.character(null_label),
      estimand_label = as.character(estimand_label)
    ),
    class = "multifer_mc"
  )
}

#' Build the cost-accounting block
#'
#' Implements Part 4 section 32.
#'
#' @param full_data_ops Integer count of `O(n*p)` operations performed.
#' @param core_updates Integer count of `O(r_x * r_y)` core-update
#'   operations.
#' @param mc_budget_spent Named integer vector, draws per unit.
#' @param cache_hits Numeric in `[0, 1]`, cache hit rate.
#' @param wall_time_phases Named numeric vector of seconds, broken
#'   down by phase (`fit`, `core`, `mc`, `lift`, ...).
#'
#' @return A list with class `multifer_cost`.
#' @export
infer_cost <- function(full_data_ops = NA_integer_,
                       core_updates = NA_integer_,
                       mc_budget_spent = integer(0),
                       cache_hits = NA_real_,
                       wall_time_phases = numeric(0)) {
  structure(
    list(
      full_data_ops = as.integer(full_data_ops),
      core_updates = as.integer(core_updates),
      mc_budget_spent = mc_budget_spent,
      cache_hits = as.numeric(cache_hits),
      wall_time_phases = wall_time_phases
    ),
    class = "multifer_cost"
  )
}

#' Build a `provenance` block
#'
#' @param adapter_id Character, the registered adapter id.
#' @param adapter_version Character.
#' @param capabilities Character vector of capability flags the adapter
#'   declared.
#' @param call The captured call object (as returned by `match.call()`).
#'
#' @return A list with class `multifer_provenance`.
#' @export
infer_provenance <- function(adapter_id = NA_character_,
                             adapter_version = NA_character_,
                             capabilities = character(0),
                             call = NULL) {
  structure(
    list(
      adapter_id = as.character(adapter_id),
      adapter_version = as.character(adapter_version),
      capabilities = as.character(capabilities),
      call = call
    ),
    class = "multifer_provenance"
  )
}

# ---- top-level constructor + validator --------------------------------------

#' Build an `infer_result`
#'
#' Composes the nine frozen sub-blocks and enforces foreign-key
#' consistency: every `unit_id` appearing in a downstream table must
#' exist in `units$unit_id`.
#'
#' @param units An `infer_units()` object.
#' @param component_tests An `infer_component_tests()` object.
#' @param variable_stability An `infer_variable_stability()` object.
#' @param score_stability An `infer_score_stability()` object.
#' @param subspace_stability An `infer_subspace_stability()` object.
#' @param assumptions An `infer_assumptions()` list.
#' @param mc An `infer_mc()` list.
#' @param cost An `infer_cost()` list.
#' @param provenance An `infer_provenance()` list.
#'
#' @return An object with class `infer_result`.
#' @export
infer_result <- function(units,
                         component_tests = infer_component_tests(),
                         variable_stability = infer_variable_stability(),
                         score_stability = infer_score_stability(),
                         subspace_stability = infer_subspace_stability(),
                         assumptions = infer_assumptions(),
                         mc = infer_mc(),
                         cost = infer_cost(),
                         provenance = infer_provenance()) {

  if (!inherits(units, "multifer_units")) {
    stop("`units` must be built with infer_units().", call. = FALSE)
  }
  if (!inherits(component_tests, "multifer_component_tests")) {
    stop("`component_tests` must be built with infer_component_tests().",
         call. = FALSE)
  }
  if (!inherits(variable_stability, "multifer_variable_stability")) {
    stop("`variable_stability` must be built with infer_variable_stability().",
         call. = FALSE)
  }
  if (!inherits(score_stability, "multifer_score_stability")) {
    stop("`score_stability` must be built with infer_score_stability().",
         call. = FALSE)
  }
  if (!inherits(subspace_stability, "multifer_subspace_stability")) {
    stop("`subspace_stability` must be built with infer_subspace_stability().",
         call. = FALSE)
  }
  if (!inherits(assumptions, "multifer_assumptions")) {
    stop("`assumptions` must be built with infer_assumptions().",
         call. = FALSE)
  }
  if (!inherits(mc, "multifer_mc")) {
    stop("`mc` must be built with infer_mc().", call. = FALSE)
  }
  if (!inherits(cost, "multifer_cost")) {
    stop("`cost` must be built with infer_cost().", call. = FALSE)
  }
  if (!inherits(provenance, "multifer_provenance")) {
    stop("`provenance` must be built with infer_provenance().",
         call. = FALSE)
  }

  valid_ids <- units$unit_id
  check_fk <- function(tbl, name) {
    if (nrow(tbl) == 0L) return(invisible())
    bad <- setdiff(tbl$unit_id, valid_ids)
    if (length(bad) > 0L) {
      stop(sprintf(
        "`%s` references unknown unit_id(s): %s",
        name, paste(bad, collapse = ", ")
      ), call. = FALSE)
    }
  }
  check_fk(component_tests, "component_tests")
  check_fk(variable_stability, "variable_stability")
  check_fk(score_stability, "score_stability")
  check_fk(subspace_stability, "subspace_stability")

  structure(
    list(
      units = units,
      component_tests = component_tests,
      variable_stability = variable_stability,
      score_stability = score_stability,
      subspace_stability = subspace_stability,
      assumptions = assumptions,
      mc = mc,
      cost = cost,
      provenance = provenance
    ),
    class = "infer_result"
  )
}

#' An empty `infer_result`
#'
#' Convenience constructor for tests and stub adapters. Returns an
#' `infer_result` with zero units and empty tables, but with valid
#' (zero-length) shapes everywhere.
#'
#' @return An `infer_result` with zero units.
#' @export
empty_infer_result <- function() {
  infer_result(units = infer_units())
}

#' Test whether an object is an `infer_result`
#' @param x Object to test.
#' @export
is_infer_result <- function(x) inherits(x, "infer_result")

# ---- print / summary --------------------------------------------------------

#' @export
print.infer_result <- function(x, ...) {
  n_units <- nrow(x$units)
  n_sig <- if (n_units > 0L) sum(x$units$selected) else 0L
  n_subspace <- if (n_units > 0L) sum(x$units$unit_type == "subspace") else 0L

  estimand_lab <- x$mc$estimand_label
  statistic_lab <- x$mc$statistic_label
  null_lab <- x$mc$null_label

  cat("<infer_result>\n")
  if (length(estimand_lab) == 1L && !is.na(estimand_lab) && nzchar(estimand_lab)) {
    cat("  estimand:         ", estimand_lab, "\n", sep = "")
  }
  if (length(statistic_lab) == 1L && !is.na(statistic_lab) && nzchar(statistic_lab)) {
    cat("  statistic:        ", statistic_lab, "\n", sep = "")
  }
  if (length(null_lab) == 1L && !is.na(null_lab) && nzchar(null_lab)) {
    cat("  null:             ", null_lab, "\n", sep = "")
  }
  cat("  units:            ", n_units,
      " (", n_subspace, " subspace, ",
      max(n_units - n_subspace, 0L), " component)\n", sep = "")
  cat("  selected:         ", n_sig, "\n", sep = "")
  cat("  component_tests:  ", nrow(x$component_tests), " rows\n", sep = "")
  cat("  variable_stab:    ", nrow(x$variable_stability), " rows\n", sep = "")
  cat("  score_stab:       ", nrow(x$score_stability), " rows\n", sep = "")
  cat("  subspace_stab:    ", nrow(x$subspace_stability), " rows\n", sep = "")
  cat("  adapter:          ", x$provenance$adapter_id, "\n", sep = "")

  unit_section <- .format_infer_units_for_display(x)
  if (length(unit_section) > 0L) {
    cat("\n")
    cat(unit_section, sep = "\n")
    cat("\n")
  }

  cat("  NOTE: variable_stability is a STABILITY measure, not a p-value.\n")
  cat("        Variable significance is deferred to Phase 3 (Part 5 section 38).\n")
  invisible(x)
}

#' @export
summary.infer_result <- function(object, ...) {
  print(object, ...)
  invisible(object)
}

.format_infer_units_for_display <- function(x) {
  if (nrow(x$units) == 0L) {
    return(character(0))
  }

  tests_by_unit <- if (nrow(x$component_tests) > 0L) {
    split(x$component_tests, x$component_tests$unit_id)
  } else {
    list()
  }
  subspace_by_unit <- if (nrow(x$subspace_stability) > 0L) {
    split(x$subspace_stability, x$subspace_stability$unit_id)
  } else {
    list()
  }

  sub_idx <- which(x$units$unit_type == "subspace")
  comp_idx <- which(x$units$unit_type != "subspace")
  lines <- character(0)

  if (length(sub_idx) > 0L) {
    lines <- c(lines, "  Subspace inference (near-tied roots):")
    lines <- c(lines, vapply(sub_idx, function(i) {
      .format_infer_unit_line(x, i, tests_by_unit, subspace_by_unit)
    }, character(1L)))
  }

  if (length(comp_idx) > 0L) {
    if (length(lines) > 0L) {
      lines <- c(lines, "")
    }
    lines <- c(lines, "  Singleton units:")
    lines <- c(lines, vapply(comp_idx, function(i) {
      .format_infer_unit_line(x, i, tests_by_unit, subspace_by_unit)
    }, character(1L)))
  }

  lines
}

.format_infer_unit_line <- function(x, i, tests_by_unit, subspace_by_unit) {
  unit_id <- x$units$unit_id[i]
  members <- x$units$members[[i]]
  members_txt <- paste(members, collapse = ",")
  selected_txt <- if (isTRUE(x$units$selected[i])) "selected" else "not selected"

  test_row <- tests_by_unit[[unit_id]]
  p_txt <- "NA"
  validity_txt <- NULL
  if (!is.null(test_row) && nrow(test_row) > 0L) {
    p_value <- test_row$p_value[1L]
    if (is.finite(p_value)) {
      p_txt <- sprintf("%.3g", p_value)
    }
    validity_txt <- test_row$validity_level[1L]
  }

  if (identical(x$units$unit_type[i], "subspace")) {
    sub_row <- subspace_by_unit[[unit_id]]
    mean_angle <- NA_real_
    max_angle <- NA_real_
    stability_label <- NULL
    if (!is.null(sub_row) && nrow(sub_row) > 0L) {
      mean_angle <- sub_row$principal_angle_mean[1L]
      max_angle <- sub_row$principal_angle_max[1L]
      stability_label <- sub_row$stability_label[1L]
    }
    return(sprintf(
      "    %s [members %s; %s; p = %s%s]: mean angle = %s, max angle = %s%s",
      unit_id,
      members_txt,
      selected_txt,
      p_txt,
      if (!is.null(validity_txt) && nzchar(validity_txt)) paste0(", ", validity_txt) else "",
      if (is.finite(mean_angle)) sprintf("%.3f", mean_angle) else "NA",
      if (is.finite(max_angle)) sprintf("%.3f", max_angle) else "NA",
      if (!is.null(stability_label) && nzchar(stability_label)) paste0(" (", stability_label, ")") else ""
    ))
  }

  sprintf(
    "    %s [member %s; %s; p = %s%s]",
    unit_id,
    members_txt,
    selected_txt,
    p_txt,
    if (!is.null(validity_txt) && nzchar(validity_txt)) paste0(", ", validity_txt) else ""
  )
}
