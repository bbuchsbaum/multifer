#' Check an adapter against its declared inference contract
#'
#' `check_infer_adapter()` is a smoke-test harness for adapter authors and
#' for `multifer`'s own engine refactors. It walks the adapter's declared
#' `(geometry, relation, target)` capability rows, compiles each requested
#' target as a recipe, runs the adapter's executable validity checks, and
#' optionally runs a small `infer()` call for that target.
#'
#' The function is deliberately report-oriented: by default it returns every
#' row it checked, including failures, so a package author can see exactly
#' which declared capability is not executable. Set `fail_on_error = TRUE`
#' to turn any failed row into an error, which is useful in downstream tests.
#'
#' @param adapter A `multifer_adapter` object or registered adapter id.
#' @param fixture Data fixture accepted by the adapter. This may be a single
#'   data object used for every checked capability, or a named list containing
#'   entries named `"geometry/relation"`, `"geometry"`, or `"default"`.
#' @param targets Character vector of targets to check, or `"default"` to
#'   check every target declared for each capability pair except
#'   `"variable_significance"`.
#' @param geometry Optional character vector limiting checked geometries.
#' @param relation Optional character vector limiting checked relations.
#' @param run One of `"smoke"` or `"compile"`. `"compile"` only compiles
#'   recipes and runs validity checks. `"smoke"` also runs `infer()` with
#'   small Monte Carlo and bootstrap budgets.
#' @param B Integer Monte Carlo budget passed to `infer()` in smoke mode.
#' @param R Integer bootstrap replicate count passed to `infer()` in smoke
#'   mode.
#' @param seed Optional integer seed passed to `infer()` in smoke mode.
#' @param strict Logical passed to `infer_recipe()`, `run_adapter_checks()`,
#'   and `infer()`.
#' @param fail_on_error Logical. If `TRUE`, stop when any checked row fails.
#' @param ... Additional arguments forwarded to `infer()` in smoke mode.
#'
#' @return A data frame of class `multifer_adapter_check` with one row per
#'   checked capability target.
#' @export
check_infer_adapter <- function(adapter,
                                fixture,
                                targets = "default",
                                geometry = NULL,
                                relation = NULL,
                                run = c("smoke", "compile"),
                                B = 19L,
                                R = 3L,
                                seed = 1L,
                                strict = TRUE,
                                fail_on_error = FALSE,
                                ...) {
  run <- match.arg(run)

  if (is.character(adapter)) {
    if (length(adapter) != 1L || is.na(adapter) || !nzchar(adapter)) {
      stop("`adapter` must be a non-empty registered adapter id or a multifer_adapter.",
           call. = FALSE)
    }
    adapter <- get_infer_adapter(adapter)
  }
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter or registered adapter id.",
         call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(fail_on_error) || length(fail_on_error) != 1L ||
      is.na(fail_on_error)) {
    stop("`fail_on_error` must be TRUE or FALSE.", call. = FALSE)
  }

  caps <- adapter_capabilities(adapter)
  pairs <- unique(caps[, c("geometry", "relation"), drop = FALSE])
  if (!is.null(geometry)) {
    if (!is.character(geometry) || length(geometry) == 0L || anyNA(geometry)) {
      stop("`geometry` must be NULL or a non-empty character vector.",
           call. = FALSE)
    }
    pairs <- pairs[pairs$geometry %in% geometry, , drop = FALSE]
  }
  if (!is.null(relation)) {
    if (!is.character(relation) || length(relation) == 0L || anyNA(relation)) {
      stop("`relation` must be NULL or a non-empty character vector.",
           call. = FALSE)
    }
    pairs <- pairs[pairs$relation %in% relation, , drop = FALSE]
  }
  if (nrow(pairs) == 0L) {
    stop("No adapter capability pairs matched the requested filters.",
         call. = FALSE)
  }

  rows <- vector("list", 0L)
  for (i in seq_len(nrow(pairs))) {
    geom <- pairs$geometry[[i]]
    rel <- pairs$relation[[i]]
    pair_targets <- .adapter_check_targets(caps, geom, rel, targets)

    for (target in pair_targets) {
      data_i <- .adapter_check_fixture_for(fixture, geom, rel)
      rows[[length(rows) + 1L]] <- .check_infer_adapter_one(
        adapter = adapter,
        data = data_i,
        geometry = geom,
        relation = rel,
        target = target,
        run = run,
        B = B,
        R = R,
        seed = seed,
        strict = strict,
        ...
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("multifer_adapter_check", class(out))

  if (isTRUE(fail_on_error) && any(!out$passed)) {
    failed <- out[!out$passed, , drop = FALSE]
    details <- sprintf(
      "%s/%s:%s failed at %s (%s)",
      failed$geometry, failed$relation, failed$target,
      failed$failed_stage, failed$detail
    )
    stop(
      paste0(
        "Adapter '", adapter$adapter_id,
        "' failed contract check(s): ",
        paste(details, collapse = "; ")
      ),
      call. = FALSE
    )
  }

  out
}

#' Check adapter-owned feature evidence declarations
#'
#' `check_feature_evidence_adapter()` smoke-tests an adapter's optional
#' `feature_stat_spec()` declarations. It is deliberately separate from
#' [check_infer_adapter()] because feature evidence is currently a sidecar
#' layer, not a frozen `infer_result` target.
#'
#' @param adapter A `multifer_adapter` object or registered adapter id.
#' @param fit A fitted adapter object used for statistic extraction.
#' @param data Original data passed to adapter feature-stat hooks.
#' @param units A `multifer_units` table.
#' @param domains Optional character vector of domains. When `NULL`, domains
#'   are resolved from the adapter and data.
#' @param fail_on_error Logical. If `TRUE`, stop when any declaration fails.
#'
#' @return A data frame with class `multifer_feature_evidence_adapter_check`.
#' @export
check_feature_evidence_adapter <- function(adapter,
                                           fit,
                                           data,
                                           units,
                                           domains = NULL,
                                           fail_on_error = FALSE) {
  if (is.character(adapter)) {
    if (length(adapter) != 1L || is.na(adapter) || !nzchar(adapter)) {
      stop("`adapter` must be a non-empty registered adapter id or a multifer_adapter.",
           call. = FALSE)
    }
    adapter <- get_infer_adapter(adapter)
  }
  if (!is_infer_adapter(adapter)) {
    stop("`adapter` must be a multifer_adapter or registered adapter id.",
         call. = FALSE)
  }
  if (!inherits(units, "multifer_units")) {
    stop("`units` must be a multifer_units table.", call. = FALSE)
  }
  if (!is.logical(fail_on_error) || length(fail_on_error) != 1L ||
      is.na(fail_on_error)) {
    stop("`fail_on_error` must be TRUE or FALSE.", call. = FALSE)
  }

  spec <- .feature_evidence_spec_for_check(adapter, fit, data)
  if (nrow(spec) == 0L) {
    out <- data.frame(
      adapter_id = character(0), statistic = character(0),
      scope = character(0), orientation = character(0),
      passed = logical(0), failed_stage = character(0),
      detail = character(0), stringsAsFactors = FALSE
    )
    class(out) <- c("multifer_feature_evidence_adapter_check", class(out))
    return(out)
  }

  .validate_feature_stat_spec(spec)
  if (is.null(domains)) {
    geom_kind <- adapter$shape_kinds[[1L]] %||% "adapter"
    domains <- .adapter_domains(adapter, fit = fit, data = data,
                                geom_kind = geom_kind)
  }

  rows <- lapply(seq_len(nrow(spec)), function(i) {
    .check_feature_evidence_spec_row(adapter, fit, data, units, domains,
                                     spec[i, , drop = FALSE])
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("multifer_feature_evidence_adapter_check", class(out))

  if (isTRUE(fail_on_error) && any(!out$passed)) {
    failed <- out[!out$passed, , drop = FALSE]
    details <- sprintf(
      "%s/%s failed at %s (%s)",
      failed$statistic, failed$scope, failed$failed_stage, failed$detail
    )
    stop(
      paste0(
        "Adapter '", adapter$adapter_id,
        "' failed feature evidence check(s): ",
        paste(details, collapse = "; ")
      ),
      call. = FALSE
    )
  }

  out
}

.feature_evidence_spec_for_check <- function(adapter, fit, data) {
  if (is.function(adapter$feature_stat_spec)) {
    spec <- adapter$feature_stat_spec(fit, data = data)
    if (!is.data.frame(spec)) {
      stop("`feature_stat_spec()` must return a data frame.", call. = FALSE)
    }
    return(spec)
  }
  data.frame(
    statistic = character(0), scope = character(0),
    orientation = character(0), stringsAsFactors = FALSE
  )
}

.validate_feature_stat_spec <- function(spec) {
  required <- c("statistic", "scope", "orientation")
  missing <- setdiff(required, names(spec))
  if (length(missing) > 0L) {
    stop(sprintf(
      "`feature_stat_spec()` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  bad_scope <- setdiff(spec$scope, c("unit", "aggregate", "both", "block",
                                    "score", "subspace"))
  if (length(bad_scope) > 0L) {
    stop(sprintf("unknown feature evidence scope(s): %s",
                 paste(unique(bad_scope), collapse = ", ")), call. = FALSE)
  }
  bad_orientation <- setdiff(spec$orientation,
                             c("signed", "unsigned", "subspace_norm",
                               "adapter_invariant"))
  if (length(bad_orientation) > 0L) {
    stop(sprintf("unknown feature evidence orientation(s): %s",
                 paste(unique(bad_orientation), collapse = ", ")),
         call. = FALSE)
  }
  invisible(spec)
}

.check_feature_evidence_spec_row <- function(adapter,
                                             fit,
                                             data,
                                             units,
                                             domains,
                                             spec_row) {
  row <- data.frame(
    adapter_id = adapter$adapter_id,
    statistic = as.character(spec_row$statistic[[1L]]),
    scope = as.character(spec_row$scope[[1L]]),
    orientation = as.character(spec_row$orientation[[1L]]),
    passed = FALSE,
    failed_stage = "statistic",
    detail = NA_character_,
    stringsAsFactors = FALSE
  )

  result <- tryCatch({
    if (is.function(adapter$feature_stat)) {
      members <- attr(units, "members")
      unit_row <- 1L
      stat <- adapter$feature_stat(
        fit,
        data = data,
        domain = domains[[1L]],
        unit_id = units$unit_id[[unit_row]],
        members = as.integer(members[[unit_row]]),
        statistic = row$statistic,
        orientation = row$orientation,
        scope = row$scope
      )
      .check_feature_stat_vector(stat)
      "feature_stat executable"
    } else if (is.function(adapter$feature_evidence_action)) {
      out <- feature_evidence_from_adapter(
        adapter = adapter,
        fit = fit,
        data = data,
        units = units,
        statistic = row$statistic,
        orientation = row$orientation,
        R = 2L,
        seed = 1L
      )
      if (!inherits(out, "multifer_feature_evidence")) {
        stop("native action did not return feature evidence.", call. = FALSE)
      }
      "feature_evidence_action executable"
    } else {
      stop("adapter declares feature evidence but provides no `feature_stat` or `feature_evidence_action` hook.",
           call. = FALSE)
    }
  }, error = function(e) e)

  if (inherits(result, "error")) {
    row$detail <- conditionMessage(result)
    return(row)
  }
  row$passed <- TRUE
  row$failed_stage <- NA_character_
  row$detail <- result
  row
}

.check_feature_stat_vector <- function(x) {
  if (!is.numeric(x) || is.matrix(x) || length(x) == 0L ||
      any(!is.finite(x))) {
    stop("feature statistic must return a non-empty finite numeric vector.",
         call. = FALSE)
  }
  invisible(x)
}

.adapter_check_targets <- function(caps, geometry, relation, targets) {
  supported <- setdiff(capabilities_for(caps, geometry, relation),
                       "variable_significance")
  if (length(supported) == 0L) {
    stop(sprintf(
      "Adapter declares no checkable targets for (%s, %s).",
      geometry, relation
    ), call. = FALSE)
  }

  if (identical(targets, "default") ||
      (is.character(targets) && length(targets) == 1L &&
       identical(targets, "default"))) {
    return(supported)
  }
  if (!is.character(targets) || length(targets) == 0L || anyNA(targets)) {
    stop("`targets` must be \"default\" or a non-empty character vector.",
         call. = FALSE)
  }

  selected <- intersect(unique(targets), supported)
  missing <- setdiff(unique(targets), supported)
  if (length(missing) > 0L) {
    stop(sprintf(
      "Target(s) not supported for (%s, %s): %s.",
      geometry, relation, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  selected
}

.adapter_check_fixture_for <- function(fixture, geometry, relation) {
  if (is.list(fixture) && !is.null(names(fixture))) {
    nm <- names(fixture)
    exact <- paste(geometry, relation, sep = "/")
    if (exact %in% nm) {
      return(fixture[[exact]])
    }
    if (geometry %in% nm) {
      return(fixture[[geometry]])
    }
    if ("default" %in% nm) {
      return(fixture[["default"]])
    }
  }
  fixture
}

.check_infer_adapter_one <- function(adapter,
                                     data,
                                     geometry,
                                     relation,
                                     target,
                                     run,
                                     B,
                                     R,
                                     seed,
                                     strict,
                                     ...) {
  row <- data.frame(
    adapter_id = adapter$adapter_id,
    geometry = geometry,
    relation = relation,
    target = target,
    compile_passed = FALSE,
    checks_passed = NA,
    infer_passed = if (identical(run, "smoke")) FALSE else NA,
    passed = FALSE,
    failed_stage = "compile",
    detail = NA_character_,
    stringsAsFactors = FALSE
  )

  recipe <- tryCatch(
    infer_recipe(
      geometry = geometry,
      relation = relation,
      targets = target,
      adapter = adapter,
      strict = strict
    ),
    error = function(e) e
  )
  if (inherits(recipe, "error")) {
    row$detail <- conditionMessage(recipe)
    return(row)
  }
  row$compile_passed <- TRUE

  checks <- tryCatch(
    run_adapter_checks(adapter, data, recipe = recipe, strict = strict),
    error = function(e) e
  )
  if (inherits(checks, "error")) {
    row$checks_passed <- FALSE
    row$failed_stage <- "checks"
    row$detail <- conditionMessage(checks)
    return(row)
  }
  row$checks_passed <- TRUE

  if (identical(run, "compile")) {
    row$passed <- TRUE
    row$failed_stage <- NA_character_
    row$detail <- "compiled and checked"
    return(row)
  }

  result <- tryCatch(
    infer(
      adapter = adapter,
      data = data,
      recipe = recipe,
      strict = strict,
      B = B,
      R = R,
      seed = seed,
      ...
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    row$infer_passed <- FALSE
    row$failed_stage <- "infer"
    row$detail <- conditionMessage(result)
    return(row)
  }
  if (!is_infer_result(result)) {
    row$infer_passed <- FALSE
    row$failed_stage <- "infer"
    row$detail <- "infer() did not return a multifer infer_result."
    return(row)
  }

  row$infer_passed <- TRUE
  row$passed <- TRUE
  row$failed_stage <- NA_character_
  row$detail <- "smoke inference returned infer_result"
  row
}

#' @export
print.multifer_adapter_check <- function(x, ...) {
  n <- nrow(x)
  n_pass <- sum(x$passed)
  cat("<multifer_adapter_check: ", n_pass, "/", n, " passed>\n", sep = "")
  if (n > 0L) {
    print.data.frame(x, row.names = FALSE)
  }
  invisible(x)
}
