# Executable validity contracts
#
# Each `multifer_adapter` carries a `checked_assumptions` list. Every
# entry is:
#
#   list(name   = <character>,
#        check  = function(data, ...),
#        detail = <character>)
#
# The `check` function is called before the engine runs with the raw
# `data` argument the user passed to `infer()` (a matrix for oneblock,
# a list(X, Y) for cross). It must return one of:
#
#   - TRUE  -- assumption is satisfied, silently pass
#   - FALSE -- assumption violated; if strict, infer() errors; else,
#     the violation is recorded in result$assumptions$checked
#   - a character scalar -- treated as FALSE with the string as the
#     violation detail (overrides the adapter's default detail)
#   - a structured list(passed = logical, detail = character) for full
#     control
#
# Errors thrown by the check function are captured and reported as
# violations named "<name>: <error message>".
#
# The output of `run_adapter_checks()` is consumed by `infer_assumptions()`
# as its `checked` argument, matching the schema documented in
# R/infer_result.R (each entry has `passed` and `detail`).

.normalize_check_result <- function(value, default_detail) {
  if (is.logical(value) && length(value) == 1L && !is.na(value)) {
    return(list(passed = as.logical(value), detail = default_detail))
  }
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    return(list(passed = FALSE, detail = value))
  }
  if (is.list(value) &&
      !is.null(value$passed) &&
      is.logical(value$passed) &&
      length(value$passed) == 1L) {
    return(list(
      passed = as.logical(value$passed),
      detail = if (!is.null(value$detail)) as.character(value$detail) else default_detail
    ))
  }
  # Anything else is treated as a violation with a diagnostic detail.
  list(
    passed = FALSE,
    detail = paste0(default_detail,
                    " (check returned unrecognized result of class ",
                    paste(class(value), collapse = "/"), ")")
  )
}

#' Run an adapter's checked assumptions against data.
#'
#' @param adapter A `multifer_adapter`.
#' @param data The raw `data` argument passed to `infer()`.
#' @param recipe Optional compiled `multifer_infer_recipe`. When supplied,
#'   check functions receive it via `recipe = ` and may use it for
#'   relation- or design-specific validity logic.
#' @param strict Logical. When `TRUE`, any failed check raises an error
#'   with the concatenated details; when `FALSE`, violations are
#'   returned in the results list for downstream reporting.
#'
#' @return A named list. Each element carries `name`, `passed`
#'   (logical), and `detail` (character). Names come from each check
#'   entry's `name` field.
#'
#' @keywords internal
#' @noRd
run_adapter_checks <- function(adapter, data, recipe = NULL, strict = TRUE) {
  checks <- adapter$checked_assumptions
  if (length(checks) == 0L) {
    return(list())
  }

  results <- vector("list", length(checks))
  for (i in seq_along(checks)) {
    entry <- checks[[i]]
    if (!is.list(entry) || is.null(entry$check) || !is.function(entry$check)) {
      stop(sprintf(
        "adapter %s: checked_assumptions[[%d]] is missing a $check function.",
        adapter$adapter_id, i
      ), call. = FALSE)
    }
    name <- if (!is.null(entry$name)) as.character(entry$name) else sprintf("check_%d", i)
    default_detail <- if (!is.null(entry$detail)) as.character(entry$detail) else name

    result <- tryCatch(
      .normalize_check_result(
        entry$check(data, recipe = recipe, adapter = adapter),
        default_detail
      ),
      error = function(e) {
        list(
          passed = FALSE,
          detail = paste0(default_detail, " (error: ", conditionMessage(e), ")")
        )
      }
    )
    result$name <- name
    results[[i]] <- result
  }
  names(results) <- vapply(results, `[[`, character(1L), "name")

  if (isTRUE(strict)) {
    failed <- !vapply(results, function(r) isTRUE(r$passed), logical(1L))
    if (any(failed)) {
      msg <- paste0(
        "adapter ", adapter$adapter_id,
        " failed validity check(s): ",
        paste(
          vapply(results[failed], function(r) sprintf("%s (%s)", r$name, r$detail),
                 character(1L)),
          collapse = "; "
        )
      )
      stop(msg, call. = FALSE)
    }
  }

  results
}
