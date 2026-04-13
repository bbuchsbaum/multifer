#' Sequential deflation test ladder driver
#'
#' Runs the Vitale-style nested hypothesis ladder, stopping at the
#' first non-rejection. Knows nothing about shape, statistic, or null
#' action -- those are provided by the caller.
#'
#' @param observed_stat_fn Function(step, deflated_data) -> numeric.
#'   Returns the observed test statistic at step \code{step} on the
#'   currently-deflated data.
#' @param null_stat_fn Function(step, deflated_data) -> numeric.
#'   Returns ONE null draw: applies the null action to \code{deflated_data}
#'   and computes the same statistic as \code{observed_stat_fn}.
#' @param deflate_fn Function(step, deflated_data) -> deflated_data.
#'   Returns the deflated data for the NEXT step (after removing the
#'   component just tested).
#' @param initial_data Whatever the caller's deflation operates on.
#' @param max_steps Integer, hard cap on how many ladder rungs to test.
#' @param B Per-rung cap on Monte Carlo draws. Also the default value used
#'   when `B_total` is not supplied.
#' @param B_total Optional integer. Global Monte Carlo budget shared across
#'   rungs via [mc_budget_allocator()]. Defaults to `B * max_steps` so the
#'   legacy fixed-B behavior is preserved when unset.
#' @param batch_size Positive integer, draws per Besag-Clifford batch
#'   inside a single rung. Default `32L`.
#' @param alpha Significance threshold for stop-at-first-non-rejection.
#'   Default \code{0.05}.
#' @param seed Integer or NULL. If not NULL, the RNG state is saved before
#'   the ladder begins and restored on exit, leaving the caller's stream
#'   unchanged.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{step_results}{List of length <= \code{max_steps}, one element
#'       per rung tested. Each element is a list with fields
#'       \code{step}, \code{observed_stat}, \code{p_value},
#'       \code{mc_se}, \code{r}, \code{B} (cap), \code{drawn}, \code{h},
#'       \code{stop_reason}, \code{batch_schedule}, \code{null_values},
#'       \code{selected} (logical, TRUE when p_value <= alpha).}
#'     \item{last_step_tested}{Integer. The highest rung that was evaluated.}
#'     \item{rejected_through}{Integer. The number of consecutive initial
#'       rejections, i.e., the estimated rank.}
#'     \item{deflated_data_final}{The last deflated_data produced before
#'       stopping.}
#'     \item{allocator}{The `multifer_budget_allocator` after the ladder
#'       finishes, exposing `$used_fn()`, `$remaining_fn()`, `$schedule`.}
#'     \item{total_draws_used}{Integer. Sum of draws across all rungs.}
#'     \item{batch_schedule}{Integer vector. Concatenated batch sizes
#'       across all rungs, suitable for `infer_mc()$batch_schedule`.}
#'     \item{stopping_boundary}{Character scalar. Describes the rule used,
#'       e.g. `"besag_clifford(h = <h>)"`.}
#'   }
#'
#' @export
ladder_driver <- function(observed_stat_fn,
                          null_stat_fn,
                          deflate_fn,
                          initial_data,
                          max_steps,
                          B,
                          B_total     = NULL,
                          batch_size  = 32L,
                          alpha       = 0.05,
                          seed        = NULL) {

  ## --- input validation -------------------------------------------------------

  if (!is.function(observed_stat_fn)) {
    stop("`observed_stat_fn` must be a function.", call. = FALSE)
  }
  if (!is.function(null_stat_fn)) {
    stop("`null_stat_fn` must be a function.", call. = FALSE)
  }
  if (!is.function(deflate_fn)) {
    stop("`deflate_fn` must be a function.", call. = FALSE)
  }
  if (!is.numeric(max_steps) || length(max_steps) != 1L ||
      is.na(max_steps) || max_steps < 1L || max_steps != as.integer(max_steps)) {
    stop("`max_steps` must be a positive integer scalar.", call. = FALSE)
  }
  max_steps <- as.integer(max_steps)

  if (!is.numeric(B) || length(B) != 1L || is.na(B) ||
      B != as.integer(B) || B < 1L) {
    stop("`B` must be a positive integer scalar.", call. = FALSE)
  }
  B <- as.integer(B)

  if (is.null(B_total)) {
    B_total <- B * max_steps
  }
  if (!is.numeric(B_total) || length(B_total) != 1L || is.na(B_total) ||
      B_total != as.integer(B_total) || B_total < 1L) {
    stop("`B_total` must be a positive integer scalar or NULL.", call. = FALSE)
  }
  B_total <- as.integer(B_total)

  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) {
    stop("`batch_size` must be a positive integer scalar.", call. = FALSE)
  }
  batch_size <- as.integer(batch_size)

  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single integer value or NULL.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  ## --- RNG state: save current state, restore on exit ------------------------

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  ## --- budget allocator ------------------------------------------------------

  allocator <- mc_budget_allocator(B_total = B_total, per_rung_cap = B)

  ## --- ladder loop ------------------------------------------------------------

  step_results     <- vector("list", max_steps)
  current_data     <- initial_data
  rejected_through <- 0L
  h_last           <- NA_integer_

  for (step in seq_len(max_steps)) {

    observed <- observed_stat_fn(step, current_data)

    cap <- allocator$checkout(request = B)
    if (cap < 1L) {
      # Pool exhausted: record an empty rung and stop testing.
      step_results[[step]] <- list(
        step           = step,
        observed_stat  = observed,
        p_value        = NA_real_,
        mc_se          = NA_real_,
        r              = 0L,
        B              = 0L,
        drawn          = 0L,
        h              = NA_integer_,
        stop_reason    = "budget_exhausted",
        batch_schedule = integer(0L),
        null_values    = numeric(0L),
        selected       = FALSE
      )
      break
    }

    mc_result <- mc_sequential_bc(
      observed_stat = observed,
      null_gen_fn   = function() null_stat_fn(step, current_data),
      B_max         = cap,
      alpha         = alpha,
      batch_size    = batch_size,
      alternative   = "greater",
      seed          = NULL   # driver owns the RNG stream
    )

    allocator$record_use(allocated = cap, used = mc_result$drawn)
    h_last <- mc_result$h

    selected <- mc_result$p_value <= alpha

    step_results[[step]] <- list(
      step           = step,
      observed_stat  = observed,
      p_value        = mc_result$p_value,
      mc_se          = mc_result$mc_se,
      r              = mc_result$r,
      B              = cap,
      drawn          = mc_result$drawn,
      h              = mc_result$h,
      stop_reason    = mc_result$stop_reason,
      batch_schedule = mc_result$batch_schedule,
      null_values    = mc_result$null_values,
      selected       = selected
    )

    if (!selected) {
      break
    }

    rejected_through <- step

    if (step < max_steps) {
      current_data <- deflate_fn(step, current_data)
    }
  }

  last_step <- step
  step_results <- step_results[seq_len(last_step)]

  total_draws_used <- as.integer(sum(vapply(
    step_results, function(sr) as.integer(sr$drawn %||% 0L), integer(1L)
  )))
  batch_schedule_all <- unlist(lapply(step_results, function(sr)
    if (is.null(sr$batch_schedule)) integer(0L) else sr$batch_schedule
  ))
  if (is.null(batch_schedule_all)) batch_schedule_all <- integer(0L)
  stopping_boundary <- sprintf("besag_clifford(h=%s)",
                               if (is.na(h_last)) "NA" else h_last)

  list(
    step_results        = step_results,
    last_step_tested    = last_step,
    rejected_through    = rejected_through,
    deflated_data_final = current_data,
    allocator           = allocator,
    total_draws_used    = total_draws_used,
    batch_schedule      = as.integer(batch_schedule_all),
    stopping_boundary   = stopping_boundary
  )
}

# Local %||% (base R >= 4.4 has it natively; this guards older versions).
`%||%` <- function(a, b) if (is.null(a)) b else a
