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
#' @param B Number of Monte Carlo draws per rung.
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
#'       \code{mc_se}, \code{r}, \code{B}, \code{null_values},
#'       \code{selected} (logical, TRUE when p_value <= alpha).}
#'     \item{last_step_tested}{Integer. The highest rung that was evaluated.}
#'     \item{rejected_through}{Integer. The number of consecutive initial
#'       rejections, i.e., the estimated rank. Zero when the first rung is
#'       not rejected.}
#'     \item{deflated_data_final}{The last deflated_data produced before
#'       stopping. Useful for diagnostics.}
#'   }
#'
#' @export
ladder_driver <- function(observed_stat_fn,
                          null_stat_fn,
                          deflate_fn,
                          initial_data,
                          max_steps,
                          B,
                          alpha = 0.05,
                          seed  = NULL) {

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

  ## --- ladder loop ------------------------------------------------------------

  step_results   <- vector("list", max_steps)
  current_data   <- initial_data
  rejected_through <- 0L

  for (step in seq_len(max_steps)) {

    observed <- observed_stat_fn(step, current_data)

    mc_result <- mc_p_value(
      observed_stat = observed,
      null_gen_fn   = function() null_stat_fn(step, current_data),
      B             = B,
      alternative   = "greater",
      seed          = NULL   # driver owns the RNG stream
    )

    selected <- mc_result$p_value <= alpha

    step_results[[step]] <- list(
      step          = step,
      observed_stat = observed,
      p_value       = mc_result$p_value,
      mc_se         = mc_result$mc_se,
      r             = mc_result$r,
      B             = B,
      null_values   = mc_result$null_values,
      selected      = selected
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
  # Trim unused slots.
  step_results <- step_results[seq_len(last_step)]

  list(
    step_results       = step_results,
    last_step_tested   = last_step,
    rejected_through   = rejected_through,
    deflated_data_final = current_data
  )
}
