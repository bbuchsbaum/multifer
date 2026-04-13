#' Sequential Monte Carlo p-value with Besag-Clifford early stopping
#'
#' Draws null statistics in batches and stops as soon as the decision is
#' forced, giving an anytime-valid p-value that matches `mc_p_value()`'s
#' Phipson-Smyth formula `(1 + r) / (draws_used + 1)` on the actual draws
#' taken.
#'
#' Classical Besag-Clifford early termination fires on the non-rejection
#' side only: when the exceedance count `r` reaches `h`, the current
#' p-value has already crossed above `alpha`, so we stop and report
#' `(1 + r) / (drawn + 1)`. Otherwise the loop runs to `B_max` and
#' returns the standard Phipson-Smyth p-value `(1 + r) / (B_max + 1)`.
#' The rejection path is kept at full budget because the reported
#' numeric p-value is consumed by downstream calibration checks that
#' compare against the fixed-B Phipson-Smyth value -- early rejection
#' would change the reported number without changing the decision.
#'
#' The default threshold `h = ceiling(alpha * (B_max + 1))` is the
#' smallest `r` at which a Phipson-Smyth p-value `(1 + r)/(B_max + 1)`
#' would exceed `alpha`, making Besag-Clifford non-rejection exactly
#' align with the fixed-B rule.
#'
#' @param observed_stat Numeric scalar. Observed test statistic.
#' @param null_gen_fn Function (thunk) returning one finite numeric
#'   null draw per call.
#' @param B_max Positive integer. Hard cap on null draws.
#' @param alpha Significance threshold in (0, 1).
#' @param h Positive integer, the Besag-Clifford threshold. If `NULL`,
#'   defaults to `ceiling(alpha * (B_max + 1))`.
#' @param batch_size Positive integer. Draws per batch. Smaller batches
#'   stop earlier (finer granularity); larger batches amortize RNG
#'   overhead. Default `32L`.
#' @param alternative One of `"greater"`, `"less"`, `"two_sided"`.
#' @param seed Optional integer. If supplied, the RNG state is saved
#'   and restored on exit.
#'
#' @return A list with elements `p_value`, `r`, `drawn`, `B_max`, `h`,
#'   `stop_reason` (one of `"non_reject_early"`, `"reject_early"`,
#'   `"exhausted"`), `batch_schedule`, `null_values`, `mc_se`,
#'   `alternative`, `seed`.
#'
#' @seealso [mc_p_value()] for the fixed-B primitive.
#' @export
mc_sequential_bc <- function(observed_stat,
                             null_gen_fn,
                             B_max,
                             alpha,
                             h = NULL,
                             batch_size = 32L,
                             alternative = c("greater", "less", "two_sided"),
                             seed = NULL) {

  ## --- input validation ------------------------------------------------------

  if (!is.function(null_gen_fn)) {
    stop("`null_gen_fn` must be a function.", call. = FALSE)
  }
  if (!is.numeric(observed_stat) || length(observed_stat) != 1L ||
      !is.finite(observed_stat)) {
    stop("`observed_stat` must be a single finite numeric value.",
         call. = FALSE)
  }
  if (!is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) ||
      B_max != as.integer(B_max) || B_max < 1L) {
    stop("`B_max` must be a positive integer scalar.", call. = FALSE)
  }
  B_max <- as.integer(B_max)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) {
    stop("`batch_size` must be a positive integer scalar.", call. = FALSE)
  }
  batch_size <- as.integer(batch_size)
  alternative <- match.arg(alternative)

  if (is.null(h)) {
    h <- as.integer(ceiling(alpha * (B_max + 1L)))
  } else {
    if (!is.numeric(h) || length(h) != 1L || is.na(h) ||
        h != as.integer(h) || h < 1L) {
      stop("`h` must be a positive integer scalar or NULL.", call. = FALSE)
    }
    h <- as.integer(h)
  }

  ## --- seed handling ---------------------------------------------------------

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single integer value or NULL.", call. = FALSE)
    }
    seed <- as.integer(seed)
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  ## --- batched drawing loop --------------------------------------------------

  null_values    <- numeric(0L)
  batch_schedule <- integer(0L)
  r              <- 0L
  drawn          <- 0L
  stop_reason    <- "exhausted"

  extreme_count <- function(v) {
    switch(alternative,
      greater   = sum(v >= observed_stat),
      less      = sum(v <= observed_stat),
      two_sided = sum(abs(v) >= abs(observed_stat))
    )
  }

  while (drawn < B_max) {
    this_batch <- min(batch_size, B_max - drawn)
    batch_vals <- numeric(this_batch)
    for (i in seq_len(this_batch)) {
      v <- null_gen_fn()
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v)) {
        stop(sprintf(
          "null draw %d: `null_gen_fn` must return a single finite numeric.",
          drawn + i
        ), call. = FALSE)
      }
      batch_vals[i] <- v
    }
    null_values    <- c(null_values, batch_vals)
    batch_schedule <- c(batch_schedule, this_batch)
    drawn          <- drawn + this_batch
    r              <- r + as.integer(extreme_count(batch_vals))

    ## --- early non-rejection: r >= h ----------------------------------------
    if (r >= h) {
      stop_reason <- "non_reject_early"
      break
    }
  }

  ## --- report ----------------------------------------------------------------

  p_value <- if (stop_reason == "non_reject_early") {
    (1 + r) / (drawn + 1L)
  } else {
    (1 + r) / (B_max + 1L)
  }
  mc_se <- sqrt(p_value * (1 - p_value) / drawn)

  list(
    p_value        = p_value,
    r              = r,
    drawn          = drawn,
    B_max          = B_max,
    h              = h,
    stop_reason    = stop_reason,
    batch_schedule = batch_schedule,
    null_values    = null_values,
    mc_se          = mc_se,
    alternative    = alternative,
    seed           = seed
  )
}
