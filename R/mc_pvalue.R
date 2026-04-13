#' Fixed-B Monte Carlo p-value
#'
#' Computes a Monte Carlo p-value using the Phipson-Smyth formula
#' `p = (1 + r) / (B + 1)`, where `r` is the number of null draws at least
#' as extreme as the observed statistic. The full null draw record is
#' returned so that downstream sequential-stopping logs (Phase 1.5) can
#' consume it without re-running the simulation.
#'
#' @param observed_stat Numeric scalar. The observed test statistic. Must be
#'   finite (not NA, NaN, or Inf).
#' @param null_gen_fn Function (thunk). Called `B` times with no arguments;
#'   each call must return a single finite numeric value.
#' @param B Positive integer scalar. Number of Monte Carlo null draws.
#' @param alternative Character scalar. One of `"greater"`, `"less"`, or
#'   `"two_sided"`. Default is `"greater"`. For `"two_sided"`, the count
#'   `r` is `sum(abs(null_values) >= abs(observed_stat))` -- the simpler
#'   magnitude comparison, which requires no re-centering and is
#'   appropriate when both the null distribution and the observed statistic
#'   are already on a symmetric scale.
#' @param seed Integer scalar or NULL. If not NULL, `set.seed(seed)` is
#'   called before drawing and the previous RNG state is restored on exit,
#'   leaving the caller's RNG stream unchanged.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`p_value`}{Numeric scalar. `(1 + r) / (B + 1)`.}
#'     \item{`r`}{Integer. Count of null draws as extreme as or more
#'       extreme than `observed_stat`.}
#'     \item{`B`}{Integer. Number of null draws, as passed.}
#'     \item{`observed_stat`}{Numeric. As passed.}
#'     \item{`null_values`}{Numeric vector of length `B`. All null draws
#'       in draw order. Retained for diagnostics and Phase 1.5 logging.}
#'     \item{`mc_se`}{Numeric. Monte Carlo standard error
#'       `sqrt(p_value * (1 - p_value) / B)`.}
#'     \item{`alternative`}{Character. The chosen alternative.}
#'     \item{`seed`}{Integer or NULL. As passed.}
#'   }
#'
#' @details
#' The two-sided counting rule used here counts the number of null
#' draws whose absolute value is at least as large as the absolute
#' value of the observed statistic -- that is,
#' \code{sum(abs(null_values) >= abs(observed_stat))}.
#' This avoids computing a null mean and is appropriate when the null
#' distribution is symmetric around zero. If the null is not symmetric,
#' callers should use `"greater"` or `"less"` explicitly.
#'
#' @references
#' Phipson B and Smyth GK (2010). "Permutation p-values should never be
#' zero: calculating exact p-values when permutations are randomly drawn."
#' Statistical Applications in Genetics and Molecular Biology, 9(1).
#'
#' @export
mc_p_value <- function(observed_stat,
                       null_gen_fn,
                       B,
                       alternative = c("greater", "less", "two_sided"),
                       seed = NULL) {

  ## --- input validation ---------------------------------------------------

  if (!is.function(null_gen_fn)) {
    stop("`null_gen_fn` must be a function.", call. = FALSE)
  }

  if (!is.numeric(observed_stat) || length(observed_stat) != 1L ||
      !is.finite(observed_stat)) {
    stop("`observed_stat` must be a single finite numeric value.",
         call. = FALSE)
  }

  if (!is.numeric(B) || length(B) != 1L || is.na(B) ||
      B != as.integer(B) || B < 1L) {
    stop("`B` must be a positive integer scalar.", call. = FALSE)
  }
  B <- as.integer(B)

  alternative <- match.arg(alternative)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop("`seed` must be a single integer value or NULL.", call. = FALSE)
    }
    seed <- as.integer(seed)
  }

  ## --- seed handling: save current RNG state, restore on exit -------------

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

  ## --- draw null distribution ---------------------------------------------

  null_values <- numeric(B)
  for (i in seq_len(B)) {
    val <- null_gen_fn()
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop(sprintf(
        "null draw %d: `null_gen_fn` must return a single finite numeric value.",
        i
      ), call. = FALSE)
    }
    null_values[i] <- val
  }

  ## --- compute r and p-value ----------------------------------------------

  r <- switch(alternative,
    greater   = sum(null_values >= observed_stat),
    less      = sum(null_values <= observed_stat),
    two_sided = sum(abs(null_values) >= abs(observed_stat))
  )

  p_value <- (1L + r) / (B + 1L)
  mc_se   <- sqrt(p_value * (1 - p_value) / B)

  ## --- return -------------------------------------------------------------

  list(
    p_value       = p_value,
    r             = r,
    B             = B,
    observed_stat = observed_stat,
    null_values   = null_values,
    mc_se         = mc_se,
    alternative   = alternative,
    seed          = seed
  )
}
