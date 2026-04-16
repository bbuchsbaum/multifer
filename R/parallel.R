#' Parallel execution backend for Phase 1.5
#'
#' Thin wrappers around `mirai::daemons()` / `mirai::mirai_map()` that give
#' the bootstrap loop and any other embarrassingly-parallel section of the
#' pipeline a single entry point. Daemons are spawned lazily on first use
#' and never auto-shut-down -- callers that want to tear down can call
#' `multifer_parallel_shutdown()` explicitly.
#'
#' Reproducibility: when `seeds` are supplied, the parallel map sets
#' each task's RNG to its pre-assigned seed before calling `FUN`.
#' Parallel and sequential runs over the same `(X, seeds)` are bit-
#' identical at task granularity, though they diverge from a single-
#' stream sequential loop that consumes the RNG incrementally.
#'
#' @name multifer_parallel
NULL

.multifer_parallel_state <- new.env(parent = emptyenv())
.multifer_parallel_state$daemons_up <- FALSE
.multifer_parallel_state$n_workers  <- 0L

#' Ensure a mirai daemon pool is running.
#'
#' Idempotent. First call spawns `n` daemons; subsequent calls with the
#' same `n` are no-ops. Calls with a different `n` reshape the pool.
#' `n` defaults to `max(1, parallel::detectCores() - 1)` or the value of
#' `getOption("multifer.workers")` when set.
#'
#' @param n Integer or NULL. Number of daemon workers.
#' @return Integer, the number of active daemons.
#' @export
multifer_parallel_init <- function(n = NULL) {
  if (is.null(n)) {
    n <- getOption("multifer.workers",
                   default = max(1L, parallel::detectCores() - 1L))
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != as.integer(n) || n < 1L) {
    stop("`n` must be a positive integer scalar or NULL.", call. = FALSE)
  }
  n <- as.integer(n)

  if (isTRUE(.multifer_parallel_state$daemons_up) &&
      isTRUE(.multifer_parallel_state$n_workers == n)) {
    return(invisible(n))
  }

  mirai::daemons(n)
  mirai::everywhere({ requireNamespace("multifer", quietly = TRUE) })
  .multifer_parallel_state$daemons_up <- TRUE
  .multifer_parallel_state$n_workers  <- n
  invisible(n)
}

#' Shut down the mirai daemon pool.
#'
#' Safe to call when no pool is up.
#'
#' @return `invisible(NULL)`.
#' @export
multifer_parallel_shutdown <- function() {
  if (isTRUE(.multifer_parallel_state$daemons_up)) {
    try(mirai::daemons(0L), silent = TRUE)
    .multifer_parallel_state$daemons_up <- FALSE
    .multifer_parallel_state$n_workers  <- 0L
  }
  invisible(NULL)
}

#' Deterministic per-task seed generation.
#'
#' Given a master seed and a task count `n`, returns a reproducible
#' integer vector of length `n`. Uses the caller's RNG stream when
#' `master_seed` is `NULL` (with the state saved and restored) and
#' `set.seed(master_seed)` otherwise.
#'
#' @param n Positive integer.
#' @param master_seed Integer or NULL.
#' @return Integer vector of length `n`.
#' @export
multifer_task_seeds <- function(n, master_seed = NULL) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L || n != as.integer(n)) {
    stop("`n` must be a positive integer scalar.", call. = FALSE)
  }
  n <- as.integer(n)

  cap <- .Machine$integer.max - 1L

  if (is.null(master_seed)) {
    return(sample.int(cap, n))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(master_seed))
  sample.int(cap, n)
}

#' Parallel `lapply` with optional per-task seed threading.
#'
#' @param X List or vector to iterate over.
#' @param FUN Function; receives one element of `X` per call, plus any
#'   `...` arguments.
#' @param ... Additional arguments passed unchanged to `FUN`.
#' @param seeds Optional integer vector of length `length(X)`. If
#'   supplied, `set.seed(seeds[i])` is called before the i-th task.
#' @param backend One of `"auto"`, `"mirai"`, `"sequential"`. `"auto"`
#'   uses mirai when a daemon pool is available or can be initialized,
#'   otherwise falls back to sequential.
#' @return A list of length `length(X)`.
#' @export
multifer_parallel_lapply <- function(X, FUN, ..., seeds = NULL,
                                     backend = c("auto", "mirai",
                                                 "sequential")) {
  backend <- match.arg(backend)
  n <- length(X)
  rng_kind <- RNGkind()
  if (!is.null(seeds)) {
    if (length(seeds) != n) {
      stop("`seeds` must have the same length as `X`.", call. = FALSE)
    }
    seeds <- as.integer(seeds)
  }
  extra_args <- list(...)

  if (backend == "auto") {
    backend <- if (requireNamespace("mirai", quietly = TRUE)) "mirai" else "sequential"
  }

  if (backend == "sequential") {
    return(lapply(seq_along(X), function(i) {
      if (!is.null(seeds)) {
        do.call(RNGkind, as.list(rng_kind))
        set.seed(seeds[i])
      }
      do.call(FUN, c(list(X[[i]]), extra_args))
    }))
  }

  ## --- mirai path ---------------------------------------------------------
  multifer_parallel_init()
  worker_fn <- function(i, X, FUN, extra_args, seeds, rng_kind) {
    if (!is.null(seeds)) {
      do.call(RNGkind, as.list(rng_kind))
      set.seed(seeds[i])
    }
    do.call(FUN, c(list(X[[i]]), extra_args))
  }
  mm <- mirai::mirai_map(
    .x    = seq_along(X),
    .f    = worker_fn,
    .args = list(
      X          = X,
      FUN        = FUN,
      extra_args = extra_args,
      seeds      = seeds,
      rng_kind   = rng_kind
    )
  )
  mm[]
}
