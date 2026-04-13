#' Sequential-MC budget allocator
#'
#' A lightweight state machine that holds a global Monte-Carlo budget and
#' lets the sequential-ladder driver draw per-rung caps that cannot exceed
#' the remaining pool. Budget released by rungs that decide early rolls
#' forward to the next undecided rung -- the "budget allocator" framing
#' from §28: easy components release unused draws to the boundary
#' component.
#'
#' Because Vitale-style deflation forces rungs to run sequentially (rung
#' `a` depends on deflations through rung `a-1`), the allocator's current
#' role is simpler than the fully parallel §28 sketch: it is a rolling
#' pool that each rung checks out a cap from and reports back unused
#' remainder to.
#'
#' @param B_total Positive integer. The global budget for the run.
#' @param per_rung_cap Optional positive integer. Hard cap on how many
#'   draws any single rung can check out in one call. `NULL` means no
#'   cap beyond the pool balance.
#'
#' @return A `multifer_budget_allocator` object: an environment with
#'   methods `checkout(request)`, `return_unused(n)`, `remaining()`,
#'   `total()`, `used()`, and `schedule()`.
#'
#' @export
mc_budget_allocator <- function(B_total, per_rung_cap = NULL) {
  if (!is.numeric(B_total) || length(B_total) != 1L || is.na(B_total) ||
      B_total != as.integer(B_total) || B_total < 1L) {
    stop("`B_total` must be a positive integer scalar.", call. = FALSE)
  }
  B_total <- as.integer(B_total)

  if (!is.null(per_rung_cap)) {
    if (!is.numeric(per_rung_cap) || length(per_rung_cap) != 1L ||
        is.na(per_rung_cap) || per_rung_cap < 1L ||
        per_rung_cap != as.integer(per_rung_cap)) {
      stop("`per_rung_cap` must be a positive integer scalar or NULL.",
           call. = FALSE)
    }
    per_rung_cap <- as.integer(per_rung_cap)
  }

  state <- new.env(parent = emptyenv())
  state$B_total      <- B_total
  state$remaining    <- B_total
  state$per_rung_cap <- per_rung_cap
  state$schedule     <- integer(0L)     # actual draws per rung
  state$allocations  <- integer(0L)     # caps granted per rung

  state$checkout <- function(request = NULL) {
    # Returns an integer B_max for the next rung. `request` is a preferred
    # size (e.g. "B I'd like if the pool were infinite"); if NULL, use
    # per_rung_cap or remaining.
    cap <- state$remaining
    if (!is.null(state$per_rung_cap)) {
      cap <- min(cap, state$per_rung_cap)
    }
    if (!is.null(request)) {
      if (!is.numeric(request) || length(request) != 1L || is.na(request) ||
          request < 1L || request != as.integer(request)) {
        stop("`request` must be a positive integer scalar or NULL.",
             call. = FALSE)
      }
      cap <- min(cap, as.integer(request))
    }
    state$allocations <- c(state$allocations, cap)
    cap
  }

  state$record_use <- function(allocated, used) {
    if (!is.numeric(allocated) || length(allocated) != 1L || allocated < 0L) {
      stop("`allocated` must be a non-negative integer scalar.", call. = FALSE)
    }
    if (!is.numeric(used) || length(used) != 1L || used < 0L ||
        used > allocated) {
      stop("`used` must be a non-negative integer scalar <= `allocated`.",
           call. = FALSE)
    }
    allocated <- as.integer(allocated); used <- as.integer(used)
    state$remaining <- state$remaining - used
    state$schedule  <- c(state$schedule, used)
    invisible(state)
  }

  state$remaining_fn <- function() state$remaining
  state$used_fn      <- function() state$B_total - state$remaining

  class(state) <- "multifer_budget_allocator"
  state
}

#' @export
print.multifer_budget_allocator <- function(x, ...) {
  cat("<multifer_budget_allocator>\n")
  cat(sprintf("  total:     %d\n", x$B_total))
  cat(sprintf("  used:      %d\n", x$used_fn()))
  cat(sprintf("  remaining: %d\n", x$remaining_fn()))
  cat(sprintf("  rungs served: %d\n", length(x$schedule)))
  invisible(x)
}
