#' Benchmark registry and access functions
#'
#' Functions for listing, retrieving, and running the four locked
#' benchmark suites defined in Phase 0. The registry is populated at
#' package load time by `.onLoad()` in `zzz.R`.
#'
#' @name bench
NULL

#' List all registered benchmarks
#'
#' Returns a data.frame with one row per registered benchmark suite.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{name}{Character. Suite identifier.}
#'     \item{geometry}{Character. Geometry label (`"oneblock"` or
#'       `"cross"`).}
#'     \item{purpose}{Character. One-line description.}
#'     \item{locked}{Logical. Always `TRUE` for Phase 0 suites.}
#'   }
#'
#' @export
list_benchmarks <- function() {
  nms <- ls(.multifer_bench_registry)
  if (length(nms) == 0L) {
    return(data.frame(
      name     = character(0),
      geometry = character(0),
      purpose  = character(0),
      locked   = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(nms, function(nm) {
    e <- .multifer_bench_registry[[nm]]
    data.frame(
      name     = nm,
      geometry = e$geometry,
      purpose  = e$purpose,
      locked   = e$locked,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Retrieve a benchmark suite by name
#'
#' @param suite_name Character scalar. One of the names returned by
#'   [list_benchmarks()].
#'
#' @return A list with elements:
#'   \describe{
#'     \item{generator}{The generator function for the suite.}
#'     \item{targets}{The locked targets list for the suite.}
#'   }
#'
#' @export
get_benchmark <- function(suite_name) {
  if (!is.character(suite_name) || length(suite_name) != 1L) {
    stop("`suite_name` must be a single character string.", call. = FALSE)
  }
  if (!exists(suite_name, envir = .multifer_bench_registry, inherits = FALSE)) {
    known <- paste(ls(.multifer_bench_registry), collapse = ", ")
    stop(sprintf(
      "Unknown benchmark `%s`. Known suites: %s.",
      suite_name, known
    ), call. = FALSE)
  }
  e <- .multifer_bench_registry[[suite_name]]
  list(generator = e$generator, targets = e$targets)
}

#' Run a benchmark generator by name
#'
#' A convenience thin wrapper that calls the named suite's generator
#' function with any additional arguments forwarded via `...`.
#'
#' @param suite_name Character scalar. Suite name; see [list_benchmarks()].
#' @param ... Additional arguments forwarded to the generator function.
#'
#' @return Whatever the generator returns (a list with `X` and
#'   optionally `Y`, plus `meta`).
#'
#' @export
run_benchmark_generator <- function(suite_name, ...) {
  b <- get_benchmark(suite_name)
  b$generator(...)
}
