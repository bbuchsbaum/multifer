#' Thin SVD cache
#'
#' A lightweight, fingerprint-keyed cache for thin SVDs. Used by Phase 1.5
#' to avoid recomputing the same decomposition across pipeline stages in a
#' single `infer()` call (the original fit, the "full observed roots" pass,
#' and any adapter hook that needs the bases).
#'
#' The cache is **per-context**: `svd_cache_new()` returns a fresh
#' environment that callers pass explicitly. A package-level default
#' cache is provided for callers that do not want to thread an
#' environment through.
#'
#' Fingerprints are computed via `digest::digest()` on `dim` + a
#' low-discrepancy sample of the matrix. This is cheap and correct for
#' idempotent re-use within a single `infer()` call. It is **not**
#' intended as a cryptographic hash or as a cross-session cache key.
#' Cache lifetimes are bounded by the enclosing `infer()` call.
#'
#' @name thin_svd_cache
NULL

.multifer_default_svd_cache <- new.env(parent = emptyenv())

#' Create a fresh thin-SVD cache context.
#'
#' @return A new environment with zeroed counters.
#' @export
svd_cache_new <- function() {
  e <- new.env(parent = emptyenv())
  e$store    <- list()
  e$lookups  <- 0L
  e$hits     <- 0L
  e$misses   <- 0L
  e
}

#' Reset a thin-SVD cache context in place.
#'
#' @param cache An environment produced by `svd_cache_new()`. If `NULL`
#'   the package-default cache is reset.
#' @return The cache environment, invisibly.
#' @export
svd_cache_reset <- function(cache = NULL) {
  if (is.null(cache)) cache <- .multifer_default_svd_cache
  cache$store    <- list()
  cache$lookups  <- 0L
  cache$hits     <- 0L
  cache$misses   <- 0L
  invisible(cache)
}

#' Number of cache hits.
#'
#' @param cache Cache environment. `NULL` uses the package default.
#' @return Integer.
#' @export
svd_cache_hits <- function(cache = NULL) {
  if (is.null(cache)) cache <- .multifer_default_svd_cache
  if (is.null(cache$hits)) 0L else as.integer(cache$hits)
}

#' Number of cache lookups (hits + misses).
#'
#' @param cache Cache environment. `NULL` uses the package default.
#' @return Integer.
#' @export
svd_cache_lookups <- function(cache = NULL) {
  if (is.null(cache)) cache <- .multifer_default_svd_cache
  if (is.null(cache$lookups)) 0L else as.integer(cache$lookups)
}

#' Hit rate over all lookups.
#'
#' @param cache Cache environment. `NULL` uses the package default.
#' @return Numeric in `[0, 1]`; `NA_real_` if there have been no lookups.
#' @export
svd_cache_rate <- function(cache = NULL) {
  if (is.null(cache)) cache <- .multifer_default_svd_cache
  lk <- svd_cache_lookups(cache)
  if (lk == 0L) return(NA_real_)
  svd_cache_hits(cache) / lk
}

#' Compute or retrieve a cached thin SVD.
#'
#' On a hit, the cached decomposition is returned and the hit counter
#' incremented. On a miss, `base::svd(x)` is computed, stored, and
#' returned. The returned list is shaped like `base::svd()` output:
#' `$u`, `$d`, `$v`.
#'
#' @param x A numeric matrix.
#' @param cache Cache environment. `NULL` uses the package default.
#' @param key Optional pre-computed fingerprint. When supplied, skips
#'   the call to `digest::digest()`. Advanced use only; callers must
#'   guarantee uniqueness.
#'
#' @return A list with elements `u`, `d`, `v`.
#' @export
cached_svd <- function(x, cache = NULL, key = NULL) {
  if (is.null(cache)) cache <- .multifer_default_svd_cache
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`cached_svd` requires a numeric matrix.", call. = FALSE)
  }

  if (is.null(key)) {
    key <- svd_cache_fingerprint(x)
  }

  cache$lookups <- (cache$lookups %||% 0L) + 1L

  hit <- cache$store[[key]]
  if (!is.null(hit)) {
    cache$hits <- (cache$hits %||% 0L) + 1L
    return(hit)
  }

  cache$misses <- (cache$misses %||% 0L) + 1L
  sv <- base::svd(x)
  # Store the thin decomposition; drop any large byproducts.
  slim <- list(u = sv$u, d = sv$d, v = sv$v)
  cache$store[[key]] <- slim
  slim
}

#' Fingerprint a matrix for use as a cache key.
#'
#' Hashes `dim` plus the full numeric content via xxhash64 -- fast enough
#' to call per lookup (microseconds for <1M doubles) and collision-free for
#' within-session re-use. Full hashing avoids the false-hit problem that
#' sparse sampling introduces under interior-cell perturbations.
#'
#' @param x A numeric matrix.
#' @return A character scalar.
#' @export
svd_cache_fingerprint <- function(x) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`svd_cache_fingerprint` requires a numeric matrix.", call. = FALSE)
  }
  digest::digest(list(dim = dim(x), x = x), algo = "xxhash64")
}

# Local %||% fallback (not re-exported; used only inside this file).
`%||%` <- function(a, b) if (is.null(a)) b else a
