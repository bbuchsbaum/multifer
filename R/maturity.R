# Canonical maturity vocabulary for multifer v1.
#
# Every public surface in the package -- README, adapter registry,
# wrapper roxygen, vignettes -- uses exactly one of these three
# labels when describing how supported a given inferential family
# is:
#
#   - "mature"  : paper-backed, exact, defended by executable
#                 calibration and parity evidence.  The package's
#                 first-class offering.
#   - "narrow"  : shipped today with a deliberately limited public
#                 surface.  The engine works end-to-end but the
#                 family surface is intentionally scoped (e.g.
#                 geneig via LDA only, predictive via PLSR only).
#   - "planned" : architecturally allocated but not shipped in v1.
#                 Names a slot the package reserves for a future
#                 engine or family.
#
# This vocabulary lives in one place -- this file -- so README, the
# adapter registry, and wrapper @details cannot drift away from each
# other.  The boundary drift-guard test
# (tests/testthat/test-v1-boundary.R) pins that invariant.

#' Canonical maturity vocabulary
#'
#' Returns the three maturity labels `multifer` uses throughout its
#' documentation, wrappers, and registry output. See
#' `notes/package_vision.md` for the support taxonomy these labels
#' express.
#'
#' @return Character vector of length three.
#' @export
multifer_maturity_levels <- function() {
  c("mature", "narrow", "planned")
}

# Package-internal map from adapter id to canonical maturity label.
# This is authoritative: if an adapter name is present here, it
# ships in v1; if it is not here, its maturity is "unknown" and the
# drift-guard test will catch it.
.v1_adapter_maturity <- c(
  svd_oneblock       = "mature",
  prcomp_oneblock    = "mature",
  cross_svd          = "mature",
  cancor_cross       = "mature",
  multivarious_pca   = "mature",
  multivarious_plsc  = "mature",
  multivarious_cca   = "mature",
  lda_refit          = "narrow",
  plsr_refit         = "narrow"
)

#' Look up the canonical maturity label for an adapter
#'
#' @param adapter_id Optional character scalar. If supplied, returns
#'   the maturity label for that adapter (or `"unknown"` if the
#'   adapter is not known to `multifer`'s v1 support table). If
#'   omitted, returns the full named vector of known maturity labels.
#'
#' @return Character. Either a single label (when `adapter_id` is
#'   supplied) or a named character vector of every known adapter's
#'   maturity label.
#'
#' @seealso [multifer_maturity_levels()], [list_infer_adapters()]
#' @export
multifer_adapter_maturity <- function(adapter_id = NULL) {
  if (is.null(adapter_id)) {
    return(.v1_adapter_maturity)
  }
  if (!is.character(adapter_id) || length(adapter_id) != 1L ||
      is.na(adapter_id)) {
    stop("`adapter_id` must be a single non-NA character string.",
         call. = FALSE)
  }
  label <- unname(.v1_adapter_maturity[adapter_id])
  if (is.na(label)) "unknown" else label
}
