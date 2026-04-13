#' multifer: Typed Perturbation Inference for Projector-Based Multivariate Models
#'
#' A shape-first inference layer for projector-based multivariate models.
#' Provides one coherent system for component significance, variable and
#' score stability, and subspace uncertainty across one-block, cross-block,
#' multiblock, and generalized-eigen shapes.
#'
#' The framework is built on four principles:
#'
#' 1. **Typed shapes** — every inferential object is a triple
#'    `(geometry, relation, design)` rather than a bare "shape".
#' 2. **Strict dispatch** — `infer()` errors on ambiguous adapters
#'    instead of silently guessing between e.g. PLS-like and CCA-like
#'    inference.
#' 3. **Machine-enforced validity** — adapters declare capabilities and
#'    checked assumptions; the framework refuses to produce `"exact"` or
#'    `"conditional"` significance output unless those assumptions pass.
#' 4. **Latent units** — result tables key off `unit_id`, which may be a
#'    singleton component or a subspace of tied roots, so tied-root
#'    handling is honest by default.
#'
#' @keywords internal
"_PACKAGE"
