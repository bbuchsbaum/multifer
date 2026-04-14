# Upstream audit: `multivarious` CCA surface needed by `multifer`

This note records the result of bead `multifer-2aa`.

The question was not "can `multifer` run CCA today?" It can, through the
dedicated `cancor_cross` adapter and the correlation path in
`adapter_cross_svd()`. The question was whether `multivarious` already exposes
the kind of first-class CCA fit surface that `multifer` can wrap in the same
way it already wraps `multivarious::pca()` and `multivarious::plsc()`.

The answer now is: **yes, enough of the surface exists to support a first
`multifer` adapter**.

## What exists in `multivarious` today

The local `multivarious` codebase already has a general two-block container:

- `cross_projector()` in
  `/Users/bbuchsbaum/code/multivarious/R/twoway_projector.R`
- `coef.cross_projector()`
- `reprocess.cross_projector()`
- `shape.cross_projector()`
- `project.cross_projector()`
- `transfer.cross_projector()`

`multivarious` now also exposes a first-class CCA surface in
`/Users/bbuchsbaum/code/multivarious/R/cca.R`:

- exported `cca()`
- class `"cca"` layered on top of `cross_projector`
- documented canonical correlations in `cor`
- paired scores in `sx` / `sy` with `scores.cca()`
- explicit subclass methods:
  - `coef.cca()`
  - `reprocess.cca()`
  - `transfer.cca()`
  - `truncate.cca()`

That is enough of a stable exported fit surface for `multifer` to wrap.

## What changed relative to the original concern

The previous blocker list is now largely satisfied:

1. Exported `cca()` fitter: **present**
2. Stable subclass contract: **present enough for wrapping**
3. Score accessors: **present via `scores.cca()`**
4. Safe truncation: **present via `truncate.cca()`**
5. Refit contract: **present via `cca()` itself**

What remains weaker than the PCA / PLSC story is not the existence of the
surface, but how opinionated it is:

- the regularized whitening choices are built into the fitter
- the subclass contract is newer and less battle-tested than `pca()` / `plsc()`
- `multifer` still needs to pin parity against the shipped `cancor_cross`
  adapter before making any default-path claims

## Downstream implication for `multifer`

The upstream surface is now good enough to justify
`adapter_multivarious_cca()`. The remaining question is not whether the adapter
is possible, but whether it should become the public default.

That requires explicit downstream parity tests:

- full-rank `lambda = 0` agreement with `stats::cancor()`
- regularized `p > n` smoke and stability checks
- `infer()`-level agreement with `cancor_cross` on the supported paired and
  nuisance-adjusted CCA designs

## Decision for `multifer`

Current recommendation:

- add `adapter_multivarious_cca()` as an available backend
- keep `infer_cca()` defaulting to `cancor_cross` until parity is pinned
- keep `cross_svd` as the dual-relation reference adapter

So the old blocker is closed. What remains is an adapter-validation task, not
an upstream surface gap.
