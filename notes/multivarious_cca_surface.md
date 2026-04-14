# Upstream audit: `multivarious` CCA surface needed by `multifer`

This note records the result of bead `multifer-2aa`.

The question was not "can `multifer` run CCA today?" It can, through the
dedicated `cancor_cross` adapter and the correlation path in
`adapter_cross_svd()`. The question was whether `multivarious` already exposes
the kind of first-class CCA fit surface that `multifer` can wrap in the same
way it already wraps `multivarious::pca()` and `multivarious::plsc()`.

The answer today is: **not yet**.

## What exists in `multivarious` today

The local `multivarious` codebase already has a general two-block container:

- `cross_projector()` in
  `/Users/bbuchsbaum/code/multivarious/R/twoway_projector.R`
- `coef.cross_projector()`
- `reprocess.cross_projector()`
- `shape.cross_projector()`
- `project.cross_projector()`
- `transfer.cross_projector()`

There is also a vignette-level CCA wrapper in
`/Users/bbuchsbaum/code/multivarious/vignettes/Extending.Rmd`:

- `fit_cca()` preprocesses `X` and `Y`
- runs `stats::cancor()`
- stores the result in a `cross_projector`
- attaches canonical correlations ad hoc via `attr(cp, "can_cor")`

So the raw ingredients are there. What is missing is a **stable exported fit
surface**.

## Why the current surface is not enough for a real `multifer` adapter

`multifer` needs a fit object whose contract is part of the package API, not a
vignette convention.

Compared with `multivarious::pca()` and `multivarious::plsc()`, the current CCA
story is missing several things:

1. No exported `cca()` fitter

- There is no documented constructor analogous to `pca()` or `plsc()`.
- The only CCA wrapper is vignette code.

2. No stable CCA class contract

- The vignette uses `classes = "cca_cross_projector"`, but that class is not an
  exported, documented fit class.
- Canonical correlations are attached ad hoc as an attribute rather than stored
  in a documented slot/field with package-level guarantees.

3. No generic `scores.cross_projector()`

- `multifer` needs stable access to the fitted training scores for both blocks.
- `plsc` solves this with a dedicated class method (`scores.plsc()`), but the
  generic `cross_projector` surface does not.
- The vignette stores `sx` and `sy`, but that is a convention, not a method.

4. No `truncate.cross_projector()`

- `multifer` adapters need a documented truncation path.
- `pca()` and `plsc()` fits support truncation at the model level.
- A general CCA `cross_projector` currently does not expose that as a stable
  method.

5. No documented refit contract

- `multifer` can always write its own `refit` hook, but to call the result
  "multivarious-backed" we need the fitter and its preprocessing behavior to be
  part of `multivarious`, not reconstructed in the adapter.

## Minimum downstream contract `multifer` needs

Before `multifer` should add `adapter_multivarious_cca()`, `multivarious`
should expose a first-class CCA fitter with this minimum contract:

### Recommended public entry point

- `cca(X, Y, ncomp = NULL, preproc_x = ..., preproc_y = ..., ...)`

The return value should be a documented subclass of `cross_projector`, e.g.
`"cca"` or `"cca_cross_projector"`.

### Required stable fields / semantics

- canonical roots stored in a documented field such as `cor` or `can_cor`
- X and Y coefficients/loadings accessible via `coef(x, source = "X"|"Y")`
- fitted preprocessors stored as `preproc_x` and `preproc_y`
- fitted X and Y training scores stored in documented fields or exposed through
  a documented method

### Required methods

- `scores.cca(x, source = c("X", "Y"))`
- `coef.cca(x, source = c("X", "Y"))`
- `truncate.cca(x, ncomp)`
- `reprocess.cca(x, new_data, source = c("X", "Y"))`
- `shape.cca(x, source = c("X", "Y"))`

If `scores()` is not added at the subclass level, `multifer` would have to
reach into raw slots (`sx`, `sy`) and treat them as API, which is exactly what
this follow-up is trying to avoid.

## Recommendation

The right upstream move is:

1. Add an exported `cca()` fitter to `multivarious`.
2. Return a documented `cross_projector` subclass with stable roots, scores,
   coefficients, preprocessors, and truncation behavior.
3. Keep the underlying representation compatible with `cross_projector`, but do
   not ask downstream packages to reconstruct the contract from vignette code.

The weaker alternative would be to export an `as_cross_projector.cancor()` or a
documented `cross_projector` constructor pattern for CCA. That is better than
the current situation, but still weaker than an actual fitter because `multifer`
would still have to own refitting and preprocessing conventions itself.

## Decision for `multifer`

Until `multivarious` ships the surface above:

- `infer_cca()` should continue to default to `cancor_cross`
- `cross_svd` should remain the dual-relation reference adapter
- `multifer` should **not** add `adapter_multivarious_cca()`

This keeps the public CCA path honest and avoids pretending that a
`multivarious` CCA backend already exists when it does not.
