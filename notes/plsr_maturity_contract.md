# PLSR Maturity Contract

This note defines what it would mean for `plsr_refit` to become mature,
and why the current v1 label remains **narrow**.

Executable evidence for the current narrow predictive boundary is
recorded in [`plsr_evidence.md`](./plsr_evidence.md).

## Current Public Claim

`infer_plsr()` ships one public predictive inference claim:

- geometry: `cross`
- relation: `predictive`
- adapter: `plsr_refit`
- target: `component_significance`
- estimand: successive cross-fitted incremental held-out predictive gains
- null: row permutation of `Y` relative to fixed `X`
- fitting backend: `pls::plsr()`

The result answers: does the next PLSR latent predictor improve
out-of-sample prediction of `Y` beyond the previously selected latent
predictors? It does not test cross-covariance roots, in-sample fitted
variance, individual variables, or response coefficients.

## Estimand

For rung `a`, the statistic is the cross-fitted increment:

```text
gain_a = R2_oos(a) - R2_oos(a - 1)
```

where `R2_oos(a)` is computed on held-out folds using predictions from
the first `a` PLSR components. Negative increments are clipped to zero
at the statistic layer. The ladder tests the gain sequence directly;
`roots_observed` stores those observed gains only because the shared
`infer_result` schema needs an ordered unit-level numeric vector.

This is the defining difference from PLSC. PLSC lives in
`(cross, covariance)` and tests latent cross-covariance roots. PLSR
lives in `(cross, predictive)` and tests out-of-sample prediction.

## Null And Fold Policy

The v1 null action permutes rows of `Y` while leaving rows of `X`
fixed. This destroys the supervised relation while preserving the
marginal geometry of both blocks.

The predictive engine computes observed and null statistics with the
same fold assignment. When the caller does not supply `folds`, the
engine constructs deterministic balanced folds from the supplied seed.
User-supplied fold labels are accepted by `run_predictive_ladder()` and
are interpreted only as fold membership; relabeling the folds must not
change the result.

## Executable Validity Contract

The current adapter checks:

- `X` and `Y` are numeric matrices;
- `X` and `Y` have matching row counts;
- both blocks contain finite values;
- at least two paired rows are present.

The engine adds stricter predictive execution checks:

- at least four paired rows are required for the predictive ladder;
- each predictive adapter claiming component significance must provide
  `refit`, `predict_response`, `residualize`, `null_action`, and a
  `component_stat()` hook with a `split` argument;
- `plsr_refit$component_stat()` refuses to compute without train/test
  indices.

This split requirement is the main safeguard against in-sample
pseudo-significance.

## What Is In Scope

The narrow shipped contract includes:

- PLSR component tests for held-out incremental predictive gain;
- row-permutation nulls of `Y` relative to fixed `X`;
- deterministic fold construction and user-supplied fold support in the
  lower-level predictive ladder;
- finite-sample Monte Carlo p-values for the stop-at-first-non-rejection
  ladder;
- bootstrap variable, score, and subspace stability summaries as
  descriptive sidecar outputs.

## What Is Out Of Scope

The current PLSR wrapper does not claim:

- covariance-root PLSC inference;
- reduced-rank regression, canonical ridge regression, kernel PLS, or
  other predictive-cross families as public wrappers;
- nested cross-validation for tuning PLSR hyperparameters;
- group-, block-, or subject-restricted predictive nulls through
  `infer_plsr()`;
- coefficient, VIP, or permutation-importance p-values as part of
  `infer_result`;
- variable-level significance from bootstrap stability tables.

`feature_importance_from_bootstrap()` and the feature-evidence sidecars
can summarize or calibrate feature-level quantities, but those are
separate from the component significance claim. A mature PLSR feature
surface should expose method-specific statistics such as VIP,
regression coefficients, or permutation importance through adapter
feature-stat hooks rather than treating generic squared loadings as the
only meaningful predictive importance estimand.

## Promotion Criteria

Promoting `plsr_refit` from `narrow` to `mature` requires all of the
following:

1. A locked calibration report showing that the held-out gain ladder has
   acceptable type-I behavior under the response-permutation null across
   sample sizes, response dimensions, noise levels, and weak-signal
   regimes.
2. Power evidence on planted predictive-rank fixtures, including cases
   where high cross-covariance does not imply held-out gain.
3. A documented fold policy for small samples, user-supplied folds, and
   repeated-measures or grouped predictive designs.
4. A decision on whether tuning of `ncomp`, algorithm choice, and any
   regularization happens outside `infer_plsr()` or inside a nested
   validation contract.
5. A feature-evidence contract for PLSR-specific statistics such as VIP,
   coefficients, and permutation importance.
6. Public docs and tests that keep component p-values separate from
   bootstrap stability and feature-importance sidecars.

Until those criteria are met, `plsr_refit` should stay `narrow`:
useful, shipped, and explicit about the predictive surface it owns.
