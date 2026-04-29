# PLSR Evidence Ledger

This note records executable evidence for the current `plsr_refit`
support boundary. It is the evidence companion to
[`plsr_maturity_contract.md`](./plsr_maturity_contract.md).

## Status

`plsr_refit` remains **narrow**. The shipped wrapper has useful
predictive calibration evidence, but it is not yet a mature claim for
the broader predictive-cross family.

The current public claim is:

- PLSR through `pls::plsr()`;
- `(cross, predictive, component_significance)`;
- successive cross-fitted held-out predictive gain;
- row permutation of `Y` relative to fixed `X`;
- deterministic or user-supplied folds in the predictive ladder;
- descriptive bootstrap stability sidecars, separate from p-values.

## Evidence Matrix

| Claim | Executable evidence | Interpretation |
| --- | --- | --- |
| The adapter declares only the shipped predictive surface | `tests/testthat/test-adapter-plsr.R`: `adapter_plsr_refit constructs with predictive capabilities` | The registry advertises predictive component significance plus descriptive stability, not covariance-root or feature-significance inference. |
| The backend hooks round-trip through `pls::plsr()` | `tests/testthat/test-adapter-plsr.R`: `adapter_plsr_refit refit and predictive hooks round-trip` | The adapter exposes roots, scores, loadings, residualization, prediction, and split-aware component statistics on supported data. |
| Explicit component caps are honored by the backend fit | `tests/testthat/test-adapter-plsr.R`: `adapter_plsr_refit honors explicit ncomp caps` | Adapter construction can restrict the fitted PLSR model without silently creating additional fitted components. |
| The null permutes responses against fixed predictors | `tests/testthat/test-adapter-plsr.R`: `plsr_refit null_action permutes Y relative to fixed X` | The predictive null destroys supervision while preserving the marginal predictor and response blocks. |
| Component statistics require train/test splits | `tests/testthat/test-adapter-plsr.R`: `plsr_refit component_stat is explicitly split-aware`; registration-gate tests for predictive adapters | In-sample predictive gains cannot satisfy the shipped adapter contract. |
| The lower-level predictive ladder handles folds coherently | `tests/testthat/test-engine-predictive.R`: generated-fold parity, fold-label invariance, and split-sensitivity tests | Fold labels are interpreted as memberships only, while different memberships can change the held-out-gain estimand. |
| The PLSR wrapper recovers planted multivariate predictive structure | `tests/testthat/test-engine-predictive.R`: `plsr_refit recovers planted predictive rank on a supervised fixture` | On a multivariate-response fixture with two planted latent predictors, the wrapper selects the planted predictive rank. |
| The response-permutation null is approximately calibrated | `tests/testthat/test-engine-predictive.R`: predictive-ladder shuffled-`Y` calibration and `plsr_refit rung-1 null calibration is within 2 SE at B = 1000` | Current smoke calibration keeps null first-rung rejection close to nominal in the tested fixtures. |
| Noisy high-dimensional predictors reduce the selected predictive surface | `tests/testthat/test-engine-predictive.R`: `plsr_refit does not invent extra rank in noisy high-dimensional X` | The ladder can retain the leading predictive unit without manufacturing an additional significant unit from many noisy `X` columns. |
| Descriptive stability separates signal units from noise units | `tests/testthat/test-engine-predictive.R`: `plsr_refit stability outputs separate selected signal units from the noise unit` | Bootstrap stability is useful for reading selected latent predictors, but remains distinct from component p-values. |

## What This Evidence Supports

The tests support saying:

> `infer_plsr()` is a shipped, validity-gated PLSR component-test wrapper
> for held-out incremental predictive gain, using a response-permutation
> null and split-aware component statistics.

That is a narrow support claim. It is not a mature claim for all
predictive-cross methods or all PLSR modeling workflows.

## Remaining Gaps Before Maturity

The following remain open before changing `plsr_refit` from `narrow` to
`mature`:

1. A larger calibration grid over sample size, response dimension,
   noise level, weak signal, and high-dimensional predictor regimes.
2. Power studies that separate predictive gain from high in-sample
   covariance.
3. A documented small-sample and grouped-fold policy.
4. A decision on whether `ncomp`, method choice, and regularization are
   tuned outside the wrapper or inside a nested validation contract.
5. A PLSR-specific feature-evidence surface for VIP, coefficients, or
   permutation importance.
6. Public examples that keep component p-values separate from
   bootstrap stability and feature-importance sidecars.

Until those gaps close, the current label should remain `narrow`.
