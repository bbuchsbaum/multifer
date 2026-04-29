# LDA Evidence Ledger

This note records the executable evidence for the current `lda_refit`
support boundary. It is the evidence companion to
[`lda_maturity_contract.md`](./lda_maturity_contract.md).

## Status

`lda_refit` remains **narrow**. The shipped claim is useful and
well-pinned, but it is not yet the full mature LDA surface described in
the promotion criteria.

The current public claim is:

- full-rank classical LDA via `MASS::lda()`;
- `(geneig, generalized_eigen, component_significance)`;
- label-permutation null with rows of `X` held fixed;
- B-metric generalized-eigen deflation;
- component/root significance only.

## Evidence Matrix

| Claim | Executable evidence | Interpretation |
| --- | --- | --- |
| The adapter declares only the shipped geneig significance surface | `tests/testthat/test-adapter-lda.R`: `adapter_lda_refit constructs and exposes only shipped geneig significance` | The registry cannot accidentally advertise bootstrap stability or variable significance as an LDA capability. |
| The null permutes labels, not rows | `tests/testthat/test-adapter-lda.R`: `lda_refit null_action permutes labels rather than rows` | The null destroys class assignment while preserving the observation geometry. |
| The rung statistic is the generalized-root tail ratio | `tests/testthat/test-adapter-lda.R`: `lda_refit component_stat returns the generalized-root tail ratio` | The adapter statistic matches the geneig ladder definition. |
| p > n and exact within-class collinearity fail before fitting | `tests/testthat/test-adapter-lda.R`: `lda_refit fails early for unsupported within-class rank regimes` | The public MASS-backed wrapper does not silently regularize high-dimensional or singular LDA. |
| Class imbalance is allowed, but too-small classes fail | `tests/testthat/test-adapter-lda.R`: `lda_refit distinguishes class imbalance from too-small classes` | Unequal class sizes are valid when every class has enough observations and full within-class rank; singleton classes fail the validity contract. |
| LDA rank is bounded by `K - 1` | `tests/testthat/test-adapter-lda.R`: `lda_refit + geneig engine recovers rank bounded by K - 1`; `tests/testthat/test-engine-geneig.R`: `run_geneig_ladder recovers a planted rank-2 four-class LDA pencil` | The ladder does not infer an open-ended component sequence for LDA. |
| Label-permutation null is approximately calibrated | `tests/testthat/test-adapter-lda.R`: `lda_refit label-permutation null is approximately calibrated`; `tests/testthat/test-engine-geneig.R`: `label-permutation null is approximately calibrated at rung 1` | Current smoke calibration keeps null p-values in a plausible range; this is not yet a full mature-tier calibration grid. |
| B-metric deflation is required | `tests/testthat/test-engine-geneig.R`: `run_geneig_ladder uses B-metric deflation rather than Euclidean deflation` | Euclidean deflation is shown to produce the wrong remaining generalized roots. |
| The geneig engine behaves sensibly under transformations and near-singular SPD metrics | `tests/testthat/test-engine-geneig.R`: congruence, scaling, and near-singular SPD metric tests | The engine-level algebra is not tied to a fragile coordinate representation. |

## What This Evidence Supports

The tests support saying:

> `infer_lda()` is a shipped, validity-gated LDA root-significance
> wrapper for full-rank classical LDA, using a label-permutation null and
> B-metric generalized-eigen deflation.

That is a narrow support claim. It is not a mature claim for all LDA
workflows.

## Remaining Gaps Before Maturity

The following remain open before changing `lda_refit` from `narrow` to
`mature`:

1. A larger calibration grid over null class balance, feature dimension,
   sample size, and number of classes.
2. Power studies over planted discriminant effects, including weak and
   near-tied roots.
3. A parity report between `MASS::lda()` scaling/scores and the
   generalized-eigen operator path on supported full-rank data.
4. A decision on whether high-dimensional or regularized LDA should be a
   separate adapter rather than an extension of `lda_refit`.
5. A metric-aware feature-evidence statistic for LDA variables.
6. A decision on whether geneig bootstrap stability belongs in the
   public LDA wrapper.

Until those gaps close, the current label should remain `narrow`.
