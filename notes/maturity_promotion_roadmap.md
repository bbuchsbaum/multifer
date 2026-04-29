# Maturity Promotion Roadmap

This note records the outcome of the narrow-adapter promotion review.
It is intentionally a roadmap, not a label-change checklist. In
`multifer`, maturity is a support claim backed by executable evidence,
not a presentation setting in `R/maturity.R`.

## Current Decisions

| Adapter | Family | Current label | Decision |
| --- | --- | --- | --- |
| `multivarious_cca` | `(cross, correlation)` | `mature` | Promoted as an explicit backend for the already mature CCA support matrix. |
| `lda_refit` | `(geneig, generalized_eigen)` | `narrow` | Keep narrow. Evidence supports the shipped LDA root-significance wrapper, not mature LDA broadly. |
| `plsr_refit` | `(cross, predictive)` | `narrow` | Keep narrow. Evidence supports the shipped predictive-gain wrapper, not mature predictive-cross inference. |

## Promotion Rule

An adapter is mature only when all of these are true:

1. The public surface is method-specific and bounded.
2. The adapter has executable validity checks for its support region.
3. Component-level evidence is calibrated or parity-checked on the
   support matrix.
4. Unsupported regimes fail early or are documented as out of scope.
5. Public docs, README tables, and drift tests agree with the registry
   label.
6. Stability and feature-evidence outputs are not confused with
   component p-values.

A future promotion should therefore begin with evidence and end with
the label flip. It should not begin by editing `R/maturity.R`.

## Track Status

### `multivarious_cca`

Status: **complete**.

`multivarious_cca` is now a mature explicit CCA backend because it has
scaffold-level parity with the mature CCA support matrix. The supporting
record is in [`cca_multivarious_parity.md`](./cca_multivarious_parity.md).

### `lda_refit`

Status: **not ready for promotion**.

The current public claim is a useful full-rank classical LDA
discriminant-root significance wrapper. The support boundary and
evidence are recorded in:

- [`lda_maturity_contract.md`](./lda_maturity_contract.md)
- [`lda_evidence.md`](./lda_evidence.md)

Remaining promotion work:

1. Run a larger calibration grid over class balance, sample size,
   feature dimension, and number of classes.
2. Add power studies over planted discriminant effects, including weak
   and near-tied roots.
3. Produce parity evidence between `MASS::lda()` and the
   generalized-eigen operator path on supported full-rank data.
4. Decide whether high-dimensional or regularized LDA is a separate
   adapter.
5. Add a metric-aware LDA feature-evidence statistic.
6. Decide whether geneig bootstrap stability belongs in the public LDA
   wrapper.

Until those items land, `lda_refit` should remain `narrow`.

### `plsr_refit`

Status: **not ready for promotion**.

The current public claim is a useful PLSR component-test wrapper for
successive cross-fitted held-out predictive gain. The support boundary
and evidence are recorded in:

- [`plsr_maturity_contract.md`](./plsr_maturity_contract.md)
- [`plsr_evidence.md`](./plsr_evidence.md)

Remaining promotion work:

1. Run a larger calibration grid over sample size, response dimension,
   noise level, weak signal, and high-dimensional predictor regimes.
2. Add power studies that separate predictive gain from high in-sample
   covariance.
3. Document the small-sample, user-fold, and grouped-fold policy.
4. Decide whether `ncomp`, method choice, and regularization are tuned
   outside the wrapper or inside a nested validation contract.
5. Add PLSR-specific feature-evidence statistics such as VIP,
   coefficients, or permutation importance.
6. Keep component p-values separate from bootstrap stability and
   feature-importance sidecars in public examples.

Until those items land, `plsr_refit` should remain `narrow`.

## Close-Out

The narrow-adapter review changed the roadmap rather than the labels:

- promote `multivarious_cca` because it shares the mature CCA scaffold;
- keep `lda_refit` narrow because the evidence is bounded to a
  classical LDA root-significance wrapper;
- keep `plsr_refit` narrow because predictive-gain maturity needs a
  broader null, fold, calibration, and feature-evidence story.

This is the stable v1 posture: mature where the support matrix is
paper-backed or parity-backed, narrow where the wrapper is shipped but
the family-level claim is not yet earned.
