# CCA Support and Calibration in `multifer`

This note records the current shipped support boundary for the CCA path in
`multifer`, along with the executable evidence that defends that boundary.

It is a package-scope note, not a paper claim. Paper 1 can remain narrower than
this and still be honest.

## Bottom line

`multifer` ships a real CCA path today:

- default public wrapper: `infer_cca()`
- default adapter: `cancor_cross`
- engine family: `(cross, correlation)`
- explicit alternate backend: `multivarious_cca`

The package supports **multi-root** CCA inference for a specific set of designs
where the whitening map, deflation rule, and null action stay aligned. Outside
that set, the package deliberately falls back to **first-root-only** inference
rather than making a stepwise claim it cannot yet defend.

`multivarious_cca` is now parity-pinned for scaffold-level component
inference on this support matrix and is labeled as a mature explicit backend
for that validity-gated surface. It is not the default `infer_cca()` backend.
Its `multivarious::cca()` fit surface brings adapter-owned loading and score
semantics, plus regularization defaults, so those summaries are documented as
backend-specific rather than as raw-equality targets.

## Supported multi-root designs

The current cross engine allows multi-root correlation inference when the null
action preserves within-group exchangeability and the whitened-space deflation
commutes with the allowed row action. In the shipped code, that means:

- `paired_rows()`
- `nuisance_adjusted(Z)`
- `nuisance_adjusted(Z, groups)`
- `blocked_rows(groups)`

The engine gate lives in [R/engine_cross.R](../R/engine_cross.R). The
`allow_multiroot_correlation` branch is the authoritative implementation of
that support matrix.

| Design | Multi-root? | Current stance |
| --- | --- | --- |
| `paired_rows()` | yes | shipped default path |
| `nuisance_adjusted(Z)` | yes | shipped supported path |
| `nuisance_adjusted(Z, groups)` | yes | shipped supported path |
| `blocked_rows(groups)` | yes | shipped supported path |
| `exchangeable_rows()` | no | conservative first-root cap |

## Conservative boundary

For correlation designs outside the support matrix, `run_cross_ladder()`
automatically caps `max_steps` at 1. This is intentional:

- first-root CCA inference remains valid,
- unsupported higher-root stepwise claims are refused,
- the package chooses conservatism over silent overreach.

This boundary is part of the intended user contract, not a temporary accident.
`exchangeable_rows()` is the simplest concrete example of a design that compiles
but is deliberately kept on this conservative first-root-only path.

## Executable validity checks

The CCA adapters do not rely on prose alone. They ship executable checks for:

- full numerical column rank in each block
- sufficient effective rows for whitening (`n > p + q`, or the residual-basis
  analogue after nuisance adjustment)
- full-rank nuisance design matrices
- grouped-design consistency and non-singleton exchangeable blocks

The relevant tests live in:

- [tests/testthat/test-validity-checks.R](../tests/testthat/test-validity-checks.R)
- [tests/testthat/test-validity-hardening.R](../tests/testthat/test-validity-hardening.R)

These tests are the main reason the current CCA path is now better described as
"shipped with an explicit support boundary" than as a vague qualified feature.

## Calibration evidence

The current regression/calibration suite already defends the shipped boundary in
three distinct ways.

### Null calibration

[tests/testthat/test-calibration.R](../tests/testthat/test-calibration.R)
contains an empirical-null study for `(cross, correlation)` showing that the
observed rejection rate stays within a loose binomial tolerance band around the
nominal alpha level on null synthetic data.

[tests/testthat/test-cca-support-matrix.R](../tests/testthat/test-cca-support-matrix.R)
extends that by pinning a grouped nuisance-adjusted null regime, so the shipped
support matrix has an explicit null-calibration check beyond the plain paired
design.

### Exactness and adapter agreement

[tests/testthat/test-synthetic-correctness.R](../tests/testthat/test-synthetic-correctness.R)
checks that:

- the `cross_svd` correlation path and `cancor_cross` recover the same
  canonical correlations,
- correlation-mode statistics obey the expected scale invariances,
- common row permutations preserve the fitted correlation roots.

These tests defend the claim that the shipped CCA adapters are computing the
same latent object, not two loosely related approximations.

[tests/testthat/test-cca-support-matrix.R](../tests/testthat/test-cca-support-matrix.R)
adds a direct support-matrix regression: the four shipped designs must permit
multi-root recovery on the same planted canonical-rank fixture, while an
unsupported design must cap the ladder to a single tested root.

[tests/testthat/test-cca-parity.R](../tests/testthat/test-cca-parity.R) extends
that adapter-agreement evidence. `cancor_cross` is the reference backend, and
both `cross_svd` and `multivarious_cca` must match its component statistics,
Monte Carlo p-values, unit formation, selected-unit decisions, and stopping rows
on paired rows, blocked rows, nuisance-adjusted designs, grouped
nuisance-adjusted designs, exchangeable-row first-root capping, aspect-ratio
stress fixtures, and rank-deficient validity failures. The tolerance contract is
recorded in [notes/cca_multivarious_parity.md](cca_multivarious_parity.md).

### Fast-path parity

[tests/testthat/test-core-update-cross.R](../tests/testthat/test-core-update-cross.R)
checks that the correlation-mode fast path matches refit bit-for-bit on:

- canonical correlations,
- aligned loadings,
- aligned scores,
- downstream stability summaries.

[tests/testthat/test-integration-phase15.R](../tests/testthat/test-integration-phase15.R)
then checks that `infer(..., relation = "correlation")` actually reports
core-update usage in the shipped path.

## What is still not claimed

This note does **not** claim:

- that every conceivable structured CCA design has multi-root validity,
- that variable-level significance exists for CCA loadings,
- that `multivarious_cca` is the public default CCA backend,
- that `multivarious_cca` regularized high-dimensional refits are part of the
  mature `infer()` boundary,
- that Paper 1 needs to treat CCA as part of its theorem-bearing center.

The next step for CCA is not another backend-name promotion. It is a targeted
feature-evidence contract for CCA loadings, scores, and regularized fits that
states which quantities are invariant, which are backend-owned, and which need
method-specific calibration.
