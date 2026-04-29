# `multivarious_cca` parity evidence

This note records the executable evidence for mote
`bd-01KQCP0ZG2J2RZAA3VHSHM1T4D`.

## Contract

`multivarious_cca` is an explicit CCA backend for the `(cross, correlation)`
family. The mature public CCA path remains `infer_cca(..., adapter =
"cancor_cross")`. For `multivarious_cca` to move toward the mature boundary, it
must agree with the mature path on the scaffold-level inferential outputs:

- component test statistics,
- Monte Carlo p-values under the same seed,
- unit formation and selected-unit decisions,
- multi-root stopping behavior on supported designs,
- shared rejection of invalid CCA inputs before inference.

The parity target is the `infer_result` scaffold, not byte-for-byte equality of
adapter-owned fit objects or loading bases.

## Executable coverage

The parity matrix lives in
[`tests/testthat/test-cca-parity.R`](../tests/testthat/test-cca-parity.R).

The reference backend is `cancor_cross`. The tested alternate backends are:

- `cross_svd`, the dual-relation reference adapter in correlation mode,
- `multivarious_cca`, the `multivarious::cca()` adapter.

The test matrix covers:

- `paired_rows()` on all shared CCA stress fixtures,
- `blocked_rows(groups)` on all shared CCA stress fixtures,
- `nuisance_adjusted(Z)`,
- `nuisance_adjusted(Z, groups)`,
- conservative first-root capping for `exchangeable_rows()`,
- explicit aspect-ratio fixtures with wide-X/narrow-Y and narrow-X/wide-Y
  shapes,
- exact rank deficiency, which must fail at the shared
  `cross_blocks_full_column_rank` validity gate for all three backends.

## Tolerances

The test uses deterministic seeds and `R = 0` so the comparison isolates the
component-significance ladder rather than bootstrap stability. For each design
and fixture:

- component statistics must match within `1e-10`,
- p-values must match within `1e-8`,
- units, selected-unit decisions, tested rows, and `stopped_at` values must
  match exactly.

Those tolerances are intentionally tight because all three infer calls compile
to the same correlation ladder. If this test fails, the likely problem is
dispatch, recipe compilation, validity gating, or a changed CCA support
boundary, not ordinary numerical noise.

## Interpretation

The current evidence supports this statement:

> `multivarious_cca` is scaffold-parity compatible with the mature CCA path for
> component-significance inference on the shipped support matrix.

It does not yet imply that `multivarious_cca` should become the default CCA
wrapper backend. Default-backend maturity still needs documentation and
calibration review, especially around adapter-owned loadings, score stability,
regularization defaults, and any future feature-evidence statistics.
