# Infer Result Feature-Evidence Schema Decision

Decision: keep `infer_result` frozen and keep feature-level evidence as
sidecar tables. Do not add a `variable_significance` or
`feature_evidence` sub-block to `infer_result` in the current schema.

## Rationale

`infer_result` is the unit-centered output contract for `infer()`.
Its built-in inferential claim is component or unit significance:
`component_tests` answers whether a latent unit is selected under the
adapter's null.

Feature-level evidence answers a different question. Depending on the
method and statistic, it may be:

- a descriptive bootstrap ratio for signed loadings;
- a null-calibrated p-value for unsigned squared loadings;
- an adapter-owned VIP, coefficient, block contribution, or
  permutation-importance statistic;
- an aggregate importance score over a selected unit set.

Those rows need their own labels for statistic, orientation,
calibration, interval method, null, adjustment, and validity. Forcing
them into `infer_result` now would either blur the p-value/stability
boundary or freeze a feature schema before the adapter-owned statistic
surface has settled.

## Current Contract

The frozen `infer_result` schema remains:

1. `units`
2. `component_tests`
3. `variable_stability`
4. `score_stability`
5. `subspace_stability`
6. `assumptions`
7. `mc`
8. `cost`
9. `provenance`

`variable_stability` remains a stability table, not a p-value table.
`variable_significance` remains blocked as an infer target. The
general feature-evidence surface is:

- `infer_feature_evidence()`
- `feature_evidence_from_bootstrap()`
- `feature_evidence_pvalues()`
- `feature_evidence_from_adapter()`

The older feature-importance helpers remain compatibility/convenience
sidecars for aggregate squared-loading style summaries.

## Consequences

- `infer()` can continue returning a small, stable object for every
  method family.
- Feature evidence can evolve without forcing an `infer_result` schema
  migration.
- Component p-values, bootstrap stability, and feature-level evidence
  stay visibly separate.
- Downstream adapters can expose method-specific feature statistics
  through sidecar hooks without claiming a frozen result block.

## Revisit Criteria

A future schema review may add a `feature_evidence` block to
`infer_result` only when all of the following are true:

1. `feature_stat_spec` / adapter-owned statistic declarations have
   stabilized.
2. Feature evidence has a settled target name, likely
   `feature_evidence`, not `variable_significance`.
3. Public examples show component p-values, stability summaries, and
   feature evidence without semantic confusion.
4. Drift tests pin the new block's schema and every method that can
   request it.

Until then, feature evidence is deliberately sidecar-only.
