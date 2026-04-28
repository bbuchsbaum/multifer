# multifer 1.0.0

## v1 support matrix

`multifer` 1.0.0 freezes the v1 public support matrix around the shipped
typed inference families:

- **Mature one-block variance inference**: PCA-family adapters, including
  `adapter_svd()`, `adapter_prcomp()`, `adapter_multivarious_pca()`, and
  `infer_pca()`.
- **Mature two-block latent-root inference**: covariance-mode cross-SVD and
  PLSC-family adapters, including `adapter_cross_svd()`,
  `adapter_multivarious_plsc()`, and `infer_plsc()`.
- **Mature CCA path with explicit support boundaries**: correlation-mode
  cross-SVD, `adapter_cancor()`, `adapter_multivarious_cca()`, and
  `infer_cca()` for paired-row and shipped nuisance-adjusted designs, with
  conservative fallback behavior outside the supported multi-root boundary.
- **Narrow generalized-eigen public surface**: `infer_lda()` and
  `adapter_lda_refit()` expose discriminant-root significance for the shipped
  LDA path.
- **Narrow predictive-gain public surface**: `infer_plsr()` and the PLSR
  adapters expose held-out predictive-gain inference for the shipped PLSR path.

Multiblock models, broader generalized-eigen families, and broader predictive
cross-family engines remain outside the v1 public support matrix.

## Public surfaces

The v1 user-facing API is intentionally narrow. Routine analyses should enter
through the method wrappers (`infer_pca()`, `infer_plsc()`, `infer_cca()`,
`infer_lda()`, and `infer_plsr()`) or through `infer_adapter()` for advanced
adapter-backed workflows. Stability summaries remain separate from component
significance: variable, score, subspace, and feature-importance summaries are
bootstrap reliability outputs, while component tests are permutation or Monte
Carlo significance outputs.

## Executable registration gates

Adapter registration and dispatch are guarded by executable checks rather than
documentation-only promises. The v1 gates cover adapter contracts, typed shape
metadata, capability declarations, maturity labels, exact-path regressions,
engine provenance labels, and family-specific support boundaries. These gates
are exercised by the package test suite, including the registration, capability
matrix, maturity, engine-provenance, oneblock, cross, geneig, predictive, CCA,
and v1-boundary tests.

## Calibration evidence matrix

The v1 calibration evidence matrix combines exact-path regression tests,
synthetic correctness checks, CCA support-matrix tests, Monte Carlo p-value
tests, bootstrap stability tests, feature-importance tests, and the retained
calibration notes under `notes/`. The evidence distinguishes mature
paper-backed paths from narrow shipped public surfaces and keeps unsupported
families out of the v1 claims.
