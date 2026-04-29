# LDA Maturity Contract

This note defines what it would mean for `lda_refit` to become mature,
and why the current v1 label remains **narrow**.

Executable evidence for the current narrow boundary is recorded in
[`lda_evidence.md`](./lda_evidence.md).

## Current Public Claim

`infer_lda()` ships one public inference claim:

- geometry: `geneig`
- relation: `generalized_eigen`
- adapter: `lda_refit`
- target: `component_significance`
- estimand: the ordered sequence of LDA discriminant roots
- null: class-label permutation with rows of `X` held fixed
- deflation: B-metric deflation of the within-class generalized-eigen
  pencil

The result answers: how many discriminant roots rise above the
label-permutation null? It does not answer which variables are
significant, which class contrasts are significant, or whether LDA
classification accuracy is better than chance.

## Executable Validity Contract

The wrapper is valid only when the raw data satisfy the adapter's
checked assumptions:

- `X` is a numeric matrix with finite entries.
- `labels` has one non-missing value per row of `X`.
- At least two classes are present.
- Each class has at least two observations.
- The within-class centered design has full column rank.

The final rank condition is deliberately conservative. The v1 adapter
delegates fitting to `MASS::lda()`, so p > n - K and collinear
within-class designs are rejected before the engine starts. The
lower-level scatter-pencil helper uses a ridge term internally, but that
does not make the public MASS-backed wrapper a high-dimensional or
regularized LDA method.

## What Is In Scope

The narrow shipped contract includes:

- label-permutation component tests for LDA discriminant roots;
- B-metric deflation, enforced by the geneig adapter gate;
- finite-sample Monte Carlo p-values for the sequential ladder;
- executable assumption checks that fail early in strict mode;
- provenance labels that identify the statistic, null, and estimand.

This is enough for a useful LDA root-significance wrapper, but not
enough for the `mature` tier.

## What Is Out Of Scope

The current LDA wrapper does not claim:

- generalized-eigen families beyond LDA;
- p > n or singular-within-class LDA;
- regularized, sparse, or shrinkage LDA;
- bootstrap variable, score, or subspace stability through `infer()`;
- variable-level significance inside `infer_result`;
- metric-aware LDA feature evidence.

`feature_importance_pvalues()` can be used as a sidecar refit example,
but its current LDA route is loading-based. A mature LDA feature surface
should expose a metric-aware statistic, such as B-metric discriminant
contribution or class-separation contribution, through the feature
evidence hook family.

## Promotion Criteria

Promoting `lda_refit` from `narrow` to `mature` requires all of the
following:

1. A locked calibration report showing that the label-permutation root
   ladder has acceptable type-I behavior and power across class balance,
   effect size, and sample-size regimes.
2. Parity evidence between the `MASS::lda()` fit path and the
   generalized-eigen operator path for roots, scores, and scaling on
   supported full-rank data.
3. Clear tests for unsupported rank regimes, including p > n and
   collinear within-class designs.
4. A decision on whether the mature surface remains classical full-rank
   LDA or adds a separate regularized LDA adapter.
5. A feature-evidence contract for LDA that avoids raw Euclidean
   loadings as the default variable statistic.
6. A decision on whether bootstrap stability for geneig belongs in the
   public wrapper or stays outside the LDA maturity claim.

Until those criteria are met, `lda_refit` should stay `narrow`: useful,
shipped, and honest about the exact surface it owns.
