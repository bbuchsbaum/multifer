# multifer

`multifer` is a typed perturbation **inference layer** for latent-root
multivariate models. It is not a fitting framework: you fit your model
with [multivarious](https://github.com/bbuchsbaum/multivarious) (or
any compatible package that provides the required accessors via a
thin adapter) and `multifer` provides the inference layer on top —
the exact collapsed Vitale P3 ladder under sequential deflation,
rank-matched residual randomization, the paired-bootstrap stability
consumers (variable / score / subspace), and the unit-centered
result schema.

It ships **mature support** for PCA-family one-block variance
inference and covariance-mode two-block methods (PLSC-family), both
based on the exact collapsed Vitale P3 ladder and an exact core-space
bootstrap, and **qualified support** for canonical-correlation models.

The package is **extensible by design**. A new inferential family is
added by writing an `infer_adapter()` that implements the accessor
contract over a fitted model: how to read its loadings, its scores,
its latent roots, and how to refit it on a perturbed data matrix.
The ladder driver, the null-action machinery, the bootstrap loop, and
the stability consumers all apply uniformly to whatever fitted-model
class an adapter declares. The goal is to stay maximally general
*within* the ordered-latent-root paradigm without pretending every
multivariate method reduces to the same root test. Methods whose
inferential target is not an association root (for example supervised
stagewise predictive methods) are not in the current vocabulary; they
belong to a different inferential family and would be added through
the same adapter contract when their validity story is worked out.

The base-R reference adapters (`svd_oneblock`, `prcomp_oneblock`,
`cross_svd`, `cancor_cross`) are kept as working examples of the
adapter contract — a concrete place to point to when explaining
"here is what it takes to plug a fitted-model class into `multifer`."
They are not the recommended fitting path: the recommended path is
to fit with `multivarious` and pass the projector via `model = `.

The package is built around three ideas:

- sequential deflation rather than one-shot testing,
- explicit null actions matched to the inferential target,
- strict validity contracts rather than silent guessing.

The package is unified at the level of its inferential scaffold:

- ordered latent objects,
- sequential removal of previously claimed structure,
- design-matched null actions,
- stop-at-first-non-rejection testing,
- latent-unit stability summaries.

Direct latent-root methods such as PCA, PLSC, and CCA fit that scaffold
directly today. Additional ordered-latent-root families can be brought
into the same scaffold through new `infer_adapter()` objects without
modifying the engine — that is what it means for `multifer` to be a
platform rather than a fixed set of methods.

The default execution path is exact. Approximate screening hooks are opt-in only
and are not part of the paper-faithful core.

## Current scope

`multifer` is farther along as a core inference engine than as a polished
end-user package. The current implementation is strongest for:

- one-block variance models, such as PCA-like decompositions,
- two-block covariance models, such as PLSC / cross-covariance SVD,
- bootstrap-based variable, score, and subspace stability summaries.

Current limitations are deliberate:

- cross correlation / CCA supports multi-root testing for the plain paired-row
  design, for common-`Z` nuisance-adjusted designs, and for nuisance-adjusted
  designs with within-block exchangeability,
- bootstrap stability defaults to sign alignment; Procrustes alignment is
  retained only for backward compatibility and is not an endorsed inferential
  target,
- `multiblock` and `geneig` are part of the architecture but not yet shipped as
  complete inference engines,
- variable significance is deferred; current variable-level output is stability,
  not p-values.

If an analysis is ambiguous or not supported for a declared
`(geometry, relation, design)` combination, `multifer` is designed to error
rather than guess.

## Installation

`multifer` is not on CRAN.

```r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/multifer")
```

Optional packages:

- `RSpectra` improves the partial-SVD fast path on larger problems.
- `multivarious` enables the `multivarious_pca` and `multivarious_plsc`
  adapters when installed.

## What `infer()` returns

The main entry point is `infer()`. It returns a frozen `infer_result` object
with:

- `units`: component or subspace units used as the common indexing scheme,
- `component_tests`: sequential significance results,
- `variable_stability`: variable-level stability summaries,
- `score_stability`: observation-level score stability summaries,
- `subspace_stability`: principal-angle stability summaries,
- `assumptions`, `mc`, `cost`, and `provenance` metadata blocks.

This is a deliberately unit-centered design: downstream tables refer to
`unit_id`, not to raw component positions.

## Minimal examples

The package ships thin convenience wrappers for the mature and qualified
families:

- `infer_pca()` for one-block variance models
- `infer_plsc()` for two-block covariance models
- `infer_cca()` for two-block correlation models

The lower-level `infer()` interface remains available when you want explicit
control over adapter, geometry, relation, and design.

### One-block PCA-like inference

```r
library(multifer)

dat <- bench_oneblock_shadowing(n = 150, p = 40, seed = 1)

res <- infer_pca(
  dat$X,
  B = 99,
  R = 49,
  alpha = 0.05,
  seed = 1
)

res
res$component_tests
res$variable_stability
```

### Two-block covariance inference

```r
library(multifer)

dat <- bench_cross_null(n = 150, p_x = 30, p_y = 20, seed = 1)

res <- infer_plsc(
  dat$X, dat$Y,
  B = 99,
  R = 49,
  alpha = 0.05,
  seed = 1
)

res$component_tests
```

### Two-block correlation inference

```r
library(multifer)

dat <- bench_cross_null(n = 150, p_x = 30, p_y = 20, seed = 1)

res <- infer_cca(
  dat$X, dat$Y,
  B = 99,
  R = 49,
  alpha = 0.05,
  seed = 1
)

res$component_tests
```

For correlation-mode cross problems, the strongest current paths are the plain
paired-row design, common-`Z` nuisance-adjusted designs, and nuisance-adjusted
designs with within-block exchangeability. Richer structured designs still
require more work on the exchangeability side.

## Built-in adapters and maturity

| Adapter                | Geometry  | Relation    | Maturity   |
|------------------------|-----------|-------------|------------|
| `svd_oneblock`         | oneblock  | variance    | **mature** |
| `prcomp_oneblock`      | oneblock  | variance    | **mature** |
| `cross_svd` (cov)      | cross     | covariance  | **mature** |
| `multivarious_pca`     | oneblock  | variance    | **mature** |
| `multivarious_plsc`    | cross     | covariance  | **mature** |
| `cross_svd` (cor)      | cross     | correlation | qualified  |
| `cancor_cross`         | cross     | correlation | qualified  |

*Mature* adapters are the paper-backed exact path. *Qualified* adapters
currently cover the paired-row design, common-`Z` nuisance-adjusted
designs, and nuisance-adjusted designs with within-block exchangeability;
richer structured designs on the correlation side are still being
exactified.

All mature-tier adapters now carry **executable validity contracts**:
`infer()` runs the adapter's `checked_assumptions` against the raw data
before the engine starts and errors in strict mode (the default) if any
check fails. Pass-through records land in `result$assumptions$checked`.

You can inspect the registry with:

```r
list_infer_adapters()
```

### Third minimal example: correlation-mode cross (CCA)

```r
library(multifer)

dat <- bench_cross_null(n = 150, p_x = 10, p_y = 8, seed = 1)

res <- infer(
  adapter  = "cancor_cross",
  data     = list(X = dat$X, Y = dat$Y),
  geometry = "cross",
  relation = "correlation",
  B = 99,
  R = 49,
  alpha = 0.05,
  seed = 1
)

res$component_tests
res$subspace_stability
```

## Design principles

### Typed shapes

Every inferential problem is described by a typed shape:

- `geometry`: oneblock, cross, multiblock, or geneig
- `relation`: variance, covariance, correlation, or generalized-eigen
- `design`: the null-action assumptions, such as exchangeable or paired rows

This is what lets the package distinguish, for example, PLS-like covariance
inference from CCA-like correlation inference even when both live in the same
`cross` geometry.

### Strict dispatch

`infer()` defaults to `strict = TRUE`. If an adapter supports multiple
relations for a geometry, the caller must disambiguate. This is meant to
prevent accidental use of the wrong inferential target.

### Validity first

The package favors conservative or partial validity claims over broad but weak
ones. A concrete example is current CCA support: multi-root inference is only
claimed for the supported paired-row and nuisance-adjusted designs, and richer
exchangeability structures are still treated conservatively. Likewise, bootstrap
stability defaults to sign alignment rather than Procrustes rotation, because
the package treats subspace uncertainty as the primary target when roots are
near-tied.

## Benchmarks and synthetic test generators

The package includes synthetic generators and benchmarking helpers for
correctness and performance work:

- `bench_oneblock_null()`
- `bench_oneblock_shadowing()`
- `bench_cross_null()`
- `bench_speed_agreement()`

These are useful both for development and for method-comparison work.

For a release-style guardrail pass on the mature families, run:

```sh
Rscript tools/mature_guardrails.R --run-bench=true
```

## Project status

Current Phase 1 core:

- oneblock significance and stability
- cross covariance significance and stability
- correlation-mode multi-root inference for paired rows, common-`Z`
  nuisance-adjusted designs, and nuisance-adjusted within-block designs
- sequential Monte Carlo ladder infrastructure
- exact covariance core-space null draws
- partial-SVD and bootstrap/stability performance improvements

Deferred or later-phase work:

- richer structured / hierarchical valid multi-root CCA
- `multiblock` and `geneig` engines
- variable significance
- vignettes, package website, and broader user-facing docs

## Paper and package

The intended split is:

- the methods paper focuses on rank-matched residual randomization for latent
  operator models, with PCA and PLSC as the main statistical scope,
- the package implements that core faithfully and also serves as the broader R
  framework around it.

So the package should be read as both:

- the reference implementation of the paper's core method,
- and an evolving perturbation-inference toolkit for latent multivariate
  models.

## License

MIT
