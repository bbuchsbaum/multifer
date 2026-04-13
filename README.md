# multifer

`multifer` provides perturbation-based inference for latent multivariate models.
It is the reference implementation of rank-matched residual randomization for
latent-root inference in PCA-like and covariance two-block models, and it also
provides a broader framework for stability analysis and method-specific
inference across latent model families.

The package is built around three ideas:

- sequential deflation rather than one-shot testing,
- explicit null actions matched to the inferential target,
- strict validity contracts rather than silent guessing.

## Current scope

`multifer` is farther along as a core inference engine than as a polished
end-user package. The current implementation is strongest for:

- one-block variance models, such as PCA-like decompositions,
- two-block covariance models, such as PLSC / cross-covariance SVD,
- bootstrap-based variable, score, and subspace stability summaries.

Current limitations are deliberate:

- cross correlation / CCA supports multi-root testing for the plain paired-row
  design and for common-`Z` nuisance-adjusted designs; more structured designs
  remain conservative,
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

The current API is explicit rather than wrapper-heavy. You specify the adapter,
geometry, and relation directly.

### One-block PCA-like inference

```r
library(multifer)

dat <- bench_oneblock_shadowing(n = 150, p = 40, seed = 1)

res <- infer(
  adapter = "prcomp_oneblock",
  data = dat$X,
  geometry = "oneblock",
  relation = "variance",
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

res <- infer(
  adapter = "cross_svd",
  data = list(X = dat$X, Y = dat$Y),
  geometry = "cross",
  relation = "covariance",
  B = 99,
  R = 49,
  alpha = 0.05,
  seed = 1
)

res$component_tests
```

For correlation-mode cross problems, the strongest current paths are the plain
paired-row design and the common-`Z` nuisance-adjusted design. More structured
correlation-mode designs are still conservative.

## Built-in adapters

Exported adapters currently include:

- `prcomp_oneblock`
- `svd_oneblock`
- `cross_svd`
- `cancor_cross`
- `multivarious_pca` when `multivarious` is installed
- `multivarious_plsc` when `multivarious` is installed

You can inspect the registry with:

```r
list_infer_adapters()
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
ones. A concrete example is current CCA support: the computational machinery is
there, but multi-root permutation validity is not claimed yet, so the current
implementation stops at the first root.

## Benchmarks and synthetic test generators

The package includes synthetic generators and benchmarking helpers for
correctness and performance work:

- `bench_oneblock_null()`
- `bench_oneblock_shadowing()`
- `bench_cross_null()`
- `bench_speed_agreement()`

These are useful both for development and for method-comparison work.

## Project status

Current Phase 1 core:

- oneblock significance and stability
- cross covariance significance and stability
- correlation-mode multi-root inference for paired rows and common-`Z`
  nuisance-adjusted designs
- sequential Monte Carlo ladder infrastructure
- partial-SVD and bootstrap/stability performance improvements

Deferred or later-phase work:

- structured / block-aware valid multi-root CCA
- exact cross-core fast path for null draws
- `multiblock` and `geneig` engines
- variable significance
- method-named wrappers such as `infer_pca()` and `infer_plsc()`
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
