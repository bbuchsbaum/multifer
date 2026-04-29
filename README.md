# multifer

**Inference for fitted multivariate latent-variable models.**

`multifer` adds component tests and bootstrap stability summaries to
models you have already fit. The usual path is to fit with
[multivarious](https://github.com/bbuchsbaum/multivarious) and pass the
projector to `multifer`; other fitted-model classes can be used through a
thin adapter.

The package focuses on ordered latent structures such as PCA, PLSC, CCA,
LDA, and PLSR. It provides the inference layer: sequential deflation,
design-matched null actions, rank-matched residual randomization, and
variable / score / subspace stability summaries in a shared result
schema.

Current support has three public maturity levels, plus an adapter-ready
extension surface:

- **mature**: PCA-family one-block variance inference, PLSC-family
  two-block covariance inference, and the shipped CCA path on supported
  paired and nuisance-adjusted designs,
- **narrow**: `infer_lda()` for generalized-eigen discriminant models and
  `infer_plsr()` for predictive-gain models,
- **adapter-ready**: multiblock inference through the generic
  adapter-driven ladder; no bundled multiblock method wrapper is shipped yet,
- **planned**: bundled multiblock method families, broader generalized-eigen
  families beyond LDA, and broader predictive-cross families beyond PLSR.

At its core, `multifer` is built around three ideas:

- sequential deflation rather than one-shot testing,
- explicit null actions matched to the inferential target,
- strict validity contracts rather than silent guessing.

The shared scaffold is:

- ordered latent objects,
- sequential removal of previously claimed structure,
- design-matched null actions,
- stop-at-first-non-rejection testing,
- latent-unit stability summaries.

Direct latent-root methods such as PCA, PLSC, and CCA fit that scaffold
directly today. Additional ordered-latent-root families can be brought
into the same scaffold through new `infer_adapter()` objects: an adapter
declares how to read loadings, scores, latent roots, and how to refit a
compatible model on perturbed data.

The default execution path is exact. Approximate screening hooks are opt-in only
and are not part of the paper-faithful core.

## Current scope

`multifer` is farther along as a core inference engine than as a polished
end-user package. The current implementation is strongest for:

- one-block variance models, such as PCA-like decompositions,
- two-block covariance models, such as PLSC / cross-covariance SVD,
- two-block correlation models through the shipped CCA path,
- generalized-eigen discriminant models via `infer_lda()`,
- predictive-gain models via `infer_plsr()`,
- bootstrap-based variable, score, and subspace stability summaries.

Current limitations are deliberate:

- cross correlation / CCA supports multi-root testing for the plain paired-row
  design, for common-`Z` nuisance-adjusted designs, for nuisance-adjusted
  designs with within-block exchangeability, and for `blocked_rows(groups)`;
  richer structured correlation designs still need more exchangeability work,
- the current `geneig` public surface is intentionally narrow: LDA is shipped,
  while broader metric-weighted / contrastive generalized-eigen models remain
  planned; `infer_lda()` defaults to `targets = "component_significance"`
  because broader geneig bootstrap / stability workflows are not yet part of
  the wrapper-level public surface,
- the predictive-gain public surface is intentionally narrow: PLSR is shipped,
  while broader predictive-cross models remain planned,
- bootstrap stability defaults to sign alignment; Procrustes alignment is
  retained only for backward compatibility and is not an endorsed inferential
  target,
- `multiblock` has an adapter-driven inference engine, bootstrap resampling,
  and result assembly, but no bundled public method wrapper such as
  `infer_multiblock()` yet,
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
- `multivarious` enables the multivarious-backed adapters. When it is
  installed, `infer_pca()` defaults to `multivarious_pca` and
  `infer_plsc()` defaults to `multivarious_plsc`; the
  `multivarious_cca` adapter is also registered as an alternate CCA
  backend, though `infer_cca()` still defaults to `cancor_cross`.
  Without `multivarious`, the PCA / PLSC wrappers fall back to
  `prcomp_oneblock` and `cross_svd`.
- `MASS` enables the `lda_refit` adapter and the `infer_lda()` wrapper.
- `pls` enables the `plsr_refit` adapter and the `infer_plsr()` wrapper.

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

The package ships thin convenience wrappers for the shipped families:

- `infer_pca()` for one-block variance models
- `infer_plsc()` for two-block covariance models
- `infer_cca()` for two-block correlation models
- `infer_lda()` for generalized-eigen discriminant models
- `infer_plsr()` for predictive-gain models

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

`infer_cca()` defaults to the dedicated `cancor_cross` adapter, which keeps
the public CCA path single-relation and unambiguous. For correlation-mode
cross problems, the strongest current paths are the plain paired-row design,
common-`Z` nuisance-adjusted designs, and nuisance-adjusted designs with
within-block exchangeability. Richer structured designs still require more
work on the exchangeability side.

## Built-in adapters and maturity

| Adapter                | Geometry  | Relation    | Maturity   |
|------------------------|-----------|-------------|------------|
| `svd_oneblock`         | oneblock  | variance           | **mature** |
| `prcomp_oneblock`      | oneblock  | variance           | **mature** |
| `cross_svd` (cov)      | cross     | covariance         | **mature** |
| `cross_svd` (cor)      | cross     | correlation        | **mature** |
| `cancor_cross`         | cross     | correlation        | **mature** |
| `multivarious_pca`     | oneblock  | variance           | **mature** |
| `multivarious_plsc`    | cross     | covariance         | **mature** |
| `multivarious_cca`     | cross     | correlation        | **narrow** |
| `lda_refit`            | geneig    | generalized_eigen  | **narrow** |
| `plsr_refit`           | cross     | predictive         | **narrow** |

The package uses exactly three maturity labels, in one place, for one
meaning each:

- **mature**: paper-backed exact path for a shipped latent-root
  family, defended by calibration and parity evidence. On the
  correlation side that means the plain paired-row design,
  common-`Z` nuisance-adjusted designs, nuisance-adjusted designs
  with within-block exchangeability, and `blocked_rows(groups)`.
- **narrow**: shipped today with a deliberately limited public
  surface. The engine works end-to-end but the family coverage is
  intentionally scoped — for example `lda_refit` exposes only
  discriminant-root significance for `(geneig, generalized_eigen)`,
  and `plsr_refit` exposes only held-out predictive gain for
  `(cross, predictive)`.
- **planned**: architecturally allocated but not shipped in v1 as a public
  method family — currently bundled multiblock adapters / wrappers beyond the
  downstream adapter path.

Call `list_infer_adapters(details = TRUE)` to read the same table
directly from the registry; the drift-guard test asserts it stays in
sync with this README.

For CCA specifically, `infer_cca()` still defaults to `cancor_cross`. The new
`multivarious_cca` adapter is available as an alternate backend, but it is not
the public default until parity with the shipped path is pinned more broadly.

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
claimed for the supported paired-row, `blocked_rows(groups)`, and
nuisance-adjusted designs, and richer exchangeability structures are still
treated conservatively. Likewise, bootstrap stability defaults to sign
alignment rather than Procrustes rotation, because the package treats subspace
uncertainty as the primary target when roots are near-tied.

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
  nuisance-adjusted designs, nuisance-adjusted within-block designs, and
  `blocked_rows(groups)`
- generalized-eigen component significance for LDA via label permutation and
  B-metric deflation
- sequential Monte Carlo ladder infrastructure
- exact covariance core-space null draws
- partial-SVD and bootstrap/stability performance improvements

Deferred or later-phase work:

- richer structured / hierarchical valid multi-root CCA
- broader `geneig` families beyond LDA
- bundled multiblock method adapters / wrappers
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

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The defaults are configured via `params$family` and `params$preset` (family = 'red', preset = 'homage'). The pkgdown site uses `template: { package: albersdown }` together with generated `pkgdown/extra.css` and `pkgdown/extra.js` so the theme is linked and activated on site pages.
<!-- albersdown:theme-note:end -->
