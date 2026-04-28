# Package Vision: `multifer` as Both Paper Implementation and Useful R Package

This note states the intended relationship between the methods paper and the package.

Companion note: for the explicit architecture rules that prioritize generality,
performance, modularity, and abstraction over ad hoc special casing, see
[architecture_rules.md](./architecture_rules.md).

The package should do two things at once:

1. serve as the **reference implementation** of the paper's core method,
2. serve as an **extremely useful R package** for perturbation-based inference on latent multivariate models.

These goals are compatible, but only if they are separated clearly.

## Core principle

The package should not mirror the paper too literally.

A good methods paper needs:

- one tight statistical story,
- a small number of exact results,
- a sharply defined scope,
- and careful validity claims.

A useful R package needs:

- stable user-facing workflows,
- sensible defaults,
- informative errors,
- practical summaries and stability outputs,
- and room for methods whose validity story differs from the paper's main theorem.

So the right framing is:

- the **paper** defines the core statistical object,
- the **package** implements that object faithfully,
- and the package also provides the broader user interface and roadmap around it.

## Cohesion doctrine

`multifer` should be unified at the level of its inferential scaffold, not by
pretending every multivariate method is the same root test.

The scaffold is:

- ordered latent objects,
- sequential removal of previously claimed structure,
- null actions matched to design and inferential target,
- stop-at-first-non-rejection testing,
- latent-unit stability summaries.

For direct latent-root models, this shared scaffold also supports a shared
ordered-root target:

- PCA,
- PLSC,
- CCA with the shipped support boundary,
- generalized-eigen models where the current public surface is LDA,
- later broader symmetric-definite `geneig` engines.

For supervised stagewise methods, the scaffold remains shared but the target may
change. PLS regression is the clearest example: it belongs in the framework, but
its inferential target should be successive predictive increments rather than a
literal reuse of an `X^T Y` root statistic.

## Recommended positioning

`multifer` should be presented as:

**A typed perturbation inference platform for latent multivariate models, with mature support for one-block variance and two-block latent-root families, narrower but shipped surfaces for generalized-eigen and predictive-gain inference, and an extensible adapter / typed-shape architecture that brings new ordered-latent families into the same shared scaffold without modifying it.**

That wording matters.

It says:

- some parts of the package are directly justified by the paper's main results and now have a mature shipped path,
- some parts are intentionally narrow public surfaces that already run today but are not the center of the paper's theorem story,
- other parts remain architectural slots awaiting a real engine.

This is honest and strong.

### What we do *not* claim

The tempting phrase "maximally general inference machinery" is not the target. The target is narrower and more defensible:

**Maximally reusable across a small number of mathematically coherent inferential families.**

Those families are ordered-root SVD models, canonical-correlation models, generalized-eigen / metric-weighted models, and supervised predictive-increment models. The paper's mainline theorem covers the first of those; the package exists to reuse as much of that scaffold as possible across the other three without pretending every multivariate method reduces to the same root test.

## Support maturity levels

Every inferential family the package ships is labeled internally with one of three maturity levels. This taxonomy is the **internal scope-decision vocabulary**: it is what the authors use to decide what is in scope for Paper 1, what is shipping today, and what is explicitly on the "not yet" list.

**Public framing intentionally differs.** The README and DESCRIPTION do not enumerate the planned tier by name; they describe extensibility in terms of the shipping `infer_adapter()` contract, so package copy never reads like a release promise. This note, `notes/paper1_rank_matched_residual_randomization.md`, and the paper outline are where the full taxonomy lives.

### Mature — shipped, defended, and central to the package story

- **One-block variance** (PCA-family): `adapter_svd`, `adapter_prcomp`, `multivarious_pca`. Collapsed Vitale P3 ladder is exact and tested against the original projected construction. Core-space bootstrap via the Fisher identity is exact.
- **Cross-block covariance** (PLSC-family): `adapter_cross_svd` in covariance mode, `multivarious_plsc`. Exact core-space bootstrap path via `(V_x, D_x, V_y, D_y)` identity is in place and machine-precision verified against refit across six shape/SNR regimes in `notes/cross_core_bootstrap_study.md`.
- **Cross-block correlation** (CCA-family): `adapter_cross_svd` in correlation mode and `adapter_cancor`. The current implementation supports multi-root testing for the paired-row design and for the shipped nuisance-adjusted designs, including grouped exchangeability. Designs outside that support matrix fall back conservatively to first-root inference.

These are what the README should lead with. Paper 1 may still emphasize PCA + PLSC most heavily, but the package-level CCA story is now a shipped and defended path rather than a merely aspirational extension.

### Narrow public surface — shipped today, but intentionally scoped

- **Generalized eigen** (`geneig`, `generalized_eigen`): the public wrapper is `infer_lda()` via `lda_refit`. This is a significance-first surface centered on discriminant roots. Broader metric-weighted / contrastive generalized-eigen families are not yet part of the public wrapper story.
- **Predictive gain** (`cross`, `predictive`): the public wrapper is `infer_plsr()` via `plsr_refit`. The inferential target is held-out predictive gain, not a recycled covariance root. The public predictive surface is intentionally narrow in v1: PLSR is shipped, broader predictive-cross families remain future work.

These surfaces are real and useful, but they are not yet the center of the package's high-level positioning sentence. Their immediate next step is boundary clarification, calibration evidence, and clearer author/user guidance rather than rapid expansion of scope.

### Planned — architectural slot or broader family expansion, not current package surface

- **Multi-block** geometry: declared in the vocabulary and the capability matrix; `bootstrap_fits()` explicitly refuses. Needs a real engine with the generalized core-space identity for `K = sum_{ij} X_i^T W_{ij} Y_j`.
- **Broader generalized-eigen families**: contrastive PCA, metric-weighted PCA, and other `geneig` special cases beyond LDA.
- **Broader supervised / predictive families**: reduced-rank regression and richer predictive-cross engines beyond the current PLSR wrapper.

These show up in documentation as future work, not as present-tense capability. Public framing should distinguish the shipped narrow surfaces from the still-unbuilt family expansions.

## The two-layer architecture

The package should have two explicit layers.

### 1. Paper-faithful core

This layer should implement the methods that the first paper actually claims.

Current / intended contents:

- oneblock variance operators: PCA-like models,
- cross covariance operators: PLSC / covariance SVD,
- exact P3 collapse where valid,
- explicit residual randomization / null-action machinery,
- strict validity labeling and refusal to guess when inference is ambiguous.

This layer is the scientific core. It should be exact, carefully tested, and easy to point to in the paper.

### 2. Community-facing package layer

This layer should make the package broadly useful to R users, even when they never read the paper.

Desired contents:

- method-named wrappers like `infer_pca()`, `infer_plsc()`, `infer_cca()`, `infer_lda()`, and `infer_plsr()`,
- stability summaries for variables, scores, and subspaces,
- good defaults for common workflows,
- clear messages when a requested inference is not valid,
- adapter / engine extensibility for advanced users,
- room for later engines such as `geneig` and supervised-latent / predictive engines.

This layer is how the package becomes genuinely useful rather than merely "the code for the paper."

Today, the package already ships thin wrappers for the mature and narrow public
surfaces (`infer_pca()`, `infer_plsc()`, `infer_cca()`, `infer_lda()`,
`infer_plsr()`). The remaining package work is mostly evidence hardening and
boundary clarification, not first introduction of wrapper-level usability.

## What the paper should and should not control

The paper should control:

- the statistical core,
- the main theorem statements,
- the mainline algorithms,
- the benchmark logic for the claimed scope.

The paper should **not** force:

- the exact package API shape,
- the full roadmap,
- all future engines to fit one theorem,
- every user-facing summary to be in scope for Paper 1.

That distinction is necessary if `multifer` is going to become a real package for the R community.

## What usefulness means for the R community

To be broadly useful, the package should optimize for three things.

### Trustworthiness

The package should refuse invalid analyses rather than silently produce answers.

That means:

- strict dispatch by default,
- explicit relation / design requirements,
- validity labels that are meaningful,
- and conservative fallbacks when a method is only partially validated.

### Convenience

Users should not need to think in internal geometry/relation language for routine workflows.

They should be able to ask for:

- PCA significance,
- PLSC stability,
- later LDA significance,

through direct wrappers and task-oriented docs.

### Extensibility

Advanced users should be able to add:

- adapters,
- null actions,
- engines,
- and specialized validity logic

without forking the package architecture.

## Recommended public framing

The cleanest public story is:

`multifer` is a typed perturbation inference platform for latent multivariate models. It delivers mature, paper-backed inference for PCA-family and two-block latent-root methods, narrower but shipped public surfaces for generalized-eigen discriminant inference and predictive-gain inference, and an extensible adapter / typed-shape architecture so new ordered-latent families can be brought into the same scaffold without modifying it.

That sentence is strong enough for a README, a package website, or a software paper.

## README positioning draft

If a top-level `README` is added, the opening positioning language should be close to this:

> `multifer` provides perturbation-based inference for latent-root multivariate models. It ships a typed inferential scaffold — ordered latent objects, sequential deflation, rank-matched residual randomization, stop-at-first-non-rejection testing, and latent-unit stability summaries — and applies that scaffold across a small set of mathematically coherent method families.
>
> **Mature support** (paper-backed, exact): PCA-family one-block variance inference and two-block latent-root methods, with the exact collapsed Vitale P3 ladder and exact core-space bootstrap where the mathematics is clean.
>
> **Mature CCA path with explicit support boundaries**: canonical-correlation models currently support multi-root inference for the paired-row and shipped nuisance-adjusted designs, with richer structured designs handled conservatively by first-root capping rather than by unsupported stepwise claims.
>
> **Shipped but intentionally narrow public surfaces**: `infer_lda()` exposes discriminant-root significance for the current generalized-eigen family, and `infer_plsr()` exposes held-out predictive-gain inference for the current predictive family. Broader family expansion remains future work.
>
> **Extensible by design**: new ordered-latent inferential families are added through the shared `infer_adapter()` contract, so the ladder driver, null-action machinery, bootstrap loop, and stability consumers apply uniformly to whatever family an adapter declares. The goal is to stay maximally general *within* a small set of coherent inferential families rather than to absorb every multivariate method.
>
> The package emphasizes strict validity contracts: it prefers refusing ambiguous or invalid analyses over silently guessing, and distinguishes component significance (permutation-based) from stability reliability (bootstrap-based) at the schema level.

## v1 support boundary

The current package boundary is not another method family. It is the
paper-quality v1 boundary that the notes, wrappers, and evidence state the same
way.

1. **The theorem-bearing core is frozen around latent-root SVD models.** One-block variance and cross-block covariance remain the cleanest theorem and speed center of gravity. Every broader claim is phrased relative to that center.

2. **The shipped CCA boundary is explicit.** The package ships a real CCA path with a support matrix, executable validity checks, calibration evidence, and a clear supported-design versus conservative-fallback boundary.

3. **The public `geneig` surface is intentionally narrow.** `infer_lda()` is explicit about what it does and does not expose: discriminant-root significance is public, broader metric-weighted / contrastive families are not.

4. **The public predictive surface is intentionally narrow.** `infer_plsr()` is explicit that predictive inference means held-out predictive gain and that v1 centers PLSR rather than a larger predictive-cross family.

5. **Validity contracts are executable, not just declared.** Every shipped adapter family carries concrete checks: rank, sample-size, nuisance-design, grouped-exchangeability, label-permutation admissibility, and similar support-boundary guards.

6. **Significance and stability stay separate.** `component_tests` are significance outputs; variable, score, subspace, and feature-importance summaries are bootstrap reliability outputs. That separation is prominent in docs and print paths.

7. **Subspace-first output remains the default higher-level object.** Near ties print and summarize as subspaces and principal angles before singleton signed directions.

The v1 boundary is the agreement among the notes, wrappers, and evidence on
those seven points. The package describes its present scope without mixing
current capability, paper scope, and future work.

## Paper / package split

The recommended split remains:

### Paper 1

Statistical methods paper:

- exact P3 collapse,
- latent-operator formulation,
- exact core-space resampling,
- PCA + PLSC mainline scope.

### Later software / framework paper

Package paper:

- adapters and typed shapes,
- user-facing wrappers,
- stability outputs,
- extensibility,
- roadmap across additional engines.

This keeps Paper 1 elegant and statistical while still allowing the package to become a broadly useful R tool.

## One sentence to keep in mind

**The paper should justify the core. The package should make the core usable, extensible, and safe.**
