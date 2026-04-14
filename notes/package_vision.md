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
- CCA with the right validity machinery,
- later symmetric-definite `geneig` engines.

For supervised stagewise methods, the scaffold remains shared but the target may
change. PLS regression is the clearest example: it belongs in the framework, but
its inferential target should be successive predictive increments rather than a
literal reuse of an `X^T Y` root statistic.

## Recommended positioning

`multifer` should be presented as:

**A typed perturbation inference platform for latent-root multivariate models, with mature support for PCA and covariance-mode two-block methods, qualified support for canonical-correlation models, and an extensible adapter / typed-shape architecture that brings new ordered-latent-root families into the same shared ladder engine without modifying it.**

That wording matters.

It says:

- some parts of the package are directly justified by the paper's main results ("mature"),
- some parts are in the framework and run today but have validity qualifications that have not yet been fully settled ("qualified"),
- other parts belong to the same framework but are architectural slots awaiting real engines ("planned").

This is honest and strong.

### What we do *not* claim

The tempting phrase "maximally general inference machinery" is not the target. The target is narrower and more defensible:

**Maximally reusable across a small number of mathematically coherent inferential families.**

Those families are ordered-root SVD models, canonical-correlation models, generalized-eigen / metric-weighted models, and supervised predictive-increment models. The paper's mainline theorem covers the first of those; the package exists to reuse as much of that scaffold as possible across the other three without pretending every multivariate method reduces to the same root test.

## Support maturity levels

Every inferential family the package ships is labeled with one of three maturity levels. This taxonomy is canonical and every public doc should match it.

### Mature — in Paper 1 scope, exact results, production-ready

- **One-block variance** (PCA-family): `adapter_svd`, `adapter_prcomp`, `multivarious_pca`. Collapsed Vitale P3 ladder is exact and tested against the original projected construction. Core-space bootstrap via the Fisher identity is exact.
- **Cross-block covariance** (PLSC-family): `adapter_cross_svd` in covariance mode, `multivarious_plsc`. Exact core-space bootstrap path via `(V_x, D_x, V_y, D_y)` identity is in place and machine-precision verified against refit across six shape/SNR regimes in `notes/cross_core_bootstrap_study.md`.

These are what the paper claims. These are what the README should lead with.

### Qualified — running today, but validity story is not yet fully settled

- **Cross-block correlation** (CCA-family): `adapter_cross_svd` in correlation mode, `adapter_cancor`. The architecture is ready. The reference adapter's multi-root deflation is currently a "safe simplification" that is exact for paired-rows and nuisance-adjusted designs but a conservative first-root-only cap for more complex designs. `checked_assumptions` is still partly declarative.

The package runs these, reports valid first-root p-values, and labels multi-root claims carefully. The next maturity step requires exact whitened-space stepwise deflation + design-aware null actions for paired-rows, nuisance-adjusted, and grouped-exchangeability designs, plus concrete validity checks. Until that lands, CCA is marketed as the *beginning* of a CCA family, not the finished version.

### Planned — architectural slot, no engine yet

- **Multi-block** geometry: declared in the vocabulary and the capability matrix; `bootstrap_fits()` explicitly refuses. Needs a real engine with the generalized core-space identity for `K = sum_{ij} X_i^T W_{ij} Y_j`.
- **Generalized eigenvalue** geometry (`geneig`): declared; no engine. Needs operator pair, metric-aware deflation, root statistic, and null family. LDA, contrastive PCA, and metric-weighted PCA land here as special cases.
- **Supervised / predictive** relation (PLS regression, reduced-rank regression): not in the relation vocabulary yet. Needs its own relation family whose inferential target is a *successive predictive increment* (`Yhat_a - Yhat_{a-1}` or a cross-fitted predictive-strength measure), not an association root. Reusing `covariance` for PLSR would corrupt the meaning of the covariance relation and is explicitly ruled out.

These show up in documentation as future work, not as present-tense capability. Public framing never says they work today.

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

- method-named wrappers like `infer_pca()`, `infer_plsc()`, `infer_cca()`, later `infer_lda()`,
- stability summaries for variables, scores, and subspaces,
- good defaults for common workflows,
- clear messages when a requested inference is not valid,
- adapter / engine extensibility for advanced users,
- room for later engines such as `geneig` and supervised-latent / predictive engines.

This layer is how the package becomes genuinely useful rather than merely "the code for the paper."

Today, the package already ships thin wrappers for the mature / qualified
families (`infer_pca()`, `infer_plsc()`, `infer_cca()`). The remaining package
work is polish, not first introduction of wrapper-level usability.

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

`multifer` is a typed perturbation inference platform for latent-root multivariate models. It delivers mature, paper-backed inference for PCA-family and covariance two-block methods, qualified support for canonical-correlation models, and an extensible adapter / typed-shape architecture so new ordered-latent-root families can be brought into the same ladder engine without modifying it.

That sentence is strong enough for a README, a package website, or a software paper.

## README positioning draft

If a top-level `README` is added, the opening positioning language should be close to this:

> `multifer` provides perturbation-based inference for latent-root multivariate models. It ships a typed inferential scaffold — ordered latent objects, sequential deflation, rank-matched residual randomization, stop-at-first-non-rejection testing, and latent-unit stability summaries — and applies that scaffold across a small set of mathematically coherent method families.
>
> **Mature support** (paper-backed, exact): PCA-family one-block variance inference and covariance-mode two-block methods (PLSC-family), both with the exact collapsed Vitale P3 ladder and exact core-space bootstrap.
>
> **Qualified support**: canonical-correlation models, currently valid for the leading canonical root and for paired-rows / nuisance-adjusted stepwise designs, with the full multi-root validity machinery tracked as an upcoming release.
>
> **Extensible by design**: new ordered-latent-root inferential families are added through the shared `infer_adapter()` contract, so the ladder driver, null-action machinery, bootstrap loop, and stability consumers apply uniformly to whatever family an adapter declares. The goal is to stay maximally general *within* the ordered-latent-root paradigm rather than to absorb every multivariate method.
>
> The package emphasizes strict validity contracts: it prefers refusing ambiguous or invalid analyses over silently guessing, and distinguishes component significance (permutation-based) from stability reliability (bootstrap-based) at the schema level.

## Roadmap to paper-quality v1

The concrete next milestone is a paper-quality v1 boundary: mature core, qualified CCA track, explicit planned engines. This is the seven-point program the package should execute against before Paper 1 submission.

1. **Freeze the theorem-bearing core around latent-root SVD models.** One-block variance and cross-block covariance are where the math is cleanest and the speed story is strongest. The collapsed Vitale P3 simplification and the exact core-space bootstrap identity both live here. Every future claim should point back to this core.

2. **Downgrade CCA from "general support" to "qualified support" until exactified.** Keep `adapter_cancor` and correlation-mode `cross_svd` in the shipping set, but label them as the beginning of a CCA family, not the finished version. The next step is exact whitened-space stepwise deflation plus design-aware null actions for paired rows, nuisance-adjusted rows, and grouped exchangeability.

3. **Create a new supervised relation for PLS regression / reduced-rank regression.** Do not reuse `covariance`. Add a relation whose inferential target is predictive gain — a held-out or cross-fitted statistic, or a `Yhat_a − Yhat_{a-1}` increment — rather than a cross-covariance root. This keeps the meaning of `covariance` uncorrupted and gives PLSR its own validity story.

4. **Build a real `geneig` engine.** The vocabulary already has `geneig` / `generalized_eigen`; give them actual mathematics: operator pair, metric-aware deflation, root statistic, and null family. Contrastive PCA, LDA, and metric-weighted PCA land as special cases. Until this engine exists, generalized-eigen support is declared, not present-tense.

5. **Make validity contracts executable, not just declared.** Every mature adapter should populate `checked_assumptions` with concrete functions: paired-row checks, centering checks, rank/conditioning checks, nuisance-design checks, grouped-exchangeability checks. This turns "strict dispatch" from an API idea into a statistical guarantee.

6. **Keep significance and stability separate.** The schema already does this: `component_tests` are permutation-based, `variable_stability` / `score_stability` / `subspace_stability` are bootstrap-based. Preserve that line — McIntosh-style PLS workflows have always distinguished LV significance from salience reliability, and this distinction is part of what makes the platform trustworthy.

7. **Promote subspace inference as the default higher-level object.** `form_units()` and `subspace_stability_from_bootstrap()` already treat near-tied roots as subspace bundles. Lean into that: when roots are close or ties are detected, user-facing output should show subspaces and principal angles, not signed loading vectors, as the primary inferential object.

The v1 boundary is reached when (1)–(6) are done, (7) is the default user experience, and the speed story is tight enough to make the §40 benchmarks pass with the exact engines.

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
