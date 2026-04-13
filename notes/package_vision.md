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

**The reference implementation of rank-matched residual randomization for latent-root inference, plus a general perturbation-inference framework for latent multivariate models.**

That wording matters.

It says:

- some parts of the package are directly justified by the paper's main results,
- other parts belong to the same framework but may require method-specific validity conditions or later engines.

This is honest and strong.

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

- method-named wrappers like `infer_pca()`, `infer_plsc()`, later `infer_lda()`,
- stability summaries for variables, scores, and subspaces,
- good defaults for common workflows,
- clear messages when a requested inference is not valid,
- adapter / engine extensibility for advanced users,
- room for later engines such as `geneig` and supervised-latent / predictive engines.

This layer is how the package becomes genuinely useful rather than merely "the code for the paper."

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

`multifer` is the reference implementation of the paper's method for PCA and PLSC, and a broader perturbation-inference framework for latent models whose exact validity depends on the engine and null action.

That sentence is strong enough for a README, a package website, or a software paper.

## README positioning draft

If a top-level `README` is added, the opening positioning language should be close to this:

> `multifer` provides perturbation-based inference for latent multivariate models.
> It is the reference implementation of rank-matched residual randomization for latent-root inference in PCA-like and covariance two-block models, and it also provides a broader framework for stability analysis and method-specific inference across latent model families.
> The package emphasizes strict validity contracts: it prefers refusing ambiguous or invalid analyses over silently guessing.

Then a short second paragraph:

> The core statistical engine is based on sequential deflation, rank-matched residual randomization, and stop-at-first-non-rejection testing of ordered latent structure. Around that core, `multifer` provides practical R workflows for component significance, variable stability, score stability, and subspace uncertainty.

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
