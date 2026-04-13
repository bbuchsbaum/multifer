# Architecture Rules for `multifer`

This note records the non-negotiable design rules for `multifer`.

The package's main goals are:

- generality,
- performance,
- modularity,
- abstraction.

Those goals are compatible only if the package resists ad hoc,
method-specific branching.

## Core rule

`multifer` should be unified at the level of its inferential scaffold,
not by pretending every multivariate method uses the same inner
statistic.

The shared scaffold is:

- ordered latent objects,
- sequential removal of previously claimed structure,
- null actions matched to design and inferential target,
- stop-at-first-non-rejection testing,
- latent-unit stability summaries.

This scaffold should remain common across:

- PCA,
- PLSC,
- CCA,
- later `geneig` engines,
- and feasible supervised extensions such as PLS regression.

## Rule 1: Prefer a small number of engine families

Do not add one-off engines for each named method unless the underlying
inferential target is genuinely different.

The intended engine split is:

- direct latent-root engines
  - PCA
  - PLSC
  - CCA, with the appropriate validity machinery
  - symmetric-definite `geneig`
- supervised-increment engines
  - PLS regression
  - later predictive / reduced-rank supervised variants

If a new method fits an existing engine family, implement it as an
operator/adapter variation, not as a bespoke branch.

## Rule 2: Treat exact reductions as reusable capabilities

Core-space updates, reduced-space resampling, and lift-back formulas
should be implemented as reusable engine capabilities, not as isolated
optimizations tied to a single method.

Examples:

- exact null core updates,
- exact bootstrap core updates,
- reduced-space lift-back of loadings and scores,
- reduced-space covariance summaries for standard errors or intervals.

Performance should come from exact low-dimensional reductions where
available, not from method-specific approximation hacks scattered across
the codebase.

## Rule 3: No alignment logic inside significance statistics

Significance should depend on invariant quantities:

- latent roots,
- relative root statistics,
- or other directly justified test targets.

Do not use label-fixing or Procrustes-like alignment inside the
significance statistic itself.

Alignment belongs, if anywhere, in the reliability/stability layer.

## Rule 4: Units, not components, are the stable abstraction

The package should attach downstream outputs to latent units rather than
raw component indices.

That means:

- singleton units when axes are identifiable,
- subspace units when roots are tied or near-tied.

This avoids method-specific hacks for axis swapping and keeps variable,
score, and subspace summaries coherent across methods.

## Rule 5: Allow the inner target to vary when needed

The scaffold is shared. The inner target is not always shared.

Direct latent-root methods can use ordered-root statistics directly:

- PCA
- PLSC
- CCA
- `geneig`

Supervised stagewise methods may require a different target:

- PLS regression should use predictive increments or a similarly
  justified supervised target, not a forced reuse of an `X^T Y` root
  test.

This distinction is essential for keeping the package both general and
mathematically honest.

## Rule 6: Prefer explicit abstractions over branching

When a new method appears to need "just one more special case," first
ask whether one of these should be extended instead:

- operator constructor,
- deflation rule,
- null-action class,
- engine capability,
- unit-formation rule,
- stability summarizer.

If the answer is yes, extend the abstraction rather than adding a branch
to a top-level engine.

## Rule 7: The package should refuse incoherent requests

Strict dispatch and validity labeling are part of the architecture, not
UI polish.

If a requested analysis is ambiguous or its validity is not yet
justified, `multifer` should:

- require explicit disambiguation,
- downgrade claims conservatively,
- or refuse the analysis.

Do not preserve a "simple API" by silently guessing.

## What this means for roadmap decisions

These rules imply:

- extend PLSC through exact core-space bootstrap before inventing more
  ad hoc alignment fixes,
- add task PLS and behavior PLS as operator variants if they fit the
  same covariance-style engine cleanly,
- implement PLS regression as a separate supervised engine family,
- add `geneig` as a direct latent-root engine rather than forcing it
  through existing oneblock or cross branches.

## One-line doctrine

Unify the scaffold. Reuse exact reductions. Keep engine families small.
Let the latent target change only when it must.
