# Paper 1 Positioning: Rank-Matched Residual Randomization for Latent Operator Models

This note captures the intended framing for a publishable first methods paper. It is deliberately narrower than the full `multifer` / "ShapeInfer" framework story.

Companion note: for the package-level positioning that explains how `multifer` should simultaneously serve as the paper's reference implementation and as a broadly useful R package, see [package_vision.md](./package_vision.md).

Executable manuscript outline: see [paper1_outline.md](../paper/paper1_outline.md).

## One-line pitch

Start from Vitale's three best ideas:

- sequential deflation,
- rank matching after resampling,
- the relative tail-root statistic,

then do three things:

1. make P3 exact and cheap,
2. lift it from PCA matrices to latent operators,
3. exploit exact core-space resampling so the method is fast enough to matter.

That yields a methods paper about **latent-root inference under user-specified null actions**, not "a flexible software framework."

## Recommended title

**Rank-matched residual randomization for latent operator models**

Suggested subtitle:

**Deflation, rank matching, and exact core-space resampling for latent-root inference**

## The central claim

Vitale is best understood as a **pattern**, not as a PCA-only trick.

A Vitale-type procedure extends whenever four ingredients exist:

1. an ordered sequence of nonnegative latent roots,
2. a meaningful way to remove the first `a - 1` latent directions,
3. a null action on the data that destroys only the structure being tested,
4. a rank-matching step after resampling.

The portable part is the shell:

- residualize / deflate,
- apply a null action on residual data,
- rank-match the randomized problem,
- test the next latent object,
- stop at the first non-rejection.

What does **not** carry over for free is the inferential meaning of the statistic. The algebra may port while the valid test target does not.

## The central theorem

The cleanest theorem is still the exact P3 collapse.

Let an SVD/eigen-type operator have ordered nonnegative roots

`rho_1 >= rho_2 >= ... >= rho_r >= 0`.

Define the tail-root statistic

`T_a = rho_a / sum_{q=a}^r rho_q`.

In Vitale's PCA setting, P3 computes a projected randomized residual and then takes the leading projected root. Algebraically, that statistic is exactly equal to taking the `a`-th root of the **unprojected** randomized residual and dividing by its tail sum.

Consequences:

- no explicit projected matrix is needed,
- no second decomposition is needed,
- only the top `a` roots of each randomized residual are needed.

This is the right "small exact theorem":

- it is easy to state,
- easy to remember,
- easy to trust,
- and it moves the method from "PCA matrix manipulation" to "ordered-root statistic on a latent operator."

## The publishable novelty beyond cleanup

The strongest publishable idea after the exact collapse is **exact core-space resampling**.

For two-block covariance models, if

`X = U_x D_x V_x^T`

and

`Y = U_y D_y V_y^T`,

then under a row permutation `P` of `Y`,

`X^T P Y = V_x [ D_x (U_x^T P U_y) D_y ] V_y^T`.

So the nonzero singular values of the permuted cross-operator are exactly the singular values of the small core

`D_x (U_x^T P U_y) D_y`.

This is an exact identity, not an approximation.

For CCA, after whitening, the permuted canonical correlations collapse further to the singular values of

`U_x^T P U_y`.

This is an even cleaner computational result, but the validity caveat is stronger:

- naive permutation CCA is not valid beyond the first canonical root,
- nuisance adjustment needs exchangeability-preserving basis handling,
- stepwise multi-root CCA needs the Winkler-style machinery, not just a fast implementation.

So:

- the **computational identity** for CCA is clean,
- the **inferential validity** is conditional and should be stated that way.

## The abstraction the paper should use

Do not abstract over "methods." Abstract over four objects:

- a latent operator `K(Z)`,
- ordered roots `rho_1, ..., rho_r`,
- a deflation map producing residual data `Z_a`,
- a null action `g in G` that preserves the null but destroys the target structure.

Then the generic algorithm is:

`T_a(Z) = rho_a(K(Z)) / sum_{q=a}^r rho_q(K(Z))`

and

`p_a = (1 + #{ b : T_a(g_b . Z_a) >= T_a(Z_a) }) / (B + 1)`.

Stop at the first non-rejection.

This makes the paper about:

- one statistic,
- one null-action abstraction,
- one sequential rule,
- many latent operators.

That is much cleaner than leading with a pluggable framework story.

## Theoretical boundary statement

This sentence should appear explicitly:

**The exact P3 collapse is portable across SVD/eigen latent-root models, but validity of the resulting ladder still depends on the null action and the inferential meaning of the root.**

An even sharper version for the notes:

**The exact P3 collapse is portable across SVD/eigen latent-root models, but not across all stagewise supervised methods.**

Both are true; the first is the more defensible paper sentence.

## Method scope: what extends directly, what does not

### Direct extension now

- PCA / one-block variance operators
- PLSC / covariance two-block SVD
- closely related SVD/root models where the latent object is genuinely an ordered singular-root sequence

These fit the paper's mainline story.

### Qualified extension inside the paper, broader shipped surface in the package

- CCA / canonical-correlation operators

The package now ships a defended CCA path for paired-row and shipped
nuisance-adjusted designs, including grouped exchangeability, and treats richer
structured designs conservatively. Paper 1 can remain narrower than that.

If CCA is mentioned in Paper 1, it should still be:

- carefully qualified,
- framed as a supported computational extension rather than the theorem-bearing center,
- and kept subordinate to the PCA / PLSC main line, whether in a short section,
  appendix, or explicit supported-extension discussion.

### Next engine, not just a new adapter

- broader symmetric-definite generalized eigenvalue problems,
- contrastive / metric-weighted problems,
- generalized-eigen families beyond the shipped LDA path.

These belong in the framework, but not as silent reuse of the existing oneblock
or cross engines. The package now ships LDA as a narrow public geneig surface;
that does not require Paper 1 to make generalized-eigen inference part of the
main theorem story.

### Different inferential target entirely

- classical PLS regression

PLSR should not be sold as a literal reuse of the PLSC recipe.

Reason:

- for univariate `y`, the raw cross-operator `X^T y` has rank 1,
- yet classical PLSR can still produce multiple latent components,
- therefore a framework based only on ordered roots of `X^T Y` cannot represent PLSR correctly.

The right target is not an association root. It is a **successive predictive increment**, e.g. a statistic built from

`Yhat_a - Yhat_{a-1}`

or an equivalent out-of-sample predictive-strength measure.

So the right distinction is:

- PLSC: association-root testing,
- PLSR: predictive-increment testing.

Reduced-rank regression is an easier supervised extension than classical PLSR
because its latent structure is more directly tied to operator roots.

The package now ships a narrow predictive public surface around PLSR and
held-out predictive gain. That is a package-level extension, not a reason to
broaden Paper 1 beyond its cleaner latent-root center.

## Recommended scope of Paper 1

Keep the paper narrower than the framework.

### Main paper

- one-block variance operator: PCA
- two-block covariance operator: PLSC / covariance SVD
- exact P3 collapse
- operator formulation
- exact core-space resampling identities
- permutation and sign-flip null actions inside one unified residual-randomization framework

The paper can explicitly note that the package also ships a supported CCA path
plus narrower LDA and PLSR surfaces, while keeping those outside the main
theorem-bearing narrative.

### Optional short section or appendix

- canonical correlation as a carefully bounded supported extension
- symmetric-definite generalized eigen after whitening
- LDA as a label-permutation special case

Only include these if they stay concise and do not blur the main statistical
center of gravity.

### Explicitly out of scope for Paper 1

- classical PLS regression
- variable-level significance
- sparse/tuned methods as mainline objects

Sparse / tuned variants should appear only as a brief held-out extension direction, not as the core contribution.

## Four results worth proving

Do not aim for a giant asymptotic theory package. A strong methods paper only needs a few crisp results.

### 1. Exact P3 collapse

The projected P3 statistic equals the `a`-th root tail statistic on the unprojected randomized residual.

### 2. Validity proposition

If the null action leaves the residual-data law invariant under `H_a`, then the Monte Carlo p-value for the `a`-th test is valid.

### 3. Sequential selection proposition

With stop-at-first-non-rejection,

`Pr(rhat > r_0) <= alpha`

because overestimation requires false rejection of the first true null.

### 4. Exact complexity / collapse result

- in PCA: each resample needs only a partial root computation,
- in two-block covariance models: each resample reduces to the small core.

That is enough theory for a rigorous first paper.

## Benchmark design for Paper 1

The benchmark section should matter almost as much as the theorem.

Compare against:

- Vitale P3,
- PA / Horn,
- DPA / deflated DPA,
- FlipPA in heteroskedastic settings,
- valid CCA permutation only if CCA is truly included.

Use four simulation families:

1. pure-noise null calibration,
2. shadowing with weak later roots behind strong earlier ones,
3. heteroskedastic noise,
4. cross-block association with strong within-block structure but weak or null between-block signal.

For cross-block cases, vary `n / p` aggressively because PLS / CCA overfit badly when the sample size is not comfortably larger than the variable count.

## Software framing

The framework is still important, but it should not lead the first paper.

Best split:

### Paper 1

**Rank-matched residual randomization for latent operator models**

Focus:

- statistical method,
- exact theorem,
- exact core-speed identities,
- careful scope statement,
- PCA + PLSC mainline.

### Package state today

The package is now broader than the first paper on purpose:

- `infer_cca()` is a shipped CCA surface with explicit support boundaries,
- `infer_lda()` is a shipped significance-first generalized-eigen surface,
- `infer_plsr()` is a shipped predictive-gain surface.

Paper 1 does not need to broaden to match every shipped wrapper. It only needs
to state that the package contains those extensions while the paper justifies
the narrower theorem-bearing core.

### Paper 2 or software paper

`multifer` / "ShapeInfer"

Focus:

- adapters and typed shapes,
- user-facing wrappers,
- stability outputs,
- extensibility,
- implementation architecture.

That split keeps Paper 1 elegant and statistical rather than turning it into a software architecture paper.

## Recommended language for the notes

The cleanest summary sentence is:

**One statistic, one null-action abstraction, one sequential rule, many latent operators.**

And the cleanest scope sentence is:

**Paper 1 should cover PCA and PLSC directly, treat CCA as a carefully bounded supported extension if included, and leave the shipped LDA / PLSR package surfaces outside the paper's central theorem claims.**
