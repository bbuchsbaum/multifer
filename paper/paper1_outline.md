# Paper 1 Outline

## Working title

**Rank-matched residual randomization for latent operator models**

Suggested subtitle:

**Deflation, rank matching, and exact core-space resampling for latent-root inference**

## One-sentence contribution

We turn Vitale-style sequential deflated permutation inference from a
PCA-specific procedure into a general latent-root testing framework, with an
exact algebraic collapse of P3 and exact core-space resampling identities that
make the method computationally practical for PCA and PLSC.

## Paper scope

Mainline scope:

- one-block variance operators: PCA
- two-block covariance operators: PLSC / covariance SVD
- permutation and sign-flip null actions within one residual-randomization
  framework

Qualified or secondary discussion only:

- CCA: computational identity yes, full validity only with caveats
- generalized eigen / LDA: brief corollary or future-work bridge

Explicitly out of scope:

- classical PLS regression
- sparse / tuned methods as mainline objects
- variable significance
- full software / framework architecture

## Target theorem and proposition list

### Theorem 1: Exact P3 collapse

Statement:

For an SVD/eigen latent-root model with ordered nonnegative roots
`rho_1 >= ... >= rho_r`, the projected P3 statistic at step `a` is exactly
equal to the unprojected tail-root statistic

`T_a = rho_a / sum_{q=a}^r rho_q`

computed on the randomized residual operator.

What this section must do:

- state the original projected statistic clearly
- prove equivalence in 5-10 lines
- isolate the computational consequence: no projected matrix, no second
  decomposition, only top `a` roots needed

### Proposition 2: Monte Carlo validity under a null action

Statement:

If the null action leaves the residual-data law invariant under `H_a`, then the
Monte Carlo p-value

`p_a = (1 + r_a) / (B + 1)`

is valid for the `a`th step.

What this section must do:

- define the residual data object `Z_a`
- state the null-action invariance condition carefully
- make the conditioning structure explicit

### Proposition 3: Overestimation control by stop-at-first-non-rejection

Statement:

If each stepwise test is level `alpha` and the estimated rank is the last
rejected step before the first non-rejection, then

`Pr(rhat > r_0) <= alpha`

for the true rank `r_0`.

What this section must do:

- formalize the nesting of hypotheses
- show that overestimation implies false rejection of the first true null

### Proposition 4: Exact complexity reduction

Statement:

For two-block covariance models with thin SVDs

`X = U_x D_x V_x^T`

and

`Y = U_y D_y V_y^T`,

the nonzero singular values of `X^T P Y` are exactly the singular values of

`D_x (U_x^T P U_y) D_y`.

What this section must do:

- prove the identity cleanly
- make the computational interpretation explicit
- separate exact identities from approximate accelerators

## Proposed section structure

### 1. Introduction

Purpose:

- motivate latent-rank inference as a problem broader than PCA
- explain why raw-root testing and one-shot permutation are inadequate
- position Vitale as the seed, not the endpoint

Required beats:

- sequential deflation matters
- rank matching matters
- the inferential target should be a relative latent-root statistic
- computational cost is the main reason exact randomization methods are ignored

Last paragraph of Introduction:

- state the three main contributions
  - exact P3 collapse
  - latent-operator generalization
  - exact core-space resampling

### 2. Problem setup and notation

Define:

- data object `Z`
- latent operator `K(Z)`
- ordered roots `rho_1, ..., rho_r`
- residual data sequence `Z_a`
- null action `g in G`

Introduce:

- hypotheses `H_a: no latent root beyond a - 1`
- test statistic
  - `T_a(Z) = rho_a(K(Z)) / sum_{q=a}^r rho_q(K(Z))`

Keep this section dry and short.

### 3. From projected P3 to a tail-root test

Contents:

- remind the reader what projected P3 does in Vitale
- present Theorem 1
- explain why this is the core abstraction shift

End of section:

- one short algorithm box for the collapsed PCA procedure

### 4. Rank-matched residual randomization

Contents:

- generic algorithm template
- residualize / deflate
- randomize residual data
- rank-match after randomization
- Monte Carlo p-value
- stop rule

This is where Proposition 2 and Proposition 3 live.

### 5. Exact core-space resampling

Contents:

- derive the two-block covariance identity
- show how it collapses each resample to a small core
- state the CCA analogue as a qualified extension

This is the main computational section and should feel like the paper's second
anchor after Theorem 1.

### 6. Instantiations

Subsections:

- PCA / one-block variance operator
- PLSC / covariance two-block operator

Optional short subsection:

- qualified note on CCA

This section should be concrete:

- what is the operator
- what is the deflation map
- what null action is used
- what validity claim is being made

### 7. Simulation design

Simulation families:

1. Pure-noise null calibration
2. Shadowing: weak later roots behind strong early roots
3. Heteroskedastic noise
4. Cross-block settings with strong within-block structure but weak or null
   between-block association

Factors to vary:

- `n`
- `p` and `q`
- signal strengths
- root spacing
- heteroskedasticity severity
- null action: permutation vs sign-flip where appropriate

Primary metrics:

- overestimation rate
- exact-rank recovery
- per-step rejection profile
- wall time

Comparators:

- Vitale P3
- Horn / parallel analysis
- DPA / deflated parallel analysis
- FlipPA in heteroskedastic settings

### 8. Results

Recommended order:

1. Null calibration
2. Shadowing recovery
3. Heteroskedastic robustness
4. Cross-block recovery
5. Computational gains

Figures to plan now:

- Figure 1: schematic of sequential residual randomization
- Figure 2: projected P3 versus collapsed tail-root view
- Figure 3: null calibration curves
- Figure 4: shadowing recovery heatmap
- Figure 5: cross-block recovery and runtime

Tables:

- Table 1: theorem / proposition summary
- Table 2: simulation grid
- Table 3: runtime comparison

### 9. Discussion

Core messages:

- the paper is about a reusable inferential pattern, not a software framework
- algebraic portability does not imply inferential portability
- PCA and PLSC are the justified mainline scope
- CCA is computationally aligned but inferentially qualified

Subsection: boundaries

- why PLSR is not a direct extension
- why generalized-eigen methods are the next engine, not a free corollary
- why sparse/tuned methods need held-out or split-sample logic

### 10. Software note

One short section only.

State:

- `multifer` is the reference implementation
- the package is broader than the paper's theorem scope
- current package validity depends on engine and design

Do not let this turn into a framework paper inside the methods paper.

## Abstract skeleton

Sentence 1:

- latent-root inference appears in PCA and many related multivariate methods,
  but practical randomization procedures remain method-specific and expensive

Sentence 2:

- Vitale's sequential deflation and rank-matched testing provide the right
  starting point

Sentence 3:

- we show that the projected P3 statistic admits an exact tail-root collapse

Sentence 4:

- this turns the method into an operator-level procedure and removes an entire
  projection/decomposition layer from each resample

Sentence 5:

- for two-block covariance models we derive exact core-space resampling
  identities that collapse each permutation to a small inner matrix

Sentence 6:

- we prove validity under null-action invariance and overestimation control
  under stop-at-first-non-rejection

Sentence 7:

- simulations show improved calibration/recovery over existing baselines and
  major computational gains

Sentence 8:

- the resulting procedure supports PCA directly, extends naturally to PLSC, and
  clarifies the boundary between portable latent-root models and methods that
  require different inferential targets

## Introduction checklist

Must cite:

- Vitale et al. for P3, deflation, and sequential testing
- Dobriban and Owen for deflated parallel analysis / shadowing motivation
- FlipPA for heteroskedastic-valid null generation
- Winkler et al. for CCA permutation validity caveats

Must not overclaim:

- do not say the paper solves all latent multivariate models
- do not imply CCA is fully solved unless the full nuisance-aware and
  multi-root validity machinery is included in the paper body
- do not present PLSR as an operator-root special case

## What to write first

Recommended drafting order:

1. Section 3: exact P3 collapse
2. Section 5: exact core-space resampling
3. Section 2: notation and setup
4. Section 4: validity propositions
5. Introduction
6. Discussion

This order writes the paper around the actual core contributions rather than
around generic framing.
