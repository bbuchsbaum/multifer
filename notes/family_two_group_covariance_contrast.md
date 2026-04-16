# Two-Group Covariance Contrast as the First Post-LDA `geneig` Family

This note defines the first non-LDA generalized-eigen family that should be
considered for `multifer`: a **two-group covariance contrast** family.

The purpose of the note is not to add a public adapter immediately. The purpose
is to set the **admission standard** for a family whose inner mechanics fit the
existing `geneig` ladder, but whose validity depends on a family-specific null
doctrine rather than on generalized-eigen algebra alone.

Companion notes:

- [engine_geneig_spec.md](./engine_geneig_spec.md)
- [architecture_rules.md](./architecture_rules.md)
- [package_vision.md](./package_vision.md)

## Status

This is a **design / governance artifact** only.

It does not imply:

- a shipped public wrapper,
- a registered adapter,
- or a claim that the family is already calibrated.

It exists to make those later decisions explicit and auditable.

## Family identity

- `family_id`: `two_group_covariance_contrast`
- `geometry`: `"geneig"`
- `relation`: `"generalized_eigen"`
- inferential family: direct latent-root

This family is meant to answer:

> Along which directions does group 1 exhibit systematically different variance
> structure from group 2, and how many ordered generalized roots rise above a
> label-exchangeability null?

## Raw-data contract

The family should operate on **raw observations plus a binary group label**, not
on precomputed covariance matrices.

Proposed raw-data payload:

```r
list(
  X = <numeric matrix, n x p>,
  g = <factor with exactly two levels>,
  ridge_rel = <non-negative scalar, default small>,
  center = TRUE
)
```

### Why raw data matters

The null action for this family is a claim about **group-label
exchangeability of rows**. That null acts upstream on the observations. It is
not scientifically equivalent to perturbing already-built covariance matrices.

Therefore:

- resampling must act on `(X, g)`,
- then rebuild `A` and `B`,
- then run the generalized-eigen ladder.

An adapter that silently resamples or perturbs `A` and `B` directly should be
treated as invalid unless that behavior is itself the declared null.

## Operator construction

Let:

- `X1 = X[g == level_1, , drop = FALSE]`
- `X2 = X[g == level_2, , drop = FALSE]`

Define group covariances:

```r
S1 = cov(X1)
S2 = cov(X2)
```

The generalized-eigen pencil is:

```r
A = S1
B = S2 + gamma * I
```

where `gamma` is a **relative ridge**:

```r
gamma = ridge_rel * trace(S_pool) / p
```

with `S_pool` the pooled covariance of all rows in `X`.

### Why relative ridge, not absolute ridge

If `gamma` is absolute, common scaling of the features changes the generalized
problem in a way that is partly a units artifact. A relative ridge tied to the
pooled scale preserves the natural common-scaling invariance of the family.

This family should therefore prefer a **relative regularization convention**.

## Estimand

The estimand is the ordered generalized spectrum of the pair `(A, B)`:

```r
A v = lambda B v
```

with `B` symmetric positive definite after ridge regularization.

Interpretation:

- `lambda > 1` indicates directions with larger variance in group 1 than in
  group 2 under the chosen regularization convention,
- `lambda < 1` indicates the opposite,
- the ordered sequence of generalized roots defines the latent object tested by
  the ladder.

The ladder target remains the generalized-root tail ratio already defined by the
`geneig` engine. No family-specific engine is required.

## Null doctrine

### Scientific null

Under the null, **group labels are exchangeable over pooled rows** and there is
no structured covariance contrast attributable to the two groups beyond random
allocation of observations.

This is the same kind of doctrine as LDA label permutation, but with covariance
contrast rather than discriminant scatter as the target.

### Resampling law

The null action is:

1. hold `X` fixed,
2. permute the binary group labels `g`,
3. preserve the original group counts,
4. rebuild `S1`, `S2`, `A`, and `B`,
5. rerun the same generalized-root statistic.

This is the only candidate null that should be considered for the first version
of this family.

### What this family explicitly does not claim

It does **not** claim validity for:

- precomputed covariance matrices without raw rows,
- arbitrary contrastive-PCA-like background/target problems,
- non-exchangeable sampling schemes,
- more than two groups,
- indefinite or rank-deficient metrics without ridge stabilization.

Those are separate family problems and need their own notes.

## Invariances

The family should ship only if these invariances are made explicit and tested.

### 1. Common orthogonal basis change

For any orthogonal matrix `Q`, replacing `X` by `X %*% Q` should preserve:

- the ordered generalized roots,
- the ladder statistics,
- the rejection decisions.

This is the core basis-invariance property of a covariance-contrast family.

### 2. Row-order invariance

Permuting the rows of `X` and permuting `g` in the same way should not change
the result.

### 3. Common positive scaling

Multiplying `X` by a positive scalar `c` should leave the inferential result
unchanged **when the family uses relative ridge regularization**.

Without relative ridge, this is not a valid invariance, which is another reason
the family should avoid an absolute ridge convention.

## Oracle

The family must define a direct spectral oracle independent of the ladder:

1. build `A` and `B` from the raw data,
2. whiten to `C = B^{-1/2} A B^{-1/2}`,
3. compute the ordered eigenvalues of `C`,
4. compare those to the adapter's observed generalized roots.

This oracle is the differential reference for:

- the observed spectrum,
- the rung statistic,
- and the effect of B-metric deflation.

Canonical planted pencils should also be used where the generalized spectrum is
known by construction.

## Admissibility checks

The family should refuse analysis unless all of the following hold:

- `X` is a finite numeric matrix,
- `g` is a factor with exactly two non-empty levels,
- each group has at least two rows,
- the pooled design supports group-label exchangeability,
- the relative ridge is non-negative and finite,
- the rebuilt `B` is SPD after ridge,
- no missing or non-finite entries are present,
- the requested design does not contradict pooled label exchangeability.

These checks belong in the family's `checked_assumptions` / admissibility layer,
not as soft documentation only.

## Proposed `family_spec()` shape

This is the target abstraction, not yet implemented:

```r
family_spec(
  family_id = "two_group_covariance_contrast",
  geometry = "geneig",
  relation = "generalized_eigen",
  estimand = "ordered generalized roots of (cov_g1, cov_g2 + ridge)",
  build_operators = <function(raw_data) -> list(A = ..., B = ..., metric = ...)>,
  null_hypothesis = "group labels exchangeable over pooled rows",
  null_action = <function(raw_data) -> raw_data_with_permuted_labels>,
  invariances = c(
    "common orthogonal feature transform",
    "row-order invariance",
    "common positive scaling under relative ridge"
  ),
  oracle = <function(raw_data) -> direct generalized spectrum>,
  admissibility_checks = <list of executable checks>
)
```

The key rule is:

> The family spec owns the scientific meaning of the pencil. The engine owns the
> outer generalized-root ladder only.

## Required test blocks

This family must not ship without all four blocks below.

### 1. Null calibration

Simulate data under the null with identical group covariance structure and
exchangeable labels.

Required checks:

- rung-1 type-I error near nominal alpha,
- calibration across balanced and moderately unbalanced group sizes,
- no pathological p-value concentration at 0 or 1 under the null,
- finite outputs throughout.

This is the core scientific validity test.

### 2. Direct spectral oracle

Use planted covariance pairs and direct whitening-based generalized eigenvalue
calculation to verify:

- observed roots,
- rung statistic,
- B-metric deflation behavior.

This is the algebraic correctness test.

### 3. Metamorphic invariance

At minimum:

- common orthogonal transform invariance,
- row-order invariance,
- common positive scaling invariance under relative ridge.

This is the family-meaning test.

### 4. Adversarial conditioning

Stress cases should include:

- imbalanced group sizes,
- near-singular but admissible group covariance,
- nearly tied leading roots,
- `p` close to per-group `n`,
- Monte Carlo budgets that are small but valid.

This is the numerical robustness test.

## Planned test file layout

The future executable test file should be:

```r
tests/testthat/test-family-two-group-cov-contrast.R
```

Recommended block structure:

```r
test_that("operator builder returns an admissible SPD metric after relative ridge", { ... })

test_that("group-label permutation null is approximately calibrated at rung 1", { ... })

test_that("observed generalized roots match the direct whitening oracle", { ... })

test_that("B-metric deflation preserves the planted remaining generalized roots", { ... })

test_that("family is invariant to common orthogonal feature transforms", { ... })

test_that("family is invariant to row reordering with aligned labels", { ... })

test_that("family is invariant to common positive scaling under relative ridge", { ... })

test_that("family remains finite for near-singular but admissible covariance pairs", { ... })

test_that("family refuses non-exchangeable or structurally incoherent requests", { ... })
```

The test file should be understood as part of the family's **admission
certificate**, not as ordinary coverage.

## Shipping rule

Do not ship a generic two-covariance generalized-eigen adapter.

Ship this family only if:

- the family spec is accepted,
- the null doctrine is explicit,
- the four test blocks exist and pass,
- the admissibility checks fail loudly on incoherent requests.

If those conditions are not met, the family should remain a note and not a
public surface.
