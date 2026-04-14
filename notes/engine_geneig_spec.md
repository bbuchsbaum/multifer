# Generalized-Eigen Engine Spec

This note is the design-note companion for Sprint 4 (`multifer-093`). It
defines the `geneig` geometry and the `run_geneig_ladder()` engine.

Companion notes:
- Package-level doctrine: [package_vision.md](./package_vision.md)
- Architecture rules: [architecture_rules.md](./architecture_rules.md)
- Latent-root engine reference: [paper1_rank_matched_residual_randomization.md](./paper1_rank_matched_residual_randomization.md)

`geneig` is the next latent-root engine. It extends the existing
collapsed Vitale P3 ladder to generalized eigenvalue problems of the
form

```
A v = lambda B v,      with B symmetric positive definite.
```

Unlike `predictive`, which changes the inferential target, `geneig`
keeps the latent-root estimand intact and only changes the inner
mechanics: deflation happens in the **B-metric**, not the Euclidean
metric. This is the non-negotiable math constraint that defines the
engine.

## Scope

The v1 flagship is **LDA** via `MASS::lda`, which is the cleanest
special case: `A = B_between`, `B = W_within`, labels permute under
the null, discriminant rank is bounded by `K - 1`.

Also supported as future extensions:

- **Contrastive PCA** (`A = Sigma_target`, `B = I + lambda Sigma_bg`
  or equivalent). Null is a sign flip on the background contribution,
  which destroys the contrast while preserving both marginals.
- **Metric-weighted PCA** (`A = X^T M X`, `B = M`). Null is row
  permutation or metric randomization depending on what the metric
  represents.

Out of scope for v1: non-positive-definite `B`, indefinite pencils,
Riemannian / manifold variants. The engine must refuse these at
construction time.

## Six-part admissibility

### 1. Ordered estimand

The estimand is the **ordered sequence of generalized eigenvalues**:

```
lambda_1 >= lambda_2 >= ... >= lambda_r >= 0
```

where `lambda_a` solves `A v_a = lambda_a B v_a` and the eigenvectors
are `B`-orthogonal: `v_a^T B v_b = delta_{ab}`.

Operationally the engine computes these via the whitened operator

```
C = B^{-1/2} A B^{-1/2}
```

which has the same eigenvalues as the original pencil. The
eigendecomposition happens in the whitened space; the ordered
eigenvalues are the estimand.

The rung statistic is the **relative tail-ratio**, identical in form to
the oneblock / cross engines:

```
stat_k(lambda) = lambda_k^2 / sum_{q >= k} lambda_q^2
```

This is the Vitale P3 collapse applied to generalized roots. Because
the form is the same as for PCA and PLSC, the existing
`mc_sequential.R` and `ladder_driver` machinery apply without
modification.

### 2. Deflation — the hard rule

**Deflation is B-orthogonal, not Euclidean.**

After step `k`, the first `k` generalized eigenvectors `v_1, ..., v_k`
are removed from the pencil by projecting `A` onto the B-orthogonal
complement of their span:

```
P_B = I - V_k (V_k^T B V_k)^{-1} V_k^T B = I - V_k V_k^T B
```

(since the `v_a` are B-orthonormal, `V_k^T B V_k = I`). The deflated
pencil is `(P_B^T A P_B, B)`, which has generalized eigenvalues
`lambda_{k+1}, lambda_{k+2}, ...`.

Equivalently, in the whitened space: deflate `C = B^{-1/2} A B^{-1/2}`
by Euclidean projection onto the orthogonal complement of
`B^{1/2} V_k`. This is mathematically identical and is the
implementation path the engine uses, because it lets the whitened
operator be deflated with standard Euclidean projectors after the
initial whitening.

**What goes wrong if you deflate in the Euclidean metric instead:**
the deflated `A` is no longer the correct operator for the remaining
generalized eigenvalues. The tail-ratio statistic then tests the wrong
sequence, and the test is silently invalid. This is the specific
failure mode the reviewer warned about and the reason the gate (§5)
fires at registration time.

### 3. Null action — varies by sub-family

Unlike the cross engines where the null is structurally shared
("permute rows of Y"), the generalized-eigen family has
sub-family-specific nulls:

- **LDA**: **label permutation**. Permuting class labels destroys the
  discriminant structure while preserving both the within-class
  covariance and the marginal distribution of `X`. This is the only
  v1 LDA null.
- **Contrastive PCA**: **sign flip** on the background term, or
  permutation of which samples are labelled "background" vs "target".
  Out of v1 scope but flagged.
- **Metric-weighted PCA**: **metric randomization** (permute the
  metric's row indices) or **row permutation** (standard PCA null)
  depending on which aspect of the problem is being tested.

The engine therefore **delegates null action to the adapter**, not to
the engine. `run_geneig_ladder` takes a `null_action` hook from the
recipe's adapter and calls it per null draw. This keeps the engine
sub-family agnostic while letting each adapter declare its own null.

The contract rule: a geneig adapter that cannot describe a valid null
action for its sub-family must not declare `component_significance`.
This is enforced by the existing capability gate; no new rule is
needed.

### 4. Rung statistic

Relative tail-ratio on generalized roots, as in §1. The form is
identical to the existing oneblock ladder, which lets the engine reuse
`ladder_driver` without modification.

The observed statistic at step `k` is computed on the deflated whitened
operator `C_k`; the null draws apply the adapter's `null_action`,
re-whiten, re-eigendecompose, and compute the same statistic.

### 5. Adapter contract

Geneig adds **no new hooks**, but it adds **one new gate rule** and one
shape-constructor requirement.

**Shape requirement (`multifer-093.2`):** `geometry("geneig", ...)`
takes an operator pair `(A, B)` with a declared metric. Construction
time validates:

- `A` is square symmetric (tolerance `max(|A|) * sqrt(eps)`).
- `B` is square symmetric positive definite. Checked via Cholesky with
  specific error message on failure: "B must be symmetric positive
  definite; Cholesky failed at pivot ...".
- `A` and `B` have matching dimensions.

If any check fails, construction errors immediately. `typed_shape`
composes the geneig geometry with `relation("generalized_eigen")` and
errors if any other relation is requested.

**Gate rule (`multifer-093.4`):** an adapter declaring
`(geneig, generalized_eigen, component_significance)` must either:

- provide a `residualize` hook annotated `b_metric = TRUE` via a new
  attribute on the function, **or**
- declare `residualize = NULL` and delegate to the engine's own
  B-metric deflation (the common case for refit-only adapters).

The gate refuses any adapter whose `residualize` is present, not
annotated `b_metric`, **and** claims `(geneig, generalized_eigen,
component_significance)`. The error message points at this spec note
by filename.

This is a registration-time check, not a runtime check. The point is
to make it structurally impossible for a downstream package to
accidentally ship Euclidean deflation for a geneig problem.

### 6. Fast-path plan

**Refit-only in v1 for `adapter_lda_refit`.** The ladder is naturally
short (LDA discriminant rank is bounded by `K - 1` classes), so the
speed gain from a fast path is modest and not worth the
implementation complexity for v1.

Future work: a warm-start path exploiting the B-whitening's reuse
across null draws when `B` is fixed (which it is under label
permutation — permuting labels changes `A` but not `W_within = B`).
This is the cleanest fast-path opportunity the engine has and should
be considered before any RSpectra-based approximation.

## Engine outline (`run_geneig_ladder`)

```
ladder_driver
  + mc_sequential
  + form_units
  + same Besag-Clifford batch schedule
```

Initial whitening happens once at the top of the ladder and is cached
by the `step_cache` environment pattern already used by the cross
engine. Per-step callbacks:

```r
observed_stat_fn(step, data) ->
    C_step = whiten(data$A_deflated, data$B)
    lambda = eigen(C_step, symmetric = TRUE)$values
    tail_ratio(lambda, step)

null_stat_fn(step, data) ->
    data_null = adapter$null_action(data)
    C_null = whiten(data_null$A_deflated, data_null$B)
    lambda_null = eigen(C_null, symmetric = TRUE)$values
    tail_ratio(lambda_null, step)

deflate_fn(step, data) ->
    # In whitened space, Euclidean projection onto the complement of
    # the top step eigenvector of C. Then lift back to the original
    # pencil's A_deflated while keeping B fixed.
    ...
```

The engine exposes a small helper `.geneig_whiten(A, B)` that takes
the pencil and returns the whitened operator plus the transform matrix
`B^{-1/2}`. Numerical stability uses the eigendecomposition of `B`
(not Cholesky) so that rank-deficient `B` fails loudly rather than
silently producing garbage.

## B = I regression (required test)

The engine must pass a **bit-for-bit regression test** against
`run_oneblock_ladder`:

- Take any oneblock PCA fixture `X` with covariance `A = X^T X / n`.
- Run `run_oneblock_ladder` on `X`.
- Run `run_geneig_ladder` on the pencil `(A, I)` with the same seed.
- Assert identical p-values, identical roots_observed, identical units.

This is not a "close enough" test. It must be exact. If it is not
exact, something is wrong either in the whitening or in the deflation
or in the tail-ratio computation, and the engine must not ship.

## Euclidean vs B-metric regression (required test)

A second regression test constructs a pencil where Euclidean deflation
would produce the wrong remaining eigenvalues. Specifically, pick
`B != I` with non-trivial off-diagonal structure, compute the true
generalized eigenvalues, and verify that:

- `run_geneig_ladder` recovers the true values.
- A hand-implemented Euclidean-deflation variant does **not**.

The gap between the two is what the gate is preventing. This test
proves the gate does what it claims.

## Interaction with existing engines

- Shares `ladder_driver`, `mc_sequential`, `form_units`, and the
  Besag-Clifford batch scheduler with all other engines.
- Does **not** share `cached_svd` because the inner math is
  eigendecomposition, not SVD. Adds a small `cached_eigen_sym` helper
  mirroring `cached_svd` for symmetric matrices.
- Does not touch the cross or oneblock engines. Adding geneig does not
  require modifying any existing engine file.

## What the spec explicitly does not promise

- No contrastive PCA adapter in v1.
- No metric-weighted PCA adapter in v1.
- No fast-path beyond refit in v1.
- No support for indefinite or rank-deficient `B`.
- No variable_significance.

## Open questions (decide during implementation, not here)

1. For LDA specifically, should the engine use `MASS::lda`'s native
   discriminant computation or whiten the scatter pencil directly?
   Current spec says "adapter_lda_refit delegates fitting to
   MASS::lda", which side-steps the question for the flagship but
   leaves it open for later geneig adapters.
2. Should the engine expose a `metric = "euclidean" | "b_metric"`
   knob so contrastive-PCA-style operator pairs with Euclidean-metric
   deflation can opt in? Current spec says no — the metric is fixed
   per-sub-family and lives on the adapter. Adding a knob would weaken
   the gate.
3. How should ties in the generalized spectrum be handled by
   `form_units`? The existing `group_near_ties` path applies
   unchanged; no geneig-specific work is needed.

## References

- Reviewer critique (2026-04-14 conversation): the B-metric deflation
  requirement and the registration-time gate.
- `notes/package_vision.md` §"Roadmap to paper-quality v1" point 4:
  "Build a real `geneig` engine."
- `notes/architecture_rules.md` Rule 1: "Prefer a small number of
  engine families" (geneig as direct-latent-root engine, not a
  oneblock branch).
- `notes/architecture_rules.md` Rule 6: "Prefer explicit abstractions
  over branching" (the B-metric gate rather than a method-specific
  branch in `run_oneblock_ladder`).
