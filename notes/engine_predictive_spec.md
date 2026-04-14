# Predictive-Gain Engine Spec

This note is the design-note companion for Sprint 3 (`multifer-d3o`). It
defines the new `predictive` relation and the `run_predictive_ladder()`
engine that implements it.

Companion notes:
- Package-level doctrine: [package_vision.md](./package_vision.md)
- Architecture rules: [architecture_rules.md](./architecture_rules.md)
- Latent-root engine reference: [paper1_rank_matched_residual_randomization.md](./paper1_rank_matched_residual_randomization.md)

This spec exists because PLS regression is **not** a latent-root model
in the sense that PCA, PLSC, CCA, and `geneig` are. Forcing it into the
existing `covariance` relation would corrupt the meaning of that
relation and would silently legitimize in-sample root tests as
significance. `architecture_rules.md` Rule 5 explicitly allows the inner
target to vary when it must; this is where it must.

## Scope

The engine must handle the supervised two-block case where a predictor
matrix `X` (n × p) is used to construct successive latent predictors of
a response matrix `Y` (n × q). The inferential question is not "does
the a-th cross-covariance root exceed noise" but **"does the a-th
latent predictor add held-out predictive gain beyond the first a - 1"**.

The v1 flagship adapter is `adapter_plsr_refit` wrapping `pls::plsr`.
Out of scope for v1: reduced-rank regression, canonical ridge
regression, non-linear kernelized variants. They may share the same
engine later via new adapters, but the engine must be designed so none
of them require branching in the engine.

## Six-part admissibility (the reviewer's test)

### 1. Ordered estimand

The estimand is the ordered sequence of **incremental predictive gains**
for latent predictors `1, 2, ..., a_max`:

```
g_a = f(Yhat_a) - f(Yhat_{a-1})
```

where `Yhat_a` is the held-out prediction of `Y` using the first `a`
latent components and `f(.)` is a scalar predictive fit functional.
v1 uses `f = 1 - PRESS / TSS` (cross-fitted R² against centered `Y`),
equivalently the held-out covariance of `Yhat_a` with the observed
`Y_oos` rescaled to a 0–1 scale. Negative increments are clipped at 0
at the test stage but not during fitting.

This estimand is **not** a root of `X^T Y` or any derived singular-value
expression. The difference matters under misspecification: a component
that aligns with a large-variance but noise-dominated direction of `Y`
can have a large `X^T Y` singular value and negligible held-out gain.
The predictive relation tests the thing the user cares about, not its
SVD proxy.

### 2. Deflation

Deflation removes the a-th latent predictor from **the predictor
construction**, not from `X^T Y`. Concretely, after fitting `a`
components, the (`a + 1`)-th rung is evaluated on the residuals of `Y`
against `Yhat_a` and on the `X`-deflation produced by the fitter
(NIPALS deflation in the flagship adapter).

This is expressed on the adapter side via `truncate(fit, a)` plus
`residualize(fit, a, data)`. The engine does not reimplement PLSR
deflation; it delegates to the adapter. Different predictive adapters
may use different deflation schemes (e.g. RRR's rank-constrained
regression vs NIPALS vs SIMPLS) and the engine is agnostic.

### 3. Null action

Default: **permute rows of `Y` relative to `X`**. This destroys the
supervised relationship while preserving both marginals, which is the
standard exchangeability null for paired-row regression problems.

For nuisance-adjusted designs, the null permutes residuals of `Y`
against a nuisance design `Z`, not raw `Y`. The engine routes through
`nuisance_residual_basis()` from `engine_cross.R` to reuse the same
exchangeability-preserving basis construction the CCA path already
uses; this keeps the predictive and correlation relations sharing
machinery rather than reimplementing it.

Blocked / grouped designs use within-group permutation of `Y` rows,
reusing `.restricted_row_permutation()`.

### 4. Rung statistic

**Cross-fitted held-out gain `g_a`** as defined in §1 above, computed
via K-fold cross-validation (v1 default: 5 folds, user-overridable).

The rung statistic at step `a` is the *incremental* gain, not the
total fit. That is the statistic the ladder tests against null draws.
This is the single most important design choice in the spec and the
source of the reviewer's hard rule.

> **Hard rule (registration-time gate):** no predictive adapter may
> claim `component_significance` unless `component_stat(fit, data, k,
> split)` takes a `split` argument and uses it to produce an
> out-of-sample statistic. Adapters that compute `component_stat` on
> the same data they fit on are refused at `infer_adapter()` time with
> a named reason.

In-sample PLSR significance is too easy to fool: adding a component
always increases in-sample fit monotonically, so any ladder built on an
in-sample statistic will over-reject. The gate exists to prevent that.

### 5. Adapter contract

The predictive family introduces **one new hook** on the adapter
contract:

```r
predict_response(fit, new_data, k = NULL)
```

Returns an `n × q` matrix of fitted responses for the first `k`
components (or all components if `k = NULL`). This is required when
the adapter claims `(cross, predictive, component_significance)`.

`component_stat` grows an optional fourth argument `split`:

```r
component_stat(fit, data, k, split = NULL)
```

where `split` is a list with `$train` and `$test` index vectors. The
engine fills `split` via its cross-fit driver; adapters use it to
route training / held-out data through `refit` and `predict_response`.

The capability gate extends with these rules:

- `component_significance` for `relation = "predictive"` requires
  `refit`, `truncate` (or `residualize`), `null_action`,
  `predict_response`, **and** a `component_stat` whose formals contain
  `split`. The last is the cross-fit discipline rule and is a hard
  registration error when missing.
- `variable_stability` / `score_stability` / `subspace_stability`
  inherit the existing latent-root requirements: a perturbation path
  (refit *or* core+update_core) and `loadings` / `scores` /
  `loadings`.
- `variable_significance` remains excluded per `adapter_contract.R`
  §38.

The `predict_response` addition is optional on the contract at large
(other relations do not use it) and required only when the predictive
gate fires.

### 6. Fast-path plan

**Refit-only in v1.** No `core` / `update_core` hook for predictive
adapters.

Justification: the predictive rung statistic is held-out, which means
every null draw requires a full cross-fit loop, which means the
innermost hot path already includes refits on a non-trivial fraction
of the data. The core-space fast path that works for covariance-mode
cross bootstrap relies on a small `k × k` update matrix; no equivalent
identity exists for the cross-fitted predictive-gain statistic because
the held-out fold structure breaks the block-SVD factorization the
covariance fast path exploits.

Future work may add a warm-start fast path for NIPALS specifically
(exploit the sequential deflation structure to avoid re-running
components 1…a on every fold), but that is a speed optimization, not a
new engine family. v1 ships refit-only and says so in the adapter's
capability listing.

## Engine outline (`run_predictive_ladder`)

The engine reuses the shared outer shell:

```
ladder_driver
  + mc_sequential
  + form_units
  + same Besag-Clifford batch schedule
```

What varies is the step callbacks:

```r
observed_stat_fn(step, data) ->
    g_step = cross_fit_gain(data, k = step) -
             cross_fit_gain(data, k = step - 1)

null_stat_fn(step, data) ->
    perm = sample_y_rows(data)         # or restricted perm for blocked
    g_step_null = cross_fit_gain(list(X = data$X, Y = data$Y[perm,]), k = step) -
                  cross_fit_gain(list(X = data$X, Y = data$Y[perm,]), k = step - 1)

deflate_fn(step, data) ->
    residualize(fit, step, data)        # adapter-supplied
```

`cross_fit_gain()` is a new helper in `engine_predictive.R` that takes
(data, k, folds) and returns the scalar held-out gain, delegating
model fits to the adapter's `refit` + `predict_response`. The fold
construction is deterministic given a seed so the ladder is
reproducible.

The ladder produces exactly the same `infer_result` schema as the
other engines — `component_tests`, `units`, `roots_observed` (filled
with the observed `g_a` sequence), plus any stability outputs.
`roots_observed` is not singular values here; the schema slot accepts
whatever ordered estimand the engine produces, with a label
distinguishing "held-out gain" from "latent root". That label is
exposed in `print.infer_result` so users are never misled about what
the numbers mean.

## Interaction with existing engines

- `run_predictive_ladder` shares `ladder_driver`, `mc_sequential`,
  `form_units`, and the nuisance-basis helpers with `run_cross_ladder`.
  No duplicated scaffold code.
- The engine does not use `cached_svd` or `top_singular_values` because
  it does not compute singular values. This is by design.
- Stability outputs (`variable_stability`, `score_stability`,
  `subspace_stability`) are produced by the same bootstrap machinery
  that already exists: the bootstrap loop calls `refit` on resamples,
  and the stability consumers run off `loadings` / `scores`. The
  meaning of "variable stability" for a predictive PLSR fit is the
  same as for a PLSC fit: how often does this variable appear in a
  well-supported latent direction across bootstraps. That meaning is
  still honest even though the significance target changed.

## What the spec explicitly does not promise

- No variable_significance.
- No fast path (core / update_core) in v1.
- No multiblock predictive variant.
- No kernelized / non-linear predictive variants.
- No claim that held-out gain matches a cross-covariance root in any
  particular regime. The predictive relation is a different estimand,
  full stop.

## Open questions (decide during implementation, not here)

1. Should the default cross-fit use K-fold or leave-one-out? K-fold is
   cheaper and preferred for the flagship. Leave-one-out may be needed
   for small-n fixtures.
2. Should `roots_observed` be cumulative held-out R² or incremental
   gains? Current spec says incremental because the ladder tests
   increments. The public print layer can show both.
3. Should the engine support a user-supplied `folds` vector for
   repeated measurements or longitudinal designs? Almost certainly
   yes, but the v1 default is a deterministic K-fold on row indices
   and the user-supplied path is a direct pass-through.

These are flagged here so they are not litigated at PR time.

## References

- Reviewer critique (2026-04-14 conversation): the six-part
  admissibility test and the in-sample statistic refusal rule.
- `notes/package_vision.md` §"Roadmap to paper-quality v1" point 3:
  "Create a new supervised relation for PLS regression / reduced-rank
  regression."
- `notes/architecture_rules.md` Rule 5: "Allow the inner target to
  vary when needed."
- `notes/architecture_rules.md` Rule 7: "The package should refuse
  incoherent requests."
