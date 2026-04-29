# PCR as the Consolidation Acceptance Spec

This note specifies what `adapter_pcr_refit` should look like once the
adapter / engine consolidation lands, what currently blocks it, and what
"consolidation done" therefore means in terms the codebase can verify.

It is intentionally not an implementation plan for shipping PCR. PCR is
chosen as the acceptance criterion because it is the smallest concrete
method that the current taxonomy cannot express cleanly. Whether the
package actually ships a public PCR surface is a separate roadmap call.

Companion notes: [architecture_rules.md](./architecture_rules.md) and
[package_vision.md](./package_vision.md).

## Status as of commits e3f2034 / 95ac053

A first wave of the consolidation has landed since this spec was drafted.
The acceptance criteria below are revised to reflect current state:

**Already complete:**

- `R/infer_plan.R::compile_infer_plan()` exists. It is the executable
  plan compiler the consolidation called for. The plan resolves
  `component_stat` once at compile time, currying `split = NULL` for the
  predictive relation, so engines see a uniform
  `component_stat(fit, data, k)` callback.
- `R/engine_adapter.R::run_adapter_ladder()` now consumes the plan
  exclusively (`plan$fit`, `plan$roots`, `plan$component_stat`,
  `plan$null_action`, `plan$residualize`, `plan$validate_data`). No
  `rel_kind` branches survive inside the engine.
- `R/adapter_check.R::check_infer_adapter()` exists with a 290-line
  contract harness and a 120-line test suite. This is the parity oracle
  the spec asked for first.
- The oneblock ladder has been compiled through plan callbacks
  (`compile_oneblock_ladder_plan()`, `R/infer_plan.R:77–150`).
  `R/engine_oneblock.R` shrunk by ~150 lines.

**Still pending:**

- The cross ladder is not yet plan-compiled.
- `infer()` still dispatches on `geom_kind` through a 7-way branch
  (`R/infer.R:192–294`).
- The bootstrap default still switches on `geom_kind`
  (`R/bootstrap.R::.default_bootstrap_data`).
- `.infer_original_fit()` still special-cases cross (`R/infer.R:504–512`).
- The capability matrix and recipe compiler still reject
  `(oneblock, predictive)`.
- `attr(residualize, "b_metric")` and `component_engine = "adapter"`
  are still in the public adapter contract.

PCR remains blocked by the pending items, primarily the capability /
recipe restriction on `(oneblock, predictive)` and the geometry-keyed
bootstrap default.

## Why PCR is the right forcing function

Principal Components Regression sits next to `pls::plsr()` in the same
CRAN package, but it stresses the design along an axis the existing
adapters do not:

- **Latent decomposition is one-block.** The basis is PCA on X.
- **Inferential target is predictive gain.** Significance is held-out
  incremental R² of Y given the leading k principal components of X.
- **Data shape is paired.** X and Y come together with row-aligned
  observations.

PCR therefore decomposes into three independent traits:

| Trait          | Value          |
| -------------- | -------------- |
| latent kind    | oneblock       |
| target family  | predictive gain|
| data shape     | cross-paired   |

The current `(geometry, relation)` enumeration cannot express this
combination. Today, predictive is locked to `geometry ∈ {cross, adapter}`
(see `R/capability_matrix.R:94–101` and `R/infer_recipe.R:104–110`),
which conflates *latent kind* (oneblock vs cross) with *data shape*
(numeric matrix vs paired list).

If the consolidation gets `data_shape × latent_kind × target_family`
truly orthogonal — even if computed at compile time rather than declared
— PCR drops in. If it does not, the only way to ship PCR is through the
`geometry = "adapter"` wildcard, which is exactly the escape hatch the
consolidation is meant to retire.

## Sketched adapter under the consolidated design

The target shape, in roughly the form a downstream package would write:

```r
adapter_pcr_refit <- function(adapter_id = "pcr_refit",
                              adapter_version = "0.0.1",
                              ncomp = NULL) {

  fit_pcr <- function(new_data) {
    Xc <- scale(new_data$X, center = TRUE, scale = FALSE)
    Yc <- scale(new_data$Y, center = TRUE, scale = FALSE)
    sv <- svd(Xc)
    k  <- min(if (is.null(ncomp)) min(nrow(Xc) - 1L, ncol(Xc), 50L)
              else as.integer(ncomp), length(sv$d))
    T  <- sweep(sv$u[, seq_len(k), drop = FALSE], 2L, sv$d[seq_len(k)], `*`)
    list(
      Xmeans       = attr(Xc, "scaled:center"),
      Ymeans       = attr(Yc, "scaled:center"),
      loadings     = sv$v[, seq_len(k), drop = FALSE],
      scores       = T,
      sdev         = sv$d[seq_len(k)],
      coefficients = solve(crossprod(T), crossprod(T, Yc)),
      ncomp        = k
    )
  }

  infer_adapter(
    adapter_id = adapter_id,
    adapter_version = adapter_version,
    capabilities = capability_matrix(
      list(geometry = "oneblock", relation = "predictive",
           targets  = c("component_significance",
                        "variable_stability",
                        "score_stability",
                        "subspace_stability"))
    ),

    roots       = function(x, ...) x$sdev^2,
    scores      = function(x, domain = "X", ...) x$scores,
    loadings    = function(x, domain = "X", ...) x$loadings,
    refit       = function(x, new_data, ...) fit_pcr(new_data),

    residualize = function(x, k, data, ...) {
      kk <- min(max(as.integer(k), 0L), x$ncomp)
      Xc <- sweep(data$X, 2L, x$Xmeans, "-")
      if (kk > 0L) {
        T <- x$scores[, seq_len(kk), drop = FALSE]
        P <- x$loadings[, seq_len(kk), drop = FALSE]
        Xc <- Xc - T %*% t(P)
      }
      list(X = Xc, Y = data$Y)
    },

    null_action = function(x, data, ...) {
      perm <- sample.int(nrow(data$Y))
      list(X = data$X, Y = data$Y[perm, , drop = FALSE])
    },

    component_stat = function(x, data, k, split = NULL, ...) {
      .pcr_split_incremental_gain(data, split, k, fit_pcr)
    },

    predict_response = function(x, new_data, k = NULL, ...) {
      kk <- if (is.null(k)) x$ncomp
            else min(max(as.integer(k), 0L), x$ncomp)
      Xc <- sweep(new_data$X, 2L, x$Xmeans, "-")
      if (kk == 0L) {
        return(matrix(x$Ymeans, nrow(new_data$X),
                      length(x$Ymeans), byrow = TRUE))
      }
      Tn <- Xc %*% x$loadings[, seq_len(kk), drop = FALSE]
      Tn %*% x$coefficients[seq_len(kk), , drop = FALSE] +
        matrix(x$Ymeans, nrow(Tn), length(x$Ymeans), byrow = TRUE)
    },

    validity_level       = "conditional",
    declared_assumptions = c("paired_rows", "continuous_response"),
    checked_assumptions  = .pcr_checks()
  )
}
```

Order of magnitude: ~80 lines of adapter logic plus the standard
`.pcr_checks()` and `.pcr_split_incremental_gain()` helpers, mirroring
the structure of `R/adapter_plsr.R`. No `core` / `update_core` (PCR has
no exact core-space identity for predictive gain). No `bootstrap_action`
(default paired-row bootstrap is correct).

This is what consolidation should make possible. The adapter is doing
nothing exotic: it supplies the same hooks the PLSR adapter already
supplies, on the same paired-data shape, with a one-block residualization
instead of a NIPALS one.

## What blocks this today — concrete file references

The blockers are taxonomy, not hooks. The required hooks all already
exist in the contract. Items marked ✓ have been resolved by the
consolidation work in `e3f2034` / `95ac053`.

1. **[pending] Capability triple `(oneblock, predictive)` is rejected
   at matrix-construction time.** `R/capability_matrix.R:94–101`
   enforces `relation == "predictive"` requires `geometry ∈ {"cross",
   "adapter"}`. The PCR adapter cannot even build its capability matrix.

2. **[pending] The recipe compiler enforces the same restriction in
   three more places.** `R/infer_recipe.R:106–110`, `:126–130`,
   `:249–253` each reject `predictive + non-cross/adapter`.

3. **[pending] `run_oneblock_ladder` is variance-test only.** Even
   after `95ac053` compiled it through plan callbacks
   (`compile_oneblock_ladder_plan()` in `R/infer_plan.R:77–150`), the
   plan compiler is locked to `relation == "variance"`
   (`R/infer_plan.R:88–94`). PCR needs a predictive ladder applied to a
   one-block latent decomposition, which is not yet a code path.

4. **[pending] The oneblock data validator expects a numeric matrix.**
   `R/engine_adapter.R::.validate_oneblock_data` and `R/infer.R:218`
   reject paired list data. PCR's data is paired (X, Y) — its *latent
   decomposition* is one-block, but its *data shape* is not.

5. **[pending] Bootstrap dispatch is geometry-keyed, not shape-keyed.**
   `R/bootstrap.R::.default_bootstrap_data` switches on `geom_kind`.
   With `geometry = "oneblock"`, the default expects `data[indices, ]`
   on a matrix. PCR needs paired-row resampling on a list, which today
   is reachable only by declaring `geometry = "cross"`.

6. **[✓ resolved in e3f2034] The relation branch inside
   `run_adapter_ladder` has been lifted into compile-time currying.**
   `R/infer_plan.R:23–28` resolves `component_stat` once: predictive
   relations get a closure that absorbs `split = NULL`, all others get
   the bare hook. The engine now calls `plan$component_stat(fit, data,
   k)` uniformly with no relation switch.

## Acceptance criterion for the consolidation epic

Consolidation is "done" when, with no further design work:

1. **[pending]** The capability matrix accepts the triple
   `(oneblock, predictive, {component_significance, variable_stability,
   score_stability, subspace_stability})` without special-case branches.
2. **[pending]** `infer()` dispatches `(oneblock, predictive)` data
   shaped as paired `list(X, Y)` to the predictive ladder driver, using
   the one-block latent decomposition for ordering and deflation, with
   no `geometry == "adapter"` fallback path involved.
3. **[✓ met in e3f2034]** `run_adapter_ladder` calls a uniform
   `component_stat(fit, data, k)` callback. Any `split` argument has
   been curried at compile time. The engine source contains no
   `if (rel_kind == ...)` branches.
4. **[pending]** `bootstrap_fits()` selects a default resampler from a
   `data_shape` trait (or the equivalent compiled plan field), not
   from `geom_kind`. Paired-list data round-trips correctly under
   `geometry = "oneblock"`.
5. **[pending]** The PCR sketch above compiles and registers without
   modification beyond the standard adapter helpers (`.pcr_checks`,
   `.pcr_split_incremental_gain`, `.validate_pcr_data`).
6. **[✓ met in e3f2034]** A contract test harness
   (`check_infer_adapter()`) accepts the PCR adapter and exercises
   every claimed target on a synthetic fixture. The harness shipped at
   `R/adapter_check.R`.

If any of the pending items still requires a per-method branch when
addressed, the consolidation has not finished.

## What this implies about the consolidation work

The spec backs into an ordered work plan. Items marked ✓ are complete
in `e3f2034` / `95ac053`.

1. **[✓]** Build `check_infer_adapter()` against the existing adapters
   as the parity oracle. Shipped in `R/adapter_check.R` with
   `tests/testthat/test-adapter-contract-harness.R`.
2. **[✓]** Lift the `(rel_kind == "predictive")` branch out of
   `run_adapter_ladder` into a compile-time currying step. Done in
   `R/infer_plan.R::compile_infer_plan()`; the engine is relation-blind.
3. **[partial]** Re-express `run_oneblock_ladder` and
   `run_cross_ladder` as default plan compilations that produce
   `run_adapter_ladder` calls. Oneblock landed in `95ac053`
   (`compile_oneblock_ladder_plan()`). Cross is still pending. The
   oneblock plan compiler is also locked to `relation == "variance"`
   and would need a parallel `compile_oneblock_predictive_ladder_plan()`
   or a generalization before PCR's ladder can be expressed.
4. **[pending]** Replace the `geom_kind` switch in
   `bootstrap.R::.default_bootstrap_data` and
   `infer.R::.infer_original_fit` with a `data_shape`-keyed default,
   inferred from the adapter at compile time.
5. **[pending]** Allow `(oneblock, predictive)` and
   `(multiblock, predictive)` in the capability matrix and recipe
   compiler. The wildcard `geometry = "adapter"` should become a
   vestigial declarative slot for adapters whose data really is opaque
   (kernel evaluations, structured tensors, parametric draws).
6. **[pending]** Move `attr(residualize, "b_metric")`,
   `component_engine = "adapter"`, and the `(adapter, predictive)`
   admissibility special-case into explicit adapter fields or capability
   metadata rows.

The remaining critical path for PCR is items 3 (cross + predictive
oneblock plan compilation), 4 (data_shape-keyed bootstrap default), and
5 (capability/recipe acceptance of `(oneblock, predictive)`). Items 1,
2, and the harness aspect of 6 are done.

## Out of scope for this spec

- **Whether to ship a public PCR surface in v1.** This is a roadmap
  question, not an architecture question. The spec proves the design
  admits PCR; the decision to ship `infer_pcr()` is independent.
- **Multi-block PLS, SO-PLS, Po-PLS, CPPLS, kernel-PLS variants.**
  Each is a related test case that the same consolidation unlocks
  (`(multiblock, predictive)` for the first three, an algorithmic
  variant of the existing PLSR adapter for kernel-PLS). They are not
  spec'd here. The discipline is the same: the consolidation is done
  iff their adapters fall out as direct expressions of the trait
  combination, not as wildcard escape hatches.
- **Performance work.** No `core` / `update_core` for PCR; refit-only is
  the v1 baseline.

## One-line doctrine

PCR is small enough to write in eighty lines and large enough to expose
every place the taxonomy still pretends data shape, latent kind, and
target family are the same dimension. Make it fall out cleanly and the
consolidation has earned its keep.
