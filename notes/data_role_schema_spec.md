# Data Role Schema — Generalizing Structured-Payload Adapters

This note specifies a **data role schema** for adapters whose data
payload is structurally richer than a numeric matrix or paired list.
Constrained PCA (`fastcpca`) is the acceptance test. The schema makes
the package handle row-aligned auxiliary structure (constraints,
metrics, weights, panel indicators, group labels) through a small
declarative grammar instead of forcing every such adapter to ship its
own `bootstrap_action`.

Companion notes:
[architecture_rules.md](./architecture_rules.md),
[package_vision.md](./package_vision.md),
[pcr_consolidation_spec.md](./pcr_consolidation_spec.md).

## What the schema does and does not shrink

The schema does **not** shrink the set of methods that need
`geometry = "adapter"`. CPCA is the worked example: its inferential
object is structurally `(Z, row, col, K, L, block)` — not any built-in
matrix or list shape. `geometry = "adapter"` is the correct geometry
for it before and after this work.

What the schema **does** shrink is the set of adapters that must
implement `bootstrap_action`. Today, any adapter with row-aligned
auxiliary structure has to write its own `bootstrap_action` because
the default resampler does not know what to co-subset. With a role
schema, the default resampler plans correctly for the safe cases and
refuses to plan for the unsafe ones. `bootstrap_action` becomes an
override for genuinely irregular resampling rules rather than a
must-implement for every structured payload.

The boundary phrasing carried into `package_vision.md`:

> **capability triple × data role schema × null spec, with
> `geometry = "adapter"` reserved for irreducible latent/data
> geometry, not merely missing taxonomy.**

## The schema vocabulary

A `data_role_schema` declares the named roles in an adapter's data
payload. Each role carries:

- **name** — character key matching what the adapter unpacks from `data`
- **kind** — one of: `"primary"`, `"constraint"`, `"metric"`,
  `"design_index"`, `"static"`
- **axes** — `c("row")`, `c("col")`, `c("row", "col")`, or
  `character(0)`
- **policy** — applies to `"metric"` roles where multiple structures
  are accepted at runtime (see below)

The role kinds are deliberately small:

- **primary** — the data matrix or matrices being decomposed
  (`Z` in CPCA; `X` and `Y` in cross fits).
- **constraint** — row- or column-aligned auxiliary subspace
  (`G`, `H` in CPCA; nuisance covariates in CCA).
- **metric** — row- or column-aligned definite quadratic form
  (`K`, `L` in CPCA; observation weights in weighted PCA).
- **design_index** — row- or column-aligned grouping or label vector
  used for nuisance / blocked / paired designs.
- **static** — carried but not row/col-aligned (hyperparameters,
  scalars, fitted offsets).

Anything that does not fit one of these is a candidate for a future
role kind, not for ad hoc handling. The grammar is intentionally
finite.

## Runtime resolution

The schema is **both declarative and runtime-resolved**:

- The adapter declares accepted role names, kinds, and axes.
- `compile_bootstrap_plan()` resolves the **actual structure** from
  the data the user supplied.

The same adapter can accept `K = NULL` (identity), `K = numeric vector`
(diagonal), or `K = matrix` (full SPD), and the legality of the
default bootstrap depends on which was supplied. The adapter does not
pre-commit to one form. It declares "I accept `K` as a row metric";
the runtime resolves which kind of `K` it actually got and decides
whether the requested bootstrap is plannable.

## Metric policies under bootstrap

A row metric `K` under row-bootstrap-with-replacement requires a
policy declaration. The accepted values:

- **`identity`** — `K = NULL` or `K = I`. No subsetting needed.
- **`diagonal_observation_weights`** — diagonal entries follow
  observations under resampling, so duplicates accumulate. Safe.
- **`diagonal_position_weights`** — diagonal entries follow positions,
  so the metric is fixed regardless of resampled rows. Safe.
- **`full_spd_no_replacement_only`** — full SPD `K` is not subsettable
  under replacement (duplicate indices guarantee rank deficiency, and
  Cholesky breaks; see `fastcpca/R/metric.R:54,132`). The default
  refuses to plan a with-replacement bootstrap and either errors with
  guidance or admits only without-replacement permutation.
- **`adapter_owned`** — the adapter handles `K` under
  `bootstrap_action`; the default resampler does not touch it.

The default bootstrap planner consults this policy at compile time.
If the policy is incompatible with the requested resampling rule, the
planner refuses to emit a default plan with a corrective message:

> Row metric `K` has policy `full_spd_no_replacement_only`; the
> default with-replacement bootstrap cannot subset it safely. Supply
> `bootstrap_action` or change the metric policy.

Symmetric policy declarations apply to column metrics under column
resampling, when that becomes a supported design.

## `null_spec` — descriptive, not categorical

`null_action` keeps its current operational signature. It gets a
companion declarative descriptor `null_spec`. The descriptor is
documentation and provenance metadata, not a dispatch key.

```r
null_spec = null_spec(
  kind            = "between_role_decoupling",
  randomized_role = "Z",
  reference_role  = "G",
  axis            = "row",
  preserves       = c("Z_row_marginals", "G_subspace", "K_policy")
)
```

Three null kinds cover the v1 surface:

- **`within_role_randomization`** — column permutation in PCA, row
  scrambling within a fixed design. Single primary role; randomization
  internal to it.
- **`between_role_decoupling`** — Y-against-X in PLSR (two primary
  roles), Z-against-`row` in CPCA design tests, treatment-against-
  covariate in nuisance-adjusted CCA. Names which role is randomized
  against which reference role.
- **`parametric_draw`** — bootstrap-of-fit-parameters for pure
  model-based nulls. Useful for future Bayesian-flavored adapters.

PCR's `null_spec` would read
`within_role_randomization` (row permutation of `Z`'s pairing with `Y`
is implementation, not categorization — semantically it is decoupling
between the X and Y primary roles, which is `between_role_decoupling`).
The point is that the descriptor is rich enough to disambiguate
without prescribing.

`infer_result$mc$null_label` already exists; `null_spec` is the
structured object that backs the label. It surfaces in provenance and
in the validity layer so downstream consumers can inspect what null
was run without having to reverse-engineer it from the function body.

## Acceptance test — CPCA schema

The schema is "done" when `adapter_cpca` declares this and the
framework behaves as specified.

```r
cpca_schema <- data_role_schema(
  Z   = role("primary",    axes = c("row", "col")),
  row = role("constraint", axis = "row"),
  col = role("constraint", axis = "col"),
  K   = role("metric",     axis = "row",
             policy = "identity"),   # one of:
                                     # identity,
                                     # diagonal_observation_weights,
                                     # diagonal_position_weights,
                                     # full_spd_no_replacement_only,
                                     # adapter_owned
  L   = role("metric",     axis = "col",
             policy = "identity")    # symmetric set
)
```

Required behavior:

1. **Default bootstrap plans only when every row-aligned role has a
   legal replacement rule under the requested design.**
2. For `K ∈ {identity, diagonal_observation_weights,
   diagonal_position_weights}`: default is row-resample `Z` and all
   row-aligned auxiliaries (`row`, `K`) with replacement, with `K`
   updating per its policy. `col`, `H`, `L` are fixed under row
   bootstrap.
3. For `K = full_spd_no_replacement_only`: the planner refuses with
   the message above. Without-replacement designs (e.g. permutation
   under a paired or blocked design) remain plannable.
4. For `K = adapter_owned`: planner refuses to plan and defers
   entirely to `bootstrap_action`.
5. Symmetric `L` policies apply when column resampling is requested.
6. The CPCA design test against row constraints declares
   `null_spec(kind = "between_role_decoupling", randomized_role = "Z",
   reference_role = "row", axis = "row", preserves = ...)`. The
   adapter writes the actual permutation in `null_action`; the
   `null_spec` surfaces in `infer_result$mc` and the validity block.

If any of these requires a CPCA-specific branch in
`compile_bootstrap_plan` or `infer()`, the schema has not done its
job.

## Implementation order

1. **`R/data_role_schema.R`** — constructors `data_role_schema()`,
   `role()`, `null_spec()`, plus validators. ~150 lines.
2. **Extend `infer_adapter()`** to accept and record `data_schema`
   and `null_spec`. Backwards-compatible: schema is optional;
   adapters without one keep current behavior.
3. **Extend `compile_infer_plan()`** to resolve the schema against
   supplied data and store the resolved structure on the plan.
4. **Extend `compile_bootstrap_plan()`** to consult the resolved
   schema for default planning. The `geom_kind` switch at
   `R/bootstrap.R:381,412,446` becomes schema-driven dispatch; the
   existing oneblock / cross / multiblock branches reduce to special
   schema templates registered at package load.
5. **Surface `null_spec` in `infer_result`** — extend `infer_mc()`
   provenance fields to carry the structured descriptor alongside the
   existing string label.
6. **Update the contract harness** `check_infer_adapter()` to verify
   schema declarations against the required hooks (e.g. an adapter
   declaring `K = role("metric", ..., policy = "adapter_owned")` must
   supply `bootstrap_action`).

Items 1–4 are the bulk of the work. Items 5–6 are documentation /
diagnostic upgrades that ride on top.

## What this does not do

- **It does not eliminate `geometry = "adapter"`.** CPCA, kernel
  methods with irreducible payloads, and parametric-draw fitters keep
  that geometry before and after this work. The schema makes living
  there more pleasant; it does not relocate the residents.
- **It does not replace `capability_matrix`.** The triple
  `(geometry, relation, target)` is unchanged. The schema is
  orthogonal to it.
- **It does not add new built-in geometries.**
  Multiblock-with-row-aligned-`K` is not a new geometry; it is
  `geometry = "adapter"` with a richer schema. Resist the urge to
  enumerate.

## Out of scope for this spec

- **Column-resampling designs.** The framework currently row-samples.
  Column-resampling designs (variable-importance bootstraps, frequency-
  domain methods) need a symmetric extension. Not in v1.
- **Cross-axis aligned roles.** Tensor methods with three or more axes
  need a richer schema. The current grammar covers `row`, `col`, and
  static; multi-mode tensors are a future generalization.
- **Whether to ship a public CPCA wrapper.** As with PCR, this is a
  roadmap question, not an architecture question. The schema proves
  the design admits CPCA; the decision to ship `infer_cpca()` is
  independent.

## One-line doctrine

The schema is the smallest declarative grammar that lets the engine
plan correctly for structured payloads without the engine knowing
which method authored them. CPCA earns the schema because no smaller
mechanism handles it without a special case; PCA and PLSR earn the
schema because their existing geometry switches *are* schema templates
in disguise, which the consolidation can now make explicit.
