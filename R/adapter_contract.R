#' The `multifer` adapter contract
#'
#' This page is the single-source reference for third-party package
#' authors who want to plug their own fitted-model class into
#' `multifer`'s inference layer. It describes the hook signatures,
#' return-value shapes, capability gates, and validity contracts that
#' the engine expects from an `infer_adapter()` object.
#'
#' The contract is what makes `multifer` extensible: the ladder driver,
#' null-action machinery, bootstrap loop, and stability consumers all
#' dispatch through the hooks below, so anything that honours the
#' contract becomes a first-class family.
#'
#' A practical design rule follows from the contract: if your fitter is
#' genuinely single-relation, declare only that one relation. The
#' dedicated `cancor_cross` adapter is the model to copy. The dual-relation
#' `cross_svd` adapter is kept as a reference implementation for the
#' strict-dispatch ambiguity rule, not as the default pattern third-party
#' packages should mimic.
#'
#' @section Required vs optional hooks:
#'
#' Hooks are passed by name through `...` to [infer_adapter()]. The
#' capability gate (Part 5 §36) enforces which hooks are required for
#' which targets. There is no single hook required for every adapter:
#' an adapter only has to provide the hooks needed for the targets it
#' claims.
#'
#' \describe{
#'   \item{`roots(x)`}{**Strongly recommended.** Return a numeric
#'     vector of ordered non-negative latent roots from the fit. Used
#'     by [form_units()] to build the result's unit table, including
#'     auto-subspace grouping of near-tied roots. Without it, the
#'     package has no canonical ordering of the fit's latent objects.}
#'
#'   \item{`scores(x, domain = NULL)`}{Required by
#'     `score_stability`. For oneblock shapes the `domain` argument
#'     is optional. For cross shapes it must accept `"X"` or `"Y"`.
#'     Return an `n × k` numeric matrix where `n` is the number of
#'     observations used in the fit and `k` is the number of latent
#'     components.}
#'
#'   \item{`loadings(x, domain = NULL)`}{Required by `variable_stability`
#'     and `subspace_stability`. Same `domain` convention as `scores`.
#'     Return a `p × k` numeric matrix (`p` is the number of variables
#'     in the relevant block).}
#'
#'   \item{`domains(x, data = NULL)`}{Optional. Return the character vector
#'     of domain labels accepted by `scores()` and `loadings()`. Oneblock
#'     adapters default to `"X"`, cross adapters default to `c("X", "Y")`,
#'     and multiblock adapters default to the names of the block-list data
#'     (or `"block1"`, `"block2"`, ... for unnamed lists). Provide this hook
#'     when domains are model-defined rather than data-name-defined.}
#'
#'   \item{`project_scores(x, data, domain = NULL)`}{Optional. Return scores
#'     for `data` projected through fit `x` in the adapter's native
#'     preprocessing/projection space. When supplied, `bootstrap_fits()` uses
#'     this hook to store aligned scores for the original data, and
#'     `score_stability_from_bootstrap()` consumes those stored scores. Without
#'     it, score stability falls back to centered-data times aligned-loadings.}
#'
#'   \item{`truncate(x, k)`}{Optional. Return a fit of the same class
#'     truncated to the leading `k` components. Used when downstream
#'     code needs a reduced-rank view of the original fit.}
#'
#'   \item{`residualize(x, k, data)`}{Required by
#'     `component_significance`. Return the data matrix (or
#'     oneblock-or-cross list) with the leading `k` components
#'     removed. This is the deflation step that feeds the ladder's
#'     next rung.}
#'
#'   \item{`refit(x, new_data)`}{Required for `variable_stability`,
#'     `score_stability`, and `subspace_stability` unless
#'     `core + update_core` or `bootstrap_action` are provided instead.
#'     Fit a new model of the same class on `new_data` (typically a bootstrap resample).
#'     `new_data` has the same shape as the original `data` argument
#'     to [infer()] (matrix for oneblock, list with `X` and `Y` for
#'     cross). Return shape: the same class as `x`.}
#'
#'   \item{`component_engine`}{Optional character scalar. Set to `"adapter"`
#'     to route component-significance through adapter-owned `component_stat`,
#'     `null_action`, and `residualize` hooks for a geometry that otherwise has
#'     a built-in reference engine. This is for methods whose latent-root
#'     statistic, null, or deflation lives in an adapter-defined geometry such
#'     as a weighted, kernel, smoothed, sparse, or otherwise non-Euclidean
#'     space. The default is `"default"`.}
#'
#'   \item{`bootstrap_action(x, data, design, replicate = NULL)`}{Optional.
#'     Adapter-owned perturbation hook for stability targets. Return a list
#'     with `fit` (a replicate fit) or `data` (a replicate data object to pass
#'     to `refit`), plus optional `resample_indices` and `info`. This hook is
#'     for methods whose exchangeable unit is not a row of the data matrix,
#'     such as subject blocks, studies, sites, sessions, or parametric draws.
#'     When absent, `bootstrap_fits()` uses its default row/aligned-block
#'     bootstrap.}
#'
#'   \item{`core(x, data)`}{Optional, pairs with `update_core` to
#'     provide the Phase 1.5 fast path. Return a lightweight core
#'     object that `update_core()` can repeatedly transform without
#'     re-running the full fit. For oneblock SVD and cross covariance,
#'     see the reference implementations in
#'     `R/adapter_oneblock_baser.R` and `R/core_space.R`.}
#'
#'   \item{`update_core(core_obj, indices = NULL, ...)`}{Optional,
#'     the fast-path counterpart to `core`. Given the `core_obj` and
#'     a row-bootstrap `indices` vector, return a new fit of the same
#'     class as `x` (i.e. the class `refit()` would have returned on
#'     the same bootstrap). The engine prefers `update_core` over
#'     `refit` whenever both are present; this is how the
#'     exact-core-space bootstrap lands.}
#'
#'   \item{`align(...)`}{Optional. Custom bootstrap alignment between a
#'     replicate fit and the reference fit. New adapters should accept named
#'     arguments `reference_loadings`, `replicate_loadings`,
#'     `replicate_scores`, `reference_fit`, `replicate_fit`, `domain`, and
#'     `method`, and return a list with `loadings` and optionally `scores`.
#'     This lets adapters perform component matching and sign/rotation under
#'     their own inner product. Without it, bootstrap alignment uses Euclidean
#'     matching plus sign or legacy Procrustes alignment on the loading
#'     matrices.}
#'
#'   \item{`null_action(x, data)`}{Required by
#'     `component_significance`. Return one null resample of `data`
#'     (same shape as input). For oneblock variance this is typically
#'     a column-wise permutation; for cross shapes it breaks the
#'     pairing by permuting rows of `Y`.}
#'
#'   \item{`component_stat(x, data, k, split = NULL)`}{Required by
#'     `component_significance`. Return the observed test statistic
#'     for rung `k` on `data`. Must be algebraically consistent with
#'     `null_action`: the ladder compares the observed statistic to
#'     repeated draws of `component_stat(x, null_action(x, data), k)`.
#'     For `(cross, predictive)`, the hook must accept a formal
#'     `split` argument and use it for held-out scoring rather than
#'     in-sample fit.}
#'
#'   \item{`predict_response(x, new_data, k = NULL)`}{Optional at the
#'     contract level, required for `(cross, predictive,
#'     component_significance)`. Return an `n × q` matrix of fitted
#'     responses on `new_data` using the first `k` components (or all
#'     available components when `k = NULL`).}
#'
#'   \item{`variable_stat(x, data, domain, k)`}{Optional alternative
#'     to `loadings` for `variable_stability`.}
#'
#'   \item{`score_stat(x, data, domain, k)`}{Optional alternative to
#'     `scores` for `score_stability`.}
#' }
#'
#' @section Capability gate rules:
#'
#' Registration-time enforcement prevents an adapter from claiming a
#' capability without the hooks to back it. The rules are:
#'
#' - `component_significance` requires `null_action`, `component_stat`,
#'   **and** `residualize`.
#' - `(cross, predictive, component_significance)` additionally
#'   requires `refit`, `predict_response`, and a split-aware
#'   `component_stat(..., split = NULL)`. This is the predictive
#'   cross-fit admissibility rule: in-sample predictive significance
#'   is refused at registration time.
#' - `variable_stability` requires a perturbation path — `refit`,
#'   `bootstrap_action`, or both `core` and `update_core` — **and** one of
#'   `variable_stat` or `loadings`.
#' - `score_stability` requires a perturbation path, `loadings`, and
#'   either `scores` or `project_scores`. For `geometry = "adapter"`,
#'   `project_scores` is required because multifer cannot infer scores
#'   from opaque adapter-owned payloads.
#' - `subspace_stability` requires a perturbation path **and** `loadings`.
#' - `variable_significance` is excluded from v1 and fails at
#'   registration time with a reference to Part 5 §38.
#'
#' Claiming a capability without the backing hooks is an error at
#' [infer_adapter()] time, not at inference time. That is the point
#' of the gate: strict dispatch refuses to build an adapter that
#' cannot back up its claims.
#'
#' @section Geneig rule:
#'
#' Geneig adapters carry one extra registration-time rule. If an
#' adapter claims `(geneig, generalized_eigen, component_significance)`,
#' its `residualize` hook must explicitly declare that it performs
#' B-metric deflation by carrying `attr(residualize, "b_metric") <- TRUE`,
#' or must declare `attr(residualize, "delegates_geneig_deflation") <- TRUE`
#' to signal that `run_geneig_ladder()` owns the deflation. Euclidean
#' residualization is not admissible for the geneig family; the gate
#' fails immediately with a pointer to `notes/engine_geneig_spec.md`.
#'
#' @section Executable validity contracts:
#'
#' Every adapter also carries two assumption layers:
#'
#' \describe{
#'   \item{`declared_assumptions`}{A character vector of assumptions
#'     the adapter **declares** but does not check (e.g.
#'     `"rows_exchangeable"`, `"centered_data"`). These are
#'     trust-the-user claims that appear in the result's
#'     `assumptions$declared` block for transparency.}
#'
#'   \item{`checked_assumptions`}{A list of **executable** validity
#'     checks. Each entry is a list with `name` (character), `check`
#'     (function taking the raw `data` argument and returning a
#'     logical, a character error message, or a `list(passed, detail)`
#'     structured result), and `detail` (character fallback message).
#'     [infer()] runs every check against the input data before the
#'     engine starts; strict mode errors on any failure, permissive
#'     mode records the failures in `result$assumptions$checked`.}
#' }
#'
#' The base-R reference adapters in `R/adapter_oneblock_baser.R` and
#' `R/adapter_cross_baser.R` populate `checked_assumptions` with the
#' standard set for their shape family — numeric-matrix check, finite
#' values, minimum dimensions, paired-row agreement, column variance,
#' and full column rank. Third-party adapters can reuse those helpers
#' or write their own.
#'
#' @section A minimal worked example:
#'
#' The smallest realistic one-block adapter is
#' [adapter_svd()], which wraps `base::svd` and claims all four
#' targets. See `R/adapter_oneblock_baser.R` for the full source.
#' The cross-family reference is [adapter_cross_svd()] in
#' `R/adapter_cross_baser.R`, which demonstrates the paired-row
#' contract and the relation-aware `component_stat` path.
#'
#' For a fitted-model wrapper (the recommended pattern when you
#' already have a fitting framework, e.g. multivarious), see
#' [adapter_multivarious_pca()] and [adapter_multivarious_plsc()] —
#' these are thin shims that delegate every hook to the underlying
#' multivarious generics (`components`, `coef`, `project`, etc.)
#' and cost only a few dozen lines each.
#'
#' @section Recommended pattern for new adapters:
#'
#' 1. Start from the adapter reference example closest to your
#'    shape family (`adapter_svd` for oneblock, `adapter_cross_svd`
#'    for cross).
#' 2. Replace the fitting body of `refit` with a call into your
#'    package's own fitter. Keep the return-value shape compatible
#'    with your `loadings`, `scores`, and `roots` hooks.
#' 3. Populate `checked_assumptions` with concrete checks that catch
#'    the failure modes that matter for your method. Reuse
#'    `.oneblock_baser_checks()` or `.cross_baser_checks()` when they
#'    apply.
#' 4. Add your adapter to the registry via
#'    [register_infer_adapter()] — typically from a zero-argument
#'    function called by your package's `.onLoad` hook.
#' 5. Write an integration test that calls [infer()] with your
#'    adapter id and a synthetic fixture.
#'
#' If your fitted-model class already lives in another package (e.g.
#' you are writing an adapter for someone else's fitter), you do not
#' need to vendor the fitter into your code — just import the classes
#' you need and delegate the hooks.
#'
#' For multiblock fits, declare `geometry = "multiblock"` and pass
#' list-of-matrix data with aligned rows through `refit()`, `null_action()`,
#' `component_stat()`, and `residualize()`. The generic multiblock ladder
#' preserves block-list data; it does not flatten blocks internally. That
#' keeps MFA/JIVE/AJIVE-style refit semantics under adapter control.
#'
#' @seealso [infer_adapter()], [register_infer_adapter()],
#'   [capability_matrix()], [infer()]
#' @name infer_adapter_contract
NULL
