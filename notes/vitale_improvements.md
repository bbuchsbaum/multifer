# Improving Vitale et al. (2017) Sequential Permutation PCA

**Source:** Vitale et al. 2017 — *Selecting the number of factors in principal component analysis by permutation testing: numerical and practical aspects* (PDF at `/Users/bbuchsbaum/Downloads/Vitale et al. 2017 - ...`)

**Status:** Parts 1–5 integrated. Part 5 tightens the framework in response to the Part 3 review. Naming note: the conceptual framing in Part 5 uses "typed shape inference" and cites the working name "ShapeInfer"; the **package name remains `multifer`** as committed at the end of Part 3. Treat "ShapeInfer" in Part 5 as the concept label, `multifer` as the shipping artifact.

**Companion note:** For the intended methods-paper framing and the precise theoretical boundary between direct SVD/eigen latent-root extensions and later supervised / generalized-eigen engines, see [paper1_rank_matched_residual_randomization.md](./paper1_rank_matched_residual_randomization.md).

**Scope shift (from Part 2).** The project is no longer "improve Vitale's PCA permutation test." It is: **turn sequential deflated permutation inference into a general latent-operator inference framework** covering PCA, PLSC, CCA, and later multi-block / generalized-eigen settings, with one shared inferential scaffold and method-specific latent targets where needed.

**Product shift (from Part 3).** The framework crystallizes into a concrete package: **`multifer`** — a shape-first inference layer built on top of `multivarious`'s projector geometry. Working name during design was "ShapeInfer"; final package name is **multifer**. Tagline: *projector-native perturbation inference for latent operator models.*

---

## 0. Starting point: what Vitale et al. got right

- **Keep sequential deflation.** The deflation loop is the right structural choice.
- **Keep the P3 logic, drop P1/P2.** Their own simulations show:
  - P1/P2 are **too conservative**.
  - **P3 gives much better null calibration** and better recovery than Horn / Dray on the same examples.
- So the goal is not to replace their framework — it is to sharpen P3 and its inference.

---

## 1. Exact algebraic simplification of the P3 statistic

**Setup.** Let $M = U\,\operatorname{diag}(s_1,\dots,s_Q)\,V^\top$ be an SVD of a residual matrix.

**Key fact.** The P3 projected matrix
$$(I - U_{1:a-1} U_{1:a-1}^\top)\, M$$
has singular values exactly $s_a, s_{a+1}, \dots, s_Q$ — i.e. the leading $a-1$ left singular directions are annihilated and the rest pass through unchanged.

**Consequence.** Vitale's projected P3 statistic simplifies to:
$$
T_a^{*}(M)
\;=\;\frac{s_a^2}{\sum_{q=a}^{Q} s_q^2}
\;=\;\frac{\lambda_a(M)}{\|M\|_F^2 - \sum_{q<a}\lambda_q(M)}
$$
where $\lambda_q(M) = s_q^2$.

**Implementation payoff.**
- **No need to materialize the projected matrix.**
- **No second SVD.**
- Only the **top $a$ singular values** of each randomized residual are needed.
- The paper's best variant becomes **one partial SVD per resample**, instead of two full decompositions.

---

## 2. Fix the Monte Carlo inference

**Problem with the paper as written.**
- Typical setting: $B = 300$ randomizations at $\alpha = 0.01$.
- P-value resolution is $1/(B+1) = 1/301 \approx 0.0033$ — comparable to $\alpha$ itself.
- Estimating a 1% tail from 300 draws is extremely noisy.
- Thresholding against the empirical 99th percentile compounds the problem.

**Fix 1 — use the finite-$B$ Monte Carlo p-value.**
$$
p_a \;=\; \frac{1 + r_a}{B_a + 1},
$$
where $r_a$ = number of randomized statistics $\ge$ the observed statistic. This is **exactly valid** regardless of $B$.

**Fix 2 — make $B_a$ adaptive (anytime-valid sequential MC tests).**
- Modern sequential Monte Carlo tests allow **stopping as soon as the decision is clearly above or below $\alpha$**.
- P-value remains valid at **arbitrary stopping times**.
- Result: **better calibration** and **much less computation on easy cases** (most components are easy; only the boundary component is hard).

---

## 3. Change the null generator: prefer sign-flips over column permutations

**Why column permutation is fragile.**
- Column permutation is only a good null **when rows are exchangeable** and the noise model truly matches that operation.
- In modern high-dimensional data, **heterogeneous noise is common** (per-row / per-column variances, mild dependence).
- Classical permutation Parallel Analysis can **fail badly** under heterogeneity.

**Better foundation: FlipPA.**
- **Random sign-flips** preserve the distribution of heterogeneous **symmetric** noise.
- Comes with:
  - **Nonasymptotic Type-I control.**
  - **Asymptotic consistency.**
- **BlockFlipPA** extends this under blockwise independence.
- If over-estimation of rank is the main concern, prefer the more **conservative upper-edge comparison**.

**Practical recommendation.** Use sign-flip / block sign-flip / design-matched restricted permutation rather than unrestricted column permutation.

---

## 4. Make the inferential target explicit

**Hypothesis ladder.** Test the nested hypotheses
$$
H_a : \operatorname{rank}(X) \le a - 1
$$
and **stop at the first non-rejection**.

**Overall Type-I control (automatic).**
If each componentwise test is level $\alpha$, then
$$
P(\hat r > r_0) \;\le\; \alpha,
$$
because over-estimating rank requires **falsely rejecting the first true null** in the chain. This is the clean formal statement for "overall Type-I control" — no extra correction needed.

**Strongest finite-sample version: split sample.**
- Estimate the **deflating subspace on one split of the rows**.
- **Test on another split.**
- Costs some power, but removes the main theoretical objection: **reusing the same data to estimate and to test**.

---

## 5. Computational upgrades

- The cleaned-up P3 statistic depends only on the **top $a$ singular values** → use **partial SVD**, never a full decomposition.
- **Very large matrices:** use **randomized low-rank methods** (randomized range finders / block Krylov / randomized Lanczos).
  - Build a small subspace capturing most of the matrix's action.
  - Often beat classical decompositions in speed.
  - Come with detailed, usable error analysis.
- **When exact inference matters more than speed:** deterministic partial eigensolver.
- **When scale matters more:** randomized block Krylov/Lanczos — use the **same accuracy target on both observed and randomized matrices** (critical for fair comparison).

---

## 6. The recommended procedure

For each component $a = 1, 2, \dots$:

1. **Deflate sequentially** (as in Vitale).
2. **Observed statistic:**
$$
T_a \;=\; \frac{\lambda_a(X)}{\sum_{q=a}^{Q}\lambda_q(X)}
\;=\; \frac{\lambda_1(E_a)}{\|E_a\|_F^2}
$$
where $E_a$ is the deflated residual at step $a$.
3. **Randomized residuals:** $M_b = G_b(E_a)$, where $G_b$ is a **sign-flip**, **block sign-flip**, or **design-matched restricted permutation**.
4. **Randomized statistic — partial SVD only:**
$$
T_a^{(b)} \;=\; \frac{\lambda_a(M_b)}{\|E_a\|_F^2 - \sum_{q<a}\lambda_q(M_b)}
$$
Only the top $a$ singular values of $M_b$ are needed.
5. **Adaptive Monte Carlo p-value:** update $p_a = (1+r_a)/(B_a+1)$, stop $B_a$ as soon as the anytime-valid decision crosses $\alpha$.
6. **Stop at first non-rejection.** Report $\hat r = a - 1$.

---

## 7. Where this sits relative to published theory

- **Classical permutation PA** is proved only to **consistently recover the larger components**; it can **miss smaller ones** (shadowing).
- **Dobriban & Owen's deterministic PA** is explicitly **faster** and uses **deflation to fight shadowing**.
- **FlipPA / BlockFlipPA** currently carry the **strongest published guarantees** under heterogeneous noise.

**Recommended stack.**
- **Base selector:** FlipPA or BlockFlipPA (strong theory).
- **Refinement near the cutoff:** the deflated tail-ratio test above (sharper, cheaper, sequentially controlled).

---

---

# Part 2 — Generalization to a Deflated Operator Inference Framework

## 8. Reframing: everything is an operator $K$

Define one operator per method; its **singular values are the latent roots** of interest.

$$
K \;=\;
\begin{cases}
X & \text{PCA}\\[1mm]
X^\top Y & \text{PLS / PLSC}\\[1mm]
(X^\top X)^{-1/2}X^\top Y (Y^\top Y)^{-1/2} & \text{CCA}
\end{cases}
$$

**Nested hypothesis ladder (method-agnostic).**
$$
H_a : \text{no $a$-th latent root beyond the first } a-1.
$$
Stop at the first non-rejection (same logic as §4 — automatic overall Type-I control).

**Package doctrine.** The scaffold is shared across methods:
- ordered latent objects,
- sequential removal of previously claimed structure,
- null actions matched to design and inferential target,
- stop-at-first-non-rejection testing,
- latent-unit stability summaries.

What may change across methods is the **inner target**. PCA, PLSC, CCA, and symmetric-definite generalized-eigen models fit the direct ordered-root story. Supervised methods such as PLS regression may still fit the scaffold, but they should use predictive-increment targets rather than being forced into an `X^T Y` root test.

**Test statistic choices.**
- **PCA / covariance-type PLS:** Vitale-style **tail-ratio** $T_a$ is natural.
- **CCA:** test the $a$-th canonical root itself, or a root-based multivariate statistic, **after stepwise removal of previously explained association**.
- **Structured PCA:** operator view accommodates generalized matrix decompositions with row/column quadratic norms (non-Euclidean geometry).
- **Multi-block (JIVE / AJIVE):** already lives in a low-rank operator world — same engine applies.

**Caveat (Winkler et al.).** Naïve permutation CCA **inflates error beyond the first root**. Valid CCA permutation requires:
1. **Stepwise residualization** of previously explained association.
2. An **exchangeability-preserving basis transform** when nuisance variables are present.

---

## 9. Core-space resampling: the big computational win for PLS/CCA

The operator view gives an even larger speedup for PLS/CCA than for PCA, because **the permuted problem collapses to a tiny core**.

**Thin SVDs (done once).**
$$
X = U_x D_x V_x^\top, \qquad Y = U_y D_y V_y^\top
$$
with ranks $r_x \ll p$, $r_y \ll q$.

**Permuted PLS cross-product.** For a row permutation $P$ of $Y$,
$$
X^\top P Y \;=\; V_x\big[D_x (U_x^\top P U_y) D_y\big] V_y^\top.
$$
Every permuted PLS singular value is determined by the **$r_x \times r_y$ core** $D_x (U_x^\top P U_y) D_y$, **not** by the $p \times q$ cross-covariance.

**Permuted CCA canonical correlations.** Even cleaner — they are the singular values of
$$
U_x^\top P U_y.
$$
With general row weights $W$:
$$
\sigma\!\left( A_W^{-1/2} C_W B_W^{-1/2} \right),
\quad
A_W = U_x^\top W U_x,\;
B_W = U_y^\top W U_y,\;
C_W = U_x^\top W U_y.
$$

**Bottom line.** After **one thin SVD per block**, each resample becomes an $r_x \times r_y$ problem. This is exactly where speed matters most — CCA/PLS overfit badly when $n$ is not comfortably larger than the variable count, so you need many resamples *and* many regularization paths.

---

## 10. Cross-fitted sparse / regularized PLS & CCA

**Problem.** Testing an **in-sample** sparse CCA/PLSC objective is circular — the objective has already been optimized on the same data.

**Fix (cross-fitted test).**
1. Estimate sparse weights on a **training split**.
2. Evaluate the covariance / correlation those fixed weights achieve on a **held-out split**.
3. **Randomize only the held-out split.**

**Signpost.** TOSCCA (fixes sparsity directly, uses permutation, argues explicitly for **out-of-sample** correlation as the inferential target) is the right reference point for this extension.

---

## 11. Two parallel resampling streams: significance vs. stability

**Design principle (from PLS literature — McIntosh/Krishnan).** Keep these as **separate outputs**:
- **Permutation** answers *"Is this latent variable stronger than noise?"* → component p-values.
- **Bootstrap** answers *"Which saliences are stable / reliable?"* → variable and score stability.
- **Procrustes alignment** is needed because resampling can swap axes or flip signs.

**Stream A — invariant null (for p-values).** Sign-flip / FlipPA / block sign-flip / restricted permutation, run through §9 core-space formulas.

**Stream B — weighted perturbation (for stability).** Prefer row weights $W_r$ and optional variable weights $W_x, W_y$ over literal resampling. With variable weights:
$$
R_x^\top R_x = V_x^\top W_x V_x, \qquad
R_y^\top R_y = V_y^\top W_y V_y,
$$
and the weighted PLS core becomes
$$
R_x\, D_x (U_x^\top W_r U_y)\, D_y\, R_y^\top.
$$
So observation- **and** variable-level perturbation never touches the full $p$- or $q$-dimensional matrices.

**Weight choices.**
- **Smooth weights** → stability maps.
- **Bernoulli variable thinning** → sensitivity to feature inclusion.

---

## 12. Variable significance vs. subspace stability vs. score uncertainty

Weighted bootstrap gives **robustness**, not null p-values. So add targeted layers:

**Variable p-values — jackstraw.**
- Jackstraw was built specifically to test which variables are genuinely associated with estimated PCs, avoiding the circularity of naïve tests.
- Already extended to **AJIVE loadings** (Yang et al.) — so the framework inherits variable-level inference for multi-block decompositions.
- **Use this layer**, not weighted bootstrap, when you need null-calibrated variable p-values.

**Subspace stability — principal angles.**
- Near **tied roots**, individual loading vectors are **not separately identifiable**.
- Report **principal angles** or an **ADLS-style angle distribution** instead of per-loading intervals.

**Score uncertainty — low-dim bootstrap reps.**
- Produce **score ellipses** and coordinate-wise intervals efficiently from the core-space representation.

---

## 13. Two-stage inference engine

1. **Screen** with a fast, reproducible method:
   - **DPA (Dobriban–Owen deterministic PA)** — faster, addresses shadowing via deflation.
   - **FlipPA / BlockFlipPA** — nonasymptotic Type-I control under heterogeneous (possibly block-dependent) symmetric noise.
2. **Refine near the boundary** with the sequential deflated tail-ratio test (§6), using:
   - Monte Carlo p-values $p_a = (1+r_a)/(B_a+1)$.
   - **Anytime-valid** sequential MC stopping.
   - Multiple-testing-aware sequential rules where relevant.
3. **Calibration benchmark.** For the Gaussian homoskedastic matrix model, use existing **exact selective-inference rank tests** as the ground-truth calibration check.

---

## 14. What is actually novel enough for a paper

"Vitale but faster" alone is a methods note. The **paper-worthy contribution** is a single unified package:

> **A deflated operator inference framework** covering PCA, PLS, CCA, and multi-block decompositions, with:
> - method-specific operators $K$ and null groups;
> - **stepwise** valid testing of latent roots;
> - **core-space** resampling formulas for speed;
> - **cross-fitted** sparse / regularized extensions;
> - **jackstraw-style** variable p-values;
> - **bootstrap / weight-based** stability for variables and scores;
> - one shared resampling engine producing **component p-values, variable p-values, variable stability, and score uncertainty in one report.**

**Why it's novel.** The pieces exist separately in the literature:
- **Vitale** — deflated permutation PCA.
- **Dobriban / Hong** — better rank-selection nulls (DPA).
- **Winkler et al.** — valid stepwise CCA permutation with nuisance-aware basis transform.
- **McIntosh / Krishnan** — PLS permutation + bootstrap pairing.
- **Chung / Storey** and **Yang et al.** — jackstraw for variable-level inference (including AJIVE).

No single framework currently unifies these around **one common operator + resampling engine**.

---

## 15. Concrete project scope

Deliver a package that:
1. **Screens** with DPA / FlipPA.
2. **Refines** with sequential deflated tail-ratio tests.
3. Runs a **core-space perturbation layer**.
4. Outputs **in one report**:
   - Component p-values.
   - Variable p-values (jackstraw).
   - Variable stability (weighted bootstrap + principal angles near ties).
   - Score uncertainty (low-dim bootstrap ellipses).

---

## Open questions / TODO (updated)

**From Part 1 (still open):**
- Precise choice of sequential MC test (Besag–Clifford, Gandy, Gandy–Hahn, newer anytime-valid).
- Exact form of "design-matched restricted permutation" when neither permutation nor sign-flip fits.
- Split-sample row splitting under longitudinal / blocked designs.
- Accuracy-target matching between observed and randomized partial SVDs.
- Benchmarking vs. Vitale P3, Horn, Dray, DPA, FlipPA, BlockFlipPA.

**New from Part 2:**
- Operator definition for **generalized / structured PCA** under non-Euclidean row/column metrics — exact form of $K$ and matching null group.
- For CCA under nuisance variables: which **exchangeability-preserving basis transform** (Winkler-style) to adopt as the default.
- **Cross-fit split strategy** for sparse CCA/PLS: random, blocked, or stratified; single split vs. multi-split averaging.
- **Procrustes alignment protocol** for bootstrap replicates under tied / near-tied roots (fall back to principal angles).
- **Jackstraw adaptation** to PLS/CCA/AJIVE latent variables (beyond PCA and AJIVE loadings already published).
- Multiple-testing handling across the $a$-indexed ladder when combined with variable-level jackstraw tests.
- API design so that a single call returns {component p-values, variable p-values, stability, score uncertainty} without redundant SVDs.
- Benchmarking additions: TOSCCA (sparse CCA), AJIVE, PLSC / McIntosh-Krishnan PLS.

---

# Part 3 — `multifer`: Product Requirements Document

## 16. Why `multifer` (design rationale)

`multivarious` already has the right **geometric substrate**:
- Classes: `projector`, `bi_projector`, `cross_projector`, `multiblock_projector`.
- Operations: `shape()`, `partial_projector()`, `compose_projector()`, `refit()`.
- Beginnings of an inference layer: generic `perm_test()` with `shuffle_fun`, `measure_fun`, `fit_fun` hooks; `bootstrap_pca()`, `bootstrap_plsc()`.
- Non-PCA shapes already supported: `cPCAplus()`, `geneig()`.

What's **missing** is a clean **typed inference layer** above that geometry. Vitale gives the algorithmic template (sequential deflation, relative latent-root test, stop at first non-significant component) but is still framed as PCA-specific.

**Key architectural decision — don't overload `shape()`.**
- `shape()` stays **lightweight and structural** (dimensions of the loading map).
- Add a second object, `infer_spec()` (internal) / `infer_recipe()` (public advanced layer), carrying **inferential semantics**:

```r
infer_spec(
  domains,        # X, Y, blocks, foreground/background
  sample_axis,    # shared / paired / grouped rows
  root_type,      # variance, covariance, correlation, generalized-eigen, loss
  truncate,       # keep first k factors
  residualize,    # method-specific deflation / nuisance removal
  refit,          # fallback slow path
  core,           # low-dimensional representation (optional)
  update_core,    # fast path for permutations / weights
  align,          # sign / permutation / Procrustes / subspace
  null_action,    # permutation / sign-flip / block action + assumptions
  component_stat, # rank / latent-root statistic
  variable_stat,  # variable-level statistic
  score_stat      # score-level statistic
)
```

**Why not one universal engine.** The same geometric shape can carry different inferential meanings. A `cross_projector` can be PLSC, CCA, or reduced-rank regression — but the valid null and stepwise test differ. CCA is the warning case: naïve permutation breaks beyond the first canonical root; nuisance residualization needs a Winkler-style lower-dimensional exchangeable basis.

**Solution — a small basis of inferential shapes:**
1. **One-block latent decomposition** — PCA, kernel / Nyström PCA, generalized PCA, contrastive single-block.
2. **Cross-block latent coupling** — PLSC, CCA, two-block SVD, reduced-rank regression.
3. **Multiblock consensus / union** — block PCA, JIVE / AJIVE, multiblock projectors.
4. **Generalized-eigen / contrastive** — cPCA++-type, whitened operators, metric-weighted decompositions.

---

## 17. Core computational principle: core-space resampling as a first-class hook

The fastest version is **not** "refit every time." It is:
1. Fit once.
2. Compress to a low-dimensional **core**.
3. Update only the core under permutations / weights.
4. Lift back only when needed.

`multivarious` already does this in two places:
- `bootstrap_pca()` — Fisher et al. exact-bootstrap idea; works on a low-dim resampled matrix then lifts back.
- `bootstrap_plsc()` — returns bootstrap ratios directly.

**Generalize as two hooks:**
```r
core(x, data, ...)
update_core(core_obj, obs_weights = NULL, var_weights = NULL, action = NULL, ...)
```

Once adapters expose score bases $U_b$ and loading bases $V_b$:
- **row permutations / row weights** only affect small matrices like $U_i^\top W U_j$;
- **variable weights** only affect small matrices like $V_b^\top M_b V_b$;
- **cross-block updates** live in $r_x \times r_y$ space instead of $p_x \times p_y$.

**Policy.** If `update_core()` exists, use the fast path; otherwise fall back to `refit()`. Abstraction without performance loss.

---

## 18. Rigor: typed `null_action` objects

Replace loose `shuffle_fun` with a **typed null-action object** carrying assumptions:

```r
null_action(
  kind        = "permute" | "signflip" | "block_permute" | "wild",
  sample_groups = NULL,
  synchronize = c("X", "Y"),
  nuisance    = NULL,
  exact       = TRUE,
  description = "paired-row permutation of Y within exchangeability blocks"
)
```

**Wins.** Validity checks, reportable assumptions, reuse across shapes.

**Monte Carlo layer is sequential.** Anytime-valid sequential MC tests replace fixed `nperm`. Essential once one run produces component tests + variable diagnostics + score intervals together.

---

## 19. Separation of inferential targets

One perturbation engine, **four distinct targets**:

| Target | Engine | Output type |
|---|---|---|
| Component significance | null-preserving perturbations | p-values |
| Variable stability | bootstrap / Bayesian boot / jackknife / variable dropout | stability measures |
| Score stability | same perturbations propagated to scores | intervals / ellipses |
| Variable significance | null-calibrated variable test (adapter-declared) | p-values |

**Critical distinction.** A **bootstrap ratio is a stability measure, not a p-value.** The framework must never silently relabel one as the other.

**Leverage existing hooks.** `partial_projector()` and `partial_inverse_projection()` (cross-projector version already allows ridge / direct-solve overrides) are the right extension points for variable weighting, feature subsets, and leave-many-variables-out diagnostics.

**Second perturbation object:**
```r
perturbation_family(
  obs  = "permute" | "bootstrap" | "bayes_boot" | "jackknife",
  vars = "identity" | "weights" | "dropout" | "thinning",
  paired = TRUE,
  synchronized = TRUE
)
```

---

## 20. Subspace stability is first-class

Near tied roots, the stable object is the **subspace**, not a signed loading vector.

`multivarious` already exposes `prinang`, `principal_angles`, `subspace_similarity` — **promote these from helpers to core inference outputs**. Sign / Procrustes alignment when roots are well-separated; always also report subspace stability. Principal angles make sense across PCA, PLSC, CCA, multiblock, and generalized-eigen models.

---

## 21. Extension paths (two supported)

**Path A — native S3:**
```r
infer_adapter.my_model <- function(x, ...) { ... }
```

**Path B — zero-dependency adapter constructor + registry:**
```r
adapter <- infer_adapter(
  shape, roots, scores, loadings, truncate, residualize, refit,
  core, update_core, align, null_action,
  component_stat, variable_stat, score_stat
)
register_infer_adapter("my_model", adapter)
```

**Why both.** External packages should integrate **without** inheriting from `multivarious` classes or importing its whole object system.

---

## 22. PRD — `multifer`

### 22.1 Product summary

`multifer` is a **shape-first inference layer** for `multivarious` providing one coherent system for:
- component significance,
- variable significance or stability,
- score uncertainty or stability,
- subspace stability,

across projector-based multivariate models, without hard-coding every method.

**Pipeline:** fit → represent as latent shape → compile to inference recipe → execute shared perturbation engine → return standardized result.

### 22.2 Problem statement

Current multivariate inference is fragmented: method-specific, inconsistently named, slow (refits every perturbation), and statistically muddled (significance and stability mixed). `multivarious` has the geometry; `multifer` adds the clean inference layer above it.

### 22.3 Goals

**Primary (v1):**
- Shape-based, not method-name-based.
- Simple API with strong defaults.
- Adapter-based extension.
- Explicit about null assumptions and validity.
- One engine for significance + stability + uncertainty where appropriate.
- Fast via low-dimensional core updates.

**Non-goals (v1):**
- A single theorem covering every multivariate estimator.
- Exact variable-level p-values for every rotated / sparse / regularized method.
- Automatic null discovery for arbitrary black-box models.
- Day-one support for every multivariate method.

### 22.4 Users

- **End users** (applied analysts): `infer_pca()`, `infer_pls()`, `infer_cca()`.
- **Advanced users** (methodologists): `infer_recipe(...)` overrides.
- **Package authors**: `register_infer_adapter("my_model", infer_adapter(...))`.

### 22.5 Product principles

1. **Shape first, aliases second.** Small set of shapes, not a growing list of named methods.
2. **One engine, two perturbation modes.** Null-preserving for significance; stability stream for uncertainty/robustness.
3. **Honesty about guarantees.** Every wrapper/adapter declares: exact / conditional / asymptotic / heuristic.
4. **Subspaces before vectors.** Near-tied roots → report subspace stability.

### 22.6 Scope (v1)

**Shape families built in:**
- one-block, cross-block, multiblock, generalized-eigen.

**Outputs:**
- component significance, variable stability, score stability, subspace stability, variable significance (adapter-declared only).

### 22.7 Public API

**Primary:**
```r
infer(model, data = NULL, shape = "auto", targets = "default", ...)
```
Behavior: auto-detect adapter from class, auto-select default shape, choose sensible null/stability engines, return standardized result.

**Convenience wrappers (the documented entry points):**
```r
infer_oneblock(model, X, ...)
infer_cross(model, X, Y, ...)
infer_multiblock(model, blocks, ...)
infer_geneig(model, A, B = NULL, ...)
```

**Method aliases — thin wrappers:**
```r
infer_pca(...)
infer_pls(...)
infer_cca(...)
infer_cpca(...)
```

**Advanced (public):**
```r
infer_recipe(shape, root = NULL, null = "auto", targets = "default",
             align = "auto", engine = "auto")
```

**Extension:**
```r
infer_adapter(...)
register_infer_adapter(class, adapter)
```

**Internal / expert-only:**
```r
infer_spec()
compile_infer_recipe()
```

### 22.8 Built-in shape wrappers — defaults

**`infer_oneblock()`** (PCA, rotated PCA, sparse PCA with cross-fit / held-out):
- Root: variance / singular-value.
- Significance: sequential deflated latent-root testing (Vitale-style generalization).
- Null: column permutation or sign-flip (recipe-dependent).
- Stability: bootstrap or weighted perturbation.
- Subspace: principal angles.

**`infer_cross()`** (PLS, CCA, reduced-rank):
- Root: covariance or correlation.
- Null: break row correspondence in one block.
- Significance: stepwise testing of cross-block roots.
- Stability: paired bootstrap / weighted row perturbation.
- Variable significance only when adapter supports it.

**`infer_multiblock()`** (consensus / joint / shared latent):
- Root: consensus / shared-root measure.
- Null: independent block perturbation or synchronized block-specific action.
- Outputs: joint component significance + per-block stability.

**`infer_geneig()`** (cPCA, discriminant, metric-weighted, generalized-eigen):
- Root: generalized eigenvalue.
- Null: recipe-dependent.
- Stability emphasizes subspaces and generalized loadings.

### 22.9 Statistical engine requirements

**Component significance:** nested sequential testing; stop-at-first-non-significant; Monte Carlo p-values; optional sequential / early-stopping resampling; adapter-supplied test statistics. Default one-block recipe = Vitale generalization (deflate, test relative latent root, P3-style residual-rank geometry).

**Variable / score stability:** bootstrap, Bayesian boot / smooth weights, variable weighting, variable thinning / dropout, component or subspace alignment, domain-specific summaries for X / Y / blocks. Default interpretation: bootstrap ratio / weight sensitivity = **stability, not p-value**.

**Variable significance:** opt-in per adapter. Framework **must not** silently relabel stability as significance.

**Subspace stability:** first-class. Principal-angle distributions, component-swap-aware summaries, subspace support intervals where possible.

### 22.10 Core computation model

Two execution modes:
- **Slow / universal:** perturb data → refit → extract stat (requires `refit()`).
- **Fast core-update:** fit once → compute core → update core under perturbation → lift. Supported models get the fast path automatically.

`core()` and `update_core()` are optional hooks; absent → fall back to `refit()`.

### 22.11 Validity contract (required on every built-in and adapter)

```r
validity = list(
  level       = c("exact", "conditional", "asymptotic", "heuristic"),
  assumptions = c(...),
  target      = c("components", "variables", "scores", "subspaces"),
  notes       = ...
)
```

Forces the framework to state what's actually justified instead of implying universal validity.

### 22.12 Output schema — `infer_result`

```
result$components   # per component: stat, p-value, MC uncertainty, stop flag, null label, validity
result$variables    # per variable per domain: estimate, stability, selection freq / weight sensitivity, p, q
result$scores       # per observation: estimate, interval / spread, leverage / influence, domain
result$subspaces    # per component set: principal-angle summary, alignment method, stability label
result$assumptions
result$engine
result$provenance
result$call
```

### 22.13 Performance requirements

- ≥ **5× speedup** over naive refit-per-perturbation for built-in shapes with `update_core()`.
- Parallel perturbations.
- Cache deflations and reusable intermediates.
- Avoid repeated full decompositions when partial / core suffice.

**Reproducibility:** deterministic seeds, recorded perturbation settings, recorded null/stability recipes, recorded adapter version + capability flags.

### 22.14 Compatibility

Must work with: native `multivarious` classes; external classes via adapters; future wrappers. **No requirement** that outside packages inherit from `projector` — adapter registration is sufficient.

### 22.15 Documentation tiers

1. **User docs** ("just give me answers"): PCA significance, PLS stability, CCA component testing, multiblock summary.
2. **Advanced docs** ("override defaults"): custom null, custom alignment, custom variable target.
3. **Developer docs** ("plug my package in"): minimal adapter, fast adapter with `update_core()`, validity contract, capability declarations.

### 22.16 Success metrics (v1)

- Users run `infer_pca()` / `infer_pls()` / `infer_cca()` / `infer_multiblock()` and get the **same result schema**.
- End users almost never need `infer_spec()`.
- Third-party packages add support with **one adapter object**.
- Built-in shapes are materially faster than refit-only perturbation / bootstrap.
- Docs cleanly separate significance from stability.

### 22.17 Rollout plan

**Phase 1 — core architecture.** `infer()`, `infer_recipe()`, `infer_adapter()`, `register_infer_adapter()`, `infer_result`. Shapes: one-block + cross-block. Outputs: component significance, variable stability, score stability, subspace stability.

**Phase 2 — fast engine.** `core()` / `update_core()`, caching, early-stop resampling, parallel execution.

**Phase 3 — broader shapes.** Multiblock, generalized-eigen / contrastive, more method aliases.

**Phase 4 — richer inference.** Variable significance where supported, multiple-testing, reporting / plotting, external adapter gallery.

### 22.18 One-line pitch

> `multifer` turns projector-based multivariate fits into a common inference language: one API for latent-root significance, variable and score stability, and subspace uncertainty, across one-block, cross-block, multiblock, and generalized-eigen shapes.

---

## 23. Immediate next step

Turn this PRD into a **GitHub epic with ~10 implementation issues**:
1. API surface
2. Adapter contract
3. One-block engine
4. Cross-block engine
5. Result schema
6. Fast core path
7. Docs
8. Tests
9. Benchmarks
10. External adapter examples

---

## Open questions / TODO (updated after Part 3)

**New from Part 3:**
- Exact boundary between `infer_spec()` (internal) and `infer_recipe()` (public advanced layer) — what does the public override surface actually expose?
- Auto-detection of `shape` from model class: dispatch table vs. capability probe vs. S3 method.
- Format and enforcement of the **validity contract** — machine-checkable vs. documentation-only.
- Cache invalidation policy between sequential deflation steps and core updates.
- Parallel backend (`future`, `mirai`, native `parallel`) — pick one default.
- Interaction between early-stopping sequential MC and multiple-testing correction across components + variables.
- Standard `print.infer_result()` / `summary.infer_result()` layout — must communicate significance vs. stability distinction without being read as p-values by default.
- Minimum adapter surface for a useful external integration — how small can a "day-one adapter" be and still work?
- Naming inside `multifer` vs. legacy `multivarious` (`perm_test`, `bootstrap_pca`, `bootstrap_plsc`): deprecation or wrapping strategy.

---

# Part 4 — Economy, Efficiency, and Scope as Design Spine

Optimization is **not a Phase 2 concern**. It is the design constraint that determines whether the scope of §22 can expand without the library collapsing. Part 4 promotes economy and efficiency from "performance requirements" to **the architectural spine** everything else hangs off.

## 24. The compression principle

Every supported operation in `multifer` must be expressible as one of three buckets:

1. **One-time full-data cost** — thin SVD, factorization, basis computation. Paid once per fit.
2. **Repeated core-only cost** — core size depends on $(r_x, r_y, k)$, **never** on $(n, p, q)$. This is where perturbation budgets are spent.
3. **One-time lift-back** — only when the user asks for full-dimensional output (variable tables, score ellipses in original space).

**Enforcement rule.** If an operation cannot be placed in one of these buckets, it **does not belong in the hot path**. This is the rule that prevents scope creep from degrading economy as new shapes are added.

---

## 25. Unification collapses code, not just runtime

The four shape families (one-block, cross-block, multiblock, generalized-eigen) share more than Part 3 admits. They are all **"singular values of a transformed operator"** — only the transform differs.

**Consequence.** Write the engine against a single abstract interface:

```r
linear_operator(
  apply,        # v |-> K v
  apply_t,      # v |-> K^T v
  dim,          # (nrow, ncol)
  core_update   # optional: fast perturbation hook
)
```

The four shape wrappers then become **thin** — they construct a `linear_operator` and hand it to the engine. That is **economy of scope**: fewer code paths → fewer tests → fewer bugs → easier adapter authoring. Dense matrices are just the trivial implementation of this interface.

---

## 26. Matrix-free by default

Partial SVDs via **Lanczos / LOBPCG / randomized range finders** need only `apply` and `apply_t`. Writing the engine matrix-free from day one costs little and unlocks, with zero additional engine code:

- **Kernel / Nyström PCA** — a `linear_operator` over a Nyström factorization.
- **Structured PCA** with non-Euclidean row/column metrics (the open question from §8).
- **Very wide $p$** without materializing $X^\top X$.
- **Streaming / out-of-core blocks** — data passes through `apply`; core stays small.
- **GPU backends** — one `linear_operator` implementation, no engine changes.
- **Tensor decompositions (Tucker, CP)** — unfoldings become operators; resampling inference reuses the same core-update machinery. `multifer` quietly becomes a tensor inference framework without a separate subsystem.

---

## 27. Amortize across the perturbation budget, not within a single fit

Part 2's §9 amortizes one fit. The bigger economy win is **amortizing across the full $B$-resample budget**:

- Precompute $U_x, U_y$ **once**.
- Precompute Gram matrices $U_x^\top U_x$, $U_y^\top U_y$, and any nuisance projectors **once**.
- Each resample touches only a permutation-indexed **view** into these precomputed matrices.
- **Shared randomized sketches** across resamples when the sketch is null-invariant:
  - **Safe:** sign-flip + Rademacher sketch (composes cleanly).
  - **Unsafe:** permutation + Gaussian sketch (does not).
  - Needs a compatibility table worked out per (`null_action`, sketch) pair. **New TODO item.**

---

## 28. Sequential MC as a budget allocator (not just an early stop)

Part 1 §2 framed anytime-valid MC as "stop cheap cases early." The stronger framing:

> **It is a budget allocator across components.**

Give the test engine a **total budget** $B_{\text{total}}$. Let easy components **release unused budget** to the boundary component. Compute concentrates exactly where precision matters — the component right at the decision boundary — instead of being uniformly wasted.

**Implication for the sequential MC layer.** Its API must expose a global budget, not a per-component `nperm`.

---

## 29. Scope expansion that comes free from the compression principle

Once the engine is **matrix-free + core-space native**, the following land with minimal additional code — they are not separate projects:

- Kernel / Nyström PCA
- Sparse PCA / sparse CCA via cross-fit (core update identical; only fit layer changes)
- Streaming / out-of-core decompositions
- GPU backends
- Tucker / CP tensor decompositions on unfoldings

This is the argument for building the compression principle **first**: it is a force multiplier on scope.

---

## 30. Optimization budget — concrete ranking (spend in this order)

1. **Thin-SVD cache + `linear_operator` interface.** Everything depends on this. Build first.
2. **Core-update paths for one-block and cross-block.** Largest constant-factor wins.
3. **Sequential MC as budget allocator.** Largest wall-clock win on realistic problems where most components are easy.
4. **Parallel perturbations over the core.** Near-free once the core is in hand. **Default backend: `mirai`** (lower overhead than `future`, good fit for short core-update tasks).
5. **Randomized partial SVD** as the default `apply`-based solver, with **deterministic Lanczos as fallback** for small / ill-conditioned cases.
6. **GPU / matrix-free backends.** Last — they are a payoff on earlier decisions, not an independent effort.

---

## 31. Two things to explicitly avoid

- **Premature C++ / Rcpp.** The hot path is BLAS/LAPACK calls that are already native. Rcpp buys little until R-level orchestration is proven, and it **costs adapter authorability** (external packages would need C++ toolchain buy-in).
- **Method-specific fast paths inside the engine.** Every `if (class == "pls")` branch is a future bug. Fast paths belong **in adapters**, not in the engine. The engine stays method-free.

---

## 32. Cost accounting in `infer_result`

The `infer_result` schema (§22.12) must carry a **`$cost` block** so the economy argument is testable:

```
result$cost
  $full_data_ops      # number of O(n*p) operations
  $core_updates       # number of O(r_x * r_y) operations per component
  $mc_budget_spent    # per-component MC draws used
  $cache_hits         # cache hit rate on Gram / projector caches
  $wall_time_phases   # time broken down by fit / core / mc / lift
```

Without this block, "materially faster" in §22.16 stays a claim. With it:
- Users have something to optimize against.
- Benchmarks become reproducible and interpretable.
- Adapter authors can see which hooks are actually reducing cost.

**Schema implication.** Add `$cost` to §22.12 as a required field of `infer_result`.

---

## 33. Updated TODO (new items from Part 4)

- **Null-action × sketch compatibility table.** Work out which (null kind, random sketch) pairs preserve validity when the sketch is shared across resamples.
- **`linear_operator` interface specification.** Exact minimum surface: `apply`, `apply_t`, `dim`, optional `core_update`, optional `gram`.
- **Dense adapter** as the reference trivial implementation — this is the one every test suite hits first.
- **Global MC budget API** for the sequential test layer (replaces per-component `nperm`).
- **Backend decision for parallel execution:** commit to `mirai` as default, document fallback to `future`.
- **Benchmark harness** that exercises the `$cost` block across shapes, adapters, and problem sizes. Required before v1 release so §22.16's "materially faster" has numbers.
- **Randomized-solver reproducibility policy** — log the sketch seed and sketch dimension alongside the stopping index from §22.13.
- **Lift-back policy.** When is full-dimensional output computed eagerly vs. lazily? Lazy by default — materialize only on `as.data.frame(result$variables)` or explicit request.

---

# Part 5 — Typed Shapes, Strict Dispatch, Machine-Enforced Validity, Latent Units

Part 5 is a **revision in response to the Part 3 review**. It does not change the core story (Vitale's template: sequential deflation, P3 correction, relative $F_a$ statistic) — it **sharpens** it around four ideas:

1. **Typed shapes**, not raw shapes.
2. **Strict dispatch**, not silent auto-magic.
3. **Machine-enforced** validity, not printed fields.
4. **Latent units**, not per-component rows.

The effect is a smaller, safer, more publishable v1.

---

## 34. Typed shapes

"Shape" alone is too weak. A `cross_projector` can mean PLSC, CCA, or reduced-rank regression, and those meanings imply different nulls. The real inferential object is:

$$
\text{typed shape} \;=\; (\text{geometry},\; \text{relation},\; \text{design})
$$

where:

- **geometry** $\in$ `{oneblock, cross, multiblock, geneig}`
- **relation** $\in$ `{variance, covariance, correlation, generalized-eigen, ...}`
- **design** $\in$ `{exchangeable_rows, paired_rows, blocked_rows, nuisance_adjusted, ...}`

**What this fixes.** The biggest silent-failure mode in the Part 3 design. A `cross` geometry can be auto-detected, but it **cannot silently collapse into PLS-like or CCA-like inference** unless the adapter declares that choice is unique.

**Conceptual label.** The project concept is **typed shape inference** (working title: *ShapeInfer*). The shipping package is still `multifer`.

---

## 35. Strict dispatch: ban silent ambiguity

`infer()` stays convenient, but silent auto-magic is banned.

**Dispatch rules.**
- **Auto-detect geometry** when obvious from the model class.
- **Require explicit relation** unless the adapter declares a single unambiguous default.
- **Require explicit design** for **significance** output unless a safe default exists for the combination.

**Example — allowed:**
```r
infer_cross(mod, X, Y, relation = "covariance")
infer_cross(mod, X, Y, relation = "correlation")
```

**Example — must error under strict dispatch:**
```r
infer(mod)   # adapter exposes more than one plausible recipe
```

**Why this matters.** This is the **single best protection against wrong answers**. An `infer()` call that quietly picks between PLS-like and CCA-like inference is exactly the class of bug the framework exists to prevent. `strict = TRUE` is the default.

---

## 36. Machine-enforced validity contract

The Part 3 validity contract (§22.11) was a data structure. Part 4 promoted cost accounting to a first-class field. Part 5 does the same for validity: **it becomes an executable gate, not a printed field.**

**Three enforcement points.**

1. **Registration-time.** An adapter **cannot claim** component significance unless it provides the hooks needed to compute it (e.g. `null_action`, `component_stat`, `residualize`). Claiming a capability without the hooks is a registration error.

2. **Compile-time.** A requested `(recipe, target)` combination **must be supported** by the adapter's capability matrix. Mismatches fail at `infer_recipe()` compilation, before any data is touched.

3. **Run-time.** **Checked assumptions must pass** before significance is computed. If a sign-flip null requires symmetric noise and the data fails a symmetry check, the result object emits `validity_level = "heuristic"` — or errors under `strict = TRUE`.

**Adapter declarations (required fields):**
```r
capabilities         # which (recipe, target) combinations work
checked_assumptions  # runtime-testable: exchangeability, symmetry, rank
declared_assumptions # trust-the-user: noise model, design correctness
validity_level       # exact | conditional | asymptotic | heuristic
```

**Rule.** `infer()` **refuses to return `"exact"` or `"conditional"` significance** unless the checked assumptions pass. This turns validity from documentation into policy.

---

## 37. Latent units as the central schema entity

Part 3's result schema (§22.12) is still **component-centric**: variables and scores attach to "component 3." Near tied roots, this is dishonest — component 3 may not be individually identifiable.

**Fix.** Introduce a first-class **latent unit**:

```r
units
  $unit_id
  $unit_type    # "component" (identifiable singleton) | "subspace" (tied/grouped)
  $members      # component indices inside this unit
  $identifiable # TRUE if sign/Procrustes aligns; FALSE if only subspace aligns
  $selected     # TRUE if the test ladder accepted this unit
```

**Every downstream table keys off `unit_id`**, not a raw component index.

- `variable_stability` attaches to a unit.
- `score_stability` attaches to a unit.
- `component_tests` has one row per unit.

**Consequence.** When roots are well-separated, each unit is a singleton component — behavior is identical to the Part 3 schema. When roots are tied, variables and scores attach to a **subspace unit**, not a fictitious individual component. Tied-root handling becomes honest without requiring the user to notice.

---

## 38. Cut variable significance from v1

Agreed with the Part 3 review, now codified.

**v1 output matrix:**

| Target | v1 status |
|---|---|
| Component significance | ✅ yes |
| Variable stability | ✅ yes |
| Score stability | ✅ yes |
| Subspace stability | ✅ yes |
| **Variable significance** | ❌ **no** (experimental / internal only) |

**Why.** Bootstrap ratios are stability summaries, not p-values. Jackstraw-style null-calibrated variable tests are load-bearing but the least-mature piece in the stack (Part 2 §12). Keeping that line sharp is one of the strongest parts of the whole proposal. **Do not blur it in v1.**

---

## 39. Narrower v1, with a real Phase 0

Phase 1 in §22.17 is still too large. Insert **Phase 0** explicitly:

### Phase 0 — freeze the contracts  **[COMPLETE — 2026-04-12]**
- Freeze **result schema** (unit-centric, with `$cost` block from §32).
- Freeze **benchmark suite** (§40 below).
- Write **one synthetic one-block example** and **one synthetic cross-block example** end-to-end.
- Implement **strict recipe compilation** and **adapter validation** (§36).
- **No engine code yet** — this phase is about contracts, examples, and gates.

**Delivered (Phase 0):**
- `R/geometry.R`, `R/relation.R`, `R/design.R`, `R/typed_shape.R` — typed shape vocabulary (§34).
- `R/infer_result.R` — frozen unit-centered schema with 9 sub-blocks, foreign-key validation on `unit_id`, `$cost` and `$mc` blocks in place (§32, §37, §41, §44).
- `R/capability_matrix.R` — `capability_matrix()` / `supports()` / `capabilities_for()` grid.
- `R/adapter.R` — `infer_adapter()` with registration-time validity gate (§36 point 1); registry (register / get / list / unregister / clear); `variable_significance` blocked with a Part 5 §38 reference.
- `R/infer_recipe.R` — strict dispatch (§35); errors on ambiguous relation, unsupported target, v1-blocked `variable_significance`, geometry-not-in-shape_kinds; `strict = FALSE` auto-picks and downgrades `validity_level` one step.
- `R/bench_generators.R` + `R/bench.R` + `R/zzz.R` + `inst/bench/*` — four locked benchmark suites (§40): `oneblock_null`, `oneblock_shadowing`, `cross_null`, `speed_agreement`; target numbers locked in metadata.
- `tests/testthat/helper-stub-adapters.R` — shared stub factories consumed by both canaries.
- `tests/testthat/test-oneblock-contract.R` — Phase 0 one-block canary (10 blocks).
- `tests/testthat/test-cross-contract.R` — Phase 0 cross-block canary (11 blocks), including the strict-mode PLS/CCA ambiguity negative test.
- 125 `test_that` blocks across 12 test files; `R CMD check → Status: OK` (0 errors, 0 warnings, 0 notes).
- Tracked in `bd` as epic `multifer-9qu` with children `9qu.1` … `9qu.10`; Phase 1 epic `multifer-kxn` is blocked on the Phase 0 epic.

### Phase 1 — minimal working engine  **[COMPLETE — 2026-04-13]**
- **Shapes:** `oneblock` + `cross` only.
- **Engine:** refit-first (slow but universal path from §22.10).
- **Outputs:** component significance, variable stability, score stability, subspace stability.
- **Explicitly excluded:** variable significance, multiblock, geneig.

**Delivered (Phase 1):**
- `R/alignment.R` — `match_components` (greedy permutation matching), `align_sign`, `align_procrustes`, `principal_angles` (prefers `multivarious::prinang` when installed), `align_loadings` pipeline. Permutation matching runs *before* sign / Procrustes per the Part 5 review.
- `R/unit_formation.R` — `form_units()`. **Conservative default**: every root is its own singleton component; near-tie grouping is opt-in behind `group_near_ties = TRUE`.
- `R/mc_pvalue.R` — `mc_p_value()` Phipson-Smyth `(1+r)/(B+1)`, fixed B (sequential stopping is Phase 1.5).
- `R/ladder.R` — shape-agnostic `ladder_driver()` shared scaffold.
- `R/engine_oneblock.R` — `run_oneblock_ladder()` implements the explicit Vitale P3 simplification (Part 1 §1) with no second SVD per null draw.
- `R/engine_cross.R` — `run_cross_ladder()` reuses `ladder_driver` for both `(cross, covariance)` and `(cross, correlation)` recipes; null action permutes rows of Y only.
- `R/bootstrap.R` — `bootstrap_fits()`. **Paired** row bootstrap for cross. Aligns each rep's loadings + scores via `match_components` then `align_sign` / `align_procrustes`.
- `R/variable_stability.R` — `variable_stability_from_bootstrap()` returns a STABILITY measure (not a p-value).
- `R/score_stability.R` — `score_stability_from_bootstrap()`. **Design lock**: scores are computed by projecting the ORIGINAL observations through each rep's fit, never on the bootstrap sample.
- `R/subspace_stability.R` — principal-angle summaries per unit.
- `R/infer.R` — top-level `infer()` dispatcher; populates `infer_result` with units, component_tests, variable/score/subspace stability, assumptions, mc reproducibility log, cost (wall-time per phase), provenance.
- `R/adapter_oneblock_baser.R` — `adapter_svd()` and `adapter_prcomp()` (no-dep oneblock references).
- `R/adapter_cross_baser.R` — `adapter_cross_svd()` (the *one real* cross reference declaring BOTH covariance and correlation; triggers strict-dispatch ambiguity errors when relation is omitted) + `adapter_cancor()` (correlation only).
- `R/adapter_multivarious_pca.R` — Tier-2 wrapper for `multivarious::pca`, soft dep, registers only when multivarious is installed.
- `R/adapter_multivarious_plsc.R` — Tier-2 wrapper for `multivarious::plsc`, single-relation (covariance only).
- `tests/testthat/helper-default-adapters.R` — `ensure_default_adapters()` test helper that re-registers the auto-loaded adapters after other test files clear the registry.
- `tests/testthat/test-integration-phase1.R` — end-to-end integration on the 4 frozen Phase 0 bench suites + multivarious round-trips.
- ~140 new test_that blocks across Phase 1 modules; `R CMD check → Status: OK` (0 errors / 0 warnings / 0 notes).
- Tracked in `bd` as epic `multifer-kxn` with children `kxn.1` … `kxn.16` (all closed).

### Phase 1.5 — fast path  **[COMPLETE — 2026-04-13]**
- `core()` / `update_core()` for `oneblock` and `cross`.
- Sequential Monte Carlo stopping as **global budget allocator** (§28).
- Caching, parallel execution (`mirai` default from §30).

**Delivered (Phase 1.5):**
- `R/linear_operator.R` — minimum operator surface (`apply`, `apply_t`, `dim`, optional `gram` and `core_update`). Dense reference wrapper via `as_linear_operator.matrix()`. (§30 rank 1.)
- `R/thin_svd_cache.R` — per-call xxhash64-keyed cache with `cached_svd()`, `svd_cache_new()` / `_reset()` / `_hits()` / `_rate()`. `infer()` resets the package-default cache at every call and reports the hit rate through `result$cost$cache_hits`.
- `R/adapter_oneblock_baser.R` — `adapter_svd()` and `adapter_prcomp()` now expose `core()` and `update_core()` via the Fisher-style `(U[idx] - 1·ū^T)·diag(d)·V^T` identity. Inner SVD is n × k; return shape matches `refit()`.
- `R/adapter_cross_baser.R` — `adapter_cross_svd()` (covariance branch only) exposes `core()`/`update_core()` via `X[idx]^T Y[idx] = V_x M V_y^T` with `M = D_x (U_x[idx]^T U_y[idx]) D_y`. Correlation branch falls back to `refit` (per-replicate QR rescaling does not admit a clean k × k update).
- `R/mc_sequential.R` — `mc_sequential_bc()` implements classical Besag-Clifford early non-rejection on the fixed-B Phipson-Smyth primitive. Default `h = ceiling(α·(B_max+1))`.
- `R/budget_allocator.R` — `mc_budget_allocator()` state machine: global `B_total`, optional `per_rung_cap`, rolling pool that accumulates unused budget across rungs. Exposed via `checkout() / record_use() / remaining_fn() / used_fn()`.
- `R/ladder.R` — `ladder_driver()` now builds an allocator and drives each rung through `mc_sequential_bc()` instead of `mc_p_value()`. Returns `allocator`, `total_draws_used`, `batch_schedule`, `stopping_boundary`.
- `R/parallel.R` — `multifer_parallel_init()` (lazy idempotent daemon pool), `multifer_parallel_shutdown()`, `multifer_task_seeds()` (deterministic per-task seed generator that preserves the caller's RNG stream), `multifer_parallel_lapply()` with `backend = c("auto", "mirai", "sequential")`.
- `R/bootstrap.R` — per-replicate body refactored into a pure closure `rep_fn()` that the parallel backend fans out. New knobs: `parallel`, `fast_path`, `core_rank = 50L`. `core_rank` is the truncation cap that converts the thin-SVD cache into wall-clock wins; without it the inner n × k SVD matches the refit cost.
- `R/infer.R` — consolidated knobs: `B`, `B_total`, `mc_batch_size`, `parallel`, `fast_path`. Populates `result$cost` fully (`full_data_ops`, `core_updates`, `mc_budget_spent`, `cache_hits`, `wall_time_phases`) and `result$mc` fully from the allocator output (`stopping_boundary`, `batch_schedule`, `stop_iteration`, `total_draws_used`, `exceedance_counts`).
- `DESCRIPTION` — `digest` and `mirai` promoted to `Imports`.
- `tests/testthat/test-integration-phase15.R` — decision-agreement tests across all four Phase 0 bench suites comparing `fast_path = "off"` vs `"auto"`. Plus Wave-specific tests: `test-linear-operator.R`, `test-thin-svd-cache.R`, `test-core-update-oneblock.R`, `test-core-update-cross.R`, `test-mc-sequential.R`, `test-budget-allocator.R`, `test-parallel-bootstrap.R`.
- Full test suite: 1146 passing, 0 failures. `R CMD check → Status: OK` (0 errors / 0 warnings; only the system-time note remains).

**§40 benchmark status (measured on the reference machine, 2026-04-13, after commit `8139fab`):**

| Suite                                   | Refit  | Fast  | Speedup      | §40 target |
|-----------------------------------------|--------|-------|--------------|-----------|
| oneblock medium (n=300, p=80)           |  0.60s | 0.42s | 1.41×        | ≥ 5×      |
| oneblock large  (n=500, p=300)          |  6.36s | 1.26s | **5.03×** ✓  | ≥ 5×      |
| cross cov medium (n=300, p=50, q=40)    |  0.38s | 0.42s | 0.90×        | ≥ 5×      |
| cross cov large  (n=500, p=250, q=200)  |  4.84s | 1.19s | **4.05×**    | ≥ 10×     |

Optimizations applied across Phase 1.5:
- Thin-SVD cache for the original fit (§30 rank 1, Wave 1).
- `core()` / `update_core()` fast paths for oneblock and cross-covariance (§17, Waves 2–3).
- Besag-Clifford sequential MC + global budget allocator (§28, Wave 4).
- `core_rank = 50L` truncation cap (Wave 6) that turns the thin-SVD cache into a real speedup in the bootstrap loop.
- Partial-SVD null draws via `RSpectra::svds` in `null_stat_fn` and `observed_stat_fn` (§30 rank 5, commit `91cc382`).
- Partial-SVD deflation via `top_svd(X, 1)` in both engines' `deflate_fn` (commit `34ae10c`).
- `matrixStats::rowQuantiles` in `.row_quantile_pair` (commit `34ae10c`).
- Rank-capped cross null draw at rung 1 via cached thin SVDs of both blocks (commit `8139fab`): replaces per-draw `crossprod(Xc, Yc[perm, ])` with a `k_x x k_y` inner core `Dx (Ux^T · Uy[perm, ]) Dy`, using the identity `top SV(Vx M Vy^T) = top SV(M)` for orthonormal `Vx / Vy`. Only runs on rung 1; deflation semantics at rung > 1 fall back to the exact path. Calibration holds because both `observed_stat_fn` and `null_stat_fn` use the same rank-capped core at rung 1.

**Absolute fast-path time arc for cross covariance large (n=500, p=250, q=200, R=30, B=99):**
- Wave 6 baseline (thin-SVD cache + core-update + matrixStats): 2.30s, speedup 2.73×
- + partial-SVD null draws (`91cc382`):                            2.00s, speedup 2.95×
- + matrixStats row quantiles + partial-SVD deflation (`34ae10c`): 1.97s, speedup 2.95×
- + rank-capped rung-1 cross null (`8139fab`):                     1.19s, speedup **4.05×**

**Why §40 hard targets are still not fully met:**
- **Medium cases (p ≤ 80, q ≤ 40)** are dominated by fixed costs (original fit, ladder rung setup, stability consumer overhead). At these sizes `RSpectra::svds` doesn't beat `base::svd` because `k` is not meaningfully less than `min(n, p)`, and the rank-capped rung-1 cross path is inactive because `cross_rank_cap = 50 ≥ min(p, q)` — no truncation happens. Hitting 5× here requires amortizing fixed costs, not inner-loop tuning.
- **Cross covariance large (≥10×)** improved from 2.95× to 4.05× with the rank-capped rung-1 path. The remaining bottlenecks are distributed: `score_stability_from_bootstrap` (~30%), `original_fit` / `thin_svd_cache` warmup (~15%), miscellaneous bootstrap+alignment (~20%). Closing to 10× requires either extending the rank-cap approximation to the bootstrap loop (which would invalidate the exact-decision-agreement property of the refit path) or refactoring `score_stability_from_bootstrap` to avoid per-observation sorts. Both are deferred to follow-up work.

The core-update mechanism is **correct** (decision agreement verified on all four Phase 0 bench suites). The large oneblock case **meets** the §40 medium target. Large cross covariance is now much closer (4.05× vs 10×). Medium cases and cross ≥10× remain follow-ups.

- Tracked in `bd` as epic `multifer-3c0`.

### Phase 2 — scope expansion
- `multiblock` shape.
- `geneig` shape (contrastive / discriminant / metric-weighted).
- Richer target types.

### Phase 3 — richer inference
- Variable significance **where the adapter explicitly supports it**.
- Multiple-testing correction across units + variables.
- Reporting / plotting, external adapter gallery.

**Why the renumbering.** Old §22.17 conflated "ship v1" with "build the fast engine." Splitting Phase 1 from Phase 1.5 lets v1 be correct and shippable **without** the core-update path — which is exactly when benchmarks of core-update vs. refit become meaningful.

---

## 40. Lock the benchmark plan into the PRD now

One of the gaps called out in the post-Part-3 review. Four benchmark suites, **frozen during Phase 0**:

1. **One-block null calibration.** Pure-noise matrices. Targets: over-selection rate, p-value calibration across $\alpha \in \{0.01, 0.05, 0.10\}$.
2. **One-block power / shadowing.** Weak later roots hidden behind strong earlier ones — exactly where deflation and P3 matter most in Vitale's own simulations.
3. **Cross-block null calibration.** Paired $(X, Y)$ with no cross-association but strong **within-block** structure. Catches the exact failure mode where naïve permutation CCA inflates error beyond the first root.
4. **Speed / agreement.** Naive refit vs. core-update, measuring both **wall-clock** and **decision agreement within Monte Carlo error**.

**Target numbers (locked into the PRD):**
- Core-update produces the **same decisions as refit within Monte Carlo error**.
- **≥ 5× speedup** on medium problems.
- **≥ 10× speedup** on large cross-block problems once core updates land.

---

## 41. Sequential Monte Carlo reproducibility — explicit log

Anytime-valid stopping is non-deterministic in runtime. The result object must log **everything needed to reproduce a run** even though wall time is data-dependent:

```r
result$mc
  $rng_seed
  $rng_kind
  $stopping_boundary    # the anytime-valid test's threshold schedule
  $batch_schedule       # draws per batch
  $stop_iteration
  $total_draws_used
  $exceedance_counts    # per component / per unit
```

With this, a run is **reproducible bit-for-bit** even though its cost profile is not.

---

## 42. Revised public API

Simpler, stricter surface:

```r
# Primary
infer(model, data = NULL, recipe = NULL, design = NULL, strict = TRUE, ...)

# Convenience wrappers
infer_oneblock(model, X, ...)
infer_cross(model, X, Y, relation, design = paired_rows(), ...)

# Thin method aliases
infer_pca(model, X, ...)
infer_pls(model, X, Y, ...)
infer_cca(model, X, Y, ...)

# Advanced (public)
infer_recipe("oneblock")
infer_recipe("cross", relation = "covariance")
infer_recipe("cross", relation = "correlation")

# Package-author extension
register_infer_adapter("my_model", adapter_from_projector(...))
```

**`infer_spec()` becomes internal.** End users never see it. Most package authors don't need it either — the adapter constructor is enough.

**`strict = TRUE` is the default.** Turning it off is allowed but requires a deliberate flag, and downgrades the validity level in the result object.

---

## 43. Three-tier pluggability

Make the cost of integration match the desired fidelity:

**Tier 1 — coercion only.**
```r
as_projector(x)         # or as_cross_projector(x)
```
The package gets the built-in refit engine and standard recipes. Zero adapter code.

**Tier 2 — adapter.**
Custom `roots`, `residualize`, `align`, and `null_action` selection. Works with the slow / universal engine from §22.10.

**Tier 3 — fast adapter.**
Also provides `core()` and `update_core()`. Gets the full §17 / §27 speedup.

**Why this matters.** Outside packages need a **path to "it works today"** that doesn't require learning the entire internals. Tier 1 is deliberately trivial.

---

## 44. Revised result schema — unit-centered

Replaces the flat schema in §22.12:

```r
res$units              # latent units (see §37)
res$component_tests    # one row per unit
res$variable_stability # keyed by unit_id
res$score_stability    # keyed by unit_id
res$subspace_stability # angles, alignment method, stability label per unit
res$assumptions        # declared + checked
res$mc                 # reproducibility log (§41)
res$cost               # from Part 4 §32
res$provenance         # adapter id, version, capability flags, call
```

The Part 3 schema fields are still present; they are **re-keyed** from raw component index to `unit_id`. Backward-compatible conceptually, cleaner semantically.

---

## 45. Name, tagline, bottom line

- **Concept name:** ShapeInfer (working title; "typed shape inference" is the conceptual phrase).
- **Package name:** `multifer` (unchanged).
- **Tagline:** *typed perturbation inference for projector-based multivariate models.*

The tagline says exactly what is new: **not just resampling, not just projectors, but a typed system that is strict about inferential meaning.**

**Bottom line.** The Part 3 review was right. The response is not a small edit. It is to tighten the framework around:

1. Typed shapes — not raw shapes.
2. Strict dispatch — not silent auto-magic.
3. Machine-enforced validity — not printed fields.
4. Latent units — not per-component rows.

Smaller, safer, more publishable at the same time.

---

## 46. Migration notes (Parts 3 → 5)

For future reference, the changes Part 5 makes to Part 3's PRD:

| Part 3 | Part 5 |
|---|---|
| `shape` as inferential object | **typed shape** = `(geometry, relation, design)` |
| `infer(mod)` with best-effort auto-dispatch | `infer(mod, strict = TRUE)` — errors on ambiguity |
| `validity` as a printed field | `validity` as a registration / compile / run-time **gate** |
| Result tables keyed by component index | Result tables keyed by **`unit_id`** |
| Phase 1 includes variable significance | Variable significance **deferred to Phase 3** |
| Single Phase 1 with engine + fast path | **Phase 0** (contracts) → **Phase 1** (refit-first) → **Phase 1.5** (core-update) |
| Benchmarks as Phase 4 concern | **Benchmarks frozen in Phase 0** with locked targets |

## 47. Updated TODO (Part 5)

- **Design vocabulary.** Concrete list of `design` objects: `exchangeable_rows()`, `paired_rows()`, `blocked_rows(groups)`, `nuisance_adjusted(Z)`. What's the minimum useful set for v1?
- **Capability matrix format.** Machine-readable grid of `(recipe, target) → supported?` that `infer_recipe()` compiles against.
- **Checked assumptions library.** What can actually be tested at run-time? Symmetry (for sign-flip), exchangeability (for permutation), rank sanity. Which ones have cheap tests and which require the user to declare?
- **Unit formation rule.** When do tied / near-tied roots collapse into a subspace unit? Threshold on eigengap? Confidence interval on eigengap? Adapter-declared?
- **Tier 1 coercion surface.** Exact minimum required by `as_projector()` / `as_cross_projector()` to get useful refit-based inference with zero extra adapter code.
- **Phase 0 exit criteria.** Formal gate: schema frozen + benchmarks frozen + two synthetic examples passing + adapter validation implemented. No engine code until this gate passes.
- **`strict = FALSE` semantics.** What exactly does relaxing strict dispatch do? Does it pick a default, or emit a warning and degrade the validity level?
- **Reconcile `ShapeInfer` / `multifer` naming in docs.** Decide whether "ShapeInfer" survives as a subtitle, a vignette title, or disappears entirely in favor of `multifer`.
