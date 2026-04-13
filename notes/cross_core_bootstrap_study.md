# Cross Core Bootstrap Study

This note tracks the exact agreement study for cross-covariance / PLSC bootstrap resampling.

## Purpose

The package already uses a core-space fast path for cross-covariance bootstrap refits through the identity

\[
X^\top W Y = V_x \left[ D_x (U_x^\top W U_y) D_y \right] V_y^\top.
\]

The question for this study is narrower than method comparison:

- does the core-space bootstrap path reproduce the raw refit path on the same resamples?
- do lifted loadings and scores agree after alignment?
- do downstream stability summaries agree?
- what speedup does the exact fast path deliver across shape and SNR regimes?

This is an implementation-validity study, not a new inferential method.

## Scenario grid

The study script covers:

- balanced separated roots
- balanced near-tied roots
- wide `Y` with higher `Y` noise
- wide `X` with higher `X` noise
- small-`n`, wide, near-tied blocks
- larger balanced blocks

Each scenario simulates paired blocks from a shared latent score matrix plus asymmetric block noise.

## What is compared

For each scenario, the script runs `bootstrap_fits()` twice on the same paired bootstrap draws:

- exact core-space fast path: `fast_path = "auto"`
- raw refit path: `fast_path = "off"`

Both runs use:

- `relation = "covariance"`
- `method_align = "sign"`
- identical seeds
- `core_rank = NULL` so the fast path remains exact rather than truncated

## Agreement metrics

Per replicate:

- leading singular values for the evaluated latent units
- centered cross-operator reconstruction
- centered block means
- aligned loadings for `X` and `Y` on the evaluated units
- aligned scores for `X` and `Y` on the evaluated units

Downstream summaries:

- variable stability
- score stability
- subspace stability

The study also records elapsed runtime and the resulting speedup. In wide settings, the raw refit path can retain a long numerical tail beyond the evaluated signal units; the study records that separately instead of treating it as a substantive disagreement.

## Script

Run:

```sh
Rscript tools/cross_core_bootstrap_study.R --n-boot=80 --outdir=tools/results/cross_core_bootstrap
```

The script writes:

- `cross_core_bootstrap_summary.csv`
- `cross_core_bootstrap_per_rep.csv`

## Preliminary findings

Run 1 used:

```sh
Rscript tools/cross_core_bootstrap_study.R --n-boot=24 --outdir=tools/results/cross_core_bootstrap_run1
```

Main result:

- for the leading signal units, the exact core-space bootstrap path matches the raw refit path to machine precision across all tested shape and SNR regimes

Observed maxima from run 1:

- singular values: about `2.2e-12`
- reconstructed cross operator: about `1.4e-12`
- aligned loadings: about `5.8e-14`
- aligned scores: about `1.1e-12`
- variable stability summaries: about `2.4e-14`
- score stability summaries: about `4.5e-13`
- subspace stability summaries: about `6.5e-15`

Interpretation:

- the core-space path is an exact implementation of the same bootstrap contract for the leading latent units
- the only persistent representation difference is that the raw refit path may keep a longer numerical tail in wide problems, while the core-space path naturally works in the active latent rank

Performance from run 1 was mixed:

- the exact fast path was not uniformly faster at these moderate sizes
- measured speedups ranged from about `0.60x` to `1.37x`

That is consistent with the current package benchmarks: exactness is already there, but substantial wall-clock gains need larger problems or further structural optimization beyond simply switching from raw refit to the current core-space bootstrap path.
