# Procrustes Bias Study

This note defines the simulation study used to compare legacy Procrustes
alignment against the package's preferred stability targets:

- sign alignment for well-separated components,
- subspace summaries for tied or near-tied components.

The goal is not to prove that Procrustes is always bad. The goal is to
measure when it changes the estimand in a way that depends on shape,
eigengap, and asymmetric block noise.

## Questions

1. Does McIntosh-style Procrustes rotation distort component-strength
   summaries under permutation null draws?
2. Does Procrustes shrink bootstrap uncertainty for axis-level loadings
   relative to sign-only alignment?
3. Are any such distortions amplified when:
   - the first two roots are near-tied,
   - one block is much wider than the other,
   - one block is much noisier than the other,
   - sample size is small relative to block width?

## Scenarios

Each scenario is a two-block latent-factor model with two shared
components. The study varies:

- `n`
- `p_x`, `p_y`
- root gap via the signal scales
- noise level in `X`
- noise level in `Y`

The default study grid is intentionally small but discriminating:

- balanced, separated roots
- balanced, near-tied roots
- wide-`Y`, noisy-`Y`, near-tied roots
- wide-`X`, noisy-`X`, near-tied roots
- small-`n`, wide blocks, near-tied roots
- large-`n`, balanced, separated roots

## Outputs

The tool script writes two CSV summaries.

### 1. Null distortion summary

This compares raw singular values from permuted cross-operators against
McIntosh-style rotated singular values obtained from the column norms of
`D Q`.

Key metrics:

- `lv1_mean_shift`
- `lv2_mean_shift`
- `tail_ratio_mean_shift`
- `mixing_mean`

These quantify how much Procrustes changes the latent-root statistic
itself under the null.

### 2. Bootstrap alignment summary

This compares `method_align = "sign"` and `method_align = "procrustes"`
for bootstrap replicates.

Key metrics:

- `truth_angle1_mean`, `truth_angle2_mean`
- `reference_angle1_mean`, `reference_angle2_mean`
- `offdiag_share_mean`
- `subspace_angle_max_mean`
- `loading_sd_mean`

Interpretation:

- If `reference_angle_*` improves while `truth_angle_*` does not, the
  method is locking bootstrap draws to the sample reference rather than
  recovering the generating axes.
- If `loading_sd_mean` is much smaller under Procrustes while
  `subspace_angle_max_mean` is unchanged, the method is suppressing
  axis-level uncertainty without changing subspace uncertainty.
- If `offdiag_share_mean` is large in near-tied settings, individual
  axes are not stable targets and subspace summaries are the right
  object to report.

## Running the study

From the repository root:

```sh
Rscript tools/procrustes_bias_study.R
```

Optional overrides:

```sh
Rscript tools/procrustes_bias_study.R --n-perm=80 --n-boot=80 --outdir=tools/results/procrustes_bias
```

The script is deterministic given its seeds and is meant as a methods
comparison harness, not a package test.
