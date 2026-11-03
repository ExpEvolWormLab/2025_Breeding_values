# Independent audit of 07b multivariate-G pleiotropy scaling

3 November 2026. Inputs: user-supplied `multivariate_G(1).zip`, especially `QTL_pleiotropy_multivariateG.csv`, `Scaling_exponents_multivariateG.csv`, `Scaling_G_uncertainty_sensitivity.csv`, and script `07b_Refit_pleiotropy_with_multivariate_G_FIXED_v2.R`.

## Reproducibility
- 1,089 R5 clusters; 194 without a significant association among the seven individual traits (Neff = NA), 895 with >= 1; after excluding tips, 881 retained for arm/center scaling (717 arms, 164 centers).
- Manhattan and Euclidean distances recomputed exactly from trait GWAS effect columns masked by significance indicators (maximum absolute differences < 4e-15 and 5e-16 respectively).
- Neff_RIL versus Neff_BLUE across 895 associated clusters: Pearson r = 0.9956, median absolute difference = 0, 95th percentile absolute difference = 0.107.
- Reported scaling coefficients reproduced using ordinary least squares on log distances versus log Neff with region interaction.

## Reported descriptive scaling
| Metric | Arm exponent (95% regression CI) | Center exponent (95% regression CI) | Reference null |
| --- | --- | --- | --- |
| Manhattan | 1.606 (1.513–1.699) | 1.702 (1.544–1.860) | 1 |
| Euclidean | 0.982 (0.892–1.073) | 1.105 (0.951–1.259) | 0.5 |

Center versus arm interactions: Holm-adjusted p=0.354 for both metrics (not demonstrated).

Posterior G-only re-fit across 200 posterior draws (2.5–97.5 percentiles): Manhattan arm 1.539–1.672, center 1.622–1.779; Euclidean arm 0.943–1.024, center 1.053–1.156. These **are not full regression/posterior confidence intervals**.

## Null-model concern
The reference slopes 1 and 0.5 arise from expectations based on the number of contributions (number of traits) under simplifying independent, exchangeable-effect assumptions. Here the predictor is `Neff` calculated from selected-trait correlation, and Neff is lower than the number of observed significant traits; an apparent excess over 1 or 0.5 can be generated without an interaction among the causal trait effects. The significant trait associations are also selected by GWAS, which further distorts the distribution of effect magnitudes. Manhattan and Euclidean statistics discard effect signs, so they cannot directly detect the alignment of trait effects with selection gradients or prove a fitness interaction.

We conducted a **diagnostic**, not an exact inferential null: for each trait and chromosome-region class, permute significant GWAS effect estimates across eligible clusters, keeping the binary associated-trait membership matrix fixed; refit the same exponent. This preserves per-trait distributions and trait-set membership but breaks across-trait co-localization of effect magnitudes. It does **not** reproduce LD, effect-estimate uncertainty, SNP-level GWAS discovery, locus-dependence or conditional significance constraints.

With 1,500 diagnostic permutations (see attached CSV), the arm slopes are not unusually steep relative to this conditional randomization. The full-data center slopes are above this diagnostic null, but this is sensitive to the focal regions and top effect-norm outliers.

## Focal-region sensitivity
Remove these *exploratorily defined* genomic intervals: I 5,883,574–7,996,414; IV 8,750,000–10,600,000 (covers IV:8,892,663 and IV:10,499,135); V 6,000,000–8,000,000. Excludes 27 center clusters; n=854; centers=137.

| Metric | Full arm | Full center | After excluding all three: arm | After excluding all three: center |
| --- | ---: | ---: | ---: | ---: |
| Manhattan | 1.606 | 1.702 | 1.606 | 1.606 |
| Euclidean | 0.982 | 1.105 | 0.982 | 0.977 |

In the same exploratory diagnostic randomization, after excluding all three regions, center exponents are not unusually high (one-sided conditional permutation p approximately 0.14 for Manhattan and 0.26 for Euclidean). This is not a preregistered or fully calibrated claim.

Trimming the highest 1% Manhattan and Euclidean distances (11 clusters) similarly attenuates center slopes to 1.617 (Manhattan) and 0.988 (Euclidean). Many excluded clusters lie in the three known focal regions, including IV:10499135, I:7311528, and V:6637695.

## Manuscript language
Supported: `The multivariate-G analysis confirms a pronounced super-reference scaling of the summed magnitudes of trait-associated genetic effects with effective trait dimensionality, conditional on GWAS discovery, and highlights especially large multivariate effects in a small number of center regions.`

Not yet supported: `Pleiotropic causal mutations synergistically enhance adaptation genome-wide`, or `synergy is stronger in centers than arms`, or `the exponent tests causally identify pleiotropy rather than LD`.

Recommended: derive a GWAS/LD/trait-covariance calibrated null using RIL genotypes and the actual association procedure. Model founder haplotypes and QTL dependency in the three focal regions; report sensitivity to these regions. Retain 07b as a supportive genetic-architecture result, not a primary causal explanation of the 802-cluster evolutionary prediction.
