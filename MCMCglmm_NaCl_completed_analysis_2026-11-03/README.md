# MCMCglmm NaCl among-RIL covariance analysis — 3 November 2026

## What is independently checked

- Three completed 07a final unstructured chain files, each with **2,000 retained MCMC samples** and **112 covariance/variance columns** (including duplicated symmetric entries), totaling 6,000 samples.
- The 112 parameter diagnostics: max Rhat 1.00300 and minimum effective sample size 5,383.48. These indicate good computational mixing **conditional on the specified model and priors**, not validation of biological identification.
- All 6,000 among-line G posterior draws are symmetric and positive definite. They were decoded from the supplied binary RDS checkpoints and reproduce the provided original-unit G and correlation matrix to numerical precision.
- `07c_G_correlation_heatmap.png`, `07c_G_eigenvalues.png`, and `07c_selected_chain_traces.png` provide diagnostic figures.
- The 7-trait correlation matrix has eigenvalues about 3.405, 1.268, 0.885, 0.585, 0.551, 0.245, and 0.062. The top axis accounts for 48.6% of the summed correlation eigenvalues. Posterior first-axis fraction 48.8%, 95% credible interval 44.0%–53.3%.
- Effective dimensionality used in the project's script 07b: `n - mean((lambda-1)^2)`. From mean correlation G, 5.901; across posterior draws, mean 5.879, 95% interval 5.641–6.104. This is **not** identical to the participation ratio (posterior mean 3.312, 95% interval 2.968–3.692). Do not interchange the measures.
- Eighteen of 21 between-trait correlation 95% marginal posterior intervals exclude zero (a descriptive count, not adjusted for multiplicity).

## Fitness-direction comparison — exploratory

The analysis transforms the original-unit among-line G into BLUE-standardized units using `BLUE_sd` from `selection_gradients_NaCl.csv`, then computes `G * beta` with fixed estimated multivariate fitness gradients, as a quantitative-genetic projection. Fixed-beta cosine between the selection vector and posterior predicted response averages 0.820 (95% posterior interval 0.748–0.888), roughly a 35-degree angle. This DOES NOT incorporate selection-gradient uncertainty or imply a validated evolutionary prediction for the evolved populations. See `07c_selection_response_alignment.png` and `07c_trait_covariance_and_response.csv`.

## Reconstruct the final RData file *without running MCMC again*

Unzip your original `RData_G_NaCl.zip` into the project root, leaving chain checkpoints in `RData_G_NaCl/`. Place `07c_Recover_G_posterior_no_refit.R` in the root and run:

```r
source("07c_Recover_G_posterior_no_refit.R")
```

This creates `RData_G_NaCl/07a_multivariate_G_final.RData`, with the objects expected by existing script 07b. It independently validates the saved matrices against your published CSV exports before saving.

## Complete the phenotype pleiotropy scaling (script 07b)

The available uploads do not contain the original script-05 output `pleiotropy_results/R5/GWAS_clusters_with_pleiotropy_R5.RData`, which holds the QTL significant-trait sets plus Manhattan/Euclidean effect distances and chromosome arm/center labels. To finish the exponent and arm-versus-center comparison, please upload that file or run:

```r
source("07b_Refit_pleiotropy_with_multivariate_G_FIXED.R")
```

Use the fixed script from our previous exchange (also included in this package), with the R5 QTL .RData in place. It creates `pleiotropy_results/R5/multivariate_G/Scaling_exponents_multivariateG.csv`, `Scaling_interactions_multivariateG.csv`, `Scaling_G_uncertainty_sensitivity.csv`, `QTL_pleiotropy_multivariateG.csv`, and a PDF comparison figure.

### Interpretation safeguards

1. **Among-RIL covariance, not proven additive G:** the random effect is `us(trait):pop_label`; no genomic relatedness kernel was fitted in this MCMC model.
2. **Model sensitivity not yet assessed:** one proper inverse-Wishart specification and unstructured residual model were run, despite excellent chain mixing. Inspect prior sensitivity, residual model, and possible date/session overlap before claims of strong constraints.
3. **Neff uncertainty is conditional on GWAS QTL trait membership.** It excludes uncertainty in the locus-effect estimates, significance-based trait inclusion, and phenotypic fitness gradient estimation.
4. **Do not treat Gβ as an out-of-sample evolutionary forecast.** It is an expected response to fixed, estimated phenotypic selection coefficients, and the model is for among-line covariance.
