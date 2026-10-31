# Frequency-class calibration + genomic prediction (R scripts 38–42)

## Goal

This package follows the corrected A00/Pool-seq sample alignment and the `33`–`36` frequency-stratified results. **GA1, GA2 and GA4 are the primary, surviving-evolution populations; GA3 is confirmed extinct and appears only in supplementary sensitivity analyses.**

It contains two separate, complementary analyses:

1. Site-specific neutral Wright–Fisher simulation calibration for rare alleles, plus **one disjoint ancestral-MAF class for each R5 QTL**.
2. Eight-phenotype rrBLUP GBLUP held-out prediction, with comparison to the existing **training-ranked top-75-SNP** predictor.

## Setup

Copy all `.R` scripts to the project root (the same working directory from which scripts 33–36 and 16 ran). Keep `09f_drift_functions.R` available either in that root or at `drift_nullfit_update/09f_drift_functions.R`. Also copy `33_Frequency_drift_helpers.R` from the preceding phase into the project root. Install `rrBLUP` if needed:

```r
install.packages('rrBLUP')
```

Paths are editable at the top of each script. All scripts write to new output directories and do not overwrite previous phases.

## Run order

```r
source('42_Synthetic_smoke_tests.R')
source('38_Calibrate_rare_WF_by_site.R')
source('39_Independent_cluster_frequency_prediction.R')
source('40_rrBLUP_GBLUP_CV.R')
source('41_GBLUP_comparisons_inference_figures.R')
```

The two branches, scripts 38–39 and 40–41, are independent after their initial inputs. You can run them in parallel R sessions on separate nodes. Do not simultaneously modify the underlying prepared RData files.

## 38 — What rare calibration does, and does not do

Reads the exact corrected SNP-by-pool read counts, one A00 reference, and the confirmed survivor status table. For 16 randomly selected SNPs per MAF bin in each scenario, simulates neutral discrete-generation WF allele-frequency paths with the site-specific A00 read posterior and observed sequencing depths. The simulations show how often the original approximate normal statistic crosses nominal thresholds. The independently selected calibration SNPs are **sampled without inspecting their selection p-values**.

It additionally computes finite-resolution conditional WF Monte Carlo p-values for the six most extreme existing founder-GWAS SNPs per low-frequency bin, plus any predeclared sites from `34_Predeclared_target_registry.csv` that occur in the founder-GWAS SNP universe. These candidate p-values are exploratory because the variants were searched for large observed statistics, and cannot be interpreted as genome-wide discoveries without accounting for all SNPs tested.

Default calibration uses `Ne=1000`. Candidate sensitivity runs use `Ne=500, 2000`. Change `sim_random`, `sim_candidates` and `random_sites_per_class` at the top of the script to increase precision. **No number of routine Monte Carlo runs can alone validate a genome-wide BH tail near 1e-7; calibration to that range would require a dedicated tail-probability method and careful ascertainment control.**

* `38_Sitewise_conditional_WF_results.csv`: site-specific MC p, null distribution SD, approximate-p false positive diagnostics; labels random/predeclared/post-selected.
* `38_Null_calibration_by_frequency.csv`: mean neutral false-positive fraction across independently sampled sites for each MAF bin.
* `38_Calibration_vs_nominal.pdf`: descriptive calibration comparison.
* `38_Sites_*.csv`: records site selection role and baseline MAF class.

Output: `selection_results/frequency_strata/rare_calibration/`.

## 39 — Mutually exclusive ancestral-frequency QTL classes

Every R5 cluster is assigned a **single class using the median ancestral MAF across its tested phase-mapped SNPs** (independent of their observed evolutionary responses). Its response is the median phase-aligned slope across **all** those SNPs, computed once for that cluster. No R5 cluster is double-counted across frequency classes. This differs deliberately from script 36, whose bin-specific clusters can reappear in multiple bins.

The four predictors (full index, three-trait index, SB-only, directly mapped beta) are the unchanged RIL GWAS effect predictors. Leave-one-chromosome-out models are fitted within each class. The matched full-index versus reduced/SB MSE differences use 2 Mb genomic-block bootstrap intervals. Class-specific results are suppressed below 15 clusters or 3 chromosomes. The `ALL` group is a benchmark, not a new MAF class.

* `39_Disjoint_cluster_membership_QC.csv`
* `39_Prediction_by_unique_MAF_class.csv`
* `39_Cluster_predictions_and_residuals.csv`
* `39_Paired_block_bootstrap_full_vs_comparators.csv`
* `39_Unique_frequency_class_prediction.pdf`
* `39_Frequency_class_heterogeneity_block_bootstrap.csv`: matched-block uncertainty comparing the full-index correlation in each class with the 10–25% class.

Output: `selection_results/frequency_strata/unique_cluster_classes/`.

## 40 — GBLUP from rrBLUP with exact train/test separation

Uses adjusted BLUEs from `RData/Preparatory_GLMM_NaCl_only_withBeta.RData`. Reconstructs beta **on the BLUE scale** from the pre-estimated seven fitness gradients; these weights are fixed and are not re-estimated from test-fold evolution. Genotype matrix = `genotype_pruned_matrix` by default, a phenotype-independent genome-wide LD-pruned input. Set `marker_set <- 'full'` for a more expensive marker-density sensitivity test.

The kinship matrix is constructed with **training-fold genotype means** and then used as a cross-relationship matrix for held-out predictions. For fold-specific training phenotypes only, `rrBLUP::mixed.solve` estimates the REML variance ratio. Prediction of unseen test RILs is `mu + K_test,train (K_train,train + Ve/Vu * I)^-1 (y_train - mu)`. This is an out-of-fold GBLUP, not the in-sample line-effect BLUP.

If present, `16_RIL_folds.csv` is used unchanged, allowing direct paired comparisons to the top-75-SNP results. Optionally, `19_RIL_source_populations.csv` enables validation in which entire source-population groups are held out; if any origins are missing, this scheme is deliberately skipped rather than classifying all unknown lines as a fictitious group.

Outputs in `phenotype_prediction/GBLUP_rrBLUP/`:

* `40_Out_of_fold_GBLUP_predictions.csv`: one row per RIL, fold and trait, plus fold training-mean baseline.
* `40_GBLUP_predictive_performance.csv`: held-out R2, Pearson r, MSE and improvement over fold training mean.
* `40_Foldwise_variance_components.csv`: training-fold Vu/Ve and a **pseudo-PVE of adjusted BLUE variation**, not strict additive heritability.
* `40_RIL_folds_and_origins.csv` and `40_Kinship_marker_QC.csv`.

## 41 — Manuscript uncertainty and comparisons

* RIL bootstrap confidence intervals for GBLUP held-out R2 and MSE; source-population block bootstrap when all origin IDs are available and at least eight groups are represented.
* Paired differences in held-out mean squared error: GBLUP minus fold training-mean baseline and **the fixed K=75 model** from script 16. This is **not** a comparison to an optimally tuned top-K model (which would require nested selection of K).
* Observed phenotypic BLUE correlations versus correlations of independently predicted GBLUP values. The latter are **not** estimates of genetic correlations G; the running MCMCglmm model is the relevant genetic-covariance analysis.
* Figure-ready PDF comparisons.

Outputs `41_Bootstrap_GBLUP_accuracy_CIs.csv`, `41_Paired_GBLUP_vs_baselines.csv`, `41_Observed_vs_predicted_trait_correlations.csv`, `41_GBLUP_vs_top75_marker_prediction.pdf`, and `41_Multitrait_correlation_recovery.pdf`.

## Critical interpretation notes

* Both the common and rare analyses are restricted to the **existing founder-GWAS-matched SNP universe**. The targeted IV/V variants absent from founders require the already completed raw-sync extractions and script 34 targeted WF tests; they will not be accidentally reclassified as founder-GWAS SNPs.
* Calibration is conditional on A00 read sampling with a Jeffreys prior and specified WF `Ne`; it does not by itself account for sequencing-mapping artifacts, founder ascertainment, selection on a linked background, or uncertain effective population size.
* Frequency strata are based on an *estimated* starting MAF; bin boundary uncertainty and residual LD between R5 clusters remain.
* Full-index beta is fixed from the present phenotypic fitness regression; its gradient uncertainty and possible nonlinear selection are not propagated through the GBLUP comparison.
* For the manuscript, use the **GA1/GA2/GA4** results in the main text. Keep the all-GA analysis with extinct GA3 in Supplementary Materials and explain that the primary analysis conditions on survival.
