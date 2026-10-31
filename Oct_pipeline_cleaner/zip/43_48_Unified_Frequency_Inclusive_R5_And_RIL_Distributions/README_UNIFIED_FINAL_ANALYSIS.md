# Frequency-inclusive R5 selection and prediction — manuscript workflow

## Scope and primary estimand

This package is a self-contained extension of scripts 33–39 and the established R5 GWAS. It **does not** use an ancestral-frequency cutoff. All 1,089 original GWAS LD clusters are entered into a registry, and every cluster with **adequate Pool-seq coverage and estimable evolution** enters the genome-wide GWAS–evolution prediction analysis, independent of selection significance. The number actually eligible is computed by script 43; do not substitute 1,089 for the eligible count.

- **Primary evolution:** surviving GA1, GA2, GA4 only. GA3 has been confirmed extinct.
- **Supplement:** all GA populations including extinct GA3, clearly separated.
- **Ancestral sample:** A00. Exclude A6G140, which may represent a different domestication/bottleneck history.
- **Frequency:** no restriction. Quality thresholds are based on A00 depth, available evolved pools and number of surviving replicates.
- **LD:** original fixed R5 phase-aligned groups, NOT new clustering based on evolution.
- **Hypothesis testing:** one GWAS- and depth-chosen tag (not chosen using observed evolution) per eligible R5 cluster, under a diploid, discrete-generation, Wright–Fisher neutral model and observed Pool-seq read depths.
- **Prediction:** median of ALL eligible phase-aligned SNP slopes per LD cluster, including those without selection significance. Leave-one-chromosome-out (LOCO) calibration for four independent GWAS predictors; block bootstrap uncertainty.
- **Founder-absent focal alleles:** chromosome IV/V etc. are kept in a **targeted raw-sync annex**. They cannot be silently counted as additional independent R5 tests without a unified, non-overlapping genotype/locus and GWAS effect definition.

## Running

Copy all `.R` files from this ZIP into your **project root**, alongside the earlier script outputs and their expected paths. Do not rename the scripts. Your corrected 20_GA_population_status.csv must contain `GA3,extinct` and the other three `survived`. Execute:

```r
source('48_Run_unified_main_analysis.R')
```

Alternatively, run in stages:

```r
source('43_Prepare_unified_R5_catalogue.R')
source('43b_Empirical_WF_null_validation.R')
source('44_Calibrated_WF_selection_all_R5_clusters.R')
source('45_Unified_GWAS_evolution_LOCO_prediction.R')
source('46_Unified_main_figures_and_focal_annex.R')
```

**Optional, separate RIL phenotype figure** after script 16 and rrBLUP script 40:

```r
source('47_RIL_predicted_vs_observed_distributions.R')
```

The latter uses only matching **random-RIL out-of-fold predictions**, and produces a distribution and observed-vs-predicted figure for each of seven component traits plus beta. It does not use or discuss source-population-held-out predictions. **Beta has no script-16 K75 prediction** and therefore is compared to GBLUP only.

## Required existing inputs

- `RData_GridLMM_NaCl/08_PoolSeq_prepared.RData` (13 pools with **corrected sync/sample mappings**)
- `selection_results/09_Selection_SNP_results.csv` (GWAS ALT orientation)
- `data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv` (fixed R5 SNP phase and group)
- `txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt` (all eight GWAS p-values for SNP tag selection)
- `selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv` (GWAS group effects for seven traits and beta)
- `selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv` (fixed full and reduced fitness-index weights)
- `selection_results/frequency_strata/33_All_SNPs_retained_GA_Ne1000.csv` (existing all-frequency corrected slopes)
- `selection_results/frequency_strata/33_All_SNPs_all_GA_Ne1000.csv` (supplement if available)
- `20_GA_population_status.csv` (final corrected statuses)
- Optional: `selection_results/frequency_strata/34a_Focal_raw_SYNC_ALT_counts.csv` and `34_Focal_WF_neutral_conditional_tests.csv` for the supplemental focal annex.

## Key configurable settings

Edit `43_Unified_config_and_helpers.R` before beginning. Defaults: `Ne=1000` diploid individuals; sensitivities 500 and 2,000; `min_A00_depth=30`, `min_evolved_replicates=2`, `min_pools=6`; 2-Mb genomic blocks; 4,000 bootstrap draws. The site-by-site sequential neutral simulations increase from 3,000 to as many as **3,000,000** draws per cluster and population scenario. Main and supplementary scenarios can thus take a long time. The main simulation writes a checkpoint after every 20 clusters, and resumes safely if restarted with **identical parameters and tag catalogue**. Logs and output names provide traceability.

`43b` performs a separate **fixed-true-p0 simulation** to check whether the conditional null test appears anti-conservative at 1% and 5% nominal levels. By default it tries true ALT frequencies .001, .005, .02, .07, .17, .35 and .5, with 400 outer simulated datasets per frequency and 2,500 nested simulations for each. If inflation is detected, script 44 **stops**. A pass is **not proof** that the extreme 10⁻⁵ tail is calibrated. Never report strict genome-wide FDR control without accepting and discussing the model assumptions and tail diagnostics.

The Monte Carlo neutral tail is a **baseline-posterior-integrated predictive test** (Jeffreys Beta(½,½) prior), not an unconditional exact frequentist test for every fixed unknown ancestral frequency. This is an important limitation at almost-zero baseline observations. The simulation-derived p-value bounds conservatively handle *finite Monte Carlo uncertainty* and sequential stopping; BY correction addresses dependence across clusters under a valid null, but cannot remedy **misspecified demographic, sampling, or ancestral-frequency assumptions**.

Direct binomial sequencing is assumed; finite library sampling of individual worms, effective population size uncertainty, migration, pooled sampling bias, common ancestral haplotype covariance beyond A00, and sequencing errors can affect significance. Existing paired data have no explicit read-error or individual-pool-size model. Sites with 0/100 A00 ALT reads are *not* treated as allele-free. The conditional neutral slope is not a selection coefficient `s`.

## Main result files

All produced under `selection_results/unified_frequency_inclusive_R5/`:

- `43_Unified_R5_cluster_registry.csv` — one row **per original R5 cluster**, with coverage eligibility and prespecified test tag.
- `43_All_SNP_phase_slope_and_inclusion_QC.csv` — every original R5 SNP; include/exclude reasons.
- `43_Catalogue_size_and_QC.csv` — actual number of newly eligible clusters versus the 675 common-frequency reference.
- `43b_Empirical_null_calibration_FIXED_p0.csv` — validation and guard.
- `44_primary_survivors_WF_cluster_selection.csv` — main selection test, MC confidence intervals, conservative BY and BH q-values; one test per cluster.
- `44_supp_all_GA_including_extinct_GA3_WF_cluster_selection.csv` — supplementary sensitivity analysis.
- `44_Selection_summary.csv` — multiplicity and MC budget counts.
- `45_All_R5_clusters_GWAS_evolution_predictions.csv` — every R5 group including untestable groups flagged as such; GWAS effect predictors, LOCO fits, neutral selection results.
- `45_LOCO_performance_main_and_supp.csv` — held-out chromosome prediction of evolving allele-frequency slopes.
- `45_Paired_genomic_block_bootstrap.csv` — comparison of full seven-trait index against reduced index, SB and direct-beta GWAS effect.
- `45_Block_bootstrap_correlations.csv` — CI for uncalibrated GWAS-predicted effect vs evolutionary slope correlations.
- `46_Main_unified_frequency_inclusive_R5.pdf` — manuscript-oriented primary figures.
- `46_Supp_GA3_inclusive_comparison.pdf` — primary-versus-extinct-GA3 sensitivity.
- `46_Supp_named_focal_trajectories.pdf` and `46_Supp_focal_raw_SYNC_annex_NOT_genomewide_FDR.csv` (if raw counts are available).

For the separate genomic prediction comparison, script 47 creates `phenotype_prediction/GBLUP_rrBLUP/47_Observed_vs_predicted_RIL_distributions.pdf` and a detailed accuracy table. This plot is not a test of full distributional equality; it describes shrinkage, extremes and held-out prediction.

## Integrity and interpretive safeguards

1. **Do not** select clusters for genomic prediction using significance in `44` — script 45 explicitly uses all quality-eligible clusters.
2. **Do not** multiply the number of tests by the number of linked SNPs: use one frozen GWAS-informed tag for the primary R5 family.
3. **Do not** equate focal-targeted raw-sync p-values with a corrected genome-wide q-value.
4. **Do not** interpret GA3 as a surviving population or include it in the main figures.
5. **Do not** interpret the LOCO R² as prediction of fitness in new RILs; it predicts *evolutionary frequency response of held-out chromosomes* from independently measured GWAS effects.
6. **Do not** assert that selecting a GWAS tag proves causal selection: linked haplotypes and hitchhiking remain plausible.
7. **Do not** claim this includes **all genome-wide de novo/rare variation**: the R5 catalogue was constructed from founder/GWAS markers. The focal annex is a separate ascertainment class.
8. **Do not** claim that rare alleles have larger causal effects from frequency changes alone.

## How to send results back

ZIP the entire `selection_results/unified_frequency_inclusive_R5/` directory. The large SNP-level QC file and MC checkpoints can be omitted if upload size is a problem. Include PDFs and the `43b` diagnostic. Also send the new script-47 PDF / accuracy table if generated. The GBLUP results already received are sufficient; no GBLUP rerun is necessary.

**Technical limitations:** This environment does not have R installed, so the scripts were statically checked for balanced syntax, function and column contracts, and alignment with the previous script outputs; they have **not** been executed on your large genotype and Pool-seq objects here. Any error or failed validation should be reported and addressed before interpreting significance.
