# Phase 12 — Drift-aware allele-frequency responses across ancestral-frequency classes

**Research goal.** Revisit the corrected Pool-seq salt-evolution experiment without treating the ancestral 0.10–0.90 restriction as a biological boundary. Compare common and rare GWAS SNPs in the same baseline-conditional neutral framework, then test a prespecified shortlist of historically interesting low-frequency variants using discrete Wright–Fisher simulations.

## What this **does**, and what it **does not** do

1. **Genome-wide within the existing founder-filtered GWAS SNP universe:** the 09f RData/CSV selection dataset contains SNPs selected during script 08 and matched to the GWAS, including some rare alleles. Script 33 recalculates a *baseline-conditional* neutral trajectory score and approximate normal p-value for all those SNPs, classifies them by **A00 posterior-mean minor-allele frequency**, and compares all GA with GA1/GA2/GA4 retained GA. These approximate p-values are **not yet genome-wide calibrated discoveries**.
2. **Prespecified low-frequency variants absent from the filtered SNP universe:** script 34a retrieves their GWAS-ALT counts directly from the original raw sync, using VCF REF/ALT annotations. Script 34 then calculates site-specific posterior-predictive Wright–Fisher neutral p-values, conditional on the A00 read counts, for both GA sets and three Ne scenarios. This is a **targeted Monte Carlo test**, not an all-sync-genome scan.
3. **Calibration:** script 34b samples a random set of SNPs within each frequency class, independently of their observed selection scores, and checks whether the *normal* approximation gives nominal null rejection rates. It **does not** validate the biological Ne model or calibrate the genome-wide tail needed for BH thresholds around 1e-7.
4. **Polygenic prediction:** script 36 reconstructs the original independent RIL GWAS beta predictors for **all** R5 clusters (not only the original 675 common-frequency clusters), checks agreement with script 13 where possible, uses R5 phase labels, and the reestimated responses, and checks correlation and leave-one-chromosome-out R² **within ancestral-frequency classes**. A given LD cluster may contribute to multiple frequency bins; these bins are not independent cluster samples.

Both the analytic null variance and Monte Carlo simulations integrate ancestral-frequency uncertainty using a Jeffreys Beta(1/2,1/2) posterior from **A00 only**. The observed A00 read count is **held fixed** in the simulations. Simulated evolved trajectories drift independently from a common latent ancestral p0 and produce read counts at the observed per-SNP depths and the verified sample times.

Allele-frequency slope is **not** an inferred selection coefficient. No new claim of causal selection, large causal effect, or strict infinitesimal model is justified by the scripts alone.

## Required files in the project directory

- `RData_GridLMM_NaCl/08_PoolSeq_prepared.RData` — must be the **sample-alignment-corrected** version (repair 08a or aligned 08), not the original swapped-column file.
- `selection_results/09_Selection_SNP_results.csv` — matches the 08 data and defines GWAS ALT orientation.
- `09f_drift_functions.R` (same directory as the scripts), or `drift_nullfit_update/09f_drift_functions.R` — the **original** script 09f helper.
- `20_GA_population_status.csv` — the existing population-status table; GA3 remains excluded from the retained analysis if still marked `unknown`, but **must not be described as biologically extinct until confirmed**.
- For the targeted raw-sync module: `Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync`, the verified `Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv`, the script-26 VCF founder signatures and R5 SNP tags from the *extended* or original deep-dive run, plus `26_Focal_variants.csv`.
- Script 36 also needs `data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv`, the complete `selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv`, and `selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv`. If present, the 675-cluster `13_Cluster_predictors_and_response.csv` is used as an exact reconstruction QC.

**Copy the `.R` files from this ZIP into your R project root**. Edit only the path variables at the top where your project differs; do not change the verified 13-pool ordering.

## Recommended run order

```r
# Quick implementation self-test on synthetic data (optional but recommended)
source('37_Synthetic_WF_null_selftest.R')

# Stage 1: all GWAS-filtered SNPs, regardless of ancestral MAF
source('33_Reanalyse_drift_all_frequency_classes.R')

# Stage 2: exact named variants absent from script 08 (reads the raw sync once)
source('34a_Extract_focal_ALT_counts_from_raw_SYNC.R')

# Stage 3: calibrated, targeted neutral tests (Ne=500,1000,2000)
source('34_Test_focal_variants_WF_conditional.R')

# Stage 4: check normal-tail calibration by MAF class
source('34b_Calibrate_normal_tail_by_frequency.R')

# Stage 5: connect results to polygenic prediction
source('36_GWAS_prediction_by_frequency_class.R')

# Stage 6: generate figures
source('35_Plot_frequency_stratified_results.R')
```

**If `28b_Exact_raw_SYNC_ALT_trajectories.csv` already contains all named variants**, script 34 can read it automatically and you can skip 34a; otherwise run 34a. The genome-wide analysis does not require rereading the giant sync. Script 34a reads it in streaming blocks and retains only VCF-annotated nominated sites.

`34_Additional_target_SNPs.csv` is an optional two-column CSV (`SNP_ID,region`) for extra targets chosen *before inspecting their script-34 p-values*. For example `III_6900000,extra_chrIII`. Empty template provided.

### Core output files

In `selection_results/frequency_strata/`:

| Filename | Interpretation |
| --- | --- |
| `33_Frequency_class_summary.csv` | Number of tested GWAS SNPs and approximate signals by starting-frequency class, full and retained GA |
| `33_All_SNPs_all_GA_Ne1000.csv` and `33_All_SNPs_retained_GA_Ne1000.csv` | Per-SNP starting-frequency posterior intervals, slope, approximate neutral score; *no genome-wide q-value* |
| `33_Matched_full_vs_retained_all_frequencies.csv` | Same-SNP full versus retained response slopes |
| `33_QC_09f_slope_match_*.csv` | Confirms preservation of the original corrected 09f slope geometry |
| `34a_Raw_SYNC_target_QC.csv` and `34_Target_count_coverage_QC.csv` | Ensures exact sites were found, REF/ALT matches and pool coverage is sufficient |
| `34_Focal_WF_neutral_conditional_tests.csv` | Targeted conditional Wright–Fisher Monte Carlo p-values and Ne sensitivities |
| `34b_Null_calibration_by_frequency.csv` | Null false-positive rates at several normal p cutoffs, within MAF classes |
| `36_Prediction_by_MAF_class.csv` | GWAS effect–evolution correlations and LOCO prediction by starting-frequency class |
| `figures/35_*.pdf`, `36_Quantitative_prediction_across_MAF_classes.pdf` | Frequency-stratified plots and focal trajectories |

## Interpretation constraints

- The original founder-filtered set is **not every segregating variant in the raw sync**. Chromosome IV/V low-frequency mutations excluded by script 08 enter **only** the nominated-site module. Their inclusion does not make genome-wide frequencies directly comparable across both ascertainment universes.
- A00 read counts of zero imply **undetected** ALT reads, **not** a proven true ancestral frequency of exactly zero. The posterior interval and site depth must be reported.
- At `n_sim_main=6000`, the smallest attainable site p-value is `1/6001 ≈ 0.000167`. This cannot certify genome-wide BH discoveries over ~200,000 SNPs. The separately named `p_BH_within_nominated_family_EXPLORATORY` is **not confirmatory** when candidates were noticed because they had unusual trajectories.
- Very rare ALT reads are sensitive to mapping, sequencing errors and contamination; the sync lacks read-level strand and mapping-quality diagnostics. Validate a subset of important low-frequency sites in BAMs.
- The baseline conjugate model assumes independent sequencing reads and does not model a finite number of sampled diploid individuals in each library, sequencing errors, migration, changes in Ne or shared pre-experiment drift. `Ne=500/1000/2000` is a sensitivity range, not independently measured.
- All GA versus retained GA is a **survivorship sensitivity**, not a proof that exclusion is biologically preferable. Confirm GA3 status from the experiment diary.
- The MAF-class analysis of polygenic prediction uses R5 clusters that may contain SNPs from multiple frequency classes; it does not separate causal allele-frequency effects or validate an infinitesimal distribution of effects.

## Files to send back

Please ZIP `selection_results/frequency_strata/`, or at minimum send `33_Frequency_class_summary.csv`, `34_Focal_WF_neutral_conditional_tests.csv`, `34_Target_count_coverage_QC.csv`, `34b_Null_calibration_by_frequency.csv`, `36_Prediction_by_MAF_class.csv` plus the `figures/` PDF files. Keep the two very large genome-wide SNP CSVs optional initially.

## Scientific follow-up after inspecting these results

If frequency-dependent calibration is adequate only in certain ranges, select the appropriate final reporting framework separately for each class. A **fully calibrated, unified genome-wide rare-plus-common selection map** would require a deeper model/importance-sampling analysis on the **unfiltered raw sync** (including non-founder variants and explicit allele orientation) and confirmation of the extreme neutral tail before any new q-value map is claimed. This package deliberately avoids fabricating that result.
