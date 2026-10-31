# Final targeted analyses: chromosome IV mutation effects and survival sensitivity

## Project context and biological interpretation

The chromosome IV variants **IV:8,892,663 (C>T)** and **IV:10,499,135 (A>G)** were confirmed as reference calls in the 16 founder haplotypes, but your RIL panel contains lines sampled from both ancestral and evolved material. There were ~140 generations of lab domestication between founder crossing and the salt evolution experiment. **Absence from founder haplotypes is NOT equivalent to absence from the derived RIL panel.** Script 19 directly checks segregation among phenotyped RILs and estimates genotype–trait associations and epistasis, if identifiable.

The salt-evolved GA samples include GA1, GA2, GA3, GA4. The current manifest contains generation-50 samples from GA1/GA2/GA4 but no GA3 generation-50 sample. **Missing sequencing data do not establish biological extinction.** The supplied status CSV therefore marks GA3 `unknown`. Edit it to `extinct` if independently documented (or to `survived` if it did survive). Unknowns are excluded from a conservative *complete-trajectory/survivor-only sensitivity* but **never called extinct by the code**.

## Script 19 — chromosome IV effects and potential epistasis

Run from project root after the saved script-00 preparation:

```r
source('19_Test_CHRIV_mutation_effects_epistasis.R')
```

Expected input: `RData/Preparatory_GLMM_NaCl_only_withBeta.RData` with `phen_NaCl`, `genotype_matrix`, and `genotype_pruned_matrix`. It uses the two original SNP dosage columns (not founder variants alone), computes a **leave-chromosome-IV-out** kinship matrix, fits single-marker and mutually adjusted association tests, and tests the two-locus dosage product. Reports an indicative 80% detectable effect size and expected power at candidate interaction magnitudes, in units of RIL BLUE standard deviations. This uses a profiled variance-ratio model on RIL BLUEs, so p-values are **approximate** and should be corroborated with the original replicate-level GridLMM for final association claims.

Outputs (`mutation_effects_chrIV/`):

- `19_Two_locus_genotype_counts.csv`: all nine possible genotype combinations and carrier counts.
- `19_Mutation_effects_and_epistasis.csv`: marginal and joint effects, epistasis estimability, conditional estimates and BH-adjusted p-values for seven phenotypes. Focus on **SB** first.
- `19_Epistasis_power.csv`: estimated interaction power for effects of 0.10, 0.25, 0.50 and 1.0 trait SD; approximate minimum effect detectable at 80% power.
- `19_Main_effect_power.csv`: marginal and conditional main-effect power for each mutation, with 80% minimum detectable effects.
- `19_RIL_mutation_genotypes.csv`: exact calls per phenotyped line, for source-stratified checks.
- `19_SB_by_CHRIV_genotypes.pdf`: descriptive SB boxplots.

### Very important: RIL ancestry/origin adjustment

If you have the true origin of every RIL, fill in `19_RIL_source_populations_TEMPLATE.csv`, save as **`19_RIL_source_populations.csv`** in your project root, with `pop_label` matching preparation and `source_population` being e.g. ancestral/evolved, or the actual source population. The script adjusts for this categorical origin. **Do not invent the labels** from genotype calls. If the two mutation genotypes are restricted to evolved-origin lines, genetic effects will be partially confounded with origin; full conditioning may remove power, and if the interaction is collinear with ancestry it is not identifiable. This is an observational association, not genetic causality.

If the mutations are almost perfectly associated, having many RILs is not enough to estimate their independent effects or interaction. The decisive quantity is the number of *recombinant two-locus haplotypes*, not total sample size.

## Script 20 — survivor/complete-trajectory salt-evolution sensitivity

Review `20_GA_population_status.csv` **before** running. Values: `survived`, `extinct`, `unknown`.

```r
source('20_Survivor_only_drift_and_GWAS.R')
```

Required inputs: the corrected `08_PoolSeq_prepared.RData`, existing `09_Selection_SNP_results.csv` allele-orientation metadata, `09f_drift_functions.R` (included in this package), corrected full `09f_SNP_drift_Ne1000.csv`, the original R5 phase metadata, and results of scripts 11 and 13. It **retains the one G0 A00 sample** and excludes `A6G140`, as in the current primary analysis. It keeps the identical GWAS, SNP universe, effect orientation and R5 clusters for a controlled comparison.

Outputs (`selection_results/survivors_only_Ne1000/`):

- `20_Retained_sample_manifest.csv` — auditable included pools.
- `20_SNP_drift_Ne1000_survivors.csv` — all screened SNPs, with re-estimated slopes, p and BH.
- `20_Compared_SNP_slopes.csv` — aligned full vs survivor-only SNP estimates.
- `20_SNP_phase_oriented_survivor_slopes.csv` and `20_Survivor_cluster_slope_summary.csv` — phase-aligned R5 cluster responses.
- `20_Survivor_overlap_summary.csv` — descriptive SNP/cluster overlap.
- `20_Matched_cluster_associations.csv` — full vs survivor-only Pearson correlations of all seven trait effects and full/reduced fitness indices on the **same clusters**.
- `20_Matched_LOCO_prediction.csv` — refitted leave-one-chromosome-out response calibration for each response set.
- `20_Paired_2Mb_block_bootstrap_summary.csv` — paired resampling CI for differences across the **same genomic blocks**.
- `20_Survivor_neutral_pilot.csv` — a small neutral calibration sanity check, not a genome-wide tail validation.

The *original* 09f and scripts 10, 11 and 13 results are not overwritten. The code refits the evolutionary response, overlap, phase alignment, GWAS–evolution correlations, held-out predictive comparisons and block bootstraps with the surviving GA populations; nothing needs to be rerun upstream of Pool-seq.

### Limitations to state in a manuscript

* Conditioning on survival is a post-selection decision, and may itself bias the evolutionary summary. Keep the all-populations analysis primary or present both; survivor-only is a **sensitivity**, not intrinsically a truer target.
* An incomplete final time point may reflect sampling/archiving, not extinction; update status from records.
* The SNP universe was assembled under the **original 13-pool coverage filtering**. This matched-universe sensitivity is not a full redefinition of the SNP inclusion threshold from raw sync. For a *fully re-filtered* survivor-only analysis, Pool-seq preparation must be rerun from raw counts with filtering restricted to retained pools.
* The neutral-pilot default is too small to validate extreme p-values and BH discovery counts; infer on continuous slopes and cluster-level robustness first.
* The RIL panel includes ancestral- and evolved-derived material; changing GA population survival status **does not justify removing RILs** from the genetic mapping analysis. If you want source-specific GWAS, that is a distinct analysis.

## Please share these results

1. `mutation_effects_chrIV/19_Two_locus_genotype_counts.csv`
2. `mutation_effects_chrIV/19_Mutation_effects_and_epistasis.csv`
3. `mutation_effects_chrIV/19_Epistasis_power.csv` and `19_Main_effect_power.csv`
4. `selection_results/survivors_only_Ne1000/20_Matched_cluster_associations.csv`
5. `selection_results/survivors_only_Ne1000/20_Matched_LOCO_prediction.csv`
6. `selection_results/survivors_only_Ne1000/20_Paired_2Mb_block_bootstrap_summary.csv`
7. `selection_results/survivors_only_Ne1000/20_Retained_sample_manifest.csv`

No R runtime is available here. Scripts were statically reviewed against your saved code objects but not run on the private genotypes/Pool-seq matrices.
