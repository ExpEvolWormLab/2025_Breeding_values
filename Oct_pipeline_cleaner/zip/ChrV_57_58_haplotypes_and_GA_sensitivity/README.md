# Chromosome V: local RIL haplotypes and evolutionary-replicate sensitivity

Continuation of scripts **54–56**, reusing the existing phase-aware R5 framework.

## Run order

From the **project root** (not from inside `Oct_pipeline_cleaner`):

```r
source('Oct_pipeline_cleaner/57_ChrV_RIL_haplotype_and_trait_tests.R')
source('Oct_pipeline_cleaner/58_ChrV_GA_leave_one_population_out.R')
```

The scripts are independent once scripts 43, 45, 54–56 have finished; you may run either one alone. If you want the optional founder comparison, first run 57, verify the original VCF's REF/ALT allele identities, then run:

```r
source('Oct_pipeline_cleaner/57b_ChrV_founder_signature_optional.R')
```

The `57b` optional script needs an indexed VCF, 16 founder sample IDs and `bcftools`. Edit `vcf` and `founder_file` at the top. It will **not** assign founder origins unless `confirm_identical_ALT_allele_coding <- TRUE` is explicitly set after allele verification. Even then, output is a four-marker signature match, **not local ancestry inference**.

## Inputs already used in your pipeline

### Script 57 — RIL haplotype and trait associations

- `RData/Preparatory_GLMM_NaCl_only_withBeta.RData` (`genotype_matrix`, `phen_NaCl`, optional `GRM_sub`)
- `selection_results/chrV_deep_dive/54_V3_12_cluster_prediction_audit.csv`
- `selection_results/chrV_deep_dive/55_Selected_RIL_markers.csv`
- `selection_results/chrV_deep_dive/56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv`

Four locked neutral-test tags: group 848 = `V_6322080`, 851 = `V_6637695`, 853 = `V_6826421`, 854 = `V_7014705`. **These are not necessarily causal sites.**

For comparison, the attached `EXPECTED_focal_tag_QC_from_existing_outputs.csv` records the phase and ancestral ALT frequencies calculated from your previous actual output. The ancestral *rare* allele for group 854 is REF, not ALT; for the other three it is ALT. **The script reads these frequencies from your current script-56 output rather than hardcoding them.**

Outputs under `selection_results/chrV_deep_dive/57_haplotypes/`:

- `57_Four_anchor_marker_and_ancestral_orientation_QC.csv`
- `57_RIL_four_marker_haplotype_assignments.csv` and `57_RIL_four_marker_haplotype_frequencies.csv`
- `57_Broad_local_background_representative_SNPs.csv`, `57_Broad_local_background_RIL_assignments.csv`, and `57_Background_vs_four_tag_pattern_counts.csv` (when >=8 representative markers and sufficient RILs)
- `57_Anchor_marker_pairwise_LD.csv`
- `57_RIL_haplotype_and_trait_BLUEs.csv`
- `57_Haplotype_association_seven_trait_BLUEs_EXPLORATORY.csv`
- `57_Within848_conditional_leads_EXPLORATORY.csv`
- `57_Four_tag_haplotype_frequencies.pdf`
- `57_Interpretation_and_QC.txt`

The script also groups RILs by **unsupervised local-genotype similarity** across one representative SNP per available R5 group, using a fixed 15% mean marker-disagreement cutoff. These group labels are descriptive background proxies, sensitive to marker choice and threshold, and not haplotype ancestry probabilities.

The biological object is a **four-marker multilocus RIL genotype pattern** (near-homozygous inbred lines), not a fully phased local haplotype from ancestral populations. Heterozygous or missing calls are excluded. RIL phenotype associations are on the previously adjusted trait BLUEs. The kinship-PC adjusted analysis is a sensitivity analysis, not a full mixed-model or causal fine-map. All p-values are exploratory, since the four tags were nominated after examining the genomic data.

### Script 58 — Leave-one-surviving-population-out

- The original complete, corrected Pool-seq data and orientation files from scripts 08, 09f, 34a/43x
- `Oct_pipeline_cleaner/09f_drift_functions.R`, `43_Unified_config_and_helpers.R`, `43x_Integrate_rawSYNC_R5_SNPs.R`
- `selection_results/unified_frequency_inclusive_R5/43_All_SNP_phase_slope_and_inclusion_QC.csv`
- `selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv`
- `selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv`
- `selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv`
- `Oct_pipeline_cleaner/20_GA_population_status.csv`

The script uses frozen eligible member SNPs from 43 for the primary 802 R5 clusters. It computes new 09f drift-aware GLS slopes using only **GA1+GA2**, **GA1+GA4**, or **GA2+GA4**, plus **one shared A00**. Within each cluster it takes the median of the re-estimated phase-aligned SNP slopes, as in script 43. It verifies that recalculated original slopes, aggregated full slopes and β/SB full LOCO predictions match script 45. This is a hard stop if source code/inputs no longer agree.

For each pair, it recalibrates β, SB and PC1 on five chromosomes and predicts the held-out chromosome. PCs are fitted only from training-chromosome **GWAS effects**, never the evolutionary responses of the test chromosome.

Outputs under `selection_results/chrV_deep_dive/58_replicate_sensitivity/`:

- `58_Corrected_sample_manifest.csv`
- `58_SNP_gls_phase_aligned_replicate_subset_slopes.csv`
- `58_R5_cluster_outcomes_fixed_membership.csv`
- `58_SNP_fitting_QC.csv`
- `58_LOCO_predictions_by_population_subset.csv`
- `58_LOCO_metrics_by_population_subset.csv`
- `58_Exploratory_paired_block_bootstrap.csv`
- `58_ChrV_four_focal_groups_population_omission.csv`
- `58_Immutable_input_and_reproduction_QC.csv`

**Critical:** This is sensitivity to dropping a replicate, not independent leave-one-population-out **prospective** forecast. The calibration still uses evolutionary responses from the retained populations on five training chromosomes. With only two evolved replicates in each sensitivity analysis, the observed outcome is noisier. The bootstrap is conditional on the fitted predictions and not an independent validation of selective mechanisms. GA3 is excluded from the primary comparison (confirmed extinct).

## What to send back

Send one ZIP containing the `57_haplotypes` and `58_replicate_sensitivity` output directories (CSV, TXT, and PDF files). If a script stops, copy the **full error message** and its input-check file; do not override a reproduction error simply to obtain plots.

## Testing status

The code was inspected against the uploaded scripts 09f, 43, 45, 53b and 54–56, and the supplied output table schemas. Delimiter/string structure was checked; it has **not** been executed in R here, because the original large RIL genotype and Pool-seq input RData files and R runtime are not available in this environment.
