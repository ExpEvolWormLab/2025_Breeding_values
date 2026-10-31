# Critical Pool-seq sample-column alignment correction

**Issue.** The refactored `08_Prepare_PoolSeq.R` used `poolSeq::af()` and `poolSeq::coverage()`, which *sort* extracted matrix columns by replicate and generation, then overwrote their native names with `populations$sample_id` in **unsorted sync-column order**. The following six labels were interchanged:

| Label assigned in the defective script 08 | Actual sample in that column |
|---|---|
| GA150 | GA1G36 |
| GA1G36 | GA150 |
| GA250 | GA2G36 |
| GA2G36 | GA250 |
| GA450 | GA4G36 |
| GA4G36 | GA450 |

The G0 samples, GA1/2/3/4 at generation 10, and GA3 generation 36 are not swapped. The original `14_Load_PoolSeq.R` sorted the population metadata after extraction, which is consistent with `poolSeq` behavior. The **revised** script 08 alone introduced this mismatch.

## Fast repair: do not reread the huge sync file

In the root of your R project, copy `08a_Repair_existing_PoolSeq_RData.R`, and run:

```r
source('08a_Repair_existing_PoolSeq_RData.R')
```

It requires `RData_GridLMM_NaCl/08_PoolSeq_prepared.RData` produced by the **buggy version of script 08**. It checks the expected manifest order, creates `RData_GridLMM_NaCl/08_PoolSeq_prepared_BEFORE_ALIGNMENT_FIX.RData`, corrects all four stored frequency/coverage/count matrices, and replaces the original `08_PoolSeq_prepared.RData` with the repaired version. It stamps `sample_alignment_validated=TRUE` to prevent repeating the correction and outputs `08_sample_alignment_repair_report.csv`. The original `08_coverage_QC.csv` was calculated from unfiltered SNPs and remains incorrectly labeled; the repair instead writes `08_coverage_QC_REPAIRED_FILTERED_SNPS.csv`. The SNP set and allele metadata are unchanged by this relabeling repair.

**Do not apply this repair to output from the corrected `08_Prepare_PoolSeq_ALIGNED.R` or to an original pipeline `SaveGLM.RData`.** The correct 08 script is provided as a full-reprocessing alternative and checks native column names explicitly.

## Regenerate downstream analysis

The previously generated results from scripts **09, 09c–09g, 10, 11, 12/12b/12c** are based on misassigned samples and are not valid for biological interpretation. In particular, all the old-vs-new model comparisons are compromised by the defect in the new pipeline. The old `16_Compute_s.R` output itself is not affected by this particular revised-script bug.

At minimum rerun, in this order, using the same previously downloaded scripts:

```r
source('09f_Drift_GLS_nullfit.R')
source('09g_WF_calibrate_nullfit.R')
source('10_Selection_GWAS_overlap_09f.R')
source('11_GWAS_selection_cluster_associations_09f.R')
```

**You may skip rerunning script 09 initially**: script 09f uses its saved CSV for the marker IDs and allele-orientation metadata, which are not changed by the sample-column permutation, and recomputes the slopes from the repaired 08 matrices. However, its old slopes/p-values in `09_Selection_SNP_results.csv` are no longer valid; rerun `09_Compute_selection_coefficients.R` if you want valid binomial comparisons or to rerun the full audit of historical results. Scripts 12/12b/12c should only be rerun if the corrected associations still disagree with the historical findings.

To avoid overwriting your *previous* selection outputs, copy the current `selection_results/` directory somewhere safe before running the scripts. The earlier selection results should be considered superseded, not comparison benchmarks.

## Full reproducible preprocessing alternative

`08_Prepare_PoolSeq_ALIGNED.R` rereads the full sync file and requires `poolSeq` plus `matrixStats`. It matches native columns named `F<generation>.R<replicate>.freq` and `.cov` to manifest rows before assigning sample IDs. It saves a mapping CSV `RData_GridLMM_NaCl/08_verified_sample_column_mapping.csv` and a stamped prepared-data object. Use this in the reproducible final pipeline instead of the original erroneous script 08. **Do not run the fast repair after generating data with this corrected script.**

## Methodological caution

The repaired data may substantially alter the estimated evolutionary slopes and the GWAS–evolution association. We cannot predict its magnitude or direction. It does not establish the correctness of other model assumptions, the chosen G0 sample, pooling variance, or allele orientation. Re-evaluate the revised neutral calibration and the composite-fitness association before making manuscript claims.
