# Scripts 03 and 04 — final code review and run instructions

## Supplementary Table S1

`Supplementary_Table_S1_updated.docx` and `.csv` reproduce the **new, user-supplied** phase-aligned summary counts. S1 retains the three columns from the manuscript. The table covers LD cutoffs 0.2, 0.4, 0.5, 0.6, 0.8, 0.9 and **0.99**, not 1.0. The prior manuscript caption's denominator of 33,510 significant SNPs is deliberately not reused because the current GWAS SNP count has not been supplied. **Supplementary Table S2 (per trait) still needs regeneration** from `02b_Summarise_clusters_CHR.R`, particularly `Summary_clusters_by_trait_R2.csv`.

## 03_Create_matrix_for_cluster_mapping_CHR.R (formerly 08)

**Confirmed original risks:**

- The old script filtered the GEMMA SNP genotype table using `%in%` and then assigned QTL groups positionally, rather than matching SNP IDs. The displayed `table()` merely printed a check and did not prevent a silent mismatch.
- The old script averaged allele dosages with arbitrary SNP REF/ALT orientation. Your new signed-LD validation showed that some high-LD SNPs have opposite effect signs exclusively because their dosage coding is reversed. Unoriented genotype averaging can cancel the shared haplotype signal.
- The old script filled missing cluster genotypes with **3**, which is not a valid alternate-allele count when genotypes are coded 0/1 or 0/1/2.
- The old script read a precomputed `genotypes_for_gemma.txt` from a separate directory even though script 00 already saves the exact SNP genotype matrix used for GridLMM.

**Revised behavior:**

- Reads the SNP genotype matrix from `RData/Preparatory_GLMM_NaCl_only_withBeta.RData` (the same source used in script 01).
- Matches SNPs by `rs` and genotype matrix column names, and produces one genotype per cluster and RIL.
- Recomputes SNP phase within each original LD component using signed `D` and a maximum-R² spanning tree, following **the same rule as script 02**, then orients `2 - dosage` (or `1 - dosage` for a 0/1 encoding) before averaging markers within a group. Full-component phase is essential when members of a sign-split subgroup are not directly linked.
- This phase computation is **repeated intentionally**, since the completed script 02 saved phase statistics but not per-SNP phase assignments. A future improvement would be to save an explicit phase column during script 02 and consume it here.
- Leaves all-missing group/RIL combinations as NA, and writes genotype QC plus per-SNP phase for reproducibility.
- Preserves legacy outputs: `data_BLUE_GridLMM_NaCl/data_clusters_BLUE_R5.txt` and `.../genotypes_matrix_clusters_R5.txt`, as well as R2/R4/R6/R8/R9/R10 (R10 represents 0.99).
- Also writes `SNP_phase_clusters_R5.csv` and `Cluster_genotypes_QC_R5.csv`.
- A **new cluster genotype** can differ from the old version even when QTL membership does not change, because reversed allele codings are now handled properly. This is an intended methodological correction.

**Before running:** Keep your current `data_BLUE_GridLMM_NaCl/ld_window_100.hap.ld` produced by script 02. Make sure `genotype_matrix` is coded in either 0/1 or 0/1/2, and check the printed dosage range. Run **R5 alone first** by temporarily setting `thresholds_clust <- 0.5`, check `CHECK_03_04_after_run.R`, then restore all seven cutoffs if you need the threshold-sensitivity analyses.

## 04_Mapping_clusters.R (formerly 09)

- Reads script 03's cluster matrix (cluster rows x RIL columns), explicitly matches the RILs by name, and transposes it to RIL rows x cluster columns for GridLMM.
- Reuses the SNP-level GWAS phenotype and LOCO kinship model exactly; caches one kinship per chromosome rather than recomputing it in each nested trait/threshold loop.
- Does **not** transform missing dosages into 3. By default, it excludes incomplete or monomorphic cluster markers. The optional `missing_policy <- "mean"` performs mean imputation of missing RIL dosages, but it must be reported if chosen.
- Writes results in per-threshold subdirectories, e.g. `GridLMM_clusters_NaCl/R5/gwas_results_R5_SF_CHR_I_NaCl.txt`. This is the layout **expected by original script 10**, although original script 09 wrote the files at the directory root. This resolves an existing inconsistency between those two scripts.
- Writes `GridLMM_clusters_NaCl/04_Cluster_mapping_QC.csv` so the number of tested, incomplete, and monomorphic clusters is explicit. If the number excluded is substantial, this needs attention before inference.
- `missing_policy = "exclude"` is close to the SNP GWAS's complete-marker policy; however, the old script's values of 3 were *not* true missing data, so association statistics can change after the fix. Additionally, orientation changes the clustered marker design; cluster-level GWAS beta signs are arbitrary with respect to the phase root and should not be compared directly to unaligned SNP effect signs.

## Interpretation and downstream consequences

- This is an analysis change, **not only cleanup**. The new phase-aware cluster genotype construction can change GWAS effects and their p-values. Refit downstream analyses from script 04 outputs.
- The cluster SNPs were selected based on trait associations before testing the averaged cluster genotype on the same panel. Cluster-GWAS significance is **post-selection / descriptive**, not an independent association-discovery validation. Avoid interpreting it as independent replication.
- The next unreviewed script (old `10_Analyze_permutation_Gemma_clusters_with_residuals.R`, to become **05**) expects nested `R#` mapping directories; this script now creates them. The old 10 omits the **0.5** threshold in its loop, so do not use it unchanged for your primary threshold. Its joins and file paths still need review.
- The old phenotype-prediction script (11) uses a fragile genotype-file read; update it to read `header=TRUE, row.names=1, check.names=FALSE` and match RIL identities explicitly when we reach that script.
- Statistical results have not been reproduced here: no private phenotype, genotype, GWAS or LD input files and no local R executable were available. A focused synthetic check script is supplied for validation after executing the R scripts on your machine.
