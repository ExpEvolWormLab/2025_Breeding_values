# Pipeline revision: scripts 05, 05b and 06

## What was changed, and why

### 05_Concatenate_cluster_GWAS.R (former script 10)
- Preserves the legacy `all_GWAS_clusters` columns and `.RData` filenames, including `beta_SF_NaCl`, `p_value_REML.cor_SF_NaCl`, `is_sign_SF_NaCl`, `rs`, `chr`, `ps`, `group`, `clusterID`, `CHROM`, and `POS`.
- Uses explicit `X_ID` to `clusterID` matching rather than implicit data-frame merges or positional matching. Retains **every** cluster from script 03, even if it was not tested by script 04; its missing effects/p-values remain `NA` and `is_sign=FALSE`.
- Computes BH-adjusted p-values within each trait across the cluster-GWAS tests actually performed. **These are post-selection descriptive q-values**, because initial marker selection, clustering, and mapping occur in the same RIL data. Do not claim independent replication or full selection-adjusted FDR control.
- Reads the `04_Cluster_mapping_QC.csv` file to distinguish an expected missing chromosome output (zero eligible clusters) from incomplete execution of script 04.
- Writes both the original `.RData` output used downstream and new explicit CSV/QC tables. Eliminates `rm(list=ls())`, `exists("all_GWAS")` and `save(list=ls())`.
- Important: `X_ID`, `beta.5`, and `p_value_REML` are the column names used in the original GridLMM outputs; inspect one of your new output headers if GridLMM changed its format.

### 05b_Assess_LD_thresholds.R (new; not a p-value optimizer)
- Runs **after 04 and 05 have finished the candidate thresholds**. Combines direct clustering properties (SNPs, clusters, singletons and physical spans), phase-oriented genotype QC, proportion of clusters mapped, and descriptive cluster association counts.
- Produces `threshold_review/LD_threshold_clustering_QC.csv`, `threshold_review/LD_threshold_cluster_associations.csv` and `threshold_review/LD_threshold_sensitivity.pdf`.
- Set `primary_threshold=0.5` provisionally, chosen for biological resolution, not maximizing the number of significant associations. Evaluate robustness of the central biological *selection-response* result at 0.4 and 0.6 later; the assessment cannot be finalized before the selection analyses are reconstructed.
- A threshold can be called `optimal` only with a prespecified independent criterion; generally recommend `chosen primary threshold` rather than `best threshold`.

### 06_Phenotype_reconstruction_diagnostics.R (former 11_TEST_PREDICTIONS_all_thresh.R)
- Matches RILs and cluster IDs *by their names*, not by sorting data frames independently.
- Uses the current adjusted RIL BLUEs (`phen_NaCl`) and selection gradients saved by script 00. The beta index is constructed from the **seven original BLUEs**, not treated as an eighth independent trait.
- Produces only explicitly labeled **in-sample** correlations and R-squared, including a separate panel for the derived beta index. There is no expectation that its regression slope equals one: allele dosage coding and phenotype scaling differ.
- Does not select the LD cutoff by maximizing this R-squared. Testing and estimating associations on the same RILs makes it optimistic. True predictive validation requires re-estimating associations and choosing clusters within training folds, ideally grouping related RILs; that is a separate analysis.
- Save figures/data under `phenotype_reconstruction/`.

## Run order
1. Finish `04_Mapping_clusters.R` (it is currently running). Do not concatenate partial outputs.
2. `source("05_Concatenate_cluster_GWAS.R")` (to process only R5, first set `thresholds_clust <- 0.5`).
3. `source("CHECK_05_after_run.R")` and review tested/untested counts.
4. Once script 04 and script 05 have completed **all seven thresholds**: `source("05b_Assess_LD_thresholds.R")`.
5. `source("06_Phenotype_reconstruction_diagnostics.R")` for the optional descriptive reconstruction sensitivity figure.
6. Later, after linking GWAS clusters with experimental-evolution trajectories: evaluate the main biological result at **0.4, 0.5, 0.6** before finalizing the methods and manuscript figures.

## Remaining methodological issues
- The primary biological threshold choice must be justified independently of statistical significance or a desired direction of selection-response association. Verify LD chaining and long-range spans, cluster genotype missingness, chromosome-heterogeneous recombination, and sensitivity of effect sizes / uncertainty.
- `is_sign` after the cluster-GWAS is **post-selection**, so cluster QTL counts are descriptive.
- GWAS genotype orientation is arbitrary across clusters, so signs of fitted cluster effects must be harmonized to experimental-evolution allele-frequency changes before any directional cross-data analysis. This is the next high-priority audit.
- Output filenames `R10` mean threshold **0.99**, never 1.0.
- The code is statically reviewed only: private input data and an R runtime are not present here. Validate against a completed R5 mapping before the full threshold loop.
