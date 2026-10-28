# Review notes — next pipeline scripts 07 and 08

## Why these two now
- **07** refactors the old `13_Characterize_pleiotropy_R5_NaCl_v2.R`. It is downstream of 05 and uses the current 00 BLUEs. It does **not** load the obsolete `RData/Preparatory_GLMM.RData` or external `Gmatrices_P.RData`.
- **08** refactors only the **active preparation** of old `14_Load_PoolSeq.R`. Per your earlier clarification, material after line 278 in the old file was not part of the live pipeline. This refactor does not import those exploratory plots or tests. Selection coefficients will be estimated in new script 09, adapted from old 16; selected-cluster/allele-frequency analyses follow after that.

## Key changes in 07
1. Works directly from script 05 `all_GWAS_clusters`; no implicit `merge()` that could duplicate or discard clusters.
2. Calculates regional assignment by matching the original SNP positions to the Rockman–Kruglyak interval file, retaining the historical center > arm > tip priority. Outputs `Cluster_location_QC_R5.csv` to flag boundary-spanning clusters or unmapped SNPs.
3. Uses **seven primary traits**, not the composite `beta`, for `Neff`, Manhattan, Euclidean, and number of affected traits. GWAS betas are already phenotype-SD standardized: no extra SD division.
4. The correlation matrix is computed from corrected, adjusted RIL BLUEs but is **not measurement-error corrected**. It is explicitly labelled as such. Conditional exponent tests against 1 (Manhattan) and 0.5 (Euclidean) are exploratory, not final evidence of causal synergistic pleiotropy.
5. Retains an optional descriptive figure and exports full cluster annotations. Outputs go to `pleiotropy_results/R5/`; the script 05 cluster-GWAS input is not overwritten.

## Key changes in 08
1. Sample order and generation/replicate metadata are written to a reusable `Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv`, which **must be checked against the sync-file columns**, especially the excluded first and third BAM-list entries. When the manifest exists, it is reused.
2. `poolSeq::read.sync(..., polarization='minor')` makes allele orientation explicit. `all.counts.caendr` is the **minor-allele read count**; `all.counts.alt.caendr` is its complement. These historical object names are retained for transition but do not imply REF/ALT.
3. Original relative-coverage filters (`max <=4`, median in `(0.15,5)`) and founder SNP restriction retained.
4. SNP and allele metadata are joined strictly using `CHROM.POS` identifiers, not positional `%in%` matching. Allele metadata returned by `poolSeq::alleles()` can be reordered, so it is matched explicitly.
5. Saves a well-scoped `RData_GridLMM_NaCl/08_PoolSeq_prepared.RData`, allele metadata CSV, sample-level coverage QC and variant-filter QC. Does not save arbitrary `.GlobalEnv` content or historical `res` objects.

## Next script (09) — important scientific decision
Old `16_Compute_s.R` fits `(1|replicate_fact)` to pooled allele counts over generations. Each read is not a biological replicate, so binomial likelihood-ratio p-values can be severely anticonservative when there is pool-to-pool variation, drift, or overdispersion. We should decide whether to use replicate-aware models or calibrated trajectories/permutations, and keep the `s` estimates in units of log-odds change **per generation**, not automatically equate them to classical population-genetic selection coefficients. Allele frequencies must be transformed to the GWAS effect allele before assessing directional concordance.

## Inputs / outputs
- 07 inputs: `RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R5_with_residuals.RData`, `RData/Preparatory_GLMM_NaCl_only_withBeta.RData`, `txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt`, `Caendr_filtering/rockman_kruglya_map_liftOver.txt`. Output: `pleiotropy_results/R5/`.
- 08 inputs: `Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync`, `Poolseq_GA_2026/all_bam_list.txt` (first execution only), `data/caendr_SNPs.txt`. Output: `RData_GridLMM_NaCl/08_PoolSeq_prepared.RData`, 08 QC CSVs, and sample manifest.

## Execution / verification
```r
source('07_Characterize_pleiotropy.R')
source('08_Prepare_PoolSeq.R')
```
Run 07 independently of 08. Check `table(all_GWAS_clusters$chr_type2)` and the region-QC CSV. Check that the 08 manifest matches the true sync-file sample columns before using its outputs in selection analyses. The scripts were statically reviewed but **not run on the user data**.
