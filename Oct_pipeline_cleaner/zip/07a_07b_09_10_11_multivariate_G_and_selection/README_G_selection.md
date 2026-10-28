# Multivariate RIL G and selection-analysis scripts — review notes

## Execution order

1. `07a_Estimate_multivariate_RIL_G.R` — estimates a **seven-trait among-RIL covariance** from original replicate-level NaCl phenotype data, not from BLUEs. Configure `run_mode='pilot'` first and inspect diagnostics; use `run_mode='final'` for the definitive analysis. Final requires three independent MCMC chains with adequate mixing. The model preserves the environmental covariates and date/session effects used earlier.
2. Your existing `07_Characterize_pleiotropy.R` (already run).
3. `07b_Refit_pleiotropy_with_multivariate_G.R` — recomputes Neff and the Manhattan/Euclidean scaling tests from the fitted G and quantifies sensitivity to posterior G draws. It **does not overwrite** old script 07 results.
4. `08_Prepare_PoolSeq.R` (already supplied): **verify `08_PoolSeq_sample_manifest.csv` against sync-file columns** and allele annotations first.
5. `09_Compute_selection_coefficients.R` — harmonizes Pool-seq minor allele to GWAS ALT, estimates site-level log-odds trends and overdispersion diagnostics.
6. `10_Selection_GWAS_overlap.R` — positional overlap between significant Pool-seq SNPs and the current R5 GWAS clusters.
7. `11_GWAS_selection_cluster_associations.R` — uses **the exact SNP phase assignments saved by script 03** to express each SNP's trend in the same allele orientation as its cluster genotype, then summarizes one value per GWAS cluster.
8. `CHECK_07a_09_after_run.R` — optional output QC.

## Important scientific choices

### G and Neff

- The prior diagonal-error subtraction was mathematically invalid for the supplied data. Do **not** use an arbitrary `nearPD` repair or maximal partial correction as if it were an unbiased G estimate.
- Model: `cbind(seven measured traits) ~ trait - 1 + trait:logD + trait:rel_humidity + trait:temperature`, random among-RIL **unstructured 7x7 covariance**, date/session trait-specific variance components, and residual **unstructured 7x7 covariance** for observations where traits are measured on the same replicate record. This partitions residual and between-RIL covariation directly.
- Between-line G is not unambiguously additive genetic G; it's repeatable among-RIL covariance conditional on modeled environmental effects, and is appropriate to assess trait redundancy. A genomic-relatedness-based model could target additive covariance but would add a distinct assumption about RIL genetic relationships.
- The new G is a posterior-mean covariance and is therefore positive definite (assuming valid samples). This is not proof of identifiability: inspect Rhat, ESS, between-line posterior covariances, and prior sensitivity.
- The script does **not** run a full bootstrap that redoes GWAS selection; the G draw sensitivity reflects uncertainty in G only. Report the scaling tests as conditional on the observed GWAS-selected traits and effect estimates.
- If traits were not measured on the same biological experimental unit, an unstructured *within-record* residual covariance is not interpretable: use `residual_structure='diagonal'` and clearly describe the limitation. Check `07a_cophenotyped_counts.csv`.
- Prior sensitivity: change the diagonal prior scales by ~2-fold and refit the final model if the estimated G correlations/Neff are consequential. If weakly identified, a reduced-rank/factor-analytic or pairwise model may be needed.
- Manuscript phrasing: "A seven-trait multivariate mixed model with RIL as a random effect estimated the among-line covariance of environmentally adjusted trait values; the corresponding correlation matrix was used to calculate effective pleiotropy." Do not claim it is a direct full measurement-error correction of Wagner et al.; it is a model-based alternative.

### Selection coefficients

- The old script `16_Compute_s.R` treated read counts as binomial trials in each pool. With drift, temporal autocorrelation, pooling and extra-binomial variance, those p-values can be anti-conservative.
- New script 09 retains the classic binomial log-odds slope `s` (logit change per generation) as a descriptive estimate. It also reports a quasi-binomial-adjusted p-value, and a within-evolving-replicate fixed-intercept sensitivity slope. Neither is a fully validated drift-aware SNP selection test. Do not interpret `s` as a classical population-genetic selection coefficient without the assumptions needed to equate them.
- When baseline ancestral pools do not belong to the same replicates as later samples, the within-replicate-intercept slope is estimated primarily from the evolved replicate time series; this is why pooled and within-replicate results may differ.
- SNPs with incompatible/missing GWAS REF/ALT metadata are excluded, never silently flipped. Frequency/count matrices are always explicitly matched by SNP ID and pool sample ID.
- Use `p_quasi_BH` in exploratory selection-overlap tables; avoid asserting formal genome-wide FDR if drift/calibration is not modeled. The q<0.01 and strict q<1e-5 selection screens are descriptive. The old script's `q<1e-5` is not a statistical necessity.
- The original `17_Manhattan_selected_clusters.R` also built **selection-only LD components** from a selection-selected SNP set. New script 10 does **not** recompute that distinct clustering. It instead performs direct positional overlap with GWAS QTLs to avoid adding another LD-threshold procedure before it is scientifically needed. If selection-only independent loci are a main manuscript statistic, that additional step remains to be done.
- Original 19's independent sorting, `%in%` filtering and `tapply()` row assumptions risked allele-frequency misalignment. New script 11 joins SNPs by ID and applies phase before summarizing selection slopes, using **one row per GWAS cluster** for descriptive comparisons.
- A cluster with more tested SNPs has a higher chance of containing a low q-value; counts/proportions are descriptive and **not** evidence of independent enrichment. A matched or permuted genomic null may be necessary for manuscript claims.

## What has NOT been empirically verified

- The user's datasets are not available here, and there is no R runtime in this environment. These are reviewed candidate scripts, not executed/validated outputs.
- Confirm the actual REF/ALT columns in `RData/Genotype_info_Cemee.RData`. The old pipeline placed them in the fourth and fifth metadata columns; script 09 logs that fallback.
- Confirm the `pop_label` isotype mapping and same-observation trait pairing against the current script 00. MCMCglmm is more computationally intensive than separate `lmer()` fits.
- Confirm that Quasi-binomial trend models are not used to overclaim selection significance; a drift-aware replicate-aware trajectory model is recommended if this is a central inferential test.
