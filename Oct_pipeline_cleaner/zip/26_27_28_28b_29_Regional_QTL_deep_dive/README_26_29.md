# Regional QTL deep dive: founder haplotypes, conditional effects, evolution, and new candidates

## Essential result from the supplied CSV tables (before running these scripts)

- The three historical beta-V positions are **separate current R5 clusters**: `V:6,322,080` → R5 **848** (GWAS beta q = 1.6384e-13); `V:6,637,695` → **851** (q = 6.3260e-20); `V:6,826,421` → **853** (q = 8.7851e-13). These are cluster-level q values, not the individual-SNP q values.
- Founder VCF signatures from script 24: `V:6,322,080` (T→C) has ALT in **0/16 founders**; `V:6,637,695` (G→A) and `V:6,826,421` (T→A) each have ALT in **CB4507 only**. Do not assume that absence from the sampled founders dates the first mutation without historical sequencing and read-depth verification.
- **None of these three focal positions has a defined common-frequency Pool-seq selection q-value** in the existing analysis. Within the V region, R5 **848** has minimum selected-SNP q ≈ **0.06861**; R5 **851** and **853** have *no common-frequency SNP*, and R5 **854** (different SNPs, 6.91–7.01 Mb) has minimum q ≈ **0.01014**. We also include broad background cluster **844** (beta q ≈ 0.0114; minimum common-SNP selection q ≈ 0.0300) in the conditional chromosome-V analysis; it spans much of the region and could tag a shared haplotype component. These values are not formal tests of selection at the three focal V SNPs.
- Prioritize paired interpretation of R5 **I:53/54** (PX174+RC301-shared versus RC301-specific; both overlap Pool-seq response), **IV:725/730 plus 727/731**, **III:598**, **X:1087**, and the incomplete common-frequency coverage of **II:128/317**.

## Files and setup

Put `26_Regions_and_R5_groups.csv`, `26_Focal_variants.csv`, and the four `.R` scripts into the *R project root* alongside script 23. The CSV configuration contains the current R5 clusters and the explicitly nominated sites. **Group IDs depend on the exact R5 pipeline version; stop and audit if rebuilt.** All data files must be on the same reference assembly and allele orientation; do not use approximate coordinate liftover without documenting it.

**26** — Edit the `vcf_file` path and, if needed, `chrom_prefix` at the top of `26_Regional_founder_haplotypes.R`. Use the exact 16 VCF sample IDs in `18_founder_sample_ids.txt`. Requires `bcftools` and the script-00 RIL genotype matrices. Extraction is limited to focal sites, representative R5 sites and ±100 kb informative local-marker windows, excluding the named focal and tag SNPs from ancestry proxy inference. This is **not phased haplotype reconstruction**: founder-identical tracts and inaccurate genotype coding can leave classifications ambiguous. Validate genotype dosage alignment and the local mismatch margin; inspect `26_Local_founder_panel_QC.csv`, `26_Founder_allele_signature_all_sites.csv`, and `26_RIL_founder_similarity_PROXIES.csv` before treating those assignments as useful covariates. Founder-allele sharing is not a mutation-age estimate.

**27** — Requires `AGHmatrix` and script-26 outputs. Fits all **seven component traits and beta** using chromosome-excluded RIL kinship and documented RIL-origin covariates where supplied. Within each region:

- Stage A: founder-like haplotype associations adjusted for source population.
- Stage A0: **founder-like association among RILs with neither focal allele**, when sufficient cells exist (only in regions with explicitly named variants).
- Stage B: nominated site/R5 representative associations conditional on founder proxies and other tagged sites, with joint nested-model improvement tests.
- Screens individual marker×marker and marker×founder interactions only with all four cells sufficiently populated. Reports rank-deficient terms instead of assigning arbitrary coefficients, and calculates approximate 80%-power effect thresholds.
- Decomposes `beta` into weighted component-trait effects *only when exact reconstructed phenotype-level beta weights check out*. This is an index decomposition, **not a causal fitness model**.

Results: `historical_QTL_R5/deep_dive/regional_models/27_Effects_founder_then_variants_all8.csv`, `27_Nested_model_tests_and_identifiability.csv`, `27_Pairwise_epistasis_all8.csv` (if testable), `27_Beta_weighted_trait_contributions.csv` (if valid), `27_Regional_forest_plots.pdf`. Variance-ratio GLS is an approximation, and a local ancestry proxy may not capture the specific haplotype or local LD surrounding the variant. A main mutation association conditional on a CB4856-like tag does not prove causality. Do not interpret unestimable mutation×ancestry as evidence of no epistasis.

**28** — Uses the **corrected aligned** script-08 Pool-seq RData and `09_Selection_SNP_results.csv` for GWAS ALT direction, plus existing full-GA and retained-GA selection summaries. Extracts counts for exact nominated SNPs, R5 tags and optionally **all 504 historical positions** (`include_all_historical_positions=TRUE`). Exports both G0 pools; preserves A00 as the primary baseline and *does not fit* the bottlenecked A6G140 to the selection curve. Exports `28_ALT_read_counts_all_pools.csv`, `28_Per_GA_trajectory_summary.csv`, `28_Exact_variant_GWAS_and_selection_summary.csv`, `28_R5_region_full_vs_retained_GA_overlap.csv`, and a focal-trajectory PDF. Where a SNP is not in the script-08/09 SNP universe, it cannot be interpreted and requires direct extraction from the original sync or an explicitly reconstructed count matrix. A SNP with p0 outside 0.10–0.90 **does not have a comparable common-filter q-value**; show its frequency trajectory descriptively rather than claiming it failed a test.

**29** — Uses the R5 cluster metadata from script 23, and optionally script 20 survivor cluster summaries, to propose additional case-study QTLs across **all six chromosomes**. Rank score rewards: strong beta or area.F GWAS q, multi-trait q values, full-GA drift-aware selected SNPs, selected SNPs after excluding GA3, sign stability of the cluster-level slope, GWAS-effect/response directional consistency, and at least 3 common SNPs. Excludes the historically nominated regions ±200 kb from the *novel* shortlist and thins nominees within 400 kb to avoid nearly identical loci. Output `29_Novel_QTL_shortlist_spatially_thinned.csv` gives transparent reasons. Ranking is **for hypothesis generation only**: selection q-values are provisional in the far WF tail, many SNPs share LD, and selection of interesting loci based on the same data makes p-values and effect sizes descriptive, not independent confirmation.

## Run order

```r
source('26_Regional_founder_haplotypes.R')
source('27_Regional_haplotype_mutation_GLS_all8.R')
source('28_Regional_Poolseq_trajectories_and_selection.R')
source('29_Rank_additional_QTL_case_studies.R')
```

For faster initial validation, run **26**, inspect local ancestry QC, then **28 and 29**; run the heavier **27** after checking the founder panels. For memory efficiency, use separate R sessions for scripts 26 and 27 if the genotype matrices are large. All outputs go to `historical_QTL_R5/deep_dive` and never overwrite scripts 19–25. Optional sources (`19_RIL_source_populations.csv`, survivor-selection tables) should be included if they exist. The GA3 status is unknown unless experimental records confirm extinction; `retained-GA` means GA1/2/4, **not necessarily only objectively surviving populations**.

## Upload for joint interpretation and paper-ready case studies

1. `26_Founder_allele_signature_all_sites.csv`, `26_Local_founder_panel_QC.csv`, `26_Local_founder_similarity_counts.csv`.
2. `27_Effects_founder_then_variants_all8.csv`, `27_Nested_model_tests_and_identifiability.csv`, `27_Pairwise_epistasis_all8.csv` if produced, and `27_Beta_weighted_trait_contributions.csv` if produced.
3. `28_Exact_variant_GWAS_and_selection_summary.csv`, `28_Per_GA_trajectory_summary.csv`, and `28_R5_region_full_vs_retained_GA_overlap.csv`.
4. `29_Novel_QTL_shortlist_spatially_thinned.csv`.

Or ZIP the complete `historical_QTL_R5/deep_dive/` folder. The original full genotype matrix and sync do not need to be uploaded.

## Scientific caveats

- A cluster with both GWAS and selected SNPs may tag **different causal sites** in LD. Need conditional fine-mapping, haplotype comparison and ideally experimental edits to argue for molecular pleiotropy or causality.
- QTL group sign is aligned to the phase-orientation of *cluster genotype*. An ALT-coded individual focal SNP can have an opposite sign. Compare like with like and use `R5_phase` to translate when available.
- Founder-specific ALT does not imply its first origin in the exact founder; sampling, reference coding and subsequent mutations must be considered.
- Modeling `beta` across its component traits does not demonstrate nonlinear fitness or direct epistasis; the separate script-22 selection model addresses SB×area.F, with its previously discussed power and relatedness caveats.
- Excluding incomplete populations is a **sensitivity** analysis, not automatically a less biased selection estimate. Evolution of extinct populations is biologically relevant.

**Runtime validation**: the file-generation environment has no R interpreter and no access to the private project genotype and Pool-seq RData. These scripts must be executed in your R environment and their outputs inspected before manuscript inferences are finalized.

**28b (optional but recommended for the three V SNPs)** — `28b_Exact_SNP_raw_SYNC_recovery.R` is a direct-read-count fallback for **all eight nominated exact variants and current R5 tags** using the original `Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync` and the verified script-08 sample manifest. It recovers positions that were removed by script 08's *founder-SNP intersection* or its coverage filter. Run it after 26 when the corresponding `28_SNP_coverage_status.csv` contains `FALSE`, or proactively to audit the three V variants. It outputs `28b_Exact_raw_SYNC_ALT_trajectories.csv`, `28b_Per_GA_raw_SYNC_changes.csv` and `28b_Raw_SYNC_exact_locus_QC.csv`. It **does not** assign genome-wide significance to previously untested low-frequency variants. Audit VCF/SYNC reference-base agreement before interpretation.
