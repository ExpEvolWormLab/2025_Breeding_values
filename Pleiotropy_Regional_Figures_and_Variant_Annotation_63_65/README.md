# Next manuscript scripts: regional figures and functional annotation

## Why these scripts

The three independent regional stories (I:5,883,574–7,996,414; IV around 8,892,663 / 10,499,135; V:6,046,042–7,278,103) need comparable figures. They should combine R5 GWAS effects, R5 drift-aware responses and **exact-site** Pool-seq trajectories, and should not claim that linked R5 groups are independent selected mutations. A separate transcript-consequence audit is needed before attributing biological effects to candidate genes.

## Scripts

- **63_Three_regions_manuscript_figures.R**: vector PDF + high-res PNG for each chromosome (three aligned panels), with completeness/QC CSV. Uses `45_All_R5_clusters_GWAS_evolution_predictions.csv` and `11_Clusters_selection_and_GWAS_09f.csv` and, when available, the output `28_ALT_read_counts_all_pools.csv` from the historical pipeline and `56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv` from chromosome V. The two exact-site series are not interchangeable with cluster median drift-aware slopes. If old I/IV trajectories are unavailable, panel C **explicitly shows missing data**. The script does **not** reconstruct Pool-seq counts from inferred effect sizes.
- **63b_Founder_genotypes_and_conditional_trait_effects.R**: optional supplementary ancestry/trait-effect panels from the historical script-26 founder dosage table, script-27 kinship-adjusted conditional effects, and 16 founder IDs. Needs the corresponding local outputs.
- **64_Annotate_focal_variants_WBcel235.R**: frozen exact-position candidate manifest, retrieve REF/ALT alleles from the original indexed CeMEE VCF with `bcftools`, verify against the **same WBcel235 FASTA** used to call variants, then use a matching Ensembl VEP cache (or matching indexed GFF3+FASTA) for transcript-level consequences. Includes `I:7,167,738; 7,491,656; 7,596,359`, `IV:8,892,663; 10,499,135`, and `V:6,322,080; 6,637,695; 6,826,421`. One gene can have multiple transcripts and distinct consequences. It does not infer causality or environmental fitness effect.
- **65_Focal_variant_evidence_matrix.R**: merges the frozen variant manifest, REF/ALT QC, VEP transcript effects, exact R5 group membership (if the SNP-phase file is present) and R5 evolutionary prediction summaries. Does not impute a variant's group by physical proximity.
- **make_preview_from_current_files.py** (optional, already run here): produces provisional figures from the summary CSVs in the current archive; exact I/IV trajectories were not available here.

## Local setup and execution

Copy the R scripts into `Oct_pipeline_cleaner/` (or adjust source paths). Run from your project root:

```r
source('Oct_pipeline_cleaner/63_Three_regions_manuscript_figures.R')
# Only when historical 26/27 outputs are available:
source('Oct_pipeline_cleaner/63b_Founder_genotypes_and_conditional_trait_effects.R')
```

For functional annotation, **edit `CFG$vcf` and `CFG$fasta`** near the start of script 64 first. Make sure `bcftools` is on PATH, and install a VEP cache or supply a matching bgzipped/indexed WormBase annotation GFF3. Then:

```r
source('Oct_pipeline_cleaner/64_Annotate_focal_variants_WBcel235.R')
source('Oct_pipeline_cleaner/65_Focal_variant_evidence_matrix.R')
```

## What to send next

1. `selection_results/manuscript_finalization/63_three_regions/`.
2. `selection_results/manuscript_finalization/64_variant_consequences/` and `65_focal_variant_evidence/`.
3. If chromosome I/IV panel C says no data: `historical_QTL_R5/deep_dive_extended/evolution/28_ALT_read_counts_all_pools.csv` plus `26_Founder_allele_signature_all_sites.csv`, `27_Effects_founder_then_variants_all8.csv`, and `27_Nested_model_tests_and_identifiability.csv`, where available.

## Functional interpretation priority

- **IV:10,499,135 (`tax-6`):** the earlier note placed it in an exon. Determine exactly whether this is synonymous, missense, splice-related or other *per transcript*. `tax-6` encodes calcineurin and has locomotor/neural associations, which makes it a strong hypothesis but not causal proof.
- **V:6,322,080 (`cest-24`):** the earlier note placed this intronically, near `T28C12.15`. Check transcript-specific intron boundaries and splice-region annotations. If intronic/noncoding, protein-impact tools should not be treated as evidence of deleteriousness.
- **I RC301/PX174 alleles:** determine genes overlapped by the shared and private SNP sets; the three SNPs in the manifest are illustrative anchors, not a complete functional annotation of the entire 2.1-Mb haplotype region.

## Important methodological caveat

**These scripts do not resolve the outstanding extreme-tail neutral-model calibration of the 12 BY candidate clusters.** The script-62 pilot identified potentially anti-conservative frequency classes and is too coarse at the ~10⁻⁵ level. Keep the 12 candidates conditional until a statistically validated extreme-tail approach is available. Regional figures are descriptive, not substitutes for valid p-values.

Code was inspected and the output logic was cross-checked with the current CSV schemas; the R scripts have **not** been executed in this environment. Python preview figures were generated from available summary data and were visually checked.
