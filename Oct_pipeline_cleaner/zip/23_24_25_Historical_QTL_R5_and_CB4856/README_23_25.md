# Historical QTL re-identification in the phase-aligned R5 GWAS pipeline

## Scope and scientific interpretation

The coordinate lists provided for the previous manuscript have been transcribed **verbatim** to `23_Historical_region_SNPs.csv`: 504 legacy coordinates in six named blocks (`beta_I`, `beta_V`, `area_II_1`, `area_II_2`, `area_III`, `area_X`) plus the two candidate chromosome-IV mutation positions, for 506 rows in total. The four area.F blocks are kept separately even though the historical text referred to three area QTLs: **only the new clustering evidence can justify any merge.** The coordinates and labels are project-specific; verify that the original and new GWAS files use the same genome assembly.

**R5 QTLs are not defined by the submitted intervals.** In scripts 02–05, only SNPs reaching the revised significant-SNP inclusion criterion (at least one phenotype BH q < 0.2) are LD-clustered, and clusters are split according to phase-aligned GWAS effect signs. A historical SNP may therefore be present in the genotype VCF or in the full SNP GWAS but missing from the R5 cluster table. That absence is informative about the new threshold, **not evidence that the SNP was absent or biologically inactive**. Region overlap by bounds or physical proximity must never be substituted silently for an exact SNP–cluster match.

## Run from the project root

Copy `23_Historical_region_SNPs.csv` and the R files alongside the other project scripts.

```r
source('23_Rediscover_historical_QTL_in_R5.R')
```

This requires:

- `txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt` (script 02)
- `RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R5_with_residuals.RData` (script 05)

It uses these inputs if they are present:

- `GridLMM_results_NaCl/All_Gwas.RData` (all tested SNPs, not only R5 members)
- `selection_results/drift_oneG0_nullfit/09f_SNP_drift_Ne1000.csv` (common-frequency drift selection)
- `selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv` (cluster response / SNP-selection counts)
- `data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv` (phase orientation)

Main outputs in `historical_QTL_R5/`:

- `23_Historical_SNP_status.csv`: **one row per historical coordinate**, exact R5 membership, tested-SNP status, original SNP-level p/q, nearest R5 SNP distance (clearly labeled PROXY), and optional drift selection.
- `23_Historical_region_summary.csv`: number of original SNPs, exact R5 matches, current R5 SNPs falling in the historical coordinate span, and candidate R5 QTL counts.
- `23_Historical_region_to_R5_QTL_candidates.csv`: **one row per region × current R5 cluster**, with match tier, original trait effect/q, all eight GWAS traits, and optional evolutionary selection support. Multiple rows per historical region are expected.
- `23_Candidate_R5_cluster_member_SNPs.csv`: exact current SNP membership and phase for matched clusters.
- `23_Historical_to_R5_region_maps.pdf`: region-by-region visual comparison. Teal = exact overlap; gold = interval only; grey = nearest proxy only.

The output does **not** automatically merge neighboring R5 groups. Multiple groups may reflect independent loci, subhaplotypes, trait-effect direction or resolution limits of LD-based grouping.

### Founder-allele composition, especially PX174/RC301

Run after 23, and configure `vcf_file`, `founder_ids_file`, `px_name` and `rc_name` in script 24. The indexed VCF and `bcftools` are required.

```r
source('24_Founder_alleles_historical_regions.R')
```

Output in `historical_QTL_R5/founder_alleles/`:

- `24_Historical_SNP_founder_genotypes_long.csv`: actual 16 founder genotypes at each historical position, where retrievable.
- `24_Historical_SNP_founder_signature_classes.csv`: PX174-only, RC301-only, shared ALT or neither, other-founder ALT carriers, R5 group and optional selection evidence.
- `24_Founder_signature_counts_by_region.csv`: founder-sharing summary by original region.

The classified **ALT** allele is relative to the current VCF REF/ALT; the biologically relevant *derived*, *associated*, or *causal* allele need not be ALT. These are allele-sharing summaries, not complete founder haplotype inference, and identical PX174/RC301 alleles may derive from shared ancestry.

### Chromosome IV: CB4856-like segment versus mutation

Run after 23 and the already completed script 21. You do not need to refit the existing kinship-adjusted chromosome-IV models.

```r
source('25_CHRIV_CB4856_haplotype_QTL_and_mutations.R')
```

This requires the script-21 outputs and (preferably) the phase-aligned cluster genotype matrix from script 03, plus optionally script 20 survivor-only selection results. Search windows are ±750 kb from `IV:8,892,663` and `IV:10,499,135`, but only R5 clusters with **actual constituent SNPs** in those windows are included; a broad cluster's bounding interval alone does not qualify.

Output in `historical_QTL_R5/chrIV_CB4856/`:

- `25_R5_chrIV_QTL_selection_CB4856_vs_mutations.csv`: R5 clusters near the mutations, SB/beta GWAS q-values, all-population and survivor-only selection evidence, and genotype correlations with each mutation and local CB4856-like ancestry. A cluster can be screened for correlation with CB4856-like ancestry **among RILs carrying neither mutation**. The `candidate_CB4856_background_QTL_without_mutations` flag uses an exploratory absolute correlation > 0.5; it is a *screen, not a formal QTL test*.
- `25_StageA_haplotype_StageB_conditional_mutations_all8.csv`: reorganizes *already fitted* script-21, eight-trait coefficients into two hierarchical stages: (A) CB4856 local ancestry H1/H2 without mutation covariates (which **can absorb mutation effects**); (B) mutation M1/M2 coefficients after adjustment for both local haplotype proxies and each other.
- `25_Nested_haplotype_then_mutation_tests.csv`: kinship-adjusted nested-model p-values already computed by script 21, including improvement when adding mutations to a haplotype model.
- `25_CB4856_and_two_mutation_genotype_cells.csv`: available classes; lack of mutation carriers on non-CB4856 backgrounds makes mutation × background interactions unidentifiable, particularly at mutation 2.

**Strongest discriminating contrast**: a mutation-carrier versus a noncarrier within CB4856-like RILs. The separate CB4856 effect can be estimated using RILs lacking the mutation, but statistical separation from other linked variants is not causal identification. The *stage A* model is an overall background association, not a mutation-free-only ancestry model. Stage-B mutation associations remain potentially confounded with a more specific linked background not captured by the binary ancestry proxy. Existing `21` p-values use a null-fitted kinship covariance approximation; they are not exact gene-editing contrasts.

**Selection overlap caveat**: `09f`/`20` genome-wide significance tests use a common ancestral-frequency filter. Variants appearing at low frequency after domestication (including the two focal mutations) may be filtered out. A cluster with no selected common-frequency SNPs is **not** demonstrated to lack selection; use its continuous phase-aligned slope and the individual mutation time series as separate evidence.

## What to upload for interpretation

Upload these four files first:

1. `historical_QTL_R5/23_Historical_region_summary.csv`
2. `historical_QTL_R5/23_Historical_region_to_R5_QTL_candidates.csv`
3. `historical_QTL_R5/founder_alleles/24_Historical_SNP_founder_signature_classes.csv`
4. `historical_QTL_R5/chrIV_CB4856/25_R5_chrIV_QTL_selection_CB4856_vs_mutations.csv`

Optionally add the `25_StageA_haplotype_StageB_conditional_mutations_all8.csv` and `25_Nested_haplotype_then_mutation_tests.csv` tables. These will allow us to name the *current* clusters for every historical locus, identify beta/area.F/SB shared associations, examine credible ancestry signatures, and decide which genomic regions merit finer joint modeling rather than treating the entire old interval as a single QTL.

## Validation

R is not available in the artifact-generation environment. The scripts have been checked against the existing pipeline's file/object names and undergone static delimiter checks; they still need to be run against your project data. No biological overlap results are fabricated here.
