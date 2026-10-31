# Phase 11: local-QTL competition and manuscript-example figures

These scripts are designed to answer the next question in the paper:
**within a focal region, do multiple local alleles/haplotypes appear to be favored for different trait effects?**

## Recommended workflow

1. Run `30_Build_extended_region_tables.R`
   - This keeps your historical regions.
   - It adds a first set of newly suggested QTL examples:
     - chromosome III groups 504 and 506 together (both loaded from the complete script-29 R5 ranking)
     - chromosome V group 886
     - chromosome X group 1010
     - chromosome II group 249

2. Re-run the regional deep-dive with the extended tables.
   - Use the patched script copies provided here if you want separate outputs,
     or edit your own copies of scripts 26/27/28.

3. Run `31_Local_QTL_competition_summary.R`
   - Main outputs:
     - `31_Term_trait_profile_summary.csv`
     - `31_Within_region_competition_summary.csv`
     - `31_Region_priority_ranking.csv`

4. Run `32_Local_QTL_competition_figures.R`
   - Main figure PDFs:
     - `32_Figure_heatmaps_effect_profiles.pdf`
     - `32_Figure_beta_decomposition.pdf`
     - `32_Figure_pairwise_competition_summary.pdf`
     - `32_Figure_selected_trajectories_preview.pdf`

## Interpretation logic

A pair of local terms is most interesting when several things happen together:
- they remain distinguishable after conditioning on one another;
- they have different multivariate trait-effect profiles;
- their strongest effects are on different traits, or their beta effects point in opposite directions;
- they also have evolutionary support (exact-site trajectories or selected R5-group summaries).

This will be especially informative for:
- `chrIV_CB4856`
- `beta_I`
- `beta_V`
- `beta_III_pair` (new)

## Notes

- The ancestry terms are still **proxy founder similarities**, not phased local ancestry.
- For some new regions the competition analysis may reveal that one nominal QTL is enough; for others it may indicate a real multi-allelic or multi-haplotype local system.
- I kept the code architecture close to your earlier scripts and the outputs are simple csv/pdf files.

## Fixes in the 2026-10-29 revision

- Script **30a** now maps **all eight** target traits to `pvalue_T1`–`pvalue_T8`, fixing the `argument is of length zero` error when new regions use traits such as `FB` or `FS`.
- Script **30** reads `29_All_R5_QTL_evidence_score.csv`, rather than the spatially thinned shortlist. This retains group **506** (which the thinned shortlist omits).
- Script **30** no longer treats rounded R5 centroids as actual focal variants. Every model tag comes from an observed SNP in the R5 group.
- Script **31** now orders the regional competition ranking from highest to lowest score.
- The scores describe **opposing association profiles**; they do not demonstrate within-region competition, genetic epistasis or simultaneous selection on distinct causal alleles. Interpreting evolutionary competition requires distinct allelic/haplotype trajectories and validated phase/orientation.

**Re-run script 30**, then 30a, 30b and 30c in order. 30b–30c have not changed in this correction.
