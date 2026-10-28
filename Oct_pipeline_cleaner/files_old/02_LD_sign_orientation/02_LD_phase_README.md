# LD phase alignment for script 02

Use **02_LD_sign_diagnostic.R** first, sourced after section 2 of
`02_Clustering_from_LD.R`, when `ld_snps`, `all_GWAS_sign`, and
`vect_P_traits` exist. It checks that the sign of VCFtools haplotype D
matches the signed RIL genotype correlations calculated from the exact
`genotype_matrix` supplied to GridLMM. The output is a diagnostic, not
a clustering result.

If the two signs agree closely, `02_Clustering_from_LD_phase_aligned.R`
implements the optional phase-alignment approach. It (1) validates D against
the GWAS genotype matrix on a reproducible 500-edge subsample, (2) derives
relative allele orientations within each LD connected component along the
maximum-R2 spanning tree, (3) changes **only a temporary copy** of beta_T*
used by `split_clusters()`, and (4) writes diagnostic reports of disagreeing
LD links. Raw GWAS beta values retained in outputs are unchanged.

`phase_align_effects <- TRUE` enables alignment. Set FALSE to retain the
original splitter. To test only the primary threshold, set
`thresholds_clust <- 0.5` near the beginning.

Cautions:
- This is a proposed method; no project data or R runtime was available
  for execution or numerical validation.
- Signed D is allele-coding dependent. If D and GWAS genotype correlations
  disagree, DO NOT multiply beta by sign(D): investigate VCF coding and
  sample correspondence.
- The signs inferred along a spanning tree may disagree with additional LD
  edges, especially in chained/heterogeneous components. Inspect the
  `phase_alignment_R5.csv` report. This is not statistical evidence that
  those SNPs define one causal locus.
- This changes a scientific clustering rule and may change cluster counts.
  The downstream numeric columns/output files are preserved.
