REGIONAL MANUSCRIPT REVIEW — 03 NOV 2026

Contents:
  66_Rebuild_three_region_integrated_figure.R
  67_WS283_GFF3_local_overlap_annotation.R
  66_Three_regions_integrated_PRELIMINARY.pdf
  66_Three_regions_integrated_PRELIMINARY.png
  66_focal_allele_ancestry_and_joint_tests.csv
  66_nested_regional_tests.csv

Use the PDFs as a first evidence-aware prototype, not a final submitted figure.
The R script 66 reproduces the analysis natively from script-63/63b outputs
in your local project. The prototype PDF was generated independently in Python
from the files uploaded in the current conversation; the R layout may differ.

Read files:
  63_three_regions/63_*_R5_prediction_window.csv
  63_three_regions/63b_Chr*_conditional_effects_plotted.csv
  63_three_regions/63_*_focal_rare_allele_trajectories.csv
  historical_QTL_R5/deep_dive/regional_models/26_Founder_allele_signature_all_sites.csv
  historical_QTL_R5/deep_dive/regional_models/27_Nested_model_tests_and_identifiability.csv
  historical_QTL_R5/deep_dive/regional_models/28_ALT_read_counts_all_pools.csv (optional fallback)

No exact-site chr IV records for the two focal positions in the supplied
28_ALT_read_counts_all_pools.csv. This is NOT the absence of selection or change.

Local GFF annotation, without VEP:
  Set CFG$gff to path of c_elegans.PRJNA13758.WS283.annotations.gff3.gz.
  Set CFG$variants to the script-65 CSV if needed.
  source('Oct_pipeline_cleaner/67_WS283_GFF3_local_overlap_annotation.R')
  No VEP, gzip uncompression, tabix, or online website needed.
  Use only BROAD gene/exon/CDS overlap; it does not calculate amino acid changes.

Phenotypic test note:
  chr I: joint SNP-terms after founder proxies explain SB, FS, FB, BS, BF;
         beta q=0.1508 (not significant).
  chr IV: chrIV_CB4856 F10499135 has strongest conditional SF and SB effects,
         both q_BH_all_regional_terms=0.000408.
  chr V: joint regional SNP terms explain SF, FB and beta, while individually
         resolved variant effects are less convincing after full correction.

Critical limitation: the current neutral-model selection-test p-value tail
has not been established to be calibrated at the BY-significant threshold.
These figures should not declare causal loci or independently selected mutations.
