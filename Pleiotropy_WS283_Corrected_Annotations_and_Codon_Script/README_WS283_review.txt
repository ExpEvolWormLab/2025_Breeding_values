Corrected WS283 annotation audit (2026-11-03)

1) 67b_WS283_curated_variant_summary.csv:
  Eight focal variants, using ONLY source==WormBase gene/transcript features.
  The old 67_WS283_broad_variant_annotation.csv was overly liberal:
  it counted auxiliary CDS prediction tracks and chose old genetic-mapping
  intervals (gmap:*) as nearest genes.

2) 67b_WS283_transcript_feature_summary.csv:
  Transcript-level GFF3 overlap results (not VEP consequences).

3) 67b_Correct_WS283_curated_features.R:
  Reproduce the two corrected tables on the user's project with base R,
  using the existing four script-67 output files.

4) 68_Offline_WS283_CDS_codon_consequences.R:
  Optional R-only prediction of the codon and amino-acid consequences for
  three coding SNPs: I:7167738 (K07G5.4), IV:10499135 (tax-6),
  V:6637695 (scd-2). Edit CFG$fasta, CFG$gff. Run script 67b first.
  No VEP required. It reads full CDS exon chains from the user's exact
  WS283 GFF3 and genomic FASTA, checks both the reference allele and CDS
  phases, and emits transcript-specific inferred substitutions.
  Any inconsistent phase receives 'needs_manual_review'. No extrapolation
  to protein function or biological causality.

All output is based on the uploaded 67_WS283_local_annotations.zip;
no genome FASTA was supplied in this conversation, so actual reference
codons/amino acids were NOT computed here.
