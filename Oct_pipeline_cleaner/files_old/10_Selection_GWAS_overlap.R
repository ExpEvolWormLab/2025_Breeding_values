## 10_Selection_GWAS_overlap.R
## Refactors the overlap component of old 17 and 18.
## Uses SNP-level selection trajectories without refitting GWAS or recomputing LD.
## This is positional overlap, not independent causal-locus colocalization.

threshold_clust <- 0.5
label <- paste0("R", round(10 * threshold_clust))
q_cutoff <- 0.01               # prespecified descriptive screen
strict_cutoff <- 1e-5         # legacy highly stringent screen
selection_q_column <- "p_quasi_BH"  # more cautious than raw binomial q
out_dir <- file.path("selection_results", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sel <- read.csv("selection_results/09_Selection_SNP_results.csv",
                stringsAsFactors = FALSE)
gwas <- read.table(file.path("txt_GridLMM_NaCl",
  paste0("all_sign_GWAS_with_clusters_", label, ".txt")),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (anyDuplicated(sel$ID) || anyDuplicated(gwas$rs))
  stop("Duplicate positions or SNP identifiers; overlap would be ambiguous")

## Explicit SNP ID matching; do not merge on numeric row order.
idx <- match(gwas$rs, sel$ID)
gwas$selection_q <- sel[[selection_q_column]][idx]
gwas$selection_s_ALT <- sel$s[idx]
gwas$selection_pools <- sel$n_pools[idx]
gwas$poolseq_covered <- !is.na(idx)
gwas$selected_q01 <- is.finite(gwas$selection_q) &
                       gwas$selection_q < q_cutoff
gwas$selected_strict <- is.finite(gwas$selection_q) &
                         gwas$selection_q < strict_cutoff

## All selected SNPs genome-wide (incl. those outside GWAS clusters).
selected <- sel[is.finite(sel[[selection_q_column]]) &
  sel[[selection_q_column]] < q_cutoff, , drop = FALSE]
selected_strict <- sel[is.finite(sel[[selection_q_column]]) &
  sel[[selection_q_column]] < strict_cutoff, , drop = FALSE]

## Only count one GWAS group once for the overlap summary.
all_groups <- sort(unique(as.character(gwas$group)))
covered_groups <- unique(as.character(gwas$group[gwas$poolseq_covered]))
overlap_groups <- unique(as.character(gwas$group[gwas$selected_q01]))
strict_groups <- unique(as.character(gwas$group[gwas$selected_strict]))

summary <- data.frame(
  threshold = threshold_clust, label = label,
  selected_snps_genomewide = nrow(selected),
  selected_snps_strict = nrow(selected_strict),
  gwas_snps_total = nrow(gwas),
  gwas_snps_poolseq_covered = sum(gwas$poolseq_covered),
  gwas_snps_selected = sum(gwas$selected_q01),
  gwas_clusters_total = length(all_groups),
  gwas_clusters_poolseq_covered = length(covered_groups),
  gwas_clusters_with_selected_snp = length(overlap_groups),
  gwas_clusters_with_strict_snp = length(strict_groups)
)

write.csv(selected, file.path(out_dir, "10_Selected_SNPs.csv"),
          row.names = FALSE)
write.csv(gwas, file.path(out_dir, "10_GWAS_SNPs_with_selection.csv"),
          row.names = FALSE)
write.csv(summary, file.path(out_dir, "10_Selection_GWAS_overlap_summary.csv"),
          row.names = FALSE)
write.table(overlap_groups,
   file.path(out_dir, "10_Selected_GWAS_clusters.txt"),
   quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("GWAS SNPs covered by poolseq:", sum(gwas$poolseq_covered),
    "; groups containing selected SNP:", length(overlap_groups), "\n")

