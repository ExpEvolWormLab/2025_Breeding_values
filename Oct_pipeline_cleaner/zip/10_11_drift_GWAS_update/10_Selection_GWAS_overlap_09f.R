## 10_Selection_GWAS_overlap_09f.R
## Descriptive overlap of R5 GWAS SNPs with the one-G0 neutral-null GLS
## selection screen (09f). Primary testing family: ancestral p0 in [0.1, 0.9].
## ALL BH calls are provisional until extreme-tail/Pool-seq sampling validation.
## SNP overlap is not proof of independent loci, colocalization, or causation.

threshold_clust <- 0.5
label <- paste0("R", round(10 * threshold_clust))
Ne <- 1000L
q_cutoff <- 0.01
q_cutoff_sensitivity <- 0.05

selection_file <- file.path("selection_results", "drift_oneG0_nullfit",
                             paste0("09f_SNP_drift_Ne", Ne, ".csv"))
gwas_file <- file.path("txt_GridLMM_NaCl",
                        paste0("all_sign_GWAS_with_clusters_", label, ".txt"))
out_dir <- file.path("selection_results", "drift_oneG0_nullfit", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sel <- read.csv(selection_file, stringsAsFactors = FALSE)
gwas <- read.table(gwas_file, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
needed_sel <- c("ID", "rs", "slope_freq", "p0_hat", "p_drift",
                "p_drift_BH_common", "p_drift_BH_all", "common_G0",
                "n_pools", "Ne", "retained_G0", "excluded_G0")
needed_gwas <- c("rs", "group", "CHROM", "POS")
if (!all(needed_sel %in% names(sel)) || !all(needed_gwas %in% names(gwas)))
  stop("Missing required fields in script 09f or clustered GWAS output")
if (anyDuplicated(sel$ID) || anyDuplicated(gwas$rs))
  stop("Duplicate SNP IDs prevent unambiguous GWAS/Pool-seq matching")
if (length(unique(sel$Ne)) != 1L || unique(sel$Ne) != Ne)
  stop("Selected 09f file does not match Ne")
if (!identical(unique(as.character(sel$retained_G0)), "A00") ||
    !identical(unique(as.character(sel$excluded_G0)), "A6G140"))
  warning("Check retained/excluded G0 provenance; currently: ",
          paste(unique(sel$retained_G0), collapse = ","), " / ",
          paste(unique(sel$excluded_G0), collapse = ","))

## Important: GWAS rs is CHROM_POS, while 09f rs is CHROM.POS.
## Therefore match to the exact 09f ID column, not rs or row order.
i <- match(as.character(gwas$rs), as.character(sel$ID))
gwas$poolseq_covered <- !is.na(i)
gwas$selection_tested <- is.finite(sel$p_drift[i])
gwas$selection_common_G0 <- !is.na(sel$common_G0[i]) & sel$common_G0[i]
gwas$selection_tested_common <- gwas$selection_tested & gwas$selection_common_G0
gwas$selection_p0_ALT <- sel$p0_hat[i]
gwas$selection_slope_freq_ALT <- sel$slope_freq[i]
gwas$selection_p_drift <- sel$p_drift[i]
gwas$selection_q_common <- sel$p_drift_BH_common[i]
gwas$selection_q_all_PROVISIONAL <- sel$p_drift_BH_all[i]
gwas$selection_n_pools <- sel$n_pools[i]
gwas$selected_common_q01 <- gwas$selection_tested_common &
  is.finite(gwas$selection_q_common) &
  gwas$selection_q_common < q_cutoff
gwas$selected_common_q05 <- gwas$selection_tested_common &
  is.finite(gwas$selection_q_common) &
  gwas$selection_q_common < q_cutoff_sensitivity

sel$selected_common_q01 <- !is.na(sel$common_G0) & sel$common_G0 &
  is.finite(sel$p_drift_BH_common) & sel$p_drift_BH_common < q_cutoff
sel$selected_common_q05 <- !is.na(sel$common_G0) & sel$common_G0 &
  is.finite(sel$p_drift_BH_common) &
  sel$p_drift_BH_common < q_cutoff_sensitivity

summarize_group <- function(z) {
  data.frame(group = as.character(z$group[1]),
    CHROM = as.character(z$CHROM[1]),
    n_GWAS_SNPs = nrow(z),
    n_poolseq_covered = sum(z$poolseq_covered),
    n_selection_tested = sum(z$selection_tested),
    n_tested_common = sum(z$selection_tested_common),
    n_selected_common_q01 = sum(z$selected_common_q01),
    n_selected_common_q05 = sum(z$selected_common_q05),
    any_selected_common_q01 = any(z$selected_common_q01),
    any_selected_common_q05 = any(z$selected_common_q05),
    stringsAsFactors = FALSE)
}
groups <- split(gwas, as.character(gwas$group))
group_tab <- do.call(rbind, lapply(groups, summarize_group))
rownames(group_tab) <- NULL

summary <- data.frame(
  LD_threshold = threshold_clust, label = label, Ne = Ne,
  retained_G0 = paste(unique(sel$retained_G0), collapse = ";"),
  excluded_G0 = paste(unique(sel$excluded_G0), collapse = ";"),
  SNPs_tested_genomewide = sum(is.finite(sel$p_drift)),
  SNPs_tested_common_genomewide = sum(sel$common_G0 &
                                       is.finite(sel$p_drift), na.rm = TRUE),
  SNPs_selected_common_q01_genomewide = sum(sel$selected_common_q01),
  SNPs_selected_common_q05_genomewide = sum(sel$selected_common_q05),
  GWAS_SNPs_total = nrow(gwas),
  GWAS_SNPs_poolseq_covered = sum(gwas$poolseq_covered),
  GWAS_SNPs_tested_common = sum(gwas$selection_tested_common),
  GWAS_SNPs_selected_common_q01 = sum(gwas$selected_common_q01),
  GWAS_SNPs_selected_common_q05 = sum(gwas$selected_common_q05),
  GWAS_clusters_total = nrow(group_tab),
  GWAS_clusters_tested_common = sum(group_tab$n_tested_common > 0),
  GWAS_clusters_with_selected_common_q01 = sum(group_tab$any_selected_common_q01),
  GWAS_clusters_with_selected_common_q05 = sum(group_tab$any_selected_common_q05),
  stringsAsFactors = FALSE)

write.csv(sel[sel$selected_common_q01, , drop = FALSE],
          file.path(out_dir, "10_Selected_common_SNPs_q01_PROVISIONAL.csv"),
          row.names = FALSE)
write.csv(gwas, file.path(out_dir, "10_GWAS_SNPs_with_09f.csv"),
          row.names = FALSE)
write.csv(group_tab, file.path(out_dir, "10_GWAS_cluster_overlap.csv"),
          row.names = FALSE)
write.csv(summary, file.path(out_dir, "10_Overlap_summary.csv"),
          row.names = FALSE)
write.table(group_tab$group[group_tab$any_selected_common_q01],
            file.path(out_dir, "10_Selected_GWAS_clusters_q01.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
print(summary, row.names = FALSE)
cat("Selection calls are provisional; overlap is descriptive.\n")
