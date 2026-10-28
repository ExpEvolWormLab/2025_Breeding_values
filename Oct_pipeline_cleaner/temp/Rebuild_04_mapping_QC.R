## Rebuild chromosome-level mapping QC without rerunning GridLMM.
## Run from the project root after scripts 03, 04 and 05.
## Matches script 04 when missing_policy == "exclude".

thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
chromosomes <- c("I", "II", "III", "IV", "V", "X")
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
threshold_label <- function(x) paste0("R", round(10 * x))

out_file <- "GridLMM_clusters_NaCl/04_Cluster_mapping_QC.csv"
qc_rows <- list()

for (threshold_clust in thresholds_clust) {
  label <- threshold_label(threshold_clust)
  qc_file <- file.path("data_BLUE_GridLMM_NaCl",
                       paste0("Cluster_genotypes_QC_", label, ".csv"))
  gqc <- read.csv(qc_file, stringsAsFactors = FALSE)
  required <- c("clusterID", "CHROM", "missing_RILs", "sd_genotype")
  if (!all(required %in% names(gqc)))
    stop("Unexpected genotype-QC columns at ", label)
  if (anyDuplicated(gqc$clusterID))
    stop("Duplicate cluster IDs at ", label)

  for (chrID in chromosomes) {
    z <- gqc[gqc$CHROM == chrID, , drop = FALSE]
    if (!nrow(z)) next  # Same behavior as script 04

    missing <- z$missing_RILs > 0 | is.na(z$missing_RILs)
    nonconstant <- is.finite(z$sd_genotype) & z$sd_genotype > 0
    tested <- !missing & nonconstant

    row <- data.frame(
      threshold = threshold_clust,
      chromosome = chrID,
      n_clusters = nrow(z),
      excluded_missing = sum(missing),
      excluded_no_variation = sum(!missing & !nonconstant),
      n_tested = sum(tested)
    )

    ## Verify the actual mapping outputs before trusting reconstructed QC.
    for (traitID in traits) {
      f <- file.path("GridLMM_clusters_NaCl", label,
        paste0("gwas_results_", label, "_", traitID,
               "_CHR_", chrID, "_NaCl.txt"))
      if (row$n_tested == 0L) {
        if (file.exists(f)) stop("Unexpected mapping file for zero tested clusters: ", f)
      } else {
        if (!file.exists(f)) stop("Missing cluster-GWAS file: ", f)
        lines <- readLines(f, warn = FALSE)
        if (length(lines) - 1L != row$n_tested)
          stop("Row count mismatch at ", f, ": expected ", row$n_tested,
               ", found ", length(lines) - 1L)
      }
    }
    qc_rows[[length(qc_rows) + 1L]] <- row
  }
}

mapping_qc <- do.call(rbind, qc_rows)
rownames(mapping_qc) <- NULL

## Only write after ALL thresholds and all mapping files have passed checks.
if (file.exists(out_file)) {
  backup <- paste0(out_file, ".previous")
  if (!file.copy(out_file, backup, overwrite = TRUE))
    stop("Could not back up original mapping QC")
  cat("Previous mapping QC saved to:", backup, "\n")
}
write.csv(mapping_qc, out_file, row.names = FALSE)
cat("Reconstructed mapping QC saved to:", out_file, "\n")
print(table(mapping_qc$threshold))
