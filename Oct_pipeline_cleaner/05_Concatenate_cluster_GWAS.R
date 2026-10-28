## 05_Concatenate_cluster_GWAS.R
## Formerly 10_Analyze_permutation_Gemma_clusters_with_residuals.R.
## Combine cluster-GWAS results from script 04, preserving all clusters.
## Run ONLY after 04_Mapping_clusters.R has finished the requested thresholds.

traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
chromosomes <- c("I", "II", "III", "IV", "V", "X")
thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
# For an initial run, use: thresholds_clust <- 0.5
threshold_label <- function(x) paste0("R", round(10 * x))
q_cutoff <- 0.05

meta_dir <- "data_BLUE_GridLMM_NaCl"
mapping_dir <- "GridLMM_clusters_NaCl"
output_dir <- "RData_GridLMM_NaCl"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mapping_qc_file <- file.path(mapping_dir, "04_Cluster_mapping_QC.csv")
mapping_qc <- if (file.exists(mapping_qc_file))
  read.csv(mapping_qc_file, stringsAsFactors = FALSE) else NULL

all_mapping_summary <- NULL

for (threshold_clust in thresholds_clust) {
  label <- threshold_label(threshold_clust)
  cat("\nCombining cluster GWAS:", label, "(R2 >", threshold_clust, ")\n")

  metadata_file <- file.path(meta_dir,
    paste0("data_clusters_BLUE_", label, ".txt"))
  df_clusters <- read.table(metadata_file, header = TRUE,
                            stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("clusterID", "group", "CHROM", "POS")
  if (!all(needed %in% names(df_clusters)))
    stop("Missing metadata fields at ", label)
  if (anyDuplicated(df_clusters$clusterID))
    stop("Duplicated cluster IDs at ", label)

  ## Preserve legacy column names expected by scripts 13, 19, etc.
  all_GWAS_clusters <- df_clusters
  all_GWAS_clusters$rs <- as.character(df_clusters$clusterID)
  all_GWAS_clusters$chr <- as.character(df_clusters$CHROM)
  all_GWAS_clusters$ps <- as.numeric(df_clusters$POS)
  all_GWAS_clusters <- all_GWAS_clusters[
    , c("rs", "chr", "ps", setdiff(names(all_GWAS_clusters),
                                     c("rs", "chr", "ps"))), drop = FALSE]
  all_GWAS_clusters <- all_GWAS_clusters[
    order(match(all_GWAS_clusters$chr, chromosomes), all_GWAS_clusters$ps),
    , drop = FALSE]
  rownames(all_GWAS_clusters) <- NULL

  for (traitID in traits) {
    tag <- paste0(traitID, "_NaCl")
    trait_files <- character()
    trait_chromosomes <- character()
    trait_expected <- integer()

    for (chrID in chromosomes) {
      out_file <- file.path(mapping_dir, label,
        paste0("gwas_results_", label, "_", traitID,
               "_CHR_", chrID, "_NaCl.txt"))

      ## If script 04 tested no clusters on a chromosome, no file is expected.
      q <- if (!is.null(mapping_qc)) {
        subset(mapping_qc, abs(threshold - threshold_clust) < 1e-8 &
               chromosome == chrID)
      } else NULL
      is_expected <- if (!is.null(q) && nrow(q) == 1L) {
        q$n_tested > 0
      } else {
        any(all_GWAS_clusters$chr == chrID)
      }
      if (is_expected && !file.exists(out_file))
        stop("Missing cluster-GWAS output: ", out_file,
             "\nFinish script 04 before running script 05.")
      if (!is_expected && file.exists(out_file))
        stop("Stale GWAS file for a chromosome with zero tested clusters: ", out_file)
      if (file.exists(out_file)) {
        trait_files <- c(trait_files, out_file)
        trait_chromosomes <- c(trait_chromosomes, chrID)
        trait_expected <- c(trait_expected,
                            if (!is.null(q) && nrow(q) == 1L) q$n_tested else NA_integer_)
      }
    }

    traitJ <- if (length(trait_files)) {
      do.call(rbind, lapply(seq_along(trait_files), function(k) {
        f <- trait_files[k]
        z <- read.table(f, header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE)
        fields <- c("X_ID", "beta.5", "p_value_REML")
        if (!all(fields %in% names(z))) stop("Unexpected GWAS columns in ", f)
        if (!is.na(trait_expected[k]) && nrow(z) != trait_expected[k])
          stop("GWAS rows disagree with script 04 QC: ", f)
        match_chr <- match(as.character(z$X_ID), all_GWAS_clusters$rs)
        if (anyNA(match_chr) ||
            any(all_GWAS_clusters$chr[match_chr] != trait_chromosomes[k]))
          stop("Cluster IDs do not belong to the expected chromosome: ", f)
        z[, fields, drop = FALSE]
      }))
    } else {
      data.frame(X_ID = character(), beta.5 = numeric(),
                 p_value_REML = numeric())
    }

    ids <- as.character(traitJ$X_ID)
    if (anyDuplicated(ids) || !all(ids %in% all_GWAS_clusters$rs))
      stop("Duplicate or unknown cluster IDs in ", label, " / ", traitID)
    bad_p <- !is.na(traitJ$p_value_REML) &
      (!is.finite(traitJ$p_value_REML) |
         traitJ$p_value_REML < 0 | traitJ$p_value_REML > 1)
    if (any(bad_p))
      stop("Invalid cluster GWAS p-values at ", label, " / ", traitID)

    ## Match by ID, never by file row order or implicit merge.
    row <- match(all_GWAS_clusters$rs, ids)
    beta <- traitJ$beta.5[row]
    p <- traitJ$p_value_REML[row]
    qval <- rep(NA_real_, length(p))
    tested <- which(is.finite(p))
    qval[tested] <- p.adjust(p[tested], method = "BH")
    significant <- !is.na(qval) & qval < q_cutoff

    all_GWAS_clusters[[paste0("beta_", tag)]] <- beta
    all_GWAS_clusters[[paste0("p_value_REML_", tag)]] <- p
    all_GWAS_clusters[[paste0("p_value_REML.cor_", tag)]] <- qval
    all_GWAS_clusters[[paste0("is_sign_", tag)]] <- significant

    all_mapping_summary <- rbind(all_mapping_summary,
      data.frame(threshold = threshold_clust, label = label,
                 trait = traitID, clusters_total = nrow(all_GWAS_clusters),
                 clusters_tested = length(tested),
                 clusters_untested = nrow(all_GWAS_clusters) - length(tested),
                 q_significant = sum(significant),
                 min_p_value = if (length(tested)) min(p[tested]) else NA_real_))
  }

  ## The same object name and output path consumed by legacy downstream code.
  output_file <- file.path(output_dir,
    paste0("GWAS_cluster_results_BLUE_concatenated_", label,
           "_with_residuals.RData"))
  save(all_GWAS_clusters, threshold_clust, file = output_file)
  write.csv(all_GWAS_clusters, file.path(output_dir,
            paste0("Cluster_GWAS_summary_", label, ".csv")), row.names = FALSE)
  cat("  Clusters:", nrow(all_GWAS_clusters), "\n")
}

write.csv(all_mapping_summary, file.path(output_dir,
          "05_Cluster_GWAS_QC_by_threshold_trait.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "05_sessionInfo.txt"))
cat("\nSaved cluster-GWAS summaries in", output_dir, "\n")
