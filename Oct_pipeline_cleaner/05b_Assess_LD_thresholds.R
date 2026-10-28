## 05b_Assess_LD_thresholds.R
## Run AFTER scripts 02b, 03, 04 and 05 have completed all requested cutoffs.
## Descriptive LD-resolution, genotype-QC and cluster-GWAS sensitivity summary.
## Does NOT select a cutoff by maximizing significance or in-sample R2.

thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
threshold_label <- function(x) paste0("R", round(10 * x))
primary_threshold <- 0.5  # chosen provisionally before downstream interpretation
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")

out_dir <- "threshold_review"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mapping_qc <- read.csv("GridLMM_clusters_NaCl/04_Cluster_mapping_QC.csv",
                       stringsAsFactors = FALSE)
association_qc <- read.csv(
  "RData_GridLMM_NaCl/05_Cluster_GWAS_QC_by_threshold_trait.csv",
  stringsAsFactors = FALSE)

threshold_summary <- NULL
trait_summary <- NULL

for (threshold_clust in thresholds_clust) {
  label <- threshold_label(threshold_clust)
  f <- file.path("txt_GridLMM_NaCl",
                 paste0("all_sign_GWAS_with_clusters_", label, ".txt"))
  cl <- read.table(f, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  if (anyDuplicated(cl$rs)) stop("Duplicate SNP ID in ", label)
  if (!all(c("group", "CHROM", "POS") %in% names(cl)))
    stop("Missing clustering columns in ", label)

  cluster_rows <- split(seq_len(nrow(cl)), as.character(cl$group))
  cluster_sizes <- lengths(cluster_rows)
  cluster_spans <- vapply(cluster_rows, function(ix) {
    if (length(unique(cl$CHROM[ix])) != 1L)
      stop("Cluster spans chromosomes at ", label)
    diff(range(as.numeric(cl$POS[ix])))
  }, numeric(1))

  gqc_file <- file.path("data_BLUE_GridLMM_NaCl",
                   paste0("Cluster_genotypes_QC_", label, ".csv"))
  gqc <- read.csv(gqc_file, stringsAsFactors = FALSE)
  if (nrow(gqc) != length(cluster_rows))
    stop("Mismatch between clustering and genotype QC at ", label)
  mq <- subset(mapping_qc, abs(threshold - threshold_clust) < 1e-8)
  aq <- subset(association_qc, abs(threshold - threshold_clust) < 1e-8)
  if (nrow(mq) != length(unique(cl$CHROM)) ||
      nrow(aq) != length(traits))
    stop("Missing chromosome/trait mapping summaries at ", label)
  if (length(unique(aq$clusters_total)) != 1L ||
      unique(aq$clusters_total) != length(cluster_rows))
    stop("Cluster totals disagree at ", label)

  threshold_summary <- rbind(threshold_summary, data.frame(
    threshold = threshold_clust,
    n_SNPs = nrow(cl),
    n_clusters = length(cluster_rows),
    n_singletons = sum(cluster_sizes == 1L),
    median_SNPs_per_cluster = median(cluster_sizes),
    max_SNPs_per_cluster = max(cluster_sizes),
    median_span_bp = median(cluster_spans),
    p95_span_bp = unname(quantile(cluster_spans, 0.95)),
    n_clusters_over_100kb = sum(cluster_spans > 100000),
    n_clusters_with_missing_genotypes = sum(gqc$missing_RILs > 0),
    n_clusters_monomorphic = sum(!is.na(gqc$sd_genotype) & gqc$sd_genotype == 0),
    n_clusters_tested = sum(mq$n_tested),
    n_clusters_excluded_missing = sum(mq$excluded_missing),
    n_clusters_excluded_monomorphic = sum(mq$excluded_no_variation)
  ))

  trait_summary <- rbind(trait_summary, aq[, c("threshold", "trait",
    "clusters_total", "clusters_tested", "clusters_untested", "q_significant")])
}

threshold_summary$is_primary <- abs(threshold_summary$threshold - primary_threshold) < 1e-8
write.csv(threshold_summary,
          file.path(out_dir, "LD_threshold_clustering_QC.csv"), row.names = FALSE)
write.csv(trait_summary,
          file.path(out_dir, "LD_threshold_cluster_associations.csv"), row.names = FALSE)

## Figures are descriptive; significance counts are NOT a tuning objective.
pdf(file.path(out_dir, "LD_threshold_sensitivity.pdf"), width = 8, height = 8)
oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(threshold_summary$threshold, threshold_summary$n_clusters,
     type = "b", pch = 16, xlab = "Minimum LD R2", ylab = "QTL clusters",
     main = "Cluster resolution", las = 1, bty = "n")
abline(v = primary_threshold, lty = 2)
plot(threshold_summary$threshold, threshold_summary$n_singletons /
       threshold_summary$n_clusters,
     type = "b", pch = 16, xlab = "Minimum LD R2",
     ylab = "Fraction of singleton clusters", las = 1, bty = "n")
abline(v = primary_threshold, lty = 2)
plot(threshold_summary$threshold, threshold_summary$p95_span_bp / 1e3,
     type = "b", pch = 16, xlab = "Minimum LD R2",
     ylab = "95th percentile span (kb)", las = 1, bty = "n")
abline(v = primary_threshold, lty = 2)
plot(threshold_summary$threshold, threshold_summary$n_clusters_tested /
       threshold_summary$n_clusters,
     type = "b", pch = 16, xlab = "Minimum LD R2",
     ylab = "Fraction of clusters mapped", las = 1, bty = "n",
     ylim = c(0, 1))
abline(v = primary_threshold, lty = 2)
par(oldpar)
dev.off()

print(threshold_summary)
cat("\nPrimary threshold remains R2 =", primary_threshold,
    "; use neighboring cutoffs (0.4 and 0.6) for robustness.\n")
cat("Choose a different primary cutoff only for a documented LD-resolution or\n",
    "genotype-quality reason, not because it yields stronger downstream significance.\n")
