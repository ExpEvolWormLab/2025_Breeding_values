## 02b_Summarise_clusters_CHR.R
## Summarize clusters produced by 02_Clustering_from_LD.R.
## Uses exactly the same R2 cutoffs and legacy filenames.

vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
threshold_label <- function(x) paste0("R", round(10 * x))

clusters_dir <- "txt_GridLMM_NaCl"
plot_dir <- "plots_BLUEs_GridLMM_NaCl"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

summary_clust <- NULL
summary_traits <- NULL

for (threshold_clust in thresholds_clust) {
  cat("Summarizing R2 >", threshold_clust, "\n")
  f <- file.path(clusters_dir, paste0("all_sign_GWAS_with_clusters_",
                                      threshold_label(threshold_clust), ".txt"))
  all_cluster_merged <- read.table(f, header = TRUE, sep = "\t")
  sizes <- table(all_cluster_merged$group)

  summary_clust <- rbind(summary_clust,
    data.frame(r2 = threshold_clust,
               nb_clust = length(sizes),
               nb_singletons = sum(sizes == 1),
               nb_SNPs = nrow(all_cluster_merged)))

  for (i in seq_along(vect_P_traits)) {
    selected <- all_cluster_merged[
      all_cluster_merged[[paste0("is_sign_T", i)]] %in% TRUE, , drop = FALSE]
    trait_sizes <- table(selected$group)
    summary_traits <- rbind(summary_traits,
      data.frame(r2 = threshold_clust, trait = vect_P_traits[i],
                 nb_clust = length(trait_sizes),
                 nb_singletons = sum(trait_sizes == 1),
                 nb_SNPs = nrow(selected)))
  }
}

write.csv(summary_clust,
          file.path(clusters_dir, "Summary_clusters_by_R2.csv"), row.names = FALSE)
write.csv(summary_traits,
          file.path(clusters_dir, "Summary_clusters_by_trait_R2.csv"), row.names = FALSE)

pdf(file.path(plot_dir, "Summary_clusters_R2.pdf"), width = 6, height = 4.5)
plot(summary_clust$r2, summary_clust$nb_clust, type = "b", pch = 16,
     xlab = expression("LD threshold " * R^2),
     ylab = "Number of QTL groups", bty = "n", las = 1,
     ylim = range(c(summary_clust$nb_clust, summary_clust$nb_singletons)))
lines(summary_clust$r2, summary_clust$nb_singletons,
      type = "b", pch = 16, col = "firebrick3")
legend("topleft", legend = c("All groups", "One-SNP groups"),
       col = c("black", "firebrick3"), lty = 1, pch = 16, bty = "n")
dev.off()

print(summary_clust)
cat("Note: original R10 filename corresponds to R2 > 0.99 (not R2 = 1).\n")
