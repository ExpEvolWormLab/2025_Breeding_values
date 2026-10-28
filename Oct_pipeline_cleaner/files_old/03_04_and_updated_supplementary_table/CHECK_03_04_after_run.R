## Run AFTER script 03 (and after 04 if mapping finished).
## A focused consistency check at your chosen R2 > 0.5 threshold.

label <- "R5"
clusters <- read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_", label, ".txt"),
                       header = TRUE, sep = "\t")
metadata <- read.table(paste0("data_BLUE_GridLMM_NaCl/data_clusters_BLUE_", label, ".txt"),
                       header = TRUE)
phase_table <- read.csv(paste0("data_BLUE_GridLMM_NaCl/SNP_phase_clusters_", label, ".csv"))
G <- as.matrix(read.table(paste0("data_BLUE_GridLMM_NaCl/genotypes_matrix_clusters_", label, ".txt"),
                          header = TRUE, row.names = 1, check.names = FALSE))

cat("Clusters:", nrow(metadata), "\n")
cat("Singletons:", sum(table(clusters$group) == 1), "\n")
cat("Cluster genotypes with missing RILs:", sum(rowSums(is.na(G)) > 0), "\n")
stopifnot(nrow(metadata) == 1089, sum(table(clusters$group) == 1) == 392,
          setequal(metadata$clusterID, rownames(G)),
          setequal(clusters$rs, phase_table$rs),
          !anyDuplicated(clusters$rs))

## Recompute a few averaged genotypes independently from saved phase labels.
prep <- new.env(parent = globalenv())
load("RData/Preparatory_GLMM_NaCl_only_withBeta.RData", envir = prep)
X <- prep$genotype_matrix
observed <- unique(as.vector(X)); observed <- observed[!is.na(observed)]
dosage_max <- if (2 %in% observed) 2 else 1

set.seed(24)
idx <- sample(seq_len(nrow(metadata)), min(20L, nrow(metadata)))
for (j in idx) {
  gp <- metadata$group[j]
  ids <- clusters$rs[clusters$group == gp]
  ph <- phase_table$phase[match(ids, phase_table$rs)]
  dosage <- t(X[, ids, drop = FALSE])
  flipped <- ph == -1
  if (any(flipped)) dosage[flipped, ] <- dosage_max - dosage[flipped, , drop = FALSE]
  expected <- colMeans(dosage, na.rm = TRUE)
  expected[is.nan(expected)] <- NA_real_
  actual <- G[metadata$clusterID[j], rownames(X)]
  stopifnot(isTRUE(all.equal(unname(actual), unname(expected),
                             tolerance = 1e-9, check.attributes = FALSE)))
}
cat("20 randomly sampled cluster genotypes match phase-aligned SNP averages.\n")

if (file.exists("GridLMM_clusters_NaCl/04_Cluster_mapping_QC.csv")) {
  qc <- read.csv("GridLMM_clusters_NaCl/04_Cluster_mapping_QC.csv")
  print(subset(qc, threshold == 0.5))
}
