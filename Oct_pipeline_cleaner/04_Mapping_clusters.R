## 04_Mapping_clusters.R  (formerly 09_Mapping_clusters.R)
## Re-run GridLMM GWAS on phase-oriented, averaged SNP genotypes per cluster.
## Keep SNP-level GWAS phenotype/model settings and LOCO kinship from script 01.
## Run after 03_Create_matrix_for_cluster_mapping_CHR.R.

library(AGHmatrix)
library(GridLMM)

vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
chromosomes <- c("I", "II", "III", "IV", "V", "X")
#thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
thresholds_clust=0.9
threshold_label <- function(x) paste0("R", round(10 * x))
mc_cores <- 10L
overwrite_results <- TRUE

## Like script 01, discard cluster markers with any missing RIL dosage.
## Change to "mean" only if you want an explicitly documented imputation.
missing_policy <- "exclude"  # either "exclude" or "mean"

input_file <- "RData/Preparatory_GLMM_NaCl_only_withBeta.RData"
cluster_dir <- "data_BLUE_GridLMM_NaCl"
output_dir <- "GridLMM_clusters_NaCl"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prep <- new.env(parent = globalenv())
load(input_file, envir = prep)
data_for_GWAS_NaCl_genotyped_NaCl <- prep$data_for_GWAS_NaCl_genotyped_NaCl
genotype_pruned_matrix <- prep$genotype_pruned_matrix
ril_ids <- rownames(prep$genotype_matrix)
rm(prep)

if (anyDuplicated(ril_ids) || anyDuplicated(rownames(genotype_pruned_matrix)) ||
    !setequal(ril_ids, rownames(genotype_pruned_matrix)))
  stop("Pruned and full genotype matrices have different RIL identifiers")
genotype_pruned_matrix <- genotype_pruned_matrix[ril_ids, , drop = FALSE]

## Calculate each leave-one-chromosome-out kinship matrix once across
## both traits AND LD thresholds, rather than many times in nested loops.
chr_pruned <- sub("_.*$", "", colnames(genotype_pruned_matrix))
GRM <- setNames(lapply(chromosomes, function(chr) {
  X <- genotype_pruned_matrix[, chr_pruned != chr, drop = FALSE]
  K <- AGHmatrix::Gmatrix(X, method = "VanRaden")
  K[ril_ids, ril_ids, drop = FALSE]
}), chromosomes)

mapping_qc <- NULL
for (threshold_clust in thresholds_clust) {
  label <- threshold_label(threshold_clust)
  cat("\nCluster GWAS for LD R2 >", threshold_clust, "\n")

  geno_file <- file.path(cluster_dir,
                          paste0("genotypes_matrix_clusters_", label, ".txt"))
  meta_file <- file.path(cluster_dir,
                          paste0("data_clusters_BLUE_", label, ".txt"))

  ## 03 writes rownames as cluster IDs and columns as RIL IDs.
  genotyped_matrix <- as.matrix(read.table(geno_file, header = TRUE,
                                           row.names = 1, check.names = FALSE))
  storage.mode(genotyped_matrix) <- "numeric"
  df_clusters <- read.table(meta_file, header = TRUE, stringsAsFactors = FALSE)
  if (anyDuplicated(rownames(genotyped_matrix)) ||
      !setequal(rownames(genotyped_matrix), df_clusters$clusterID) ||
      !setequal(colnames(genotyped_matrix), ril_ids))
    stop("Cluster IDs or RIL IDs do not match at ", label)
  genotyped_matrix <- genotyped_matrix[df_clusters$clusterID, ril_ids, drop = FALSE]

  threshold_out <- file.path(output_dir, label)
  dir.create(threshold_out, recursive = TRUE, showWarnings = FALSE)

  for (chrID in chromosomes) {
    cat("  Chromosome", chrID, "\n")
    ix <- which(df_clusters$CHROM == chrID)
    if (!length(ix)) next
    X_CHR <- t(genotyped_matrix[ix, , drop = FALSE])

    ## Apply the missing-data rule to the cluster genotype, never to a
    ## fabricated genotype value of 3 (the old script's sentinel).
    n_missing <- colSums(!is.finite(X_CHR))
    if (missing_policy == "mean") {
      for (j in which(n_missing > 0)) {
        valid <- is.finite(X_CHR[, j])
        if (any(valid)) X_CHR[!valid, j] <- mean(X_CHR[valid, j])
      }
    } else if (missing_policy != "exclude") {
      stop("Choose missing_policy = 'exclude' or 'mean'")
    }

    complete <- colSums(!is.finite(X_CHR)) == 0
    nonconstant <- vapply(seq_len(ncol(X_CHR)), function(j) {
      x <- X_CHR[, j]
      if (!all(is.finite(x))) return(FALSE)
      length(unique(x)) > 1L
    }, logical(1))
    keep <- complete & nonconstant

    mapping_qc <- rbind(mapping_qc,
      data.frame(threshold = threshold_clust, chromosome = chrID,
                 n_clusters = ncol(X_CHR),
                 excluded_missing = sum(!complete),
                 excluded_no_variation = sum(complete & !nonconstant),
                 n_tested = sum(keep)))

    if (!any(keep)) {
      warning("No complete polymorphic clusters on ", chrID, " at ", label)
      next
    }
    X_CHR <- X_CHR[, keep, drop = FALSE]
    stopifnot(identical(rownames(X_CHR), rownames(GRM[[chrID]])))

    for (traitID in vect_P_traits) {
      cat("    Trait", traitID, "\n")
      data_for_GWAS_NaCl_genotyped_NaCl$y <- as.numeric(
        scale(data_for_GWAS_NaCl_genotyped_NaCl[[traitID]]))

      out_file <- file.path(threshold_out,
        paste0("gwas_results_", label, "_", traitID, "_CHR_", chrID, "_NaCl.txt"))
      if (!overwrite_results && file.exists(out_file)) next

      fit <- GridLMM::GridLMM_GWAS(
        formula = y ~ 1 + logD + rel_humidity + temperature +
          (1 | date_str) + (1 | pop_label) + (1 | session_id),
        test_formula = ~1, reduced_formula = ~0,
        data = data_for_GWAS_NaCl_genotyped_NaCl,
        X = X_CHR, X_ID = "pop_label",
        relmat = list(pop_label = GRM[[chrID]]),
        mc.cores = mc_cores, method = "REML"
      )
      write.table(fit$results, out_file, quote = FALSE, row.names = FALSE)
    }
  }
}

write.csv(mapping_qc, file.path(output_dir, "04_Cluster_mapping_QC.csv"),
          row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "04_sessionInfo.txt"))
cat("Finished cluster GWAS; outputs in", output_dir, "\n")
