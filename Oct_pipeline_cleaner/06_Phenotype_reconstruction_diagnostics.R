## 06_Phenotype_reconstruction_diagnostics.R
## Replaces exploratory 11_TEST_PREDICTIONS_all_thresh.R.
## Reconstruct adjusted RIL phenotypes using significant cluster GWAS effects.
## IMPORTANT: discovery, effect estimation, and reconstruction use the same RILs.
## These are IN-SAMPLE descriptive diagnostics, NOT predictive validation.
## Run AFTER scripts 00, 03, 04 and 05 at the requested thresholds.

traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
source_traits <- traits[1:7]
thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
threshold_label <- function(x) paste0("R", round(10 * x))

out_dir <- "phenotype_reconstruction"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## The BLUEs come from the current script 00, not the historical Gmatrices.RData.
prep <- new.env(parent = globalenv())
load("RData/Preparatory_GLMM_NaCl_only_withBeta.RData", envir = prep)
phen <- prep$phen_NaCl
phen$pop_label <- as.character(phen$pop_label)
if (anyDuplicated(phen$pop_label)) stop("Duplicate RIL in BLUE table")

## The composite response is constructed from the seven BLUEs and the
## selection gradients already estimated in script 00. It is not an eighth
## independent phenotype; report its results separately from the seven traits.
Z <- sapply(source_traits, function(trait) {
  (phen[[trait]] - prep$BLUE_mean[trait]) / prep$BLUE_sd[trait]
})
phen$beta <- drop(Z %*% prep$beta_NaCl[source_traits])
rm(prep, Z)

reconstruction_results <- NULL

for (threshold_clust in thresholds_clust) {
  label <- threshold_label(threshold_clust)
  cat("\nIn-sample phenotype reconstruction:", label, "\n")

  result_file <- file.path("RData_GridLMM_NaCl",
    paste0("GWAS_cluster_results_BLUE_concatenated_", label,
           "_with_residuals.RData"))
  results_env <- new.env(parent = globalenv())
  load(result_file, envir = results_env)
  a <- results_env$all_GWAS_clusters

  genotype_file <- file.path("data_BLUE_GridLMM_NaCl",
                      paste0("genotypes_matrix_clusters_", label, ".txt"))
  G <- as.matrix(read.table(genotype_file, header = TRUE, row.names = 1,
                             check.names = FALSE))
  storage.mode(G) <- "numeric"

  if (anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)) ||
      anyDuplicated(a$rs) || !setequal(rownames(G), a$rs))
    stop("Cluster identifiers do not match at ", label)
  G <- G[as.character(a$rs), , drop = FALSE]

  common_rils <- intersect(colnames(G), phen$pop_label)
  if (length(common_rils) < 15L) stop("Too few matched RILs at ", label)
  G <- G[, common_rils, drop = FALSE]
  observed <- phen[match(common_rils, phen$pop_label), , drop = FALSE]

  for (trait in traits) {
    suffix <- paste0(trait, "_NaCl")
    sig <- a[[paste0("is_sign_", suffix)]] %in% TRUE
    beta <- a[[paste0("beta_", suffix)]]
    selected <- which(sig & is.finite(beta))

    ## Center cluster dosages: the arbitrary reference allele then affects
    ## the intercept, not the reconstruction correlation.
    if (length(selected)) {
      X <- G[selected, , drop = FALSE]
      X <- sweep(X, 1, rowMeans(X, na.rm = TRUE), "-")
      predicted <- colSums(sweep(X, 1, beta[selected], "*"))
    } else {
      predicted <- rep(NA_real_, length(common_rils))
    }

    y <- observed[[trait]]
    ok <- is.finite(y) & is.finite(predicted)
    if (sum(ok) >= 15L && sd(y[ok]) > 0 && sd(predicted[ok]) > 0) {
      correlation <- cor(y[ok], predicted[ok])
      r2 <- correlation^2
    } else {
      correlation <- NA_real_
      r2 <- NA_real_
    }

    reconstruction_results <- rbind(reconstruction_results,
      data.frame(threshold = threshold_clust, trait = trait,
                 n_RIL = sum(ok), n_significant_clusters = length(selected),
                 correlation = correlation, R2_in_sample = r2,
                 kind = if (trait == "beta") "derived_selection_index" else "original_trait"))
  }
}

write.csv(reconstruction_results,
          file.path(out_dir, "06_In_sample_reconstruction_by_threshold.csv"),
          row.names = FALSE)

## Separate primary measured traits from the constructed selection index.
real_traits <- subset(reconstruction_results, kind == "original_trait")
summary_results <- aggregate(R2_in_sample ~ threshold, data = real_traits,
                             FUN = function(x) median(x, na.rm = TRUE))
write.csv(summary_results,
          file.path(out_dir, "06_Median_in_sample_R2_original_traits.csv"),
          row.names = FALSE)

pdf(file.path(out_dir, "06_In_sample_reconstruction_sensitivity.pdf"),
    width = 8, height = 4.2)
oldpar <- par(mfrow = c(1, 2), mar = c(4.4, 4.5, 2, 1))
plot(summary_results$threshold, summary_results$R2_in_sample,
     type = "b", pch = 16, ylim = c(0, 1), bty = "n", las = 1,
     xlab = "Minimum LD R2", ylab = "Median in-sample R2 (7 traits)",
     main = "Original phenotypes")
abline(v = 0.5, lty = 2)
idx <- subset(reconstruction_results, trait == "beta")
plot(idx$threshold, idx$R2_in_sample,
     type = "b", pch = 16, ylim = c(0, 1), bty = "n", las = 1,
     xlab = "Minimum LD R2", ylab = "In-sample R2",
     main = "Derived selection phenotype")
abline(v = 0.5, lty = 2)
par(oldpar)
dev.off()

cat("\nThese R2 values are in-sample. They cannot establish the optimal\n",
    "LD threshold or the out-of-sample predictive accuracy of your GWAS.\n")
