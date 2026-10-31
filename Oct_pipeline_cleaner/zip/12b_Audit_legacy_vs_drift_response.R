## 12b_Audit_legacy_vs_drift_response.R
## Isolate response-scale and estimator discrepancies on IDENTICAL SNPs and
## clusters, without rerunning GWAS, Pool-seq models, or simulations.
## Readable alongside 12_Comparison_associations.csv. No package required.

snp_file <- "selection_results/drift_oneG0_nullfit/R5/11_SNP_phase_aligned_09f.csv"
qtl_file <- "selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv"
legacy_file <- "RData_GridLMM_NaCl/Estimates_s.txt"
out_dir <- "selection_results/drift_oneG0_nullfit/R5/discrepancy_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

snp <- read.csv(snp_file, stringsAsFactors = FALSE, check.names = FALSE)
qtl <- read.csv(qtl_file, stringsAsFactors = FALSE, check.names = FALSE)
old <- read.table(legacy_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
required_snp <- c("rs", "group", "phase", "tested_common",
                  "p0_cluster", "slope_freq_cluster")
required_qtl <- c("group", "beta_beta_NaCl", "p_value_REML.cor_beta_NaCl")
if (!all(required_snp %in% names(snp)) ||
    !all(required_qtl %in% names(qtl)) ||
    !all(c("rs", "s") %in% names(old)))
  stop("Missing expected columns in 11 SNP/cluster outputs or legacy Estimates_s.txt")
if (anyDuplicated(snp$rs) || anyDuplicated(qtl$group))
  stop("Duplicate SNP or group keys; cannot match safely")

old$ID <- sub("\\.", "_", as.character(old$rs))
if (anyDuplicated(old$ID)) stop("Duplicate legacy SNP positions")
ix <- match(as.character(snp$rs), old$ID)
snp$legacy_s_input <- suppressWarnings(as.numeric(old$s[ix]))
## This assumes old 16_Compute_s.R success counts correspond to GWAS ALT,
## as intended by the old polarization code. VERIFY at individual loci.
snp$old_logit_oriented <- snp$phase * snp$legacy_s_input
snp$new_freq_oriented <- snp$slope_freq_cluster
p0 <- as.numeric(snp$p0_cluster)
snp$old_freq_local <- snp$old_logit_oriented * p0 * (1 - p0)
## A logistic-trajectory conversion over 50 generations (not a fitted slope).
## It is a scale check, NOT a reanalysis of the old Pool-seq likelihood.
snp$old_freq_logistic50 <- (plogis(qlogis(p0) + 50 * snp$old_logit_oriented) - p0) / 50
snp$new_logit_local <- snp$new_freq_oriented / (p0 * (1 - p0))
snp$new_logit_local[!is.finite(snp$new_logit_local)] <- NA_real_
keep <- !is.na(snp$tested_common) & snp$tested_common &
  is.finite(snp$new_freq_oriented) & is.finite(snp$old_logit_oriented) &
  is.finite(p0) & p0 > 0 & p0 < 1
snp <- snp[keep, , drop = FALSE]
if (!nrow(snp)) stop("No common-frequency SNPs with both old and new estimates")

safe_median <- function(x) { x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
safe_mean <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
cols <- c("old_logit_oriented", "old_freq_local", "old_freq_logistic50",
          "new_freq_oriented", "new_logit_local")
cluster_split <- split(snp, as.character(snp$group))
cl <- do.call(rbind, lapply(cluster_split, function(d) {
  x <- data.frame(group = as.character(d$group[1]), n_matched = nrow(d),
                  mean_p0 = mean(d$p0_cluster, na.rm = TRUE),
                  stringsAsFactors = FALSE)
  for (v in cols) {
    x[[paste0(v, "_median")]] <- safe_median(d[[v]])
    x[[paste0(v, "_mean")]] <- safe_mean(d[[v]])
  }
  x
}))
rownames(cl) <- NULL
cl <- merge(cl, qtl[, required_qtl], by = "group", all.x = TRUE, sort = FALSE)
cl$composite_GWAS_q05 <- is.finite(cl$p_value_REML.cor_beta_NaCl) &
  cl$p_value_REML.cor_beta_NaCl < .05

## Direct correlations compare estimators on the exact SAME SNPs or clusters.
cor_row <- function(a, b, label, cohort) {
  good <- is.finite(a) & is.finite(b)
  data.frame(level = label, cohort = cohort, n = sum(good),
    Pearson_r = if (sum(good) > 3) cor(a[good], b[good]) else NA_real_,
    Spearman_r = if (sum(good) > 3) cor(a[good], b[good], method = "spearman") else NA_real_,
    sign_agreement = if (sum(good) > 0) mean(sign(a[good]) == sign(b[good])) else NA_real_)
}
compare <- list(
  cor_row(snp$old_logit_oriented, snp$new_freq_oriented,
          "SNP_old_logit_vs_new_freq", "common"),
  cor_row(snp$old_freq_local, snp$new_freq_oriented,
          "SNP_old_local_freq_vs_new_freq", "common"),
  cor_row(snp$old_freq_logistic50, snp$new_freq_oriented,
          "SNP_old_logistic50_freq_vs_new_freq", "common")
)
for (cohort in c("all", "composite_GWAS_q05")) {
  d <- if (cohort == "all") cl else cl[cl$composite_GWAS_q05, , drop = FALSE]
  for (suffix in c("median", "mean")) {
    compare[[length(compare) + 1L]] <-
      cor_row(d[[paste0("old_logit_oriented_", suffix)]],
              d[[paste0("new_freq_oriented_", suffix)]],
              paste0("cluster_old_logit_vs_new_freq_", suffix), cohort)
    compare[[length(compare) + 1L]] <-
      cor_row(d[[paste0("old_freq_local_", suffix)]],
              d[[paste0("new_freq_oriented_", suffix)]],
              paste0("cluster_old_local_freq_vs_new_freq_", suffix), cohort)
    compare[[length(compare) + 1L]] <-
      cor_row(d[[paste0("old_freq_logistic50_", suffix)]],
              d[[paste0("new_freq_oriented_", suffix)]],
              paste0("cluster_old_logistic50_freq_vs_new_freq_", suffix), cohort)
  }
}
comparisons <- do.call(rbind, compare)

fit_row <- function(d, xvar, response, cohort) {
  ok <- is.finite(d$beta_beta_NaCl) & is.finite(d[[xvar]])
  d <- d[ok, , drop = FALSE]
  out <- data.frame(cohort = cohort, response = response, n = nrow(d),
                    pearson_r = NA_real_, spearman_r = NA_real_,
                    slope = NA_real_, p_descriptive = NA_real_, R2 = NA_real_)
  if (nrow(d) < 5L || sd(d$beta_beta_NaCl) == 0 || sd(d[[xvar]]) == 0) return(out)
  fit <- lm(d[[xvar]] ~ d$beta_beta_NaCl)
  out$pearson_r <- cor(d$beta_beta_NaCl, d[[xvar]])
  out$spearman_r <- cor(d$beta_beta_NaCl, d[[xvar]], method = "spearman")
  out$slope <- unname(coef(fit)[2]); out$p_descriptive <- coef(summary(fit))[2, 4]
  out$R2 <- summary(fit)$r.squared
  out
}
fits <- list()
for (cohort in c("all", "composite_GWAS_q05")) {
  d <- if (cohort == "all") cl else cl[cl$composite_GWAS_q05, , drop = FALSE]
  for (suffix in c("median", "mean")) for (v in cols) {
    var <- paste0(v, "_", suffix)
    fits[[length(fits) + 1L]] <- fit_row(d, var, var, cohort)
  }
}
fits <- do.call(rbind, fits)

write.csv(comparisons, file.path(out_dir, "12b_Old_new_response_correlations.csv"), row.names = FALSE)
write.csv(fits, file.path(out_dir, "12b_Matched_beta_associations.csv"), row.names = FALSE)
write.csv(cl, file.path(out_dir, "12b_Matched_cluster_slopes.csv"), row.names = FALSE)
cat("\nDirect old/new estimator comparison (SNPs and clusters):\n")
print(comparisons, row.names = FALSE, digits = 4)
cat("\nComposite GWAS effects vs each response on the matched SNP set:\n")
print(fits, row.names = FALSE, digits = 4)
cat("\nSaved files in:", out_dir, "\n")
cat("NOTE: Logistic conversions check SCALE, not the effects of different G0 samples,\n")
cat("      random effects, or covariance estimation. Validate legacy allele orientation.\n")
