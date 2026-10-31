## 12_Compare_old_new_beta_associations.R
## Diagnose changes in the composite-fitness GWAS/evolution relationship.
## Run from the project root AFTER scripts 09f, 10 and 11.
## Does not modify any upstream results or declare SNPs selected.

cluster_file <- "selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv"
snp_file <- "selection_results/drift_oneG0_nullfit/R5/11_SNP_phase_aligned_09f.csv"
legacy_file <- "RData_GridLMM_NaCl/Estimates_s.txt"
out_dir <- "selection_results/drift_oneG0_nullfit/R5/discrepancy_audit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qtl <- read.csv(cluster_file, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("group", "beta_beta_NaCl", "p_value_REML.cor_beta_NaCl",
            "slope_median_common", "slope_mean_common", "slope_median_all_sensitivity",
            "n_slope_common", "n_slope_all")
if (!all(needed %in% names(qtl)))
  stop("Missing script-11 cluster fields: ", paste(setdiff(needed, names(qtl)), collapse = ", "))
if (anyDuplicated(as.character(qtl$group))) stop("Duplicate GWAS cluster IDs")

## Trait-specific GWAS significance is a prespecified *GWAS-only* filter.
## Note: selection on |GWAS beta| can also change apparent association strength.
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")

summarise_fit <- function(d, x_col, y_col, cohort, response) {
  ok <- is.finite(d[[x_col]]) & is.finite(d[[y_col]])
  d <- d[ok, , drop = FALSE]
  ans <- data.frame(cohort = cohort, response = response, n = nrow(d),
                    pearson_r = NA_real_, spearman_r = NA_real_,
                    slope = NA_real_, p_descriptive = NA_real_, R2 = NA_real_,
                    sign_agreement = NA_real_, stringsAsFactors = FALSE)
  if (nrow(d) < 4L || sd(d[[x_col]]) == 0 || sd(d[[y_col]]) == 0) return(ans)
  xx <- d[[x_col]]; yy <- d[[y_col]]
  mod <- lm(yy ~ xx)
  ans$pearson_r <- cor(xx, yy)
  ans$spearman_r <- cor(xx, yy, method = "spearman")
  ans$slope <- unname(coef(mod)[2])
  ans$p_descriptive <- coef(summary(mod))[2, "Pr(>|t|)"]
  ans$R2 <- summary(mod)$r.squared
  nonzero <- is.finite(xx) & is.finite(yy) & xx != 0 & yy != 0
  ans$sign_agreement <- if (any(nonzero)) mean(sign(xx[nonzero]) == sign(yy[nonzero])) else NA_real_
  ans
}

compare <- list(); k <- 0L
append_fit <- function(d, x, y, cohort, response) {
  k <<- k + 1L
  compare[[k]] <<- summarise_fit(d, x, y, cohort, response)
}

## Same GWAS beta effect, varied ascertainment and cluster response.
beta_q <- "p_value_REML.cor_beta_NaCl"
beta_col <- "beta_beta_NaCl"
cohorts <- list(
  all_clusters = rep(TRUE, nrow(qtl)),
  composite_GWAS_q05 = is.finite(qtl[[beta_q]]) & qtl[[beta_q]] < .05,
  composite_GWAS_q01 = is.finite(qtl[[beta_q]]) & qtl[[beta_q]] < .01,
  multiSNP_common = is.finite(qtl$n_slope_common) & qtl$n_slope_common >= 2
)
responses <- c(median_common = "slope_median_common",
               mean_common = "slope_mean_common",
               median_all = "slope_median_all_sensitivity")
for (cohort in names(cohorts)) {
  dd <- qtl[cohorts[[cohort]], , drop = FALSE]
  for (r in names(responses))
    append_fit(dd, beta_col, responses[[r]], cohort, r)
}

## Direct comparison between component traits and the composite on the
## IDENTICAL composite-GWAS-significant cluster set.
for (tr in traits) {
  x <- paste0("beta_", tr, "_NaCl")
  if (!x %in% names(qtl)) next
  dd <- qtl[cohorts$composite_GWAS_q05, , drop = FALSE]
  append_fit(dd, x, "slope_median_common", "composite_GWAS_q05", paste0("trait_", tr))
}

## Optional: compare original 16_Compute_s.R response on EXACT same SNP keys
## and GWAS cluster definitions, keeping the latest phase and GWAS betas fixed.
## This isolates response estimation from cluster construction/phenotype changes.
if (file.exists(legacy_file) && file.exists(snp_file)) {
  snps <- read.csv(snp_file, stringsAsFactors = FALSE, check.names = FALSE)
  old <- read.table(legacy_file, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c("rs", "s") %in% names(old)) &&
      all(c("rs", "group", "phase", "tested_common", "slope_freq_cluster") %in% names(snps))) {
    old$ID <- sub("\\.", "_", as.character(old$rs))
    if (!anyDuplicated(old$ID)) {
      mi <- match(as.character(snps$rs), old$ID)
      snps$legacy_s_ALT <- suppressWarnings(as.numeric(old$s[mi]))
      ## Original cbind() count order appears to use GWAS ALT as the first allele.
      ## Confirm polarization with individual loci before causal interpretation.
      snps$legacy_s_cluster <- snps$phase * snps$legacy_s_ALT
      snps$legacy_s_cluster_common <- ifelse(snps$tested_common,
                                               snps$legacy_s_cluster, NA_real_)
      agg <- lapply(split(snps, as.character(snps$group)), function(d) {
        z <- d$legacy_s_cluster_common
        z <- z[is.finite(z)]
        data.frame(group = as.character(d$group[1]),
                   n_legacy = length(z),
                   legacy_s_median_common = if (length(z)) median(z) else NA_real_,
                   legacy_s_mean_common = if (length(z)) mean(z) else NA_real_)
      })
      legacy_qtl <- do.call(rbind, agg)
      rownames(legacy_qtl) <- NULL
      ix <- match(as.character(qtl$group), legacy_qtl$group)
      qtl$n_legacy <- legacy_qtl$n_legacy[ix]
      qtl$legacy_s_median_common <- legacy_qtl$legacy_s_median_common[ix]
      qtl$legacy_s_mean_common <- legacy_qtl$legacy_s_mean_common[ix]
      for (cohort in c("all_clusters", "composite_GWAS_q05")) {
        dd <- qtl[cohorts[[cohort]], , drop = FALSE]
        append_fit(dd, beta_col, "legacy_s_median_common", cohort, "legacy_logit_median_common")
        append_fit(dd, beta_col, "legacy_s_mean_common", cohort, "legacy_logit_mean_common")
      }
      write.csv(qtl[, c("group", "beta_beta_NaCl", "p_value_REML.cor_beta_NaCl",
                        "slope_median_common", "slope_mean_common",
                        "legacy_s_median_common", "legacy_s_mean_common",
                        "n_slope_common", "n_legacy")],
                file.path(out_dir, "12_Cluster_old_vs_new_response.csv"), row.names = FALSE)
      valid <- is.finite(qtl$slope_median_common) & is.finite(qtl$legacy_s_median_common)
      cat("Legacy and new cluster responses on common-frequency SNPs:",
          sum(valid), "clusters; Pearson r =",
          if (sum(valid) > 2L) cor(qtl$slope_median_common[valid],
                                    qtl$legacy_s_median_common[valid]) else NA_real_, "\n")
    } else warning("Original selection results have duplicate locus IDs; skipped legacy comparison")
  } else warning("Original selection result fields differ; skipped legacy comparison")
} else cat("Optional original Estimates_s.txt not found; running new-pipeline cohort comparison only.\n")

ans <- do.call(rbind, compare)
write.csv(ans, file.path(out_dir, "12_Comparison_associations.csv"), row.names = FALSE)
write.csv(data.frame(group = qtl$group,
                     beta_GWAS_q = qtl[[beta_q]],
                     selected_beta_GWAS_q05 = cohorts$composite_GWAS_q05,
                     n_slope_common = qtl$n_slope_common,
                     n_slope_all = qtl$n_slope_all),
          file.path(out_dir, "12_GWAS_composite_cluster_counts.csv"), row.names = FALSE)
print(ans, row.names = FALSE, digits = 4)
cat("\nSaved discrepancy audit in:", out_dir, "\n")
cat("Interpret regressions as descriptive, not independent-locus inferential tests.\n")
