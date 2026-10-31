## 09b_Diagnose_empty_selection.R
## Investigate why no SNP meets the selection threshold in scripts 09 / 10.
## This is diagnostic only: it does not redefine selected SNPs.

results_file <- "selection_results/09_Selection_SNP_results.csv"
out_dir <- "selection_results"
sel <- read.csv(results_file, stringsAsFactors = FALSE)

required <- c("P", "P_BH", "p_quasi", "p_quasi_BH",
              "P_within", "P_within_BH", "p_within_quasi",
              "p_within_quasi_BH", "dispersion", "n_pools")
missing <- setdiff(required, names(sel))
if (length(missing)) stop("Columns missing: ", paste(missing, collapse = ", "))

summarize_test <- function(name) {
  x <- sel[[name]]
  xf <- x[is.finite(x)]
  data.frame(test = name,
             total_SNPs = length(x),
             finite_p = length(xf),
             min_p = if (length(xf)) min(xf) else NA_real_,
             median_p = if (length(xf)) median(xf) else NA_real_,
             n_below_0.05 = sum(xf < 0.05),
             n_below_0.01 = sum(xf < 0.01),
             n_below_1e_5 = sum(xf < 1e-5))
}

cols <- c("P", "P_BH", "p_quasi", "p_quasi_BH",
          "P_within", "P_within_BH", "p_within_quasi", "p_within_quasi_BH")
summary_tests <- do.call(rbind, lapply(cols, summarize_test))
print(summary_tests, row.names = FALSE, digits = 5)
write.csv(summary_tests, file.path(out_dir, "09b_selection_test_summary.csv"),
          row.names = FALSE)

cat("\nNumber of pools contributing to each SNP:\n")
print(table(sel$n_pools, useNA = "ifany"))

cat("\nDispersion for tested SNPs:\n")
print(summary(sel$dispersion[is.finite(sel$dispersion)]))
cat("Fraction with dispersion > 5:",
    mean(sel$dispersion[is.finite(sel$dispersion)] > 5), "\n")

cat("\nTen SNPs with smallest quasi-binomial p-value:\n")
ord <- order(sel$p_quasi, na.last = TRUE)
if (length(ord)) {
  keep_cols <- intersect(c("ID", "CHROM", "POS", "n_pools",
                            "s", "p_quasi", "p_quasi_BH",
                            "P", "P_BH", "dispersion"), names(sel))
  print(utils::head(sel[ord, keep_cols, drop = FALSE], 10), row.names = FALSE)
}

if (file.exists("selection_results/09_allele_orientation_QC.csv")) {
  meta <- read.csv("selection_results/09_allele_orientation_QC.csv",
                   stringsAsFactors = FALSE)
  cat("\nPool-seq/GWAS allele orientation:\n")
  print(table(meta$orientation, useNA = "ifany"))
}

cat("\nDiagnostic complete. No thresholds or scientific results were changed.\n")
