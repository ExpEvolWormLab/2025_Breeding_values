## CHECK_07a_09_after_run.R
## Run the relevant section only after creating its files.
if (file.exists("RData_G_NaCl/07a_multivariate_G_final.RData")) {
  g_env <- new.env(parent = globalenv())
  load("RData_G_NaCl/07a_multivariate_G_final.RData", envir = g_env)
  cat("Eigenvalues of among-RIL covariance G:\n")
  print(eigen(g_env$G_RIL, symmetric = TRUE)$values)
  cat("Genetic correlations:\n")
  print(round(g_env$R_traits, 2))
  cat("MCMC diagnostics: maximum Rhat / minimum ESS:\n")
  print(c(max(g_env$diagnostics$Rhat, na.rm = TRUE),
          min(g_env$diagnostics$ESS, na.rm = TRUE)))
  stopifnot(min(eigen(g_env$G_RIL, symmetric = TRUE)$values) > 0,
            all(abs(diag(g_env$R_traits) - 1) < 1e-8))
}
if (file.exists("selection_results/09_Selection_SNP_results.csv")) {
  sel <- read.csv("selection_results/09_Selection_SNP_results.csv")
  stopifnot(!anyDuplicated(sel$ID),
            all(sel$orientation %in% c("minor_is_alt", "minor_is_ref")),
            all(sel$P[is.finite(sel$P)] >= 0 &
                sel$P[is.finite(sel$P)] <= 1),
            all(sel$p_quasi[is.finite(sel$p_quasi)] >= 0 &
                sel$p_quasi[is.finite(sel$p_quasi)] <= 1))
  cat("SNPs with selection trajectory fitted:", sum(is.finite(sel$s)), "\n")
  cat("Median overdispersion:", median(sel$dispersion, na.rm = TRUE), "\n")
  cat("Quasi-binomial q<0.01:",
      sum(sel$p_quasi_BH < 0.01, na.rm = TRUE), "\n")
}
