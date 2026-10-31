## 09c_Drift_GLS_genomewide.R
## Fast, approximate neutral drift-aware trajectory test, ALL SNPs.
## Retains ONE biologically justified generation-zero sample and all four evolved populations.
## This is NOT a classical estimate of s and does NOT replace WF simulation.
## Run BEFORE 09d. Uses outputs from scripts 08 and 09.

library(parallel)
source("Oct_pipeline_cleaner/09c_drift_functions.R")

pool_file <- "RData_GridLMM_NaCl/08_PoolSeq_prepared.RData"
selection_file <- "selection_results/09_Selection_SNP_results.csv"
out_dir <- "selection_results/drift_oneG0"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Confirm which G0 underwent the independently documented bottleneck.
# Assumes A6G140 was bottlenecked, consistent with earlier exploratory exclusion.
exclude_G0 <- "A6G140"
Ne_values <- c(1000)     # Main run. Optional sensitivity: c(500, 1000, 2000)
pool_individuals <- NA_real_  # Set TRUE number in each DNA pool if known.
mc_cores <- 6L
chunk_size <- 500L

p <- read_drift_inputs(pool_file, selection_file, exclude_G0 = exclude_G0)
sel <- p$sel
pop <- p$pop
k_alt <- p$k_alt
depth <- p$depth

cat("Excluded bottlenecked G0:", exclude_G0, "\n")
cat("Retained samples: ", paste(pop$sample_id, collapse = ", "), "\n", sep = "")
write.csv(pop, file.path(out_dir, "09c_retained_pool_manifest.csv"), row.names = FALSE)
print(pop[, c("sample_id", "generations", "replicate")])
cat("GWAS-ALT matched SNPs:", nrow(k_alt), "\n")
if (.Platform$OS.type == "windows") mc_cores <- 1L

for (Ne in Ne_values) {
  cat("\nAnalyzing drift covariance at Ne =", Ne, "\n")
  idx <- seq_len(nrow(k_alt))
  chunks <- split(idx, ceiling(idx / chunk_size))

  ans <- mclapply(chunks, function(ii) {
    t(vapply(ii, function(i) {
      z <- fit_drift_site(k_alt[i, ], depth[i, ], pop, Ne,
                          pool_individuals = pool_individuals)
      c(p0 = z$p0, slope = z$slope, SE = z$SE, z = z$z, p = z$p,
        n_pools = z$n_pools, n_G0 = z$n_G0,
        n_evolved_reps = z$n_evolved_reps,
        baseline_diff = z$baseline_diff, baseline_Z = z$baseline_Z)
    }, FUN.VALUE = c(p0 = 0, slope = 0, SE = 0, z = 0, p = 0,
                     n_pools = 0, n_G0 = 0, n_evolved_reps = 0,
                     baseline_diff = 0, baseline_Z = 0)))
  }, mc.cores = mc_cores, mc.preschedule = TRUE)
  ans <- do.call(rbind, ans)
  if (nrow(ans) != nrow(sel)) stop("SNP outputs changed order")

  result <- data.frame(
    rs = sel$rs,
    ID = sel$ID,
    CHROM = sel$CHROM,
    POS = sel$POS,
    orientation = sel$orientation,
    p0_hat = ans[, "p0"],
    slope_freq = ans[, "slope"],
    SE_freq = ans[, "SE"],
    Z_drift = ans[, "z"],
    p_drift = ans[, "p"],
    n_pools = as.integer(ans[, "n_pools"]),
    n_G0 = as.integer(ans[, "n_G0"]),
    n_evolved_reps = as.integer(ans[, "n_evolved_reps"]),
    baseline_freq_difference = ans[, "baseline_diff"],
    baseline_Z = ans[, "baseline_Z"],
    s_logit_original_twoG0 = sel$s,
    P_original_twoG0 = sel$P,
    p_quasi_original_twoG0 = sel$p_quasi,
    stringsAsFactors = FALSE
  )
  result$p_drift_BH <- p.adjust(result$p_drift, "BH")
  result$p_drift_BY <- p.adjust(result$p_drift, "BY")
  result$boundary_p0 <- is.finite(result$p0_hat) &
    (result$p0_hat < .02 | result$p0_hat > .98)
  result$Ne <- Ne
  result$excluded_G0 <- exclude_G0
  result$retained_G0 <- pop$sample_id[pop$generations == 0L]
  result$test <- "WF_covariance_GLS_approximation"

  fn <- file.path(out_dir, paste0("09c_SNP_drift_Ne", Ne, ".csv"))
  write.csv(result, fn, row.names = FALSE)
  cat("Wrote", fn, "\n")
  cat("Baseline_Z is NA with a single G0; this is expected.\n")
  print(data.frame(
    Ne = Ne, tested = sum(is.finite(result$p_drift)),
    q01_BH = sum(result$p_drift_BH < .01, na.rm = TRUE),
    q05_BH = sum(result$p_drift_BH < .05, na.rm = TRUE),
    q01_BY = sum(result$p_drift_BY < .01, na.rm = TRUE),
    boundary_p0 = sum(result$boundary_p0, na.rm = TRUE),
    n_G0 = unique(result$n_G0[is.finite(result$n_G0)])[1],
    min_raw_p = min(result$p_drift, na.rm = TRUE)))
}

capture.output(sessionInfo(), file = file.path(out_dir, "09c_sessionInfo.txt"))
cat("NOTE: Test calibration must be assessed with script 09d before\n",
    "interpreting BH/BY discoveries as evidence of selection.\n")
