## Drift-aware, neutral-null-frequency GLS screening test for Pool-seq.
## Compares to prior one-G0 approach without overwriting it.
## Strongest validity for common ancestral alleles; calibration follows in 09g.

library(parallel)
source("09f_drift_functions.R")

pool_file <- "RData_GridLMM_NaCl/08_PoolSeq_prepared.RData"
selection_file <- "selection_results/09_Selection_SNP_results.csv"
out_dir <- "selection_results/drift_oneG0_nullfit"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## EDIT only if a different G0 sample is confirmed bottlenecked.
exclude_G0 <- "A6G140"
Ne_values <- c(1000)
pool_individuals <- NA_real_
mc_cores <- 6L
chunk_size <- 500L
## Prespecified common-allele sensitivity set; boundary estimates unreliable.
common_range <- c(0.10, 0.90)

p <- read_drift_inputs(pool_file, selection_file, exclude_G0 = exclude_G0)
sel <- p$sel
pop <- p$pop
K <- p$k_alt
D <- p$depth
write.csv(pop, file.path(out_dir, "09f_retained_pool_manifest.csv"), row.names = FALSE)
cat("Retained G0:", pop$sample_id[pop$generations == 0], "\n")
if (.Platform$OS.type == "windows") mc_cores <- 1L

for (Ne in Ne_values) {
  idx <- seq_len(nrow(K))
  chunks <- split(idx, ceiling(idx / chunk_size))
  ans <- mclapply(chunks, function(ii) {
    t(vapply(ii, function(i) {
      z <- fit_drift_site_null(K[i, ], D[i, ], pop, Ne,
                               pool_individuals = pool_individuals)
      c(p0 = z$p0, p0_null = z$p0_null, slope = z$slope,
        SE = z$SE, z = z$z, p = z$p,
        n_pools = z$n_pools, n_G0 = z$n_G0,
        n_evolved_reps = z$n_evolved_reps)
    }, FUN.VALUE = c(p0 = 0, p0_null = 0, slope = 0,
                     SE = 0, z = 0, p = 0,
                     n_pools = 0, n_G0 = 0, n_evolved_reps = 0)))
  }, mc.cores = mc_cores, mc.preschedule = TRUE)
  ans <- do.call(rbind, ans)
  if (nrow(ans) != nrow(sel)) stop("Output rows changed order")
  result <- data.frame(
    rs = sel$rs, ID = sel$ID, CHROM = sel$CHROM, POS = sel$POS,
    orientation = sel$orientation,
    p0_hat = ans[, "p0"],
    p0_null = ans[, "p0_null"],
    slope_freq = ans[, "slope"],
    SE_freq = ans[, "SE"], Z_drift = ans[, "z"],
    p_drift = ans[, "p"],
    n_pools = as.integer(ans[, "n_pools"]),
    n_G0 = as.integer(ans[, "n_G0"]),
    n_evolved_reps = as.integer(ans[, "n_evolved_reps"]),
    stringsAsFactors = FALSE
  )
  result$common_G0 <- is.finite(result$p0_hat) &
    result$p0_hat >= common_range[1] & result$p0_hat <= common_range[2]
  result$p_drift_BH_all <- p.adjust(result$p_drift, "BH")
  result$p_drift_BY_all <- p.adjust(result$p_drift, "BY")
  result$p_drift_BH_common <- NA_real_
  ix <- which(result$common_G0 & is.finite(result$p_drift))
  result$p_drift_BH_common[ix] <- p.adjust(result$p_drift[ix], "BH")
  result$Ne <- Ne
  result$excluded_G0 <- exclude_G0
  result$retained_G0 <- pop$sample_id[pop$generations == 0]
  result$test <- "WF_GLS_null_intercept_variance_approximation"
  fn <- file.path(out_dir, paste0("09f_SNP_drift_Ne", Ne, ".csv"))
  write.csv(result, fn, row.names = FALSE)
  cat("Wrote", fn, "\n")
  print(data.frame(
    Ne = Ne, tested = sum(is.finite(result$p_drift)),
    n_common = length(ix),
    BH01_all_PROVISIONAL = sum(result$p_drift_BH_all < .01, na.rm = TRUE),
    BH01_common_PROVISIONAL = sum(result$p_drift_BH_common < .01, na.rm = TRUE),
    min_BH_all = min(result$p_drift_BH_all, na.rm = TRUE)))
}
capture.output(sessionInfo(), file = file.path(out_dir, "09f_sessionInfo.txt"))
cat("All discovery counts remain PROVISIONAL until neutral tail calibration.\n")
