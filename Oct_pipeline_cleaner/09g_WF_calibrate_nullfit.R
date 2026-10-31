## Calibrate 09f GLS test under exact discrete-generation Wright-Fisher drift.
## By-frequency Monte Carlo checks and independent loci, not genome-wide MC FDR.
source("Oct_pipeline_cleaner/09f_drift_functions.R")

pool_file <- "RData_GridLMM_NaCl/08_PoolSeq_prepared.RData"
selection_file <- "selection_results/09_Selection_SNP_results.csv"
out_dir <- "selection_results/drift_oneG0_nullfit"
exclude_G0 <- "A6G140"   # Must match 09f.

Ne_values <- c(500, 1000, 2000)
n_sim <- 500L             # Start modestly, increase for tail assessment.
random_per_freq_bin <- 12L # 8 bins x 12 loci x 500 x 3 Ne (approx).
pool_individuals <- NA_real_
set.seed(20261028)

p <- read_drift_inputs(pool_file, selection_file, exclude_G0 = exclude_G0)
sel <- p$sel
pop <- p$pop
K <- p$k_alt
D <- p$depth
res <- read.csv(file.path(out_dir, "09f_SNP_drift_Ne1000.csv"),
                stringsAsFactors = FALSE)
if (!identical(res$rs, sel$rs)) stop("Script 09f and 09 SNP orders differ")
valid <- which(is.finite(res$p0_hat) & res$p0_hat > 0 & res$p0_hat < 1 &
               is.finite(res$p_drift))
bins <- c(0, .02, .05, .10, .25, .50, .75, .90, .98, 1)
frequency_bin <- cut(res$p0_hat, bins, include.lowest = TRUE)
ix_per_bin <- split(valid, frequency_bin[valid], drop = TRUE)
selected <- unlist(lapply(ix_per_bin, function(ii)
  sample(ii, size = min(length(ii), random_per_freq_bin))), use.names = FALSE)
cat("Calibrating", length(selected), "loci x", n_sim,
    "simulations x", length(Ne_values), "Ne scenarios\n")

loci <- vector("list", length(selected) * length(Ne_values))
by_bin <- list()
q <- 0L
for (Ne in Ne_values) {
  for (j in seq_along(selected)) {
    i <- selected[j]
    fit <- fit_drift_site_null(K[i, ], D[i, ], pop, Ne,
                               pool_individuals = pool_individuals)
    if (!is.finite(fit$z)) next
    sim <- simulate_wf_site(fit$p0, D[i, ], pop, Ne = Ne,
                            n_sim = n_sim, pool_individuals = pool_individuals)
    ## Exact algebraic refit for every simulation, using one fixed GLS
    ## geometry (depths, generations, Ne). Much faster than repeated chol().
    ## Both the slope and the restricted null-frequency estimate are refit.
    slope_sim <- drop(sim %*% fit$weights)
    p_null_sim <- drop(sim %*% fit$null_weights)
    good <- is.finite(K[i, ]) & is.finite(D[i, ]) &
      D[i, ] > 0 & K[i, ] >= 0 & K[i, ] <= D[i, ]
    lower <- 0.5 / (sum(D[i, good]) + 1)
    p_null_sim <- pmax(lower, pmin(1 - lower, p_null_sim))
    z <- slope_sim / (fit$slope_scale * sqrt(p_null_sim * (1 - p_null_sim)))
    z <- z[is.finite(z)]
    if (!length(z)) next
    pp <- 2 * pnorm(abs(z), lower.tail = FALSE)
    q <- q + 1L
    loci[[q]] <- data.frame(
      Ne = Ne, ID = sel$ID[i], freq_bin = as.character(frequency_bin[i]),
      p0_hat = fit$p0, p0_null = fit$p0_null, n_sim = length(z),
      meanZ = mean(z), sdZ = sd(z),
      fp05 = mean(pp < .05), fp01 = mean(pp < .01),
      fp001 = mean(pp < .001), fp0001 = mean(pp < .0001),
      stringsAsFactors = FALSE)
    by_bin[[paste(Ne, i)]] <- data.frame(Ne = Ne,
         freq_bin = as.character(frequency_bin[i]), z = z)
    if (j %% 10L == 0L)
      cat("Ne", Ne, "loci", j, "/", length(selected), "\n")
  }
}
loci <- do.call(rbind, loci[seq_len(q)])
write.csv(loci, file.path(out_dir, "09g_WF_locus_calibration.csv"), row.names = FALSE)

zz <- do.call(rbind, by_bin)
zz$p <- 2 * pnorm(abs(zz$z), lower.tail = FALSE)
## Aggregate all null draws, with equal simulations per sampled locus.
aggregate_fn <- function(x) data.frame(
    n = length(x), meanZ = mean(x), sdZ = sd(x),
    fp05 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .05),
    fp01 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .01),
    fp001 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .001),
    fp0001 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .0001))
## Group without parsing labels (which contain decimal points).
groups <- split(seq_len(nrow(zz)), interaction(zz$Ne, zz$freq_bin, drop = TRUE))
summary_by_bin <- do.call(rbind, lapply(groups, function(ii) {
  a <- aggregate_fn(zz$z[ii]);
  cbind(Ne = zz$Ne[ii[1]], freq_bin = zz$freq_bin[ii[1]], a)
}))
row.names(summary_by_bin) <- NULL
write.csv(summary_by_bin,
          file.path(out_dir, "09g_WF_calibration_by_frequency.csv"), row.names = FALSE)
groups_Ne <- split(zz$z, zz$Ne)
summary_all <- do.call(rbind, lapply(names(groups_Ne), function(N) {
  cbind(Ne = as.integer(N), aggregate_fn(groups_Ne[[N]]))
}))
row.names(summary_all) <- NULL
write.csv(summary_all,
          file.path(out_dir, "09g_WF_null_calibration_summary.csv"), row.names = FALSE)
print(summary_all, row.names = FALSE)
cat("\nBy frequency bin:\n")
print(summary_by_bin, row.names = FALSE)
cat("Neutral calibration is scenario-conditional and does not account\n",
    "for unmodeled biological Pool-seq sampling or selection heterogeneity.\n")
