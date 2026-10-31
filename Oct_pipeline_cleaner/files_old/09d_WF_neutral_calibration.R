## 09d_WF_neutral_calibration.R
## EXPLICIT Wright-Fisher neutral simulations for a manageable SNP subset.
## Use to check whether 09c's approximate p-values are calibrated.
## Monte Carlo p-values from this subset MUST NOT be BH-adjusted genome-wide.
## Run after 09c at Ne=1000. Keeps BOTH generation-zero samples.

source("09c_drift_functions.R")

pool_file <- "RData_GridLMM_NaCl/08_PoolSeq_prepared.RData"
selection_file <- "selection_results/09_Selection_SNP_results.csv"
drfile <- "selection_results/drift/09c_SNP_drift_Ne1000.csv"
out_dir <- "selection_results/drift"

## Pilot settings: modest computational cost.
## For final sensitivity: n_sim = 2000 or more; increase random_per_stratum.
Ne_values <- c(500, 1000, 2000)
n_sim <- 500L
random_per_stratum <- 3L
n_top <- 20L
pool_individuals <- 10000  # True number of pooled individuals, if known
set.seed(20261028)

p <- read_drift_inputs(pool_file, selection_file)
sel <- p$sel
pop <- p$pop
K <- p$k_alt
D <- p$depth
res <- read.csv(drfile, stringsAsFactors = FALSE)
if (!identical(res$rs, sel$rs))
  stop("Script 09c rows do not match the script 09 SNP set")

valid <- which(is.finite(res$Z_drift) & is.finite(res$p0_hat) &
                 res$p0_hat >= .02 & res$p0_hat <= .98)
if (!length(valid)) stop("No valid non-boundary SNPs to calibrate")

## Stratify across ancestral p0 and mean sequencing depth.
p0_breaks <- c(.02, .1, .25, .5, .75, .9, .98)
p0_bin <- cut(res$p0_hat, p0_breaks, include.lowest = TRUE)
mean_depth <- rowMeans(D, na.rm = TRUE)
depth_breaks <- unique(as.numeric(quantile(mean_depth[valid],
                                          c(0, .33, .67, 1), na.rm = TRUE)))
if (length(depth_breaks) < 3L) {
  depth_breaks <- c(min(mean_depth[valid]) - 1,
                    max(mean_depth[valid]) + 1)
}
depth_bin <- cut(mean_depth, depth_breaks, include.lowest = TRUE)
stratum <- interaction(p0_bin, depth_bin, drop = TRUE)

split_ix <- split(valid, stratum[valid], drop = TRUE)
random_ix <- unlist(lapply(split_ix, function(ii)
  sample(ii, min(length(ii), random_per_stratum))), use.names = FALSE)
top_ix <- head(valid[order(res$p_drift[valid])], n_top)
selected_ix <- unique(c(random_ix, top_ix))
role <- ifelse(selected_ix %in% top_ix, "top_candidate", "stratified_random")

cat("Neutral WF calibration: ", length(selected_ix), " SNPs x ",
    n_sim, " simulations x ", length(Ne_values), " Ne scenarios\n", sep = "")
cat("Each simulation follows four independent populations for up to ",
    max(pop$generations), " generations, plus both G0 libraries.\n", sep = "")

out <- vector("list", length(selected_ix) * length(Ne_values))
ix_out <- 0L
calibration_z <- vector("list", length(Ne_values))
names(calibration_z) <- as.character(Ne_values)

for (Ne in Ne_values) {
  null_z <- numeric(0)
  null_z_random <- numeric(0)
  for (a in seq_along(selected_ix)) {
    i <- selected_ix[a]
    fit <- fit_drift_site(K[i, ], D[i, ], pop, Ne,
                          pool_individuals = pool_individuals)
    if (!is.finite(fit$z)) next

    ## Simulate each independent evolutionary replicate with discrete WF drift.
    sim <- simulate_wf_site(fit$p0, D[i, ], pop,
                             Ne = Ne, n_sim = n_sim,
                             pool_individuals = pool_individuals)
    Z <- as.numeric(sim %*% fit$weights) / fit$SE
    Z <- Z[is.finite(Z)]
    if (!length(Z)) next
    pmc <- (1 + sum(abs(Z) >= abs(fit$z))) / (length(Z) + 1)

    ix_out <- ix_out + 1L
    out[[ix_out]] <- data.frame(
      Ne = Ne, ID = sel$ID[i], rs = sel$rs[i],
      role = role[a], stratum = as.character(stratum[i]),
      initial_freq = fit$p0, mean_depth = mean_depth[i],
      z_observed = fit$z, p_analytic = fit$p,
      p_WF_MC = pmc, n_sim = length(Z),
      null_mean_Z = mean(Z), null_sd_Z = sd(Z),
      null_tail_005 = mean(abs(Z) > qnorm(.975)),
      null_tail_001 = mean(abs(Z) > qnorm(.995)),
      stringsAsFactors = FALSE
    )
    null_z <- c(null_z, Z)
    if (role[a] == "stratified_random")
      null_z_random <- c(null_z_random, Z)
    if (a %% 10L == 0L)
      cat("  Ne=", Ne, ": ", a, "/", length(selected_ix),
          " SNPs checked\n", sep = "")
  }
  calibration_z[[as.character(Ne)]] <- null_z_random
}

out <- do.call(rbind, out[seq_len(ix_out)])
write.csv(out, file.path(out_dir, "09d_WF_locus_calibration.csv"),
          row.names = FALSE)

summ <- do.call(rbind, lapply(Ne_values, function(Ne) {
  x <- calibration_z[[as.character(Ne)]]
  data.frame(Ne = Ne, n_simulated_test_statistics = length(x),
             mean_Z = mean(x), sd_Z = sd(x),
             frac_p_less_005 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .05),
             frac_p_less_001 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .01),
             frac_p_less_0001 = mean(2 * pnorm(abs(x), lower.tail = FALSE) < .001))
}))
write.csv(summ, file.path(out_dir, "09d_WF_null_calibration_summary.csv"),
          row.names = FALSE)
print(summ, row.names = FALSE)

pdf(file.path(out_dir, "09d_WF_null_QQ.pdf"), width = 9, height = 3.5)
old <- par(mfrow = c(1, length(Ne_values)), mar = c(4, 4, 3, 1))
for (Ne in Ne_values) {
  zz <- calibration_z[[as.character(Ne)]]
  if (length(zz) >= 5L) {
    qqnorm(zz, main = paste("Ne", Ne), pch = 16,
           cex = .3, xlab = "Normal reference", ylab = "Neutral simulated Z")
    qqline(zz, col = "firebrick")
  } else plot.new()
}
par(old); dev.off()

cat("Saved WF calibration. DO NOT adjust p_WF_MC over all genome-wide SNPs:\n",
    "Monte Carlo p-values apply only to the selected SNP subset,\n",
    "and their minimum resolution is 1 / (n_sim + 1).\n")
