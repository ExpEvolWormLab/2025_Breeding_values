## 15_Final_prediction_chrIV_analysis.R
## Final, focused checks for the multivariate GWAS -> evolutionary-response result.
## Run from the project root AFTER corrected scripts 08a, 09f, 11 and 13.
## NO refitting of GWAS, Pool-seq drift models, or MCMCglmm.
## The main analysis uses ancestral polymorphisms (G0 ALT p in [0.1,0.9]),
## phase-aligned cluster slopes, and the same 675 R5 LD clusters as script 13.
##
## Outputs are descriptive and block-bootstrap sensitivities, not causal tests.
## Candidate chromosome-IV mutations are entered separately: the founder-SNP
## coverage filtering in script 08 generally excludes genuinely novel alleles.

## ------------------------------------------------------------
## 0. Settings
## ------------------------------------------------------------
label <- "R5"
base_dir <- file.path("selection_results", "drift_oneG0_nullfit", label,
                      "multivariate_prediction")
out_dir <- file.path(base_dir, "final_chrIV")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20261028L)

block_bp <- 2e6
n_boot <- 2000L
n_placebo <- 2000L
windows_bp <- c(250000, 500000, 1000000, 2000000)
chromosomes <- c("I", "II", "III", "IV", "V", "X")
mutation_file <- "15_chrIV_mutations.csv"       # user-supplied positions
mutation_counts_file <- "15_chrIV_mutation_counts.csv"  # optional, independently measured
manifest_file <- "Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv"

read_required <- function(fname, fields) {
  path <- file.path(base_dir, fname)
  if (!file.exists(path)) stop("Missing script-13 output: ", path)
  z <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(fields, colnames(z))
  if (length(missing)) stop(fname, " lacks: ", paste(missing, collapse = ", "))
  z
}
finite_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}
pred_R2 <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  if (sum(ok) < 5L) return(NA_real_)
  sst <- sum((y[ok] - mean(y[ok]))^2)
  if (sst <= 0) return(NA_real_)
  1 - sum((y[ok] - yhat[ok])^2) / sst
}
err_pair <- function(d, preferable, comparator) {
  mean((d$response - d[[comparator]])^2 -
       (d$response - d[[preferable]])^2)
}
models <- c("full_reconstructed", "full_mapped", "three_trait",
            "SB", "area.F")
labels <- c(full_reconstructed = "Full index, reconstructed",
            full_mapped = "Mapped composite beta",
            three_trait = "SF + SB + area.F",
            SB = "SB alone", area.F = "area.F alone")

## ------------------------------------------------------------
## 1. Load MATCHED cluster table and held-out predictions
## ------------------------------------------------------------
qtl <- read_required("13_Cluster_predictors_and_response.csv",
                     c("group", "chromosome", "position", "response",
                       "pred_full_reconstructed", "pred_three_trait",
                       "pred_full_GWAS", "beta_area.F_NaCl", "beta_SB_NaCl"))
cv <- read_required("13_Leave_chromosome_out_fitted.csv",
                    c("group", "chromosome", "response", "full_reconstructed",
                      "full_mapped", "three_trait", "SB", "area.F"))
if (anyDuplicated(qtl$group) || anyDuplicated(cv$group))
  stop("Cluster IDs are not unique")
ix <- match(as.character(cv$group), as.character(qtl$group))
if (anyNA(ix) || nrow(qtl) != nrow(cv))
  stop("Script 13's cluster and held-out tables do not have the same R5 groups")
q <- qtl[ix, , drop = FALSE]
if (!identical(as.character(cv$chromosome), as.character(q$chromosome)) ||
    max(abs(cv$response - q$response), na.rm = TRUE) > 1e-10)
  stop("Cluster chromosome or evolutionary-response mismatch between files")
cv$position <- as.numeric(q$position)
cv$block <- paste(cv$chromosome, floor(cv$position / block_bp), sep = ":")
if (any(!is.finite(cv$position)) || any(!is.finite(cv$response)) ||
    !all(cv$chromosome %in% chromosomes)) stop("Invalid cluster positions/response")
if (any(!is.finite(as.matrix(cv[, models, drop = FALSE]))))
  stop("Some script-13 held-out predictions are nonfinite")
cv <- cv[order(match(cv$chromosome, chromosomes), cv$position), , drop = FALSE]
q <- q[match(as.character(cv$group), as.character(q$group)), , drop = FALSE]
rownames(cv) <- NULL

## ------------------------------------------------------------
## 2. Fit and evaluate a truly held-out area.F-only benchmark
## ------------------------------------------------------------
## The area.F effect is fitted to the evolutionary response on 5 chromosomes;
## the sixth is used ONLY for evaluating the fitted regression. This makes
## the area.F-only benchmark directly comparable to script 13's LOCO scores.
cv$area.F_LOCO <- NA_real_
for (chr in chromosomes) {
  train <- cv$chromosome != chr
  test <- !train
  x_train <- as.numeric(q$`beta_area.F_NaCl`[train])
  x_test <- as.numeric(q$`beta_area.F_NaCl`[test])
  y_train <- cv$response[train]
  if (sd(x_train) > 0) {
    fit <- lm(y_train ~ x_train)
    cv$area.F_LOCO[test] <- unname(coef(fit)[1] + coef(fit)[2] * x_test)
  } else {
    cv$area.F_LOCO[test] <- mean(y_train)
  }
}
if (any(!is.finite(cv$area.F_LOCO))) stop("area.F benchmark prediction failed")
models <- c(models, "area.F_LOCO")
labels <- c(labels, area.F_LOCO = "area.F alone (new LOCO fit)")
## Keep the newly refitted area.F model as the ONLY area.F evaluation column.
## Script 13 already generated area.F LOCO predictions; compare as a strict QC.
max_difference <- max(abs(cv$area.F_LOCO - cv$area.F))
cat("Maximum discrepancy in area.F held-out predictions:", max_difference, "\n")
if (max_difference > 1e-7)
  warning("area.F re-fit and script-13 area.F predictions differ. Check GWAS effect columns.")

## ------------------------------------------------------------
## 3. Matched summary: all clusters, non-IV, and IV only
## ------------------------------------------------------------
make_summary <- function(which_rows, cohort) {
  dd <- cv[which_rows, , drop = FALSE]
  qq <- q[match(dd$group, q$group), , drop = FALSE]
  cols <- c(full_reconstructed = "pred_full_reconstructed",
            full_mapped = "pred_full_GWAS",
            three_trait = "pred_three_trait",
            SB = "beta_SB_NaCl", area.F_LOCO = "beta_area.F_NaCl")
  do.call(rbind, lapply(names(cols), function(m) {
    x <- qq[[cols[[m]]]]
    yhat <- dd[[m]]
    data.frame(cohort = cohort, model = m, n_clusters = nrow(dd),
      in_sample_r = finite_cor(x, dd$response),
      heldout_R2 = pred_R2(dd$response, yhat),
      heldout_MSE = mean((dd$response - yhat)^2),
      stringsAsFactors = FALSE)
  }))
}
cohorts <- list(all_675 = rep(TRUE, nrow(cv)),
                without_IV = cv$chromosome != "IV",
                IV_only = cv$chromosome == "IV")
summary_result <- do.call(rbind, Map(make_summary, cohorts, names(cohorts)))
write.csv(summary_result, file.path(out_dir, "15_Model_comparison.csv"), row.names = FALSE)

## To avoid conflating model comparison with regression re-fitting, compute
## chromosome-stratified differences in SQUARED ERRORS using exactly the
## same fitted predictions (no training fold is re-estimated here).
comparisons <- list(c("full_reconstructed", "three_trait"),
                    c("full_mapped", "three_trait"),
                    c("full_reconstructed", "SB"),
                    c("full_reconstructed", "area.F_LOCO"),
                    c("full_reconstructed", "full_mapped"))
contrast_one <- function(d, chr, a, b) {
  delta <- (d$response - d[[b]])^2 - (d$response - d[[a]])^2
  data.frame(chromosome = chr, preferred = a, comparator = b,
             n_clusters = nrow(d), delta_MSE = mean(delta),
             n_positive_cluster_errors = sum(delta > 0),
             frac_clusters_favor_preferred = mean(delta > 0),
             stringsAsFactors = FALSE)
}
by_chr <- do.call(rbind, lapply(chromosomes, function(chr) {
  d <- cv[cv$chromosome == chr, , drop = FALSE]
  do.call(rbind, lapply(comparisons, function(pair)
    contrast_one(d, chr, pair[1], pair[2])))
}))
write.csv(by_chr, file.path(out_dir, "15_By_chromosome_paired_errors.csv"),
          row.names = FALSE)

## ------------------------------------------------------------
## 4. Bootstrap prediction-error differences by 2-Mb genomic blocks
## ------------------------------------------------------------
## Predictions remain FIXED, as already trained on other chromosomes.
## This bootstraps the distribution of held-out prediction errors over
## genomic blocks; it is not nested model refitting and not an independent
## biological replication of the selection experiment.
boot_groups <- split(seq_len(nrow(cv)), cv$block)
all_blocks <- unname(boot_groups)
block_boot <- do.call(rbind, lapply(comparisons, function(pair) {
  a <- pair[1]; b <- pair[2]
  loss_diff <- (cv$response - cv[[b]])^2 - (cv$response - cv[[a]])^2
  vals <- replicate(n_boot, {
    ii <- unlist(all_blocks[sample.int(length(all_blocks),
                                        length(all_blocks), replace = TRUE)],
                 use.names = FALSE)
    mean(loss_diff[ii])
  })
  data.frame(preferred = a, comparator = b, n_blocks = length(all_blocks),
             n_boot = n_boot, delta_MSE = mean(loss_diff),
             lower95 = unname(quantile(vals, .025)),
             upper95 = unname(quantile(vals, .975)),
             fraction_bootstrap_positive = mean(vals > 0),
             stringsAsFactors = FALSE)
}))
write.csv(block_boot, file.path(out_dir, "15_Heldout_block_bootstrap.csv"),
          row.names = FALSE)

## ------------------------------------------------------------
## 5. Chromosome IV: where do the prediction errors occur?
## ------------------------------------------------------------
iv <- cv[cv$chromosome == "IV", , drop = FALSE]
qi <- q[match(iv$group, q$group), , drop = FALSE]
iv$GWAS_SB <- as.numeric(qi$beta_SB_NaCl)
iv$GWAS_beta <- as.numeric(qi$pred_full_GWAS)
iv$delta_sqerr_full_vs_three <-
  (iv$response - iv$three_trait)^2 -
  (iv$response - iv$full_reconstructed)^2
iv$delta_sqerr_beta_vs_three <-
  (iv$response - iv$three_trait)^2 -
  (iv$response - iv$full_mapped)^2
iv$delta_sqerr_full_vs_SB <-
  (iv$response - iv$SB)^2 -
  (iv$response - iv$full_reconstructed)^2
iv$abs_error_full <- abs(iv$response - iv$full_reconstructed)
iv$abs_error_three <- abs(iv$response - iv$three_trait)
write.csv(iv, file.path(out_dir, "15_CHRIV_cluster_prediction_errors.csv"),
          row.names = FALSE)

## ------------------------------------------------------------
## 6. Optional mutation metadata: MUST be independently documented
## ------------------------------------------------------------
## These are not ordinary founder SNPs. Do NOT assume the filtered
## 08_PoolSeq_prepared.RData contains them. Record chromosome-IV coordinate
## on the SAME genome assembly as qtl$position and explicitly verify origin.
if (!file.exists(mutation_file)) {
  template <- data.frame(mutation_id = c("CHRIV_mutation_1", "CHRIV_mutation_2"),
                         chromosome = c("IV", "IV"),
                         position_bp = c(NA_real_, NA_real_),
                         sb_association_note = c("", ""),
                         origin_evidence = c("", ""),
                         stringsAsFactors = FALSE)
  write.csv(template, mutation_file, row.names = FALSE, na = "")
  message("Created editable mutation position template: ", mutation_file)
}
mut <- read.csv(mutation_file, stringsAsFactors = FALSE)
need_mut <- c("mutation_id", "chromosome", "position_bp")
if (!all(need_mut %in% names(mut)))
  stop("Mutation positions CSV requires: ", paste(need_mut, collapse = ", "))
if (anyDuplicated(mut$mutation_id)) stop("Duplicate mutation IDs")
mut$position_bp <- suppressWarnings(as.numeric(mut$position_bp))
if (any(is.finite(mut$position_bp) & mut$chromosome != "IV"))
  stop("Positions supplied in the chrIV file must be on chromosome IV")
has_mut_positions <- all(is.finite(mut$position_bp)) && nrow(mut) >= 1L
if (!has_mut_positions)
  message("Mutation coordinates are not populated; spatial-locality analysis skipped.")

## Regional contrast: compare removing candidate neighborhoods with removing
## neighborhoods around random observed chrIV cluster positions. This is
## EXPLORATORY and conditional on the chosen windows, not a test of causality.
region_result <- data.frame()
placebo_result <- data.frame()
if (has_mut_positions) {
  positions <- mut$position_bp
  iv$nearest_mutation_bp <- vapply(iv$position, function(x)
    min(abs(x - positions)), numeric(1))
  write.csv(iv, file.path(out_dir, "15_CHRIV_cluster_prediction_errors.csv"),
            row.names = FALSE)
  for (win in windows_bp) {
    near <- iv$nearest_mutation_bp <= win
    if (sum(near) < 1L || sum(!near) < 8L) next
    baseline_delta <- mean(iv$delta_sqerr_full_vs_three)
    outside_delta <- mean(iv$delta_sqerr_full_vs_three[!near])
    inside_delta <- mean(iv$delta_sqerr_full_vs_three[near])
    shift <- outside_delta - baseline_delta
    ## Positions are sampled from observed IV clusters to keep the number
    ## and density of available evaluation loci geographically plausible.
    placebo <- replicate(n_placebo, {
      anchor <- sample(iv$position, length(positions), replace = FALSE)
      out <- vapply(iv$position, function(x)
        min(abs(x - anchor)) > win, logical(1))
      if (sum(out) < 8L || sum(!out) < 1L) NA_real_
      else mean(iv$delta_sqerr_full_vs_three[out]) - baseline_delta
    })
    placebo <- placebo[is.finite(placebo)]
    region_result <- rbind(region_result, data.frame(
      window_bp = win, n_near = sum(near), n_outside = sum(!near),
      delta_MSE_before_exclusion = baseline_delta,
      delta_MSE_inside = inside_delta, delta_MSE_outside = outside_delta,
      improvement_in_delta_after_exclusion = shift,
      percentile_vs_random_IV_centers = if (length(placebo))
        mean(placebo <= shift) else NA_real_,
      n_valid_placebo = length(placebo), stringsAsFactors = FALSE))
    if (length(placebo)) placebo_result <- rbind(placebo_result,
      data.frame(window_bp = win, shift = placebo))
  }
  write.csv(region_result, file.path(out_dir, "15_CHRIV_mutation_window_sensitivity.csv"),
            row.names = FALSE)
  write.csv(placebo_result, file.path(out_dir, "15_CHRIV_random_window_sensitivity.csv"),
            row.names = FALSE)
}

## ------------------------------------------------------------
## 7. Optional trajectories of the TWO mutations, from raw counts
## ------------------------------------------------------------
## Template is a different table from the SNPs filtered to founder alleles.
## Provide mutant_reads and depth from the ORIGINAL locus-specific sync or
## another independently checked read-count source, with a defined mutant
## allele. The plot can also accept independently reported frequency values
## in the allele_frequency column, but counts/depth are preferred.
if (!file.exists(mutation_counts_file)) {
  if (file.exists(manifest_file)) {
    m <- read.csv(manifest_file, stringsAsFactors = FALSE)
    required_manifest <- c("sample_id", "generations", "replicate")
    if (!all(required_manifest %in% names(m)))
      stop("Sample manifest missing expected columns")
    tmpl <- merge(data.frame(mutation_id = mut$mutation_id,
                             stringsAsFactors = FALSE),
                  m[, required_manifest], by = NULL)
    tmpl <- tmpl[, c("mutation_id", "sample_id", "generations", "replicate")]
    tmpl$mutant_reads <- NA_real_
    tmpl$depth <- NA_real_
    tmpl$allele_frequency <- NA_real_
    write.csv(tmpl, mutation_counts_file, row.names = FALSE, na = "")
    message("Created editable mutation-trajectory template: ", mutation_counts_file)
  }
}
traj <- NULL
if (file.exists(mutation_counts_file)) {
  mm <- read.csv(mutation_counts_file, stringsAsFactors = FALSE)
  if (all(c("mutation_id", "sample_id", "generations", "replicate",
            "mutant_reads", "depth", "allele_frequency") %in% names(mm))) {
    mm$mutant_reads <- suppressWarnings(as.numeric(mm$mutant_reads))
    mm$depth <- suppressWarnings(as.numeric(mm$depth))
    mm$allele_frequency <- suppressWarnings(as.numeric(mm$allele_frequency))
    count_ok <- is.finite(mm$mutant_reads) & is.finite(mm$depth) & mm$depth > 0 &
      mm$mutant_reads >= 0 & mm$mutant_reads <= mm$depth
    mm$frequency <- mm$allele_frequency
    mm$frequency[count_ok] <- mm$mutant_reads[count_ok] / mm$depth[count_ok]
    mm$generations <- as.numeric(mm$generations)
    mm$replicate <- as.character(mm$replicate)
    if (any(count_ok & is.finite(mm$allele_frequency) &
            abs(mm$mutant_reads / pmax(mm$depth, 1) - mm$allele_frequency) > 0.03))
      warning("Some supplied mutation allele frequencies disagree with counts by >0.03")
    valid_traj <- is.finite(mm$frequency) & mm$frequency >= 0 &
      mm$frequency <= 1 & is.finite(mm$generations)
    if (any(valid_traj)) {
      traj <- mm[valid_traj, , drop = FALSE]
      write.csv(traj, file.path(out_dir, "15_CHRIV_mutation_trajectories.csv"),
                row.names = FALSE)
    }
  } else warning("Mutation counts table has missing columns; skipping plot")
}

## ------------------------------------------------------------
## 8. Figures: chrIV locus diagnostic, plus optional mutation trajectories
## ------------------------------------------------------------
fig_file <- file.path(out_dir, "15_Chromosome_IV_prediction_diagnostics.pdf")
pdf(fig_file, width = 11, height = 8.2, useDingbats = FALSE)
old <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3, 1.2), oma = c(1.3, 0.4, 0.8, 0.3))
mark_positions <- function() {
  if (has_mut_positions) {
    abline(v = mut$position_bp / 1e6, col = "#BD573E", lty = 2, lwd = 1.6)
    legend("topright", legend = "Candidate mutation positions", lty = 2,
           col = "#BD573E", bty = "n", cex = .72)
  }
}
plot(iv$position / 1e6, iv$delta_sqerr_full_vs_three * 1e8,
     col = "#386b87", pch = 16, cex = .65,
     xlab = "Chromosome IV position (Mb)",
     ylab = "Reduction in squared error (x 1e-8)",
     main = "A  Paired prediction errors", bty = "n")
abline(h = 0, lty = 2, col = "grey55")
mark_positions()
block_mean <- aggregate(delta_sqerr_full_vs_three ~ block,
                        data = iv, FUN = mean)
block_mid <- aggregate(position ~ block, data = iv, FUN = median)
block_mean$pos <- block_mid$position[match(block_mean$block, block_mid$block)]
plot(block_mean$pos / 1e6, block_mean$delta_sqerr_full_vs_three * 1e8,
     col = "#386b87", pch = 19, cex = .85,
     xlab = "Chromosome IV position (Mb)",
     ylab = "Mean block reduction (x 1e-8)",
     main = "B  2-Mb genomic blocks", bty = "n")
abline(h = 0, lty = 2, col = "grey55")
mark_positions()
chr_cmp <- by_chr[by_chr$preferred == "full_reconstructed" &
                      by_chr$comparator == "three_trait", ]
chr_cmp <- chr_cmp[match(chromosomes, chr_cmp$chromosome), ]
barplot(chr_cmp$delta_MSE * 1e8, names.arg = chr_cmp$chromosome,
        col = ifelse(chr_cmp$chromosome == "IV", "#BD573E", "#386b87"),
        border = NA, xlab = "Chromosome", ylab = "Reduction in MSE (x 1e-8)",
        main = "C  Leave-one-chromosome-out", las = 1)
abline(h = 0, lty = 2, col = "grey40")
barcols <- c("#386b87", "#6e9997", "#A48252", "#87829c", "#7a848e")
sumall <- summary_result[summary_result$cohort == "all_675", ]
sumall <- sumall[match(c("full_reconstructed", "three_trait", "SB",
                         "area.F_LOCO", "full_mapped"), sumall$model), ]
barplot(sumall$heldout_R2, names.arg = c("Full 7", "Three", "SB", "area.F", "Mapped"),
        col = barcols, border = NA, ylab = "Held-out predictive R2",
        main = "D  Matched 675-cluster prediction", las = 1)
abline(h = 0, lty = 2, col = "grey50")
mtext("Positive error reduction means the full index predicts better. Mutations are candidate explanations, not tested causes.",
      side = 1, outer = TRUE, cex = .75)
par(old)
dev.off()

if (!is.null(traj)) {
  pdf(file.path(out_dir, "15_CHRIV_mutation_trajectories.pdf"),
      width = 9.5, height = 5.2, useDingbats = FALSE)
  old <- par(no.readonly = TRUE)
  unique_mut <- unique(traj$mutation_id)
  par(mfrow = c(1, length(unique_mut)), mar = c(4.3, 4.5, 3, 1))
  for (id in unique_mut) {
    dd <- traj[traj$mutation_id == id, , drop = FALSE]
    plot(NA, xlim = range(c(0, dd$generations)), ylim = c(0, max(.35, dd$frequency)),
         xlab = "Generation", ylab = "Mutant allele frequency", main = id,
         bty = "n")
    abline(h = 0, col = "grey85")
    evol <- sort(unique(dd$replicate[dd$generations > 0]))
    cols <- setNames(grDevices::hcl.colors(max(3, length(evol)), "Dark 3")
                     [seq_along(evol)], evol)
    for (r in evol) {
      sub <- dd[dd$replicate == r & dd$generations > 0, , drop = FALSE]
      sub <- sub[order(sub$generations), , drop = FALSE]
      lines(sub$generations, sub$frequency, col = cols[[r]], lwd = 2)
      points(sub$generations, sub$frequency, col = cols[[r]], pch = 19)
    }
    g0 <- dd[dd$generations == 0, , drop = FALSE]
    if (nrow(g0)) points(g0$generations, g0$frequency,
                         pch = 4, col = "black", cex = 1.3, lwd = 2)
    legend("topleft", legend = c(paste("Evolved", evol), "Generation 0"),
           col = c(unname(cols), "black"), lty = c(rep(1, length(evol)), NA),
           pch = c(rep(19, length(evol)), 4), bty = "n", cex = .76)
  }
  par(old)
  dev.off()
}

## ------------------------------------------------------------
## 9. Concise written QC; do not claim spatial localization without positions
## ------------------------------------------------------------
report <- c(
  paste0("R5 clusters: ", nrow(cv), "; chromosome IV: ", nrow(iv)),
  paste0("Overall full-vs-three held-out difference in MSE: ",
         signif(err_pair(cv, "full_reconstructed", "three_trait"), 5)),
  paste0("Chromosome IV full-vs-three held-out difference in MSE: ",
         signif(err_pair(iv, "full_reconstructed", "three_trait"), 5)),
  paste0("Non-IV full-vs-three held-out difference in MSE: ",
         signif(err_pair(cv[cv$chromosome != "IV", ],
                         "full_reconstructed", "three_trait"), 5)),
  paste0("Candidate mutation positions verified in same assembly: ",
         has_mut_positions),
  paste0("Independent candidate mutation trajectories provided: ",
         !is.null(traj)),
  "The one-chromosome-out comparison is a genomic sensitivity analysis, not independent experimental replication.",
  "Mutation-localization window and placebo analyses are exploratory; neither proves mutation causality.",
  "Do not interpret founder-SNP slope estimates as direct tests of novel mutation trajectories."
)
writeLines(report, file.path(out_dir, "15_Analysis_QC.txt"))
cat(paste(report, collapse = "\n"), "\n")
cat("Finished. Outputs: ", out_dir, "\n")
