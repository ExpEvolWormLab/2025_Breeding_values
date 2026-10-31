## 12c: controlled diagnostic of the legacy-vs-drift discrepancy.
## One high-coverage, common-frequency SNP per composite-GWAS cluster.
## Only diagnostic: NOT a replacement for the full cluster-level analysis.
## Models are refitted at the identical ALT-oriented read counts, varying
## exactly one feature at a time (two vs one G0; GLMM vs drift GLS).

library(lme4)

functions_file <- if (file.exists('Oct_pipeline_cleaner/09f_drift_functions.R'))
  '09f_drift_functions.R' else 'drift_nullfit_update/09f_drift_functions.R'
source('Oct_pipeline_cleaner/09f_drift_functions.R')

pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
selection_file <- 'selection_results/09_Selection_SNP_results.csv'
phase_file <- 'selection_results/drift_oneG0_nullfit/R5/11_SNP_phase_aligned_09f.csv'
qtl_file <- 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv'
legacy_file <- 'RData_GridLMM_NaCl/Estimates_s.txt'
out_dir <- 'selection_results/drift_oneG0_nullfit/R5/discrepancy_audit'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Must be the independently documented bottlenecked G0.
exclude_G0 <- 'A6G140'
Ne <- 1000
max_clusters <- Inf  # Set 100 for a quick pilot; Inf uses all eligible clusters.

p <- read_drift_inputs(pool_file, selection_file, exclude_G0 = character())
K <- p$k_alt
D <- p$depth
pop2 <- p$pop
use1 <- pop2$sample_id != exclude_G0
pop1 <- pop2[use1, , drop = FALSE]

phase <- read.csv(phase_file, stringsAsFactors = FALSE)
qtl <- read.csv(qtl_file, stringsAsFactors = FALSE)
old <- read.table(legacy_file, header = TRUE, sep = '\t',
                  stringsAsFactors = FALSE, check.names = FALSE)
old$ID <- sub('\\.', '_', as.character(old$rs))
if (anyDuplicated(old$ID)) stop('Duplicate legacy SNP IDs')

## Only clusters associated with the composite phenotype. The candidate
## SNP is selected by read depth, NEVER by slope or significance.
qtl_keep <- as.character(qtl$group[
  is.finite(qtl$p_value_REML.cor_beta_NaCl) &
    qtl$p_value_REML.cor_beta_NaCl < .05])
phase <- phase[!is.na(phase$tested_common) & phase$tested_common &
                 as.character(phase$group) %in% qtl_keep &
                 as.character(phase$rs) %in% old$ID, , drop = FALSE]
ix <- match(as.character(phase$rs), as.character(p$sel$ID))
phase <- phase[!is.na(ix), , drop = FALSE]
ix <- ix[!is.na(ix)]
phase$sel_row <- ix
phase$median_depth <- apply(D[ix, , drop = FALSE], 1, median, na.rm = TRUE)
cluster_sizes <- table(as.character(phase$group))
phase <- phase[order(as.character(phase$group), -phase$median_depth), , drop = FALSE]
phase <- phase[!duplicated(as.character(phase$group)), , drop = FALSE]
phase <- phase[seq_len(min(nrow(phase), max_clusters)), , drop = FALSE]
cat('Representative composite-QTL SNPs:', nrow(phase), '\n')

fit_glm <- function(k, d, pop) {
  good <- is.finite(k) & is.finite(d) & d > 0 & k >= 0 & k <= d
  if (sum(good) < 5 || length(unique(pop$generations[good])) < 3) return(NA_real_)
  dat <- data.frame(k = round(k[good]), other = round(d[good] - k[good]),
                    generations = pop$generations[good])
  z <- try(suppressWarnings(glm(cbind(k, other) ~ generations,
                                 family = binomial(), data = dat)), silent = TRUE)
  if (inherits(z, 'try-error')) return(NA_real_)
  unname(coef(z)['generations'])
}

fit_glmm <- function(k, d, pop) {
  good <- is.finite(k) & is.finite(d) & d > 0 & k >= 0 & k <= d
  if (sum(good) < 5 || length(unique(pop$generations[good])) < 3) return(NA_real_)
  dat <- data.frame(k = round(k[good]), other = round(d[good] - k[good]),
                    generations = pop$generations[good],
                    replicate_fact = factor(pop$replicate[good]))
  z <- try(suppressWarnings(glmer(cbind(k, other) ~ generations +
                    (1 | replicate_fact), family = binomial(), data = dat)),
           silent = TRUE)
  if (inherits(z, 'try-error')) return(NA_real_)
  unname(fixef(z)['generations'])
}

rows <- vector('list', nrow(phase))
for (i in seq_len(nrow(phase))) {
  snp_row <- phase$sel_row[i]
  k2 <- as.numeric(K[snp_row, ])
  d2 <- as.numeric(D[snp_row, ])
  k1 <- k2[use1]
  d1 <- d2[use1]
  z1 <- fit_drift_site_null(k1, d1, pop1, Ne = Ne)
  z2 <- fit_drift_site_null(k2, d2, pop2, Ne = Ne)
  o <- match(as.character(phase$rs[i]), old$ID)
  grp <- match(as.character(phase$group[i]), as.character(qtl$group))
  orient <- as.numeric(phase$phase[i])
  p0 <- as.numeric(phase$p0_cluster[i])
  rows[[i]] <- data.frame(
    group = as.character(phase$group[i]), SNP_ID = as.character(phase$rs[i]),
    n_SNPs_in_cluster = as.integer(cluster_sizes[as.character(phase$group[i])]),
    beta = qtl$beta_beta_NaCl[grp],
    phase = orient, p0_cluster = p0,
    median_depth = phase$median_depth[i],
    legacy_s_raw = as.numeric(old$s[o]),
    glm_alt_2G0 = fit_glm(k2, d2, pop2),
    glm_alt_1G0 = fit_glm(k1, d1, pop1),
    glmm_alt_2G0 = fit_glmm(k2, d2, pop2),
    glmm_alt_1G0 = fit_glmm(k1, d1, pop1),
    gls_alt_2G0 = z2$slope,
    gls_alt_1G0 = z1$slope,
    stringsAsFactors = FALSE
  )
  if (i %% 25L == 0L) cat('Completed', i, 'of', nrow(phase), '\n')
}
d <- do.call(rbind, rows)
write.csv(d, file.path(out_dir, '12c_Controlled_locus_slopes.csv'), row.names = FALSE)

## All GWAS-comparisons must be phase-aligned (script 03 convention).
## Binomial logit slopes are converted locally to frequency slopes for a
## descriptive *scale check*; p0 may differ between the two G0 scenarios.
for (v in c('glm_alt_2G0', 'glm_alt_1G0',
            'glmm_alt_2G0', 'glmm_alt_1G0')) {
  d[[paste0(v, '_freq_approx')]] <-
    d[[v]] * d$phase * d$p0_cluster * (1 - d$p0_cluster)
}
for (v in c('gls_alt_2G0', 'gls_alt_1G0'))
  d[[paste0(v, '_cluster')]] <- d[[v]] * d$phase

## Do not ASSUME the legacy fitted count was ALT rather than REF: show both.
d$legacy_as_ALT_freq <- d$legacy_s_raw * d$phase *
  d$p0_cluster * (1 - d$p0_cluster)
d$legacy_as_REF_freq <- -d$legacy_as_ALT_freq

metric <- function(x, y, label) {
  ok <- is.finite(x) & is.finite(y)
  data.frame(comparison = label, n = sum(ok),
    pearson_r = if (sum(ok) >= 5) cor(x[ok], y[ok]) else NA_real_,
    spearman_r = if (sum(ok) >= 5) cor(x[ok], y[ok], method = 'spearman') else NA_real_,
    sign_agreement = if (sum(ok)) mean(sign(x[ok]) == sign(y[ok])) else NA_real_)
}

## Differences isolate one modeling choice on a common set of loci.
comparisons <- do.call(rbind, list(
  metric(d$legacy_s_raw, d$glmm_alt_2G0, 'Legacy file vs refitted ALT GLMM, 2 G0 (orientation audit)'),
  metric(d$glmm_alt_2G0, d$glmm_alt_1G0, 'ALT GLMM, 2 G0 vs 1 G0'),
  metric(d$gls_alt_2G0, d$gls_alt_1G0, 'Drift GLS, 2 G0 vs 1 G0'),
  metric(d$glmm_alt_1G0_freq_approx, d$gls_alt_1G0_cluster,
         'Matched ALT GLMM vs drift GLS, 1 G0 (frequency units)'),
  metric(d$glmm_alt_2G0_freq_approx, d$gls_alt_2G0_cluster,
         'Matched ALT GLMM vs drift GLS, 2 G0 (frequency units)')
))

response_cols <- c('legacy_as_ALT_freq', 'legacy_as_REF_freq',
                   'glm_alt_2G0_freq_approx', 'glm_alt_1G0_freq_approx',
                   'glmm_alt_2G0_freq_approx', 'glmm_alt_1G0_freq_approx',
                   'gls_alt_2G0_cluster', 'gls_alt_1G0_cluster')
associations <- do.call(rbind, lapply(response_cols, function(v) {
  z <- d[[v]]
  ok <- is.finite(z) & is.finite(d$beta)
  data.frame(response = v, n = sum(ok),
             r_beta = if (sum(ok) >= 5) cor(d$beta[ok], z[ok]) else NA_real_,
             p_descriptive = if (sum(ok) >= 5)
               cor.test(d$beta[ok], z[ok])$p.value else NA_real_)
}))
write.csv(comparisons, file.path(out_dir, '12c_Model_pairwise_comparisons.csv'), row.names = FALSE)
write.csv(associations, file.path(out_dir, '12c_Beta_by_model.csv'), row.names = FALSE)
cat('\nModel comparison:\n'); print(comparisons, row.names = FALSE)
cat('\nBeta associations:\n'); print(associations, row.names = FALSE)
cat('\nWARNING: Representative-locus pilot only. These descriptive p-values\n',
    'ignore LD and QTL ascertainment. Do not use as manuscript tests.\n')
