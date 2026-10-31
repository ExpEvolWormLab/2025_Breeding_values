## 13_Test_multivariate_selection_prediction.R
## Does a biologically weighted combination of seven component-trait GWAS
## effects predict the direction of allele-frequency change during evolution?
##
## Run after the CORRECTED Pool-seq sample alignment and scripts 09f, 11.
## Uses R5 LD clusters and phase-oriented 09f frequency slopes (per generation).
## Primary: all clusters with complete GWAS effects and common-frequency slopes.
## Secondary: beta-GWAS q < .05 and q < .01 on matched clusters.
##
## This script DOES NOT refit GWAS/Pool-seq models or modify their outputs.
## Resampling/shift and weight-comparison statistics are robustness diagnostics;
## none fully model GWAS effect uncertainty, ascertainment or residual LD.

## ------------------------------------------------------------
## 0. Settings
## ------------------------------------------------------------
threshold_clust <- 0.5
label <- paste0("R", round(10 * threshold_clust))
Ne <- 1000L
block_bp <- 2e6
n_boot <- 1500L
n_shift <- 1999L
n_random_weights <- 4000L
n_gradient_draws <- 1500L
seed <- 20261028L
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")
direct_traits <- c("SF", "SB", "area.F")
chromosomes <- c("I", "II", "III", "IV", "V", "X")

cluster_file <- file.path("selection_results", "drift_oneG0_nullfit", label,
                          "11_Clusters_selection_and_GWAS_09f.csv")
gradient_file <- "RData/selection_gradients_NaCl.csv"
preparation_file <- "RData/Preparatory_GLMM_NaCl_only_withBeta.RData"
output_dir <- file.path("selection_results", "drift_oneG0_nullfit",
                        label, "multivariate_prediction")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

## ------------------------------------------------------------
## 1. Align the QTL table, current selection gradients and GWAS scales
## ------------------------------------------------------------
qtl <- read.csv(cluster_file, stringsAsFactors = FALSE, check.names = FALSE)
grad <- read.csv(gradient_file, stringsAsFactors = FALSE, check.names = FALSE)
needed_qtl <- c("group", "chr", "ps", "beta_beta_NaCl",
                "p_value_REML.cor_beta_NaCl", "slope_median_common",
                "slope_mean_common", "n_slope_common",
                paste0("beta_", traits, "_NaCl"))
needed_grad <- c("trait", "beta", "SE", "beta_for_mapping")
if (!all(needed_qtl %in% names(qtl)))
  stop("Missing script-11 columns: ",
       paste(setdiff(needed_qtl, names(qtl)), collapse = ", "))
if (!all(needed_grad %in% names(grad)))
  stop("Missing selection-gradient columns: ",
       paste(setdiff(needed_grad, names(grad)), collapse = ", "))
if (anyDuplicated(qtl$group) || anyDuplicated(grad$trait))
  stop("Duplicate QTL or selection-gradient identifiers")
if (!all(traits %in% grad$trait))
  stop("Selection-gradient file does not contain all seven traits")

## Use the exact weight vector that produced the CURRENT mapped composite.
## Its numerical values need not equal coefficients printed in the eLife paper.
grad <- grad[match(traits, grad$trait), , drop = FALSE]
w_map <- setNames(as.numeric(grad$beta_for_mapping), traits)
if (any(!is.finite(w_map))) stop("Invalid current composite-index weights")

## Trait GWAS coefficients are on individually standardized response scales
## (script 04 calls scale(y)). The composite was built on standardized
## replicate traits in script 00 and then standardized once more for mapping.
## Thus a component-GWAS coefficient contributes w_i * SD(X_i) / SD(index)
## to the reconstructed full-index GWAS coefficient.
if (!file.exists(preparation_file)) {
  alternate <- paste0(preparation_file, ".RData")
  if (file.exists(alternate)) preparation_file <- alternate
}
prep <- new.env(parent = globalenv())
load(preparation_file, envir = prep)
if (exists("data_for_GWAS_NaCl_genotyped_NaCl", envir = prep,
           inherits = FALSE)) {
  ph <- prep$data_for_GWAS_NaCl_genotyped_NaCl
} else if (exists("all_BLUEs_genotyped_NaCl", envir = prep,
                  inherits = FALSE)) {
  ph <- prep$all_BLUEs_genotyped_NaCl
} else {
  stop("Preparatory file lacks the mapping phenotype table")
}
if (!all(c(traits, "beta") %in% names(ph)))
  stop("Mapping phenotype table lacks seven traits or beta")
X_ph <- as.matrix(ph[, traits, drop = FALSE])
storage.mode(X_ph) <- "double"

complete_index <- rowSums(is.finite(X_ph)) == length(traits) &
  is.finite(ph$beta)

index_check <- drop(
  X_ph[complete_index, , drop = FALSE] %*% w_map
)

index_agreement <- cor(
  index_check,
  ph$beta[complete_index]
)

index_max_error <- max(
  abs(index_check - ph$beta[complete_index])
)

trait_SD_GWAS <- apply(X_ph, 2, sd, na.rm = TRUE)
index_SD <- sd(ph$beta, na.rm = TRUE)
w_full <- w_map * trait_SD_GWAS / index_SD
w_3map <- w_map
w_3map[!names(w_3map) %in% direct_traits] <- 0
index_3 <- drop(X_ph %*% w_3map)
#index_3_SD <- sd(index_3)
#if (!is.finite(index_3_SD) || index_3_SD <= 0)
#  stop("Three-trait index has no variance")
#w_three <- w_3map * trait_SD_GWAS / index_3_SD
w_three <- w_3map * trait_SD_GWAS / index_SD
expected_sign <- c(SF = -1, SB = 1, area.F = 1)
actual_sign <- sign(w_map[names(expected_sign)])
if (any(actual_sign != expected_sign)) {
  warning("Current pipeline weights differ in sign from the previously published ",
          "directions (SF negative, SB and area.F positive). Check gradient ",
          "provenance before biological interpretation.")
}
weight_qc <- data.frame(
  trait = traits, beta_BLUE = as.numeric(grad$beta),
  beta_SE = as.numeric(grad$SE),
  beta_CI_low = if ("lower95" %in% names(grad)) grad$lower95 else NA_real_,
  beta_CI_high = if ("upper95" %in% names(grad)) grad$upper95 else NA_real_,
  p_gradient = if ("p_value" %in% names(grad)) grad$p_value else NA_real_,
  beta_for_mapping = as.numeric(w_map),
  GWAS_pheno_SD = as.numeric(trait_SD_GWAS),
  weight_full_GWAS_scale = as.numeric(w_full),
  weight_three_GWAS_scale = as.numeric(w_three),
  included_in_three = traits %in% direct_traits,
  stringsAsFactors = FALSE)
write.csv(weight_qc, file.path(output_dir, "13_Gradient_weight_QC.csv"),
          row.names = FALSE)
rm(ph, X_ph, index_check, index_3)
gc(verbose = FALSE)

## ------------------------------------------------------------
## 2. One independent observation PER LD cluster (not per SNP)
## ------------------------------------------------------------
effect_cols <- paste0("beta_", traits, "_NaCl")
B <- as.matrix(qtl[, effect_cols, drop = FALSE])
storage.mode(B) <- "double"
colnames(B) <- traits
qtl$pred_full_GWAS <- as.numeric(qtl$beta_beta_NaCl)
qtl$pred_full_reconstructed <- drop(B %*% w_full)
qtl$pred_three_trait <- drop(B %*% w_three)
qtl$pred_SF <- B[, "SF"]
qtl$pred_SB <- B[, "SB"]
qtl$pred_area.F <- B[, "area.F"]
qtl$response <- as.numeric(qtl$slope_median_common)
qtl$chromosome <- as.character(qtl$chr)
qtl$position <- as.numeric(qtl$ps)
qtl$block <- paste(qtl$chromosome,
                   floor(qtl$position / block_bp), sep = ":")
qtl$beta_GWAS_q <- as.numeric(qtl$p_value_REML.cor_beta_NaCl)
qtl$cohort_beta_q05 <- is.finite(qtl$beta_GWAS_q) & qtl$beta_GWAS_q < .05
qtl$cohort_beta_q01 <- is.finite(qtl$beta_GWAS_q) & qtl$beta_GWAS_q < .01

## A common complete-case set is essential for comparing correlations.
needed_common <- c(effect_cols, "pred_full_GWAS", "pred_full_reconstructed",
                   "pred_three_trait", "response", "position")
valid <- rowSums(!is.finite(as.matrix(qtl[, needed_common,
                                      drop = FALSE]))) == 0 &
  qtl$chromosome %in% chromosomes &
  !is.na(qtl$group) & qtl$n_slope_common > 0
valid[is.na(valid)] <- FALSE
D <- qtl[valid, , drop = FALSE]
if (nrow(D) < 30L) stop("Too few common-frequency QTLs after matched filtering")
if (length(unique(D$chromosome)) < 3L)
  stop("Fewer than three chromosomes are represented in the primary set")
D <- D[order(match(D$chromosome, chromosomes), D$position), , drop = FALSE]
rownames(D) <- NULL
write.csv(D[, c("group", "chromosome", "position", "block", "n_slope_common",
                "beta_GWAS_q", "response", "slope_mean_common",
                "pred_full_GWAS", "pred_full_reconstructed", "pred_three_trait",
                effect_cols), drop = FALSE],
          file.path(output_dir, "13_Cluster_predictors_and_response.csv"),
          row.names = FALSE)

## ------------------------------------------------------------
## 3. Observed associations on matched QTL sets
## ------------------------------------------------------------
predictors <- c(full_mapped = "pred_full_GWAS",
                full_reconstructed = "pred_full_reconstructed",
                three_trait = "pred_three_trait",
                SF = "pred_SF", SB = "pred_SB", area.F = "pred_area.F")
cohorts <- list(all_common = rep(TRUE, nrow(D)),
                beta_GWAS_q05 = D$cohort_beta_q05,
                beta_GWAS_q01 = D$cohort_beta_q01)

safe_cor <- function(x, y, method = "pearson") {
  good <- is.finite(x) & is.finite(y)
  if (sum(good) < 4L || sd(x[good]) == 0 || sd(y[good]) == 0)
    return(NA_real_)
  cor(x[good], y[good], method = method)
}

assoc <- function(d, xcol, cohort, predictor) {
  x <- d[[xcol]]; y <- d$response
  fit <- if (nrow(d) >= 5L && sd(x) > 0 && sd(y) > 0)
    lm(y ~ x) else NULL
  data.frame(cohort = cohort, predictor = predictor,
    n = nrow(d),
    r = safe_cor(x, y),
    spearman_r = safe_cor(x, y, "spearman"),
    R2_in_sample = safe_cor(x, y)^2,
    slope = if (is.null(fit)) NA_real_ else unname(coef(fit)[2]),
    p_lm_descriptive = if (is.null(fit)) NA_real_ else
      unname(coef(summary(fit))[2, "Pr(>|t|)"]),
    sign_agreement = mean(sign(x) == sign(y) & x != 0 & y != 0),
    stringsAsFactors = FALSE)
}
association_list <- list()
for (cohort in names(cohorts)) {
  d <- D[cohorts[[cohort]], , drop = FALSE]
  for (pred in names(predictors)) {
    association_list[[length(association_list) + 1L]] <-
      assoc(d, predictors[[pred]], cohort, pred)
  }
}
association_summary <- do.call(rbind, association_list)
write.csv(association_summary, file.path(output_dir, "13_Association_summary.csv"),
          row.names = FALSE)

## Mean-vs-median aggregation and SNP-coverage sensitivity, keeping the same
## GWAS phenotypes and frequency-filtered cluster universe.
coverage_sensitivity <- do.call(rbind, lapply(c(1L, 2L, 10L, 30L), function(min_snp) {
  dd <- D[D$n_slope_common >= min_snp, , drop = FALSE]
  do.call(rbind, lapply(c("response", "slope_mean_common"), function(yname) {
    data.frame(min_common_SNPs = min_snp, response = yname, n_clusters = nrow(dd),
      r_full_mapped = safe_cor(dd$pred_full_GWAS, dd[[yname]]),
      r_full_reconstructed = safe_cor(dd$pred_full_reconstructed, dd[[yname]]),
      r_three = safe_cor(dd$pred_three_trait, dd[[yname]]),
      stringsAsFactors = FALSE)
  }))
}))
write.csv(coverage_sensitivity, file.path(output_dir,
          "13_Aggregation_coverage_sensitivity.csv"), row.names = FALSE)

full_cor <- safe_cor(D$pred_full_GWAS, D$pred_full_reconstructed)
scale_QC <- data.frame(n_clusters = nrow(D),
  correlation_mapped_vs_reconstructed_full = full_cor,
  index_reconstruction_max_abs_error = index_max_error,
  index_SD = index_SD, three_index_SD = index_3_SD,
  n_genomic_blocks = length(unique(D$block)),
  n_chromosomes = length(unique(D$chromosome)),
  stringsAsFactors = FALSE)
write.csv(scale_QC, file.path(output_dir, "13_Index_reconstruction_QC.csv"),
          row.names = FALSE)
if (is.finite(full_cor) && full_cor < .95)
  warning("Mapped beta GWAS effects and reconstructed full effects correlate <0.95. ",
          "Outcome-specific GridLMM fitting may matter; compare both analyses.")

## ------------------------------------------------------------
## 4. Block bootstrap of effect sizes AND correlated differences
## ------------------------------------------------------------
## A chromosome:2-Mb block is resampled as an intact unit. This reduces,
## but does not eliminate, residual LD or ascertainment dependence.
bootstrap_one <- function(d, cohort, n_boot) {
  split_idx <- split(seq_len(nrow(d)), as.character(d$block))
  n_blocks <- length(split_idx)
  if (n_blocks < 6L) return(NULL)
  observed <- c(full = safe_cor(d$pred_full_GWAS, d$response),
                reconstructed = safe_cor(d$pred_full_reconstructed, d$response),
                three = safe_cor(d$pred_three_trait, d$response))
  storage <- matrix(NA_real_, nrow = n_boot, ncol = 5L,
                    dimnames = list(NULL, c("full", "reconstructed", "three",
                                            "delta_full_minus_three",
                                            "delta_reconstructed_minus_three")))
  for (b in seq_len(n_boot)) {
    chosen <- sample.int(n_blocks, n_blocks, replace = TRUE)
    ii <- unlist(split_idx[chosen], use.names = FALSE)
    dd <- d[ii, , drop = FALSE]
    rr <- c(safe_cor(dd$pred_full_GWAS, dd$response),
            safe_cor(dd$pred_full_reconstructed, dd$response),
            safe_cor(dd$pred_three_trait, dd$response))
    storage[b, ] <- c(rr, rr[1] - rr[3], rr[2] - rr[3])
  }
  est <- c(observed, delta_full_minus_three = observed[1] - observed[3],
           delta_reconstructed_minus_three = observed[2] - observed[3])
  est <- as.numeric(est)
  data.frame(cohort = cohort, metric = colnames(storage),
             n_clusters = nrow(d), n_blocks = n_blocks,
             estimate = est,
             lower95 = apply(storage, 2, function(z)
               quantile(z, .025, na.rm = TRUE, names = FALSE)),
             upper95 = apply(storage, 2, function(z)
               quantile(z, .975, na.rm = TRUE, names = FALSE)),
             frac_bootstrap_gt_zero = colMeans(storage > 0, na.rm = TRUE),
             stringsAsFactors = FALSE)
}
boot_result <- do.call(rbind, Filter(Negate(is.null),
  lapply(c("all_common", "beta_GWAS_q05"), function(cohort)
    bootstrap_one(D[cohorts[[cohort]], , drop = FALSE], cohort, n_boot))))
write.csv(boot_result, file.path(output_dir,
          "13_Genomic_block_bootstrap.csv"), row.names = FALSE)

## ------------------------------------------------------------
## 5. Are the observed gradient directions special?
## ------------------------------------------------------------
## Two descriptive controls:
##   (a) ALL 2^3 sign configurations of magnitudes for SF/SB/area.F
##   (b) Random directions in 3- or 7-dimensional trait-effect space.
## These are WEIGHT-COMPARISON BENCHMARKS, NOT RANDOMIZATION p-values.
B_D <- as.matrix(D[, effect_cols, drop = FALSE]); colnames(B_D) <- traits

sign_configs <- expand.grid(SF = c(-1, 1), SB = c(-1, 1),
                            area.F = c(-1, 1))
sign_result <- do.call(rbind, lapply(seq_len(nrow(sign_configs)), function(k) {
  s <- unlist(sign_configs[k, ], use.names = TRUE)
  ww <- setNames(rep(0, length(traits)), traits)
  ww[direct_traits] <- abs(w_three[direct_traits]) * s[direct_traits]
  score <- drop(B_D %*% ww)
  data.frame(sign_SF = s["SF"], sign_SB = s["SB"],
             sign_area.F = s["area.F"],
             matches_observed_gradient_signs =
               all(s == sign(w_three[direct_traits])),
             pearson_r = safe_cor(score, D$response),
             stringsAsFactors = FALSE)
}))
write.csv(sign_result, file.path(output_dir,
          "13_Three_trait_sign_combinations.csv"), row.names = FALSE)

random_comparison <- list()
random_detail <- list()
## Keep random directions uniformly distributed on the unit sphere. No
## independent shuffling of SNP effects, preserving empirical pleiotropy.
for (space in c("three_trait", "seven_trait")) {
  selected <- if (space == "three_trait") direct_traits else traits
  w_observed <- if (space == "three_trait") w_three else w_full
  obs <- safe_cor(drop(B_D %*% w_observed), D$response)
  random_r <- numeric(n_random_weights)
  for (k in seq_len(n_random_weights)) {
    w_random <- rnorm(length(selected))
    w_random <- w_random / sqrt(sum(w_random^2))
    random_r[k] <- safe_cor(drop(B_D[, selected, drop = FALSE] %*% w_random),
                            D$response)
  }
  random_comparison[[space]] <- data.frame(
    weight_space = space, n_random = n_random_weights,
    observed_r = obs,
    fraction_random_r_ge_observed = mean(random_r >= obs),
    fraction_random_abs_r_ge_observed = mean(abs(random_r) >= abs(obs)),
    median_random_r = median(random_r), q025_random_r = quantile(random_r, .025),
    q975_random_r = quantile(random_r, .975), stringsAsFactors = FALSE)
  random_detail[[space]] <- data.frame(weight_space = space,
                                      draw = seq_along(random_r), r = random_r)
}
random_summary <- do.call(rbind, random_comparison)
write.csv(random_summary, file.path(output_dir,
          "13_Alternative_weight_benchmarks.csv"), row.names = FALSE)
write.csv(do.call(rbind, random_detail), file.path(output_dir,
          "13_Alternative_weight_draws.csv"), row.names = FALSE)

## ------------------------------------------------------------
## 5b. Selection-gradient estimation uncertainty (conditional sensitivity)
## ------------------------------------------------------------
## Use the covariance of the fitted seven-trait phenotypic selection model.
## This propagates only variation in the gradient estimate. It does NOT
## account for GWAS coefficient uncertainty, linkage or evolutionary drift.
## The mapping-scale conversion is beta_BLUE * replicate_sd / BLUE_sd.
## For a correlation the positive composite-SD normalization cancels.
gradient_draw_summary <- NULL
if (exists("model_beta", envir = prep, inherits = FALSE)) {
  gradient_cov <- try(vcov(prep$model_beta)[traits, traits, drop = FALSE],
                      silent = TRUE)
  if (!inherits(gradient_cov, "try-error") &&
      all(is.finite(gradient_cov)) &&
      !inherits(try(chol(gradient_cov), silent = TRUE), "try-error")) {
    gradient_factor <- chol(gradient_cov)
    draws <- matrix(rnorm(n_gradient_draws * length(traits)),
                    ncol = length(traits)) %*% gradient_factor
    draws <- sweep(draws, 2, as.numeric(grad$beta), "+")
    colnames(draws) <- traits
    convert <- as.numeric(grad$replicate_sd / grad$BLUE_sd) *
      as.numeric(trait_SD_GWAS)
    if (length(convert) == length(traits) && all(is.finite(convert))) {
      obs_sets <- list(all_common = rep(TRUE, nrow(D)),
                       beta_GWAS_q05 = D$cohort_beta_q05)
      gradient_draw_summary <- do.call(rbind,
        lapply(names(obs_sets), function(cohort) {
          idx <- which(obs_sets[[cohort]])
          y <- D$response[idx]
          bx <- B_D[idx, , drop = FALSE]
          do.call(rbind, lapply(c("full", "three"), function(model) {
            w <- draws
            if (model == "three") w[, !traits %in% direct_traits] <- 0
            w <- sweep(w, 2, convert, "*")
            score <- bx %*% t(w)
            rr <- apply(score, 2, function(z) safe_cor(z, y))
            data.frame(cohort = cohort, index = model,
              n_gradient_draws = length(rr),
              median_r = median(rr, na.rm = TRUE),
              lower95 = quantile(rr, .025, na.rm = TRUE, names = FALSE),
              upper95 = quantile(rr, .975, na.rm = TRUE, names = FALSE),
              fraction_positive = mean(rr > 0, na.rm = TRUE),
              stringsAsFactors = FALSE)
          }))
        }))
      write.csv(gradient_draw_summary,
        file.path(output_dir, "13_Selection_gradient_uncertainty.csv"),
        row.names = FALSE)
    }
  }
}
if (is.null(gradient_draw_summary))
  warning("Could not propagate selection-gradient covariance; inspect ",
          "the saved model_beta and gradient scale columns.")

## ------------------------------------------------------------
## 6. Chromosome-shift sensitivity (spatially restricted, not exact null)
## ------------------------------------------------------------
chrom_shift <- function(d, xname, B = n_shift) {
  x <- d[[xname]]; y <- d$response
  obs <- safe_cor(x, y)
  if (!is.finite(obs)) return(NA_real_)
  by_chr <- split(seq_len(nrow(d)), factor(d$chromosome,
                                           levels = chromosomes), drop = TRUE)
  extreme <- 0L
  for (b in seq_len(B)) {
    yp <- y
    for (ix in by_chr) {
      n <- length(ix)
      if (n < 2L) next
      k <- sample.int(n, 1L) - 1L
      if (k > 0L) yp[ix] <- y[ix[c((k + 1L):n, seq_len(k))]]
    }
    r <- safe_cor(x, yp)
    if (is.finite(r) && abs(r) >= abs(obs)) extreme <- extreme + 1L
  }
  (extreme + 1) / (B + 1)
}
shift_result <- do.call(rbind, lapply(c("all_common", "beta_GWAS_q05"), function(cohort) {
  d <- D[cohorts[[cohort]], , drop = FALSE]
  do.call(rbind, lapply(c("full_mapped", "full_reconstructed", "three_trait"), function(pred) {
    data.frame(cohort = cohort, predictor = pred, n = nrow(d),
      r = safe_cor(d[[predictors[[pred]]]], d$response),
      p_chrom_shift_exploratory = chrom_shift(d, predictors[[pred]]),
      stringsAsFactors = FALSE)
  }))
}))
write.csv(shift_result, file.path(output_dir,
          "13_Chromosome_shift_sensitivity.csv"), row.names = FALSE)

## ------------------------------------------------------------
## 7. Conditional pleiotropy: seven effects in one joint model
## ------------------------------------------------------------
## Descriptive conditional model; coefficient signs need not equal phenotypic
## selection gradients due to LD, pleiotropy, unmeasured traits and scaling.
Xs <- scale(B_D)
ys <- as.numeric(scale(D$response))
Xs <- as.matrix(Xs)
colnames(Xs) <- traits
joint <- data.frame(y = ys, as.data.frame(Xs, check.names = FALSE),
                    check.names = FALSE)
form <- reformulate(paste0("`", traits, "`"), response = "y")
fit_joint <- lm(form, data = joint)
coeff <- coef(summary(fit_joint))
ols_results <- do.call(rbind, lapply(traits, function(tr) {
  if (!tr %in% rownames(coeff)) return(data.frame(
    trait = tr, beta_conditional = NA_real_, SE = NA_real_,
    p_descriptive = NA_real_, lower95 = NA_real_, upper95 = NA_real_))
  b <- coeff[tr, "Estimate"]; se <- coeff[tr, "Std. Error"]
  data.frame(trait = tr, beta_conditional = b, SE = se,
             p_descriptive = coeff[tr, "Pr(>|t|)"],
             lower95 = b - qt(.975, df.residual(fit_joint)) * se,
             upper95 = b + qt(.975, df.residual(fit_joint)) * se)
}))
ols_results$direct_gradient_expected_sign <- c(-1, 1, NA, NA, NA, NA, 1)
ols_results$gradient_weight_GWAS_scale <- as.numeric(w_full)
ols_results$Pearson_marginal <- vapply(traits, function(tr)
  safe_cor(D[[paste0("beta_", tr, "_NaCl")]], D$response), numeric(1))
write.csv(ols_results, file.path(output_dir,
          "13_Conditional_pleiotropy_OLS.csv"), row.names = FALSE)

VIF <- vapply(seq_along(traits), function(i) {
  r2 <- summary(lm(Xs[, i] ~ Xs[, -i, drop = FALSE]))$r.squared
  if (r2 >= 1) Inf else 1 / (1 - r2)
}, numeric(1))
write.csv(data.frame(trait = traits, VIF = VIF),
          file.path(output_dir, "13_Pleiotropy_VIF.csv"), row.names = FALSE)

## Ridge coefficients in standardized predictor and response units; fixed
## penalties for sensitivity, not tuned using the evaluation chromosomes.
ridge_coef <- function(X, y, lambda) {
  Z <- scale(X)
  zz <- as.numeric(scale(y))
  as.numeric(solve(crossprod(Z) + diag(lambda, ncol(Z)),
                   crossprod(Z, zz)))
}
ridge_result <- do.call(rbind, lapply(c(1, 10, 50), function(lambda) {
  data.frame(trait = traits, lambda = lambda,
             standardized_beta = ridge_coef(B_D, D$response, lambda))
}))
write.csv(ridge_result, file.path(output_dir,
          "13_Conditional_pleiotropy_ridge.csv"), row.names = FALSE)

## ------------------------------------------------------------
## 8. Leave-one-chromosome-out (LOCO) predictive sensitivity
## ------------------------------------------------------------
## Each response-model calibration is trained on 5 chromosomes and evaluated
## on the held-out chromosome. This is NOT independently replicated biology.
predict_simple <- function(train, test, xcol) {
  x <- train[[xcol]]; y <- train$response
  if (sd(x) == 0) return(rep(mean(y), nrow(test)))
  fit <- lm(y ~ x)
  as.numeric(coef(fit)[1] + coef(fit)[2] * test[[xcol]])
}
predict_ridge <- function(Xtr, ytr, Xte, lambda) {
  mu <- colMeans(Xtr)
  sdev <- apply(Xtr, 2, sd)
  sdev[sdev == 0 | !is.finite(sdev)] <- 1
  Ztr <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdev, "/")
  Zte <- sweep(sweep(Xte, 2, mu, "-"), 2, sdev, "/")
  y0 <- mean(ytr)
  b <- solve(crossprod(Ztr) + diag(lambda, ncol(Ztr)),
             crossprod(Ztr, ytr - y0))
  as.numeric(y0 + Zte %*% b)
}

cv_models <- c("intercept", "full_mapped", "full_reconstructed", "three_trait",
               "SF", "SB", "area.F", "joint_7_OLS", "joint_7_ridge10")
CV <- matrix(NA_real_, nrow = nrow(D), ncol = length(cv_models),
             dimnames = list(NULL, cv_models))
for (chr in sort(unique(D$chromosome))) {
  is_test <- D$chromosome == chr
  tr <- D[!is_test, , drop = FALSE]
  te <- D[is_test, , drop = FALSE]
  if (nrow(tr) < 20 || nrow(te) < 2) next
  CV[is_test, "intercept"] <- mean(tr$response)
  for (m in c("full_mapped", "full_reconstructed", "three_trait",
              "SF", "SB", "area.F"))
    CV[is_test, m] <- predict_simple(tr, te, predictors[[m]])
  A <- as.matrix(tr[, effect_cols, drop = FALSE]); Z <- as.matrix(te[, effect_cols, drop = FALSE])
  CV[is_test, "joint_7_OLS"] <- predict_ridge(A, tr$response, Z, 1e-8)
  CV[is_test, "joint_7_ridge10"] <- predict_ridge(A, tr$response, Z, 10)
}
cv_summary <- do.call(rbind, lapply(colnames(CV), function(m) {
  yy <- D$response; pp <- CV[, m]
  ix <- is.finite(yy) & is.finite(pp)
  data.frame(model = m, n = sum(ix),
    RMSE = sqrt(mean((yy[ix] - pp[ix])^2)),
    R2_LOCO = 1 - sum((yy[ix] - pp[ix])^2) / sum((yy[ix] - mean(yy[ix]))^2),
    Pearson_LOCO = safe_cor(yy[ix], pp[ix]),
    stringsAsFactors = FALSE)
}))
write.csv(cv_summary, file.path(output_dir,
          "13_Leave_chromosome_out_prediction.csv"), row.names = FALSE)
write.csv(data.frame(group = D$group, chromosome = D$chromosome,
                     response = D$response, as.data.frame(CV, check.names = FALSE)),
          file.path(output_dir, "13_Leave_chromosome_out_fitted.csv"),
          row.names = FALSE)

leave_chr <- do.call(rbind, lapply(unique(D$chromosome), function(chr) {
  d <- D[D$chromosome != chr, , drop = FALSE]
  data.frame(removed_chromosome = chr, n = nrow(d),
    r_full = safe_cor(d$pred_full_GWAS, d$response),
    r_full_reconstructed = safe_cor(d$pred_full_reconstructed, d$response),
    r_three = safe_cor(d$pred_three_trait, d$response))
}))
write.csv(leave_chr, file.path(output_dir,
          "13_Leave_one_chromosome_correlations.csv"), row.names = FALSE)

rm(prep)
gc(verbose = FALSE)

## ------------------------------------------------------------
## 9. Figures and summary
## ------------------------------------------------------------
pdf(file.path(output_dir, "13_Multivariate_selection_prediction.pdf"),
    width = 12, height = 8)
oldpar <- par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3, 1))
plot(seq_along(traits), grad$beta, pch = 16, cex = 1.2,
     xaxt = "n", xlab = "Trait", ylab = "Multivariate selection gradient",
     main = "A. Gradients used in current GWAS index", bty = "n",
     ylim = range(c(0, grad$beta - 2 * grad$SE,
                    grad$beta + 2 * grad$SE)))
axis(1, at = seq_along(traits), labels = traits)
arrows(seq_along(traits), grad$beta - 1.96 * grad$SE,
       seq_along(traits), grad$beta + 1.96 * grad$SE,
       code = 3, angle = 90, length = .035)
abline(h = 0, lty = 2, col = "grey60")

plot(D$pred_full_GWAS, D$response, pch = 16, cex = .6,
     col = adjustcolor("grey35", alpha.f = .45),
     xlab = "Composite GWAS effect (standardized)",
     ylab = "Phase-aligned frequency slope / generation",
     main = sprintf("B. Composite GWAS vs evolution (r = %.3f)",
                    safe_cor(D$pred_full_GWAS, D$response)), bty = "n")
abline(lm(response ~ pred_full_GWAS, data = D), col = "steelblue4", lwd = 2)
abline(h = 0, v = 0, col = "grey70", lty = 2)

sel_assoc <- subset(association_summary, cohort == "all_common" &
  predictor %in% c("full_mapped", "full_reconstructed", "three_trait", "SF", "SB", "area.F"))
barplot(sel_assoc$r, names.arg = sel_assoc$predictor, las = 2,
        ylab = "Pearson r (same clusters)",
        main = "C. Full vs reduced fitness prediction", bty = "n",
        col = "grey65")
abline(h = 0, lty = 2)

plot(seq_along(traits), ols_results$beta_conditional,
     ylim = range(c(0, ols_results$lower95, ols_results$upper95), na.rm = TRUE),
     pch = 16, xaxt = "n", xlab = "Trait",
     ylab = "Conditional slope (standardized)",
     main = "D. Joint model of seven GWAS effects", bty = "n")
axis(1, at = seq_along(traits), labels = traits)
arrows(seq_along(traits), ols_results$lower95,
       seq_along(traits), ols_results$upper95,
       code = 3, angle = 90, length = .035)
abline(h = 0, lty = 2, col = "grey60")
par(oldpar)
dev.off()

pdf(file.path(output_dir, "13_Weight_direction_benchmarks.pdf"),
    width = 10, height = 4)
oldpar <- par(mfrow = c(1, 2))
for (space in names(random_detail)) {
  hist(random_detail[[space]]$r, breaks = 45, col = "grey75", border = "white",
       xlab = "Correlation with evolutionary response",
       main = paste("Random", space, "weight directions"))
  abline(v = random_summary$observed_r[random_summary$weight_space == space],
         col = "firebrick3", lwd = 2)
}
par(oldpar)
dev.off()

print(scale_QC, row.names = FALSE)
print(subset(association_summary, cohort == "all_common"),
      row.names = FALSE, digits = 3)
cat("\nBlock bootstrap differences (full - three-trait):\n")
print(subset(boot_result, grepl("delta", metric)),
      row.names = FALSE, digits = 3)
cat("\nLeave-one-chromosome-out response prediction:\n")
print(cv_summary, row.names = FALSE, digits = 3)
cat("\nGradient direction QC:\n")
print(weight_qc[, c("trait", "beta_BLUE", "beta_for_mapping",
                    "included_in_three")], row.names = FALSE)
cat("\nSaved outputs in: ", output_dir, "\n", sep = "")
cat("IMPORTANT: spatial resampling, CV and alternative-weight rankings are ",
    "sensitivity analyses, not independent causal or exact selection tests.\n", sep = "")
