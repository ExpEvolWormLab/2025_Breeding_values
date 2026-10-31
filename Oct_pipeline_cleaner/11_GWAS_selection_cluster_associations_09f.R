## 11_GWAS_selection_cluster_associations_09f.R
## Compare phase-aligned, cluster-level GWAS effects with continuous allele-
## frequency slopes from the drift-aware model 09f (frequency units/generation).
## Primary descriptive association uses common ancestral SNPs [0.1, 0.9].
## Does not treat linked SNPs as independent selection events.
## Regression p-values and chromosome-shift p-values are EXPLORATORY;
## the shift null approximately preserves within-chromosome spatial structure.

threshold_clust <- 0.5
label <- paste0("R", round(10 * threshold_clust))
Ne <- 1000L
q_cutoff <- .01
out_dir <- file.path("selection_results", "drift_oneG0_nullfit", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

selection_file <- file.path("selection_results", "drift_oneG0_nullfit",
                             paste0("09f_SNP_drift_Ne", Ne, ".csv"))
phase_file <- file.path("data_BLUE_GridLMM_NaCl",
                         paste0("SNP_phase_clusters_", label, ".csv"))
cluster_file <- file.path("RData_GridLMM_NaCl",
                           paste0("GWAS_cluster_results_BLUE_concatenated_",
                                  label, "_with_residuals.RData"))
sel <- read.csv(selection_file, stringsAsFactors = FALSE)
phase <- read.csv(phase_file, stringsAsFactors = FALSE)
if (!all(c("rs", "group", "phase") %in% names(phase)) ||
    !all(c("ID", "slope_freq", "common_G0", "p0_hat",
           "p_drift", "p_drift_BH_common") %in% names(sel)))
  stop("Missing phase or drift-selection fields")
if (anyDuplicated(phase$rs) || anyDuplicated(sel$ID))
  stop("Duplicate SNPs prevent safe phase-to-allele matching")
if (anyNA(phase$phase) || any(!phase$phase %in% c(-1, 1)))
  stop("Phase must be precisely +1 or -1 for every clustered SNP")

## Genotype in script 03: phase +1 -> ALT dosage, -1 -> 2-ALT dosage.
## Thus d(f_cluster)/d(generation) = phase * d(f_ALT)/d(generation).
## NOTE: 09f slopes are ALLELE FREQUENCY per generation, not log odds or s.
i <- match(as.character(phase$rs), as.character(sel$ID))
phase$poolseq_covered <- !is.na(i)
phase$selection_tested <- is.finite(sel$p_drift[i])
phase$common_G0 <- !is.na(sel$common_G0[i]) & sel$common_G0[i]
phase$tested_common <- phase$selection_tested & phase$common_G0
phase$p0_ALT <- sel$p0_hat[i]
phase$p0_cluster <- ifelse(phase$phase == 1, phase$p0_ALT,
                           1 - phase$p0_ALT)
phase$slope_freq_ALT <- sel$slope_freq[i]
phase$slope_freq_cluster <- phase$phase * phase$slope_freq_ALT
phase$slope_freq_cluster_common <- ifelse(phase$tested_common,
                                           phase$slope_freq_cluster, NA_real_)
phase$p_drift <- sel$p_drift[i]
phase$q_common <- sel$p_drift_BH_common[i]
phase$selected_common_q01 <- phase$tested_common &
  is.finite(phase$q_common) & phase$q_common < q_cutoff

cluster_env <- new.env(parent = globalenv())
load(cluster_file, envir = cluster_env)
qtl <- cluster_env$all_GWAS_clusters
if (!all(c("group", "chr", "ps") %in% names(qtl)) ||
    anyDuplicated(as.character(qtl$group)))
  stop("Unexpected GWAS cluster metadata or duplicate cluster groups")

safemed <- function(x) { x <- x[is.finite(x)]; if (length(x)) median(x) else NA_real_ }
safemean <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
safesd <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1) sd(x) else NA_real_ }
safemin <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x) else NA_real_ }

cluster_one <- function(z) {
  y <- z$slope_freq_cluster_common[is.finite(z$slope_freq_cluster_common)]
  all <- z$slope_freq_cluster[is.finite(z$slope_freq_cluster)]
  m <- if (length(y)) median(y) else NA_real_
  agreement <- if (length(y) >= 2L && is.finite(m) && m != 0)
    mean(sign(y) == sign(m)) else NA_real_
  data.frame(group = as.character(z$group[1]),
    n_SNPs = nrow(z), n_poolseq_covered = sum(z$poolseq_covered),
    n_tested_common = sum(z$tested_common), n_slope_all = length(all),
    n_slope_common = length(y),
    slope_median_common = m,
    slope_mean_common = safemean(y),
    slope_SD_common = safesd(y),
    slope_median_all_sensitivity = safemed(all),
    common_sign_agreement = agreement,
    min_q_common = safemin(z$q_common),
    n_selected_common_q01 = sum(z$selected_common_q01),
    any_selected_common_q01 = any(z$selected_common_q01),
    stringsAsFactors = FALSE)
}
qc <- do.call(rbind, lapply(split(phase, as.character(phase$group)), cluster_one))
rownames(qc) <- NULL
j <- match(as.character(qtl$group), qc$group)
if (anyNA(j) || length(j) != nrow(qc))
  stop("Phase-group and mapped GWAS-group universes do not agree")
qtl <- cbind(qtl, qc[j, setdiff(names(qc), "group"), drop = FALSE])

write.csv(phase, file.path(out_dir, "11_SNP_phase_aligned_09f.csv"),
          row.names = FALSE)
write.csv(qtl, file.path(out_dir, "11_Clusters_selection_and_GWAS_09f.csv"),
          row.names = FALSE)

## Exploratory spatially restricted null: rotate the ordered cluster-level
## response within each chromosome, preserving the ordering and autocorrelation
## of each chromosome's response vector. This is NOT an exact GWAS null.
## A random chromosome-wide circular shift does not account for all
## ascertainment or chromosome-arm effects; treat p-values as sensitivity only.
set.seed(20261028)
n_shift <- 2000L
shift_p <- function(x, y, chr, pos, B = n_shift) {
  if (length(x) < 10L || sd(x) == 0 || sd(y) == 0 ||
      length(unique(chr)) < 2L) return(NA_real_)
  ord <- order(match(chr, c("I", "II", "III", "IV", "V", "X")), pos)
  x <- x[ord]; y <- y[ord]; chr <- chr[ord]
  idx_by_chr <- split(seq_along(chr), chr)
  obs <- abs(cor(x, y))
  extreme <- 0L
  for (b in seq_len(B)) {
    yy <- y
    for (ix in idx_by_chr) {
      n <- length(ix)
      if (n < 2L) next
      k <- sample.int(n, 1L) - 1L
      if (k > 0L) yy[ix] <- y[ix[c((k + 1L):n, seq_len(k))]]
    }
    if (abs(cor(x, yy)) >= obs) extreme <- extreme + 1L
  }
  (extreme + 1) / (B + 1)
}

traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
trait_cols <- paste0("beta_", traits, "_NaCl")
trait_cols <- trait_cols[trait_cols %in% names(qtl)]
if (!length(trait_cols)) stop("No standardized GWAS-effect columns found")

association_one <- function(col, y_col) {
  keep <- is.finite(qtl[[col]]) & is.finite(qtl[[y_col]]) &
    is.finite(as.numeric(qtl$ps)) & !is.na(qtl$chr)
  d <- qtl[keep, , drop = FALSE]
  if (nrow(d) < 10L || sd(d[[col]]) == 0 || sd(d[[y_col]]) == 0) {
    return(data.frame(trait = col, response = y_col, n_clusters = nrow(d),
      slope_lm = NA_real_, SE_lm = NA_real_, p_lm_descriptive = NA_real_,
      Pearson_r = NA_real_, Spearman_r = NA_real_,
      p_chrom_shift_exploratory = NA_real_))
  }
  mod <- lm(d[[y_col]] ~ d[[col]])
  co <- coef(summary(mod))[2, ]
  data.frame(trait = col, response = y_col, n_clusters = nrow(d),
    slope_lm = unname(co["Estimate"]),
    SE_lm = unname(co["Std. Error"]),
    p_lm_descriptive = unname(co["Pr(>|t|)"]),
    Pearson_r = cor(d[[col]], d[[y_col]]),
    Spearman_r = cor(d[[col]], d[[y_col]], method = "spearman"),
    p_chrom_shift_exploratory = shift_p(d[[col]], d[[y_col]],
                                         d$chr, as.numeric(d$ps)),
    stringsAsFactors = FALSE)
}

primary <- do.call(rbind, lapply(trait_cols, association_one,
                                 y_col = "slope_median_common"))
primary$p_lm_Holm_descriptive <- p.adjust(primary$p_lm_descriptive, "holm")
primary$p_chrom_shift_Holm_exploratory <-
  p.adjust(primary$p_chrom_shift_exploratory, "holm")
all_sensitivity <- do.call(rbind, lapply(trait_cols, association_one,
                            y_col = "slope_median_all_sensitivity"))
all_sensitivity$p_lm_Holm_descriptive <-
  p.adjust(all_sensitivity$p_lm_descriptive, "holm")
all_sensitivity$p_chrom_shift_Holm_exploratory <-
  p.adjust(all_sensitivity$p_chrom_shift_exploratory, "holm")

write.csv(primary,
  file.path(out_dir, "11_GWAS_selection_associations_COMMON_exploratory.csv"),
  row.names = FALSE)
write.csv(all_sensitivity,
  file.path(out_dir, "11_GWAS_selection_associations_ALL_sensitivity.csv"),
  row.names = FALSE)

pdf(file.path(out_dir, "11_GWAS_vs_drift_frequency_slopes_COMMON.pdf"),
    width = 12, height = 6)
oldpar <- par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (col in trait_cols) {
  selected <- qtl$any_selected_common_q01
  selected[is.na(selected)] <- FALSE
  cols <- adjustcolor(ifelse(selected, "firebrick3", "grey40"), alpha.f = .45)
  plot(qtl[[col]], qtl$slope_median_common, pch = 16, cex = .65,
       col = cols, xlab = paste("GWAS effect:", sub("^beta_|_NaCl$", "", col)),
       ylab = "Cluster-oriented frequency slope/generation",
       main = col, bty = "n")
  abline(h = 0, v = 0, col = "grey60", lty = 2)
}
par(oldpar)
dev.off()

qc_overall <- data.frame(
  LD_threshold = threshold_clust, Ne = Ne,
  total_clusters = nrow(qtl),
  clusters_with_common_frequency_slopes = sum(qtl$n_slope_common > 0),
  clusters_with_any_frequency_slope = sum(qtl$n_slope_all > 0),
  clusters_with_common_q01_SNP = sum(qtl$any_selected_common_q01),
  n_shift = n_shift,
  stringsAsFactors = FALSE)
write.csv(qc_overall, file.path(out_dir, "11_Cluster_selection_QC.csv"),
          row.names = FALSE)
print(qc_overall, row.names = FALSE)
print(primary[, c("trait", "n_clusters", "Pearson_r", "Spearman_r",
                  "p_chrom_shift_exploratory")], row.names = FALSE)
cat("Finished. Cluster correlations and association p-values are exploratory,",
    "not independent-locus causal evidence.\n")
