## 11_GWAS_selection_cluster_associations.R
## Refactors central cluster association in old 19.
## Critical: s and allele frequencies are phase-aligned to cluster genotype
## orientation, just like script 03's averaged SNP dosages.
## Association is descriptive: SNPs within LD clusters are not independent.

threshold_clust <- 0.5
label <- paste0("R", round(10 * threshold_clust))
out_dir <- file.path("selection_results", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sel <- read.csv("selection_results/09_Selection_SNP_results.csv",
                stringsAsFactors = FALSE)
phase <- read.csv(file.path("data_BLUE_GridLMM_NaCl",
  paste0("SNP_phase_clusters_", label, ".csv")), stringsAsFactors = FALSE)
if (anyDuplicated(phase$rs) || anyDuplicated(sel$ID))
  stop("Non-unique SNP ID in phase or selection file")

## Here s is on GWAS ALT allele log-odds, from script 09.
## A phase=-1 SNP has cluster-oriented genotype = 2 - ALT dosage;
## therefore its cluster-oriented frequency is 1 - pALT and slope = -sALT.
i <- match(phase$rs, sel$ID)
phase$s_ALT <- sel$s[i]
phase$s_cluster <- phase$phase * phase$s_ALT
phase$p_quasi_BH <- sel$p_quasi_BH[i]
phase$has_poolseq <- !is.na(i)

## Match cluster association metadata explicitly.
cluster_env <- new.env(parent = globalenv())
load(file.path("RData_GridLMM_NaCl",
  paste0("GWAS_cluster_results_BLUE_concatenated_", label,
         "_with_residuals.RData")), envir = cluster_env)
qtl <- cluster_env$all_GWAS_clusters

split_rows <- split(seq_len(nrow(phase)), as.character(phase$group))
qc <- do.call(rbind, lapply(names(split_rows), function(group) {
  z <- phase[split_rows[[group]], , drop = FALSE]
  s <- z$s_cluster[is.finite(z$s_cluster)]
  q <- z$p_quasi_BH[is.finite(z$p_quasi_BH)]
  data.frame(group = as.character(group), n_SNPs = nrow(z),
             n_poolseq = sum(z$has_poolseq), n_s_estimated = length(s),
             s_mean = if (length(s)) mean(s) else NA_real_,
             s_median = if (length(s)) median(s) else NA_real_,
             s_sd = if (length(s) > 1) sd(s) else NA_real_,
             min_SNP_q = if (length(q)) min(q) else NA_real_,
             n_selected_SNPs_q01 = sum(q < .01),
             stringsAsFactors = FALSE)
}))

j <- match(as.character(qtl$group), as.character(qc$group))
qtl <- cbind(qtl, qc[j, setdiff(names(qc), "group"), drop = FALSE])
qtl$selected_any_SNP <- !is.na(qtl$n_selected_SNPs_q01) &
                         qtl$n_selected_SNPs_q01 > 0
write.csv(qtl, file.path(out_dir, "11_Clusters_selection_and_GWAS.csv"),
          row.names = FALSE)
write.csv(phase, file.path(out_dir, "11_SNP_selection_phase_aligned.csv"),
          row.names = FALSE)

## Descriptive relationships, avoiding pseudoreplication of SNPs within QTLs
trait_cols <- paste0("beta_", c("SF", "SB", "FS", "FB", "BS", "BF",
                              "area.F", "beta"), "_NaCl")
trait_cols <- trait_cols[trait_cols %in% names(qtl)]
model_summary <- do.call(rbind, lapply(trait_cols, function(col) {
  d <- qtl[is.finite(qtl[[col]]) & is.finite(qtl$s_median), , drop = FALSE]
  if (nrow(d) < 5 || sd(d[[col]]) == 0 || sd(d$s_median) == 0)
    return(data.frame(trait = col, n = nrow(d), slope = NA_real_,
                      SE = NA_real_, p = NA_real_, r = NA_real_))
  mod <- lm(d$s_median ~ d[[col]])
  co <- coef(summary(mod))[2, ]
  data.frame(trait = col, n = nrow(d), slope = co["Estimate"],
             SE = co["Std. Error"], p = co["Pr(>|t|)"],
             r = cor(d[[col]], d$s_median), row.names = NULL)
}))
model_summary$p_Holm <- p.adjust(model_summary$p, "holm")
write.csv(model_summary, file.path(out_dir, "11_Selection_GWAS_trait_associations.csv"),
          row.names = FALSE)

pdf(file.path(out_dir, "11_Selection_GWAS_trait_associations.pdf"),
    width = 10, height = 6)
oldpar <- par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (col in trait_cols) {
  plot(qtl[[col]], qtl$s_median, pch = 21, cex = .7,
       bg = adjustcolor(ifelse(qtl$selected_any_SNP,
                    "firebrick3", "grey50"), alpha.f = .5),
       xlab = paste("GWAS effect:", col),
       ylab = "Median selection log-odds slope",
       main = col, bty = "n")
  abline(h = 0, v = 0, lty = 2, col = "grey50")
}
par(oldpar); dev.off()
cat("Finished cluster selection/GWAS description in", out_dir, "\n")
