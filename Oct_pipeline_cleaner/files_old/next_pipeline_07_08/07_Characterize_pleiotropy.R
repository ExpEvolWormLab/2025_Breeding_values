## 07_Characterize_pleiotropy.R
## Formerly 13_Characterize_pleiotropy_R5_NaCl_v2.R.
## Descriptive architecture and conditional Wagner/Frachon scaling.
## Run after scripts 02, 03 and 05. Uses seven original phenotypes for Neff.
## The composite 'beta' phenotype is NOT counted as an independent trait.

threshold_clust <- 0.5
label <- paste0("R", round(10 * threshold_clust))
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")

out_dir <- file.path("pleiotropy_results", label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## 1. Load clusters and RIL BLUEs from the current pipeline
## ------------------------------------------------------------
clust_env <- new.env(parent = globalenv())
load(file.path("RData_GridLMM_NaCl",
  paste0("GWAS_cluster_results_BLUE_concatenated_", label,
         "_with_residuals.RData")), envir = clust_env)
all_GWAS_clusters <- clust_env$all_GWAS_clusters

prep <- new.env(parent = globalenv())
load("RData/Preparatory_GLMM_NaCl_only_withBeta.RData", envir = prep)
phen_NaCl <- prep$phen_NaCl

snp_clusters <- read.table(file.path("txt_GridLMM_NaCl",
  paste0("all_sign_GWAS_with_clusters_", label, ".txt")),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## A short ID check is important: the SNP table is NOT ordered like cluster table.
if (anyDuplicated(snp_clusters$rs) || anyDuplicated(all_GWAS_clusters$group))
  stop("Duplicate SNP or QTL-group identifiers")
if (!setequal(snp_clusters$group, all_GWAS_clusters$group))
  stop("SNP groups differ from mapped QTL groups")

## ------------------------------------------------------------
## 2. Genomic location (retain original priority center > arm > tip)
## ------------------------------------------------------------
rkmap <- read.table("Caendr_filtering/rockman_kruglya_map_liftOver.txt",
                    header = TRUE, stringsAsFactors = FALSE)

## A SNP can overlap more than one boundary interval; use explicit priority.
location_priority <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  if (any(grepl("center", x, ignore.case = TRUE))) return("center")
  if (any(grepl("arm", x, ignore.case = TRUE))) return("arm")
  if (any(grepl("tip", x, ignore.case = TRUE))) return("tip")
  "unclassified"
}

snp_clusters$chr_type <- NA_character_
for (i in seq_len(nrow(rkmap))) {
  idx <- as.character(snp_clusters$CHROM) == as.character(rkmap$CHROM[i]) &
    as.numeric(snp_clusters$POS) >= as.numeric(rkmap$Start[i]) &
    as.numeric(snp_clusters$POS) <= as.numeric(rkmap$End[i])
  ## Preserve the first assigned interval, rather than overwriting silently.
  snp_clusters$chr_type[idx & is.na(snp_clusters$chr_type)] <- as.character(rkmap$type[i])
}

## Retain the historical tip boundary on V, where the map is incomplete.
idxV <- as.character(snp_clusters$CHROM) == "V" &
  as.numeric(snp_clusters$POS) >= 20341612 &
  is.na(snp_clusters$chr_type)
snp_clusters$chr_type[idxV] <- "right_tip"

rows_by_group <- split(seq_len(nrow(snp_clusters)),
                       as.character(snp_clusters$group))
summary_groups <- do.call(rbind, lapply(names(rows_by_group), function(g) {
  z <- snp_clusters[rows_by_group[[g]], , drop = FALSE]
  types <- unique(z$chr_type[!is.na(z$chr_type)])
  data.frame(
    group = g,
    nb_SNPs = nrow(z),
    clust_length = diff(range(as.numeric(z$POS))),
    chr_type2 = location_priority(types),
    mixed_location = length(types) > 1,
    unclassified_SNPs = sum(is.na(z$chr_type)),
    stringsAsFactors = FALSE
  )
}))

j <- match(as.character(all_GWAS_clusters$group), summary_groups$group)
all_GWAS_clusters$nb_SNPs <- summary_groups$nb_SNPs[j]
all_GWAS_clusters$clust_length <- summary_groups$clust_length[j]
all_GWAS_clusters$chr_type2 <- summary_groups$chr_type2[j]
all_GWAS_clusters$mixed_location <- summary_groups$mixed_location[j]

## ------------------------------------------------------------
## 3. Pleiotropy of seven original traits, and effect magnitudes
## ------------------------------------------------------------
beta_cols <- paste0("beta_", traits, "_NaCl")
sig_cols <- paste0("is_sign_", traits, "_NaCl")
if (!all(c(beta_cols, sig_cols) %in% names(all_GWAS_clusters)))
  stop("Missing seven-trait GWAS effect or significance columns")

B <- as.matrix(all_GWAS_clusters[, beta_cols, drop = FALSE])
storage.mode(B) <- "numeric"
S <- as.matrix(all_GWAS_clusters[, sig_cols, drop = FALSE])
S[is.na(S)] <- FALSE
storage.mode(S) <- "logical"

## Effects are already standardized by phenotypic SD in script 01.
## No additional trait-specific SD scaling is applied.
if (any(!is.finite(B[S])))
  stop("Nonfinite significant effects prevent Manhattan/Euclidean calculation")
B[!S] <- 0
all_GWAS_clusters$Ntraits <- rowSums(S)
all_GWAS_clusters$Manhattan <- rowSums(abs(B))
all_GWAS_clusters$Euclidean <- sqrt(rowSums(B^2))

## Correlation of adjusted RIL means: useful baseline but NOT a full
## measurement-error-corrected multivariate genetic covariance matrix.
P <- phen_NaCl[, traits, drop = FALSE]
P <- P[complete.cases(P), , drop = FALSE]
R_traits <- cor(P)
colnames(R_traits) <- rownames(R_traits) <- traits

neff_fun <- function(affected, R) {
  n <- length(affected)
  if (n == 0L) return(NA_real_)
  if (n == 1L) return(1)
  lambda <- eigen(R[affected, affected, drop = FALSE],
                  symmetric = TRUE, only.values = TRUE)$values
  ## Population variance of eigenvalues, denominator n (not var(lambda)).
  neff <- n - mean((lambda - 1)^2)
  max(1, min(n, neff))
}

all_GWAS_clusters$Neff <- vapply(seq_len(nrow(all_GWAS_clusters)), function(i) {
  neff_fun(traits[S[i, ]], R_traits)
}, numeric(1))

all_GWAS_clusters$sign_overall <- all_GWAS_clusters$Ntraits > 0
qtl <- subset(all_GWAS_clusters, sign_overall)

write.csv(all_GWAS_clusters,
  file.path(out_dir, paste0("All_clusters_characterized_", label, ".csv")),
  row.names = FALSE)
write.csv(qtl, file.path(out_dir, paste0("QTL_pleiotropy_", label, ".csv")),
          row.names = FALSE)
write.csv(summary_groups,
          file.path(out_dir, paste0("Cluster_location_QC_", label, ".csv")),
          row.names = FALSE)
write.csv(R_traits, file.path(out_dir, "RIL_BLUE_trait_correlations.csv"))

## Avoid silently changing the original script 05 object.
save(all_GWAS_clusters, threshold_clust,
  file = file.path(out_dir,
    paste0("GWAS_clusters_with_pleiotropy_", label, ".RData")))

## ------------------------------------------------------------
## 4. Descriptive figures
## ------------------------------------------------------------
pdf(file.path(out_dir, paste0("QTL_description_", label, ".pdf")),
    width = 9, height = 7)
oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(all_GWAS_clusters$nb_SNPs, all_GWAS_clusters$clust_length + 1,
     log = "xy", pch = 16, cex = .55, col = adjustcolor("grey30", .5),
     xlab = "SNPs per QTL", ylab = "Cluster span + 1 bp", bty = "n")
barplot(table(factor(all_GWAS_clusters$chr_type2,
                     levels = c("center", "arm", "tip", "unclassified"))),
        ylab = "QTL clusters", main = "Genomic region", las = 1)
barplot(table(factor(all_GWAS_clusters$Ntraits, levels = 0:length(traits))),
        xlab = "Significant phenotypes", ylab = "QTL clusters",
        main = "Seven-trait pleiotropy", las = 1)
plot(qtl$Neff, qtl$Manhattan, pch = 16, cex = .6,
     col = adjustcolor("grey30", .5),
     xlab = "Neff (uncorrected RIL BLUE correlations)",
     ylab = "Manhattan total effect", bty = "n")
par(oldpar)
dev.off()

## ------------------------------------------------------------
## 5. Conditional synergistic-pleiotropy tests by genomic region
## ------------------------------------------------------------
## Power law: T = a_g * Neff ^ c_g. Log-linear errors are assumed.
## These p-values condition on GWAS significance and on estimated Neff.
## They are exploratory until ascertainment and correlation uncertainty are addressed.

test_scaling <- function(dat, response, c0) {
  d <- subset(dat, chr_type2 %in% c("arm", "center") &
                    Ntraits > 0 & Neff > 0 &
                    is.finite(dat[[response]]) & dat[[response]] > 0)
  d$region <- factor(d$chr_type2, levels = c("arm", "center"))
  d$logN <- log(d$Neff)
  d$z <- log(d[[response]]) - c0 * d$logN

  if (any(table(d$region) < 4) ||
      any(tapply(d$Neff, d$region, function(z) length(unique(z))) < 2)) {
    warning("Not enough Neff variation for ", response)
    return(NULL)
  }

  model_null <- lm(z ~ region, data = d)
  model_common <- lm(z ~ region + logN, data = d)
  model_full <- lm(z ~ region * logN, data = d)
  b <- coef(model_full)
  V <- vcov(model_full)
  df <- df.residual(model_full)

  contrast <- function(L) {
    est <- sum(L * b)
    se <- sqrt(drop(t(L) %*% V %*% L))
    c(exponent = c0 + est,
      SE = se,
      lower95 = c0 + est - qt(.95, df) * se,
      p_one_sided = pt(est / se, df, lower.tail = FALSE))
  }
  Larm <- setNames(rep(0, length(b)), names(b))
  Larm["logN"] <- 1
  Lcenter <- Larm
  Lcenter["regioncenter:logN"] <- 1
  tests <- rbind(arm = contrast(Larm), center = contrast(Lcenter))
  data.frame_tests <- data.frame(distance = response, region = rownames(tests),
                                 tests, row.names = NULL)
  difference <- data.frame(distance = response,
    center_minus_arm = unname(b["regioncenter:logN"]),
    p_interaction = anova(model_common, model_full)[["Pr(>F)"]][2])
  list(data = d, null = model_null, common = model_common,
       full = model_full, tests = data.frame_tests, difference = difference)
}

M <- test_scaling(qtl, "Manhattan", 1)
E <- test_scaling(qtl, "Euclidean", 0.5)
models <- Filter(Negate(is.null), list(M, E))
if (length(models)) {
  slopes <- do.call(rbind, lapply(models, `[[`, "tests"))
  slopes$p_Holm <- p.adjust(slopes$p_one_sided, "holm")
  differences <- do.call(rbind, lapply(models, `[[`, "difference"))
  differences$p_Holm <- p.adjust(differences$p_interaction, "holm")
  write.csv(slopes, file.path(out_dir, "Scaling_exponent_tests.csv"),
            row.names = FALSE)
  write.csv(differences, file.path(out_dir, "Scaling_region_interactions.csv"),
            row.names = FALSE)
  print(slopes)
  print(differences)

  pdf(file.path(out_dir, "Pleiotropy_scaling.pdf"), width = 9, height = 4.5)
  oldpar <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  cols <- c(arm = "forestgreen", center = "deepskyblue3")
  for (response in c("Manhattan", "Euclidean")) {
    fit <- if (response == "Manhattan") M else E
    c0 <- if (response == "Manhattan") 1 else 0.5
    if (is.null(fit)) {
      plot.new(); title(main = paste(response, "(not estimable)"))
      next
    }
    d <- fit$data
    plot(d$Neff, d[[response]], pch = 21, cex = .8,
         bg = adjustcolor(cols[as.character(d$region)], .55),
         xlab = "Neff", ylab = paste(response, "effect"),
         main = paste(response, ": null c =", c0), las = 1, bty = "n")
    xx <- seq(min(d$Neff), max(d$Neff), length.out = 150)
    for (g in levels(d$region)) {
      nd <- data.frame(logN = log(xx),
        region = factor(rep(g, length(xx)), levels = levels(d$region)))
      yy <- exp(predict(fit$full, newdata = nd) + c0 * log(xx))
      yy0 <- exp(predict(fit$null, newdata = nd) + c0 * log(xx))
      lines(xx, yy, lwd = 3, col = cols[g])
      lines(xx, yy0, lwd = 2, col = cols[g], lty = 2)
    }
    legend("topleft", legend = c("arm", "center", "no synergy"),
           col = c(cols, "grey30"), lty = c(1, 1, 2),
           lwd = c(3, 3, 2), bty = "n", cex = .8)
  }
  par(oldpar)
  dev.off()
}

capture.output(sessionInfo(), file = file.path(out_dir, "07_sessionInfo.txt"))
cat("QTL characterization complete:", out_dir, "\n")
