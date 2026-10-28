## Check whether VCFtools D describes the allele phase of GridLMM SNP dosages.
## Run AFTER section 2 of 02_Clustering_from_LD.R (after ld_snps is created),
## and before changing the sign-splitting rule.

prep <- new.env(parent = globalenv())
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData', envir = prep)
X <- prep$genotype_matrix  # SAME SNP dosage matrix used by GridLMM

n_pairs_to_check <- 2000L
set.seed(123)

if (!'D' %in% names(ld_snps))
  stop('The VCFtools LD file has no D column')

E <- ld_snps[is.finite(ld_snps$R.2) & ld_snps$R.2 > 0.5 &
               ld_snps$ID1 %in% colnames(X) &
               ld_snps$ID2 %in% colnames(X), , drop = FALSE]
E$D <- as.numeric(as.character(E$D))
E <- E[is.finite(E$D) & E$D != 0, , drop = FALSE]
if (!nrow(E)) stop('No overlapping LD edges with finite nonzero D')

E <- E[sample.int(nrow(E), min(n_pairs_to_check, nrow(E))), , drop = FALSE]
E$r_genotype <- vapply(seq_len(nrow(E)), function(k) {
  a <- X[, E$ID1[k]]
  b <- X[, E$ID2[k]]
  suppressWarnings(cor(a, b, use = 'pairwise.complete.obs'))
}, numeric(1))
E <- E[is.finite(E$r_genotype), , drop = FALSE]

cat('LD pairs evaluated:', nrow(E), '\n')
cat('Agreement: sign(D) versus signed genotype correlation:',
    mean(sign(E$D) == sign(E$r_genotype)), '\n')
cat('Median |VCF haplotype R2 - squared genotype correlation|:',
    median(abs(E$R.2 - E$r_genotype^2)), '\n')

# Compare opposite raw beta signs with phase-aware incompatibility.
# This is diagnostic, NOT a change to QTL clustering.
beta_diagnostic <- do.call(rbind, lapply(seq_along(vect_P_traits), function(j) {
  a <- match(E$ID1, all_GWAS_sign$rs)
  b <- match(E$ID2, all_GWAS_sign$rs)
  valid <- !is.na(a) & !is.na(b)
  p1 <- all_GWAS_sign[[paste0('beta_T', j)]][a[valid]]
  p2 <- all_GWAS_sign[[paste0('beta_T', j)]][b[valid]]
  s1 <- all_GWAS_sign[[paste0('is_sign_T', j)]][a[valid]]
  s2 <- all_GWAS_sign[[paste0('is_sign_T', j)]][b[valid]]
  r  <- E$r_genotype[valid]
  both <- (s1 %in% TRUE) & (s2 %in% TRUE) & is.finite(p1) & is.finite(p2)
  raw_opposite <- both & (sign(p1) != sign(p2))
  phase_incompatible <- both & (sign(p1) != sign(r) * sign(p2))
  data.frame(trait = vect_P_traits[j],
             n_pairs_both_significant = sum(both),
             opposite_raw = sum(raw_opposite),
             opposite_after_phase = sum(phase_incompatible),
             rescued_by_phase = sum(raw_opposite & !phase_incompatible))
}))
print(beta_diagnostic, row.names = FALSE)

write.csv(beta_diagnostic, 'txt_GridLMM_NaCl/LD_phase_beta_diagnostic.csv',
          row.names = FALSE)
