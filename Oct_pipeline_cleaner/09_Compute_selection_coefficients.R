## 09_Compute_selection_coefficients.R
## Formerly 16_Compute_s.R. Reads output from script 08.
## Alleles are explicitly harmonized to GWAS ALT dosage before estimation.
## s_logit = slope in log-odds / generation; NOT automatically a classical s.
## Binomial results retained for comparison; quasi-binomial p-values are a
## more cautious screen but do not fully model drift or temporal correlation.

library(parallel)

pool_file <- "RData_GridLMM_NaCl/08_PoolSeq_prepared.RData"
snp_file <- "RData/Genotype_info_Cemee.RData"
out_dir <- "selection_results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mc_cores <- 8L
chunk_size <- 500L
also_fit_within_replicate <- TRUE

p_env <- new.env(parent = globalenv())
load(pool_file, envir = p_env)
g_env <- new.env(parent = globalenv())
load(snp_file, envir = g_env)
SNP <- as.data.frame(g_env$snp_positions)

## ------------------------------------------------------------
## 1. Match alleles by CHROM.POS and the exact REF/ALT alleles
## ------------------------------------------------------------
if (!all(c("CHROM", "POS", "ID") %in% names(SNP)))
  stop("Missing GWAS SNP CHROM, POS, ID metadata")
ref_col <- grep("^REF$|^ref[._]?GWAS$", names(SNP), ignore.case = TRUE)
alt_col <- grep("^ALT$|^alt[._]?GWAS$", names(SNP), ignore.case = TRUE)
if (length(ref_col) != 1 || length(alt_col) != 1) {
  ## This was the layout of the original Genotype_info_Cemee.RData
  if (ncol(SNP) < 5) stop("Cannot identify reference/alternative allele columns")
  ref_col <- 4L; alt_col <- 5L
  message("Using SNP metadata columns 4 and 5 as REF and ALT; verify their names.")
}
SNP$ref.GWAS <- toupper(as.character(SNP[[ref_col]]))
SNP$alt.GWAS <- toupper(as.character(SNP[[alt_col]]))
SNP$ID_sync <- paste(as.character(SNP$CHROM), as.numeric(as.character(SNP$POS)), sep = ".")
if (anyDuplicated(SNP$ID_sync))
  stop("GWAS SNP table has repeated chromosome-position identifiers")

meta <- p_env$allele_metadata
k <- match(meta$ID, SNP$ID_sync)
meta$ref.GWAS <- SNP$ref.GWAS[k]
meta$alt.GWAS <- SNP$alt.GWAS[k]
minor <- toupper(as.character(meta$minor.SYNC))
major <- toupper(as.character(meta$major.SYNC))
ref <- toupper(as.character(meta$ref.GWAS))
alt <- toupper(as.character(meta$alt.GWAS))
match_same <- minor == alt & major == ref
match_reverse <- minor == ref & major == alt
match_same[is.na(match_same)] <- FALSE
match_reverse[is.na(match_reverse)] <- FALSE
meta$orientation <- ifelse(match_same, "minor_is_alt",
                           ifelse(match_reverse, "minor_is_ref", "incompatible"))
write.csv(meta, file.path(out_dir, "09_allele_orientation_QC.csv"),
          row.names = FALSE)
cat("GWAS ALT allele matching:", table(meta$orientation), "\n")

keep <- meta$orientation != "incompatible"
if (!any(keep)) stop("No compatible GWAS/Pool-seq alleles")
meta <- meta[keep, , drop = FALSE]
minor_count <- p_env$all.counts.caendr[keep, , drop = FALSE]
depth <- p_env$all.cov.caendr[keep, , drop = FALSE]
if (!identical(rownames(minor_count), meta$ID) ||
    !identical(colnames(minor_count), p_env$populations$sample_id))
  stop("Minor-allele read counts differ in order from SNP/sample metadata")

## All k counts now correspond to the GWAS ALT allele used in GridLMM.
k_alt <- minor_count
rev <- meta$orientation == "minor_is_ref"
k_alt[rev, ] <- depth[rev, , drop = FALSE] - minor_count[rev, , drop = FALSE]

pop <- p_env$populations
pop$generations <- as.numeric(pop$generations)
pop$replicate_fact <- factor(pop$replicate)
X_global <- model.matrix(~ generations, data = pop)
X_within <- model.matrix(~ generations + replicate_fact, data = pop)

## ------------------------------------------------------------
## 2. Fit log-odds trajectories and flag overdispersion
## ------------------------------------------------------------
fit_one <- function(i) {
  successes <- as.numeric(k_alt[i, ])
  coverage <- as.numeric(depth[i, ])
  good <- is.finite(successes) & is.finite(coverage) &
          coverage > 0 & successes >= 0 & successes <= coverage
  missing_result <- c(P = NA_real_, s_logit = NA_real_,
     p_quasi = NA_real_, dispersion = NA_real_,
     P_within = NA_real_, s_within = NA_real_,
     p_within_quasi = NA_real_, n_pools = sum(good))
  if (sum(good) < 6L || length(unique(pop$generations[good])) < 2L)
    return(missing_result)
  if (sum(successes[good]) == 0 ||
      sum(coverage[good] - successes[good]) == 0)
    return(missing_result)

  fit_matrix <- function(X) {
    X <- X[good, , drop = FALSE]
    if (qr(X)$rank < ncol(X) || nrow(X) <= ncol(X) + 1L)
      return(rep(NA_real_, 4))
    fit <- suppressWarnings(try(stats::glm.fit(
      x = X, y = successes[good] / coverage[good],
      weights = coverage[good], family = binomial(),
      control = glm.control(maxit = 50)), silent = TRUE))
    if (inherits(fit, "try-error") || !isTRUE(fit$converged) ||
        any(!is.finite(fit$coefficients))) return(rep(NA_real_, 4))
    j <- match("generations", colnames(X))
    s <- unname(fit$coefficients[j])
    if (!is.finite(s) || abs(s) > 3) return(rep(NA_real_, 4))
    p_fit <- pmin(pmax(fit$fitted.values, 1e-10), 1 - 1e-10)
    df <- fit$df.residual
    phi <- max(1, sum(coverage[good] *
       (successes[good] / coverage[good] - p_fit)^2 /
       (p_fit * (1 - p_fit))) / df)
    Xw <- X * sqrt(fit$weights)
    V <- try(solve(crossprod(Xw)), silent = TRUE)
    if (inherits(V, "try-error") || !is.finite(V[j, j]) || V[j, j] <= 0)
      return(rep(NA_real_, 4))
    z <- s / sqrt(V[j, j])
    c(s = s,
      P = 2 * pnorm(abs(z), lower.tail = FALSE),
      quasi = 2 * pt(abs(z / sqrt(phi)), df, lower.tail = FALSE),
      phi = phi)
  }
  a <- fit_matrix(X_global)
  if (also_fit_within_replicate) b <- fit_matrix(X_within) else
    b <- rep(NA_real_, 4)
  c(P = unname(a["P"]), s_logit = unname(a["s"]),
    p_quasi = unname(a["quasi"]), dispersion = unname(a["phi"]),
    P_within = unname(b["P"]), s_within = unname(b["s"]),
    p_within_quasi = unname(b["quasi"]), n_pools = sum(good))
}

## Parallelize across chunks, not one operating-system task per SNP.
idx <- seq_len(nrow(k_alt))
chunks <- split(idx, ceiling(idx / chunk_size))
cat("Fitting", length(idx), "SNPs in", length(chunks), "chunks\n")
if (.Platform$OS.type == "windows") mc_cores <- 1L
res <- parallel::mclapply(chunks, function(ii)
  lapply(ii, fit_one), mc.cores = mc_cores, mc.preschedule = TRUE)
res <- do.call(rbind, unlist(res, recursive = FALSE))

results <- data.frame(
  rs = meta$ID,
  ID = paste(meta$CHROM, meta$POS, sep = "_"),
  CHROM = meta$CHROM, POS = meta$POS,
  ref.GWAS = meta$ref.GWAS, alt.GWAS = meta$alt.GWAS,
  orientation = meta$orientation,
  P = as.numeric(res[, "P"]),
  s = as.numeric(res[, "s_logit"]),
  p_quasi = as.numeric(res[, "p_quasi"]),
  dispersion = as.numeric(res[, "dispersion"]),
  P_within = as.numeric(res[, "P_within"]),
  s_within = as.numeric(res[, "s_within"]),
  p_within_quasi = as.numeric(res[, "p_within_quasi"]),
  n_pools = as.integer(res[, "n_pools"])
)
for (v in c("P", "p_quasi", "P_within", "p_within_quasi")) {
  out <- paste0(v, "_BH")
  results[[out]] <- p.adjust(results[[v]], method = "BH")
}
write.csv(results, file.path(out_dir, "09_Selection_SNP_results.csv"),
          row.names = FALSE)
write.csv(data.frame(
  n_input = nrow(p_env$allele_metadata),
  n_matched = nrow(results), n_tested = sum(is.finite(results$P)),
  n_quasi_q01 = sum(results$p_quasi_BH < .01, na.rm = TRUE),
  n_binomial_q01 = sum(results$P_BH < .01, na.rm = TRUE)),
  file.path(out_dir, "09_selection_QC.csv"), row.names = FALSE)
save(results, pop, meta, file = file.path(out_dir, "09_Selection_results.RData"),
     compress = "gzip")
capture.output(sessionInfo(), file = file.path(out_dir, "09_sessionInfo.txt"))
cat("Finished selection-trajectory summary:", nrow(results), "matched SNPs\n")
