## 07c_Recover_G_posterior_no_refit.R
## Run from the project root. Does NOT execute MCMCglmm.
## Reconstructs the complete posterior G using three finished chain checkpoints,
## and recreates the 07a_multivariate_G_final.RData required by script 07b.

output_dir <- "RData_G_NaCl"
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")
files <- file.path(output_dir, sprintf("07a_final_unstructured_chain%02d.rds", 1:3))
if (!all(file.exists(files))) {
  stop("Missing chain files: ", paste(files[!file.exists(files)], collapse = ", "))
}
ck <- lapply(files, readRDS)
if (!all(vapply(ck, function(x) is.list(x) && is.matrix(x$VCV) &&
                 !is.null(colnames(x$VCV)) && is.list(x$meta), logical(1))))
  stop("A checkpoint is invalid. Expected list(VCV=..., meta=...).")
if (!all(vapply(ck, function(x) identical(x$meta, ck[[1]]$meta), logical(1))))
  stop("Checkpoint metadata differs between chains; do not combine models.")
if (!all(vapply(ck, function(x) identical(colnames(x$VCV),
                   colnames(ck[[1]]$VCV)), logical(1))))
  stop("The chains have different parameter columns.")
meta <- ck[[1]]$meta
if (!identical(as.character(meta$traits), traits)) stop("Trait order changed")
if (!identical(meta$run_mode, "final") ||
    !identical(meta$residual_structure, "unstructured"))
  stop("This recovery script expects the final unstructured model.")
if (!all(vapply(ck, function(x) all(is.finite(x$VCV)), logical(1))))
  stop("At least one chain has nonfinite values")
cat("Recovered", length(ck), "chains with", sum(vapply(ck, function(x)
    nrow(x$VCV), integer(1))), "total retained draws.\n")

parse_block <- function(vcv, suffix, traits) {
  ii <- grep(paste0("\\.", suffix, "$"), colnames(vcv))
  nt <- length(traits)
  if (length(ii) != nt * nt)
    stop("Expected ", nt * nt, " covariance columns in ", suffix,
         "; found ", length(ii))
  x <- array(NA_real_, c(nt, nt, nrow(vcv)),
             dimnames = list(traits, traits, NULL))
  for (j in ii) {
    name <- sub(paste0("\\.", suffix, "$"), "", colnames(vcv)[j])
    ab <- sub("^trait", "", strsplit(name, ":", fixed = TRUE)[[1]])
    if (length(ab) != 2L || !all(ab %in% traits))
      stop("Unexpected covariance column: ", colnames(vcv)[j])
    if (all(is.finite(x[ab[1], ab[2], ])))
      stop("Duplicate covariance element: ", name)
    x[ab[1], ab[2], ] <- vcv[, j]
  }
  if (anyNA(x)) stop("Missing covariance elements in ", suffix)
  if (max(abs(x - aperm(x, c(2, 1, 3)))) > 1e-8)
    stop("Posterior G is not symmetric for ", suffix)
  x
}

combine_arrays <- function(arr) {
  a <- arr[[1]]
  n <- sum(vapply(arr, function(z) dim(z)[3], integer(1)))
  out <- array(NA_real_, c(dim(a)[1], dim(a)[2], n),
               dimnames = list(dimnames(a)[[1]], dimnames(a)[[2]], NULL))
  start <- 1L
  for (z in arr) {
    ix <- seq.int(start, length.out = dim(z)[3])
    out[, , ix] <- z
    start <- max(ix) + 1L
  }
  out
}

G_posterior <- combine_arrays(lapply(ck, function(x)
  parse_block(x$VCV, "pop_label", traits)))
R_posterior <- combine_arrays(lapply(ck, function(x)
  parse_block(x$VCV, "units", traits)))
G_RIL_scaled <- apply(G_posterior, c(1, 2), mean)
trait_mean <- meta$trait_mean[traits]
trait_sd <- meta$trait_sd[traits]
if (any(!is.finite(trait_sd) | trait_sd <= 0)) stop("Bad trait scaling")
D <- diag(as.numeric(trait_sd), length(traits))
G_RIL <- D %*% G_RIL_scaled %*% D
rownames(G_RIL) <- colnames(G_RIL) <- traits
R_traits <- stats::cov2cor(G_RIL_scaled)
pair_n <- meta$pair_n
settings <- meta$settings
run_mode <- meta$run_mode
residual_structure <- meta$residual_structure

diagnostics <- read.csv(file.path(output_dir, "07a_MCMC_diagnostics.csv"))
variance_summary <- read.csv(file.path(output_dir,
 "07a_G_variance_credible_intervals.csv"))
corr_summary <- read.csv(file.path(output_dir,
 "07a_G_correlation_credible_intervals.csv"))
if (!all(is.finite(diagnostics$Rhat) & is.finite(diagnostics$ESS)) ||
    max(diagnostics$Rhat) >= 1.05 || min(diagnostics$ESS) < 200)
  stop("Diagnostics fail minimum convergence check")

## Compare reconstructed means against existing files before exporting.
g_old <- as.matrix(read.csv(file.path(output_dir, "07a_G_among_RILs.csv"),
                             row.names = 1, check.names = FALSE))[traits, traits]
r_old <- as.matrix(read.csv(file.path(output_dir, "07a_G_genetic_correlation.csv"),
                             row.names = 1, check.names = FALSE))[traits, traits]
if (max(abs(g_old - G_RIL)) > 1e-8 || max(abs(r_old - R_traits)) > 1e-8)
  stop("Recovered G does not match the previous 07a exported tables.")

## Conditional-on-model posterior summaries of genetic correlation eigenstructure.
eigdraw <- vapply(seq_len(dim(G_posterior)[3]), function(k) {
  base::eigen(stats::cov2cor(G_posterior[, , k]),
              symmetric = TRUE, only.values = TRUE)$values
}, numeric(length(traits)))
stopifnot(is.matrix(eigdraw))
Neff <- length(traits) - colMeans((eigdraw - 1)^2)
part_ratio <- length(traits)^2 / colSums(eigdraw^2)
summ <- function(x) c(mean = mean(x), lower95 = unname(quantile(x, .025)),
                      upper95 = unname(quantile(x, .975)))
write.csv(data.frame(measure = c("Neff", "participation_ratio", "PC1_fraction"),
  rbind(summ(Neff), summ(part_ratio), summ(eigdraw[1, ] / length(traits))),
  row.names = NULL), file.path(output_dir, "07c_posterior_eigensummary.csv"),
  row.names = FALSE)

target <- file.path(output_dir, "07a_multivariate_G_final.RData")
save(G_RIL, G_RIL_scaled, R_traits, G_posterior, R_posterior,
     traits, diagnostics, variance_summary, corr_summary, trait_mean,
     trait_sd, pair_n, settings, run_mode, residual_structure,
     file = target, compress = "gzip")
cat("Recovered G object saved: ", target, "\n", sep = "")
cat("Maximum Rhat = ", max(diagnostics$Rhat),
    "; minimum ESS = ", min(diagnostics$ESS), "\n", sep = "")
cat("Neff posterior 95% interval: ", paste(signif(quantile(Neff,
    c(.025, .5, .975)), 4), collapse = ", "), "\n", sep = "")
cat("Next: run 07b_Refit_pleiotropy_with_multivariate_G_FIXED.R",
    " if the original R5 QTL .RData is available.\n")
