## 07a_Estimate_multivariate_RIL_G.R
## Direct estimate of the seven-trait among-RIL covariance from replicated
## NaCl phenotypes. This is NOT the covariance of estimated BLUEs.
## Requires MCMCglmm and coda. Run from the project root.

library(MCMCglmm)
library(coda)

traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")
output_dir <- "RData_G_NaCl"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## Use 'pilot' to check estimability, then 'final' for results.
run_mode <- "final"                  # 'pilot' or 'final'
residual_structure <- "unstructured"  # 'unstructured' or 'diagonal'
n_chains <- 3L
settings <- if (run_mode == "final")
  list(nitt = 300000L, burnin = 50000L, thin = 125L) else
  list(nitt = 40000L, burnin = 10000L, thin = 30L)

## ------------------------------------------------------------
## 1. Match the phenotypic filtering and RIL labels used in script 00
## ------------------------------------------------------------
phen <- read.delim("data/all_phenotypes_from_movies_fixed.txt",
                   stringsAsFactors = FALSE, check.names = FALSE)
phen$pop_label <- as.character(phen$pop_label)
phen$block_id <- as.character(phen$block_id)
phen$date_str <- as.character(phen$date_str)
phen$session_id <- substr(phen$date_str, 1, 4)
phen$session_id[!is.na(phen$session_id) & phen$session_id == "2014"] <- "2013"
phen$session_id[!is.na(phen$block_id) &
                substr(phen$block_id, 1, 1) == "F"] <- "F"
phen$session_id[!is.na(phen$block_id) &
                substr(phen$block_id, 1, 2) == "SB"] <- "SB"
phen <- subset(phen, !is.na(pop_label) & nzchar(pop_label) &
  !is.na(block_id) & !(block_id %in% paste0("B", 300:350)) &
  !(pop_label %in% c("A00", "CX12311L0", "LSJ2L0")))

env_covariates <- c("temperature", "rel_humidity", "logD")
for (v in env_covariates)
  phen[[v]] <- as.numeric(as.character(phen[[v]]))
phen <- phen[complete.cases(phen[, env_covariates]), , drop = FALSE]
for (v in env_covariates) phen[[v]] <- as.numeric(scale(phen[[v]]))
phen <- phen[phen$env_label == "NaCl", , drop = FALSE]

removed <- read.table("data/removed_isotypes_Cemee_2025.txt",
                      stringsAsFactors = FALSE)[[1]]
correspondence <- read.csv("data/correspondance_isotypes.csv",
                           stringsAsFactors = FALSE)
clean_iso <- function(x) gsub("_sorted", "", gsub("CeMee", "", x, fixed = TRUE),
                              fixed = TRUE)
removed <- unique(clean_iso(removed))
correspondence$Line <- clean_iso(correspondence$Line)
correspondence$Representant_line <- clean_iso(correspondence$Representant_line)
correspondence <- unique(correspondence[, c("Line", "Representant_line")])
replaced_iso <- subset(correspondence,
                       Line %in% removed & Line != Representant_line)
if (anyDuplicated(replaced_iso$Line))
  stop("Ambiguous RIL isotype correspondence; inspect the mapping")
k <- match(phen$pop_label, replaced_iso$Line)
phen$pop_label[!is.na(k)] <- replaced_iso$Representant_line[k[!is.na(k)]]

required <- c(traits, "pop_label", "date_str", "session_id", env_covariates)
phen <- phen[, required, drop = FALSE]
for (v in traits) phen[[v]] <- as.numeric(as.character(phen[[v]]))
phen$pop_label <- factor(phen$pop_label)
phen$date_str <- factor(phen$date_str)
phen$session_id <- factor(phen$session_id)
phen <- droplevels(phen)

## Scale each response for numerical stability. Correlations are invariant
## to these positive rescalings. Keep scaling parameters to export raw G.
trait_mean <- vapply(phen[, traits, drop = FALSE], mean,
                     numeric(1), na.rm = TRUE)
trait_sd <- vapply(phen[, traits, drop = FALSE], sd,
                   numeric(1), na.rm = TRUE)
if (any(!is.finite(trait_sd) | trait_sd <= 0))
  stop("Invalid trait SD")
for (v in traits) phen[[v]] <- (phen[[v]] - trait_mean[v]) / trait_sd[v]

n_t <- length(traits)
pair_n <- outer(traits, traits, Vectorize(function(a, b)
  sum(is.finite(phen[[a]]) & is.finite(phen[[b]]))))
dimnames(pair_n) <- list(traits, traits)
write.csv(pair_n, file.path(output_dir, "07a_cophenotyped_counts.csv"))
if (residual_structure == "unstructured" &&
    any(pair_n[upper.tri(pair_n)] < 15))
  stop("Some trait pairs have <15 co-measured records. Check pairing or try residual_structure='diagonal'.")
cat("NaCl records:", nrow(phen), "; RILs:", nlevels(phen$pop_label),
    "; dates:", nlevels(phen$date_str),
    "; sessions:", nlevels(phen$session_id), "\n")

## ------------------------------------------------------------
## 2. Multivariate Gaussian random-RIL model
## ------------------------------------------------------------
## RIL covariance = G. Date/session are nuisance random effects with
## independent variances by trait; replicate-level residual covariance is
## unstructured when all phenotypes refer to the same observation unit.
response <- paste0("cbind(", paste(traits, collapse = ", "), ")")
fixed_form <- as.formula(paste(response,
  "~ trait - 1 + trait:logD + trait:rel_humidity + trait:temperature"))
random_form <- ~ us(trait):pop_label + idh(trait):date_str + idh(trait):session_id
rcov_form <- if (residual_structure == "unstructured")
  ~ us(trait):units else ~ idh(trait):units

## Inverse-Wishart priors are regularizing; evaluate sensitivity, especially
## if RIL genetic variances are small. n_t+1 is proper and moderately weak.
prior <- list(
  G = list(
    G1 = list(V = diag(n_t) * 0.25, nu = n_t + 1),
    G2 = list(V = diag(n_t) * 0.08, nu = n_t + 1),
    G3 = list(V = diag(n_t) * 0.08, nu = n_t + 1)),
  R = list(V = diag(n_t) * 0.65, nu = n_t + 1)
)

chains_VCV <- vector("list", n_chains)
for (chain in seq_len(n_chains)) {
  cat("Fitting MCMC chain", chain, "of", n_chains, "\n")
  set.seed(20261028L + chain)
  fit <- MCMCglmm::MCMCglmm(
    fixed = fixed_form, random = random_form, rcov = rcov_form,
    family = rep("gaussian", n_t), data = phen, prior = prior,
    nitt = settings$nitt, burnin = settings$burnin,
    thin = settings$thin, verbose = FALSE, pr = FALSE
  )
  chains_VCV[[chain]] <- as.matrix(fit$VCV)
  rm(fit)
  gc(verbose = FALSE)
}

## ------------------------------------------------------------
## 3. Reconstruct G (and residual R) from named posterior columns
## ------------------------------------------------------------
parse_block <- function(vcv, suffix, traits) {
  # Expected MCMCglmm label: traitSF:traitSB.pop_label
  nms <- colnames(vcv)
  cols <- grep(paste0("\\.", suffix, "$"), nms)
  expected <- length(traits) * (length(traits) + 1) / 2
  if (length(cols) != expected)
    stop("Expected ", expected, " covariance terms for '", suffix,
         "', found ", length(cols), ". Available names: ",
         paste(head(nms, 25), collapse = ", "))
  arr <- array(NA_real_, c(length(traits), length(traits), nrow(vcv)),
               dimnames = list(traits, traits, NULL))
  for (j in cols) {
    nm <- sub(paste0("\\.", suffix, "$"), "", nms[j])
    pair <- strsplit(nm, ":", fixed = TRUE)[[1]]
    pair <- sub("^trait", "", pair)
    if (length(pair) != 2 || !all(pair %in% traits))
      stop("Unexpected MCMCglmm term name: ", nms[j])
    arr[pair[1], pair[2], ] <- vcv[, j]
    arr[pair[2], pair[1], ] <- vcv[, j]
  }
  if (anyNA(arr)) stop("Incomplete covariance reconstruction for ", suffix)
  arr
}

concat_draws <- function(arrays) {
  d <- dim(arrays[[1]])
  out <- array(NA_real_, c(d[1], d[2], sum(vapply(arrays, function(a)
    dim(a)[3], integer(1)))), dimnames = list(dimnames(arrays[[1]])[[1]],
      dimnames(arrays[[1]])[[2]], NULL))
  pos <- 1L
  for (a in arrays) {
    ids <- seq.int(pos, length.out = dim(a)[3])
    out[, , ids] <- a
    pos <- max(ids) + 1L
  }
  out
}
G_posterior <- concat_draws(lapply(chains_VCV, parse_block,
                                    suffix = "pop_label", traits = traits))
R_posterior <- if (residual_structure == "unstructured")
  concat_draws(lapply(chains_VCV, parse_block,
                      suffix = "units", traits = traits)) else NULL
G_RIL_scaled <- apply(G_posterior, c(1, 2), mean)
G_RIL_scaled <- (G_RIL_scaled + t(G_RIL_scaled)) / 2
D <- diag(trait_sd[traits], n_t)
G_RIL <- D %*% G_RIL_scaled %*% D
rownames(G_RIL) <- colnames(G_RIL) <- traits
R_traits <- cov2cor(G_RIL_scaled)

## Diagnostics on variance/covariance components, not only fixed effects
mcmc_chains <- do.call(coda::mcmc.list,
                      lapply(chains_VCV, coda::as.mcmc))
psrf <- coda::gelman.diag(mcmc_chains, multivariate = FALSE)$psrf
ESS <- coda::effectiveSize(mcmc_chains)
diagnostics <- data.frame(term = rownames(psrf),
                          Rhat = psrf[, 1],
                          ESS = as.numeric(ESS[rownames(psrf)]),
                          row.names = NULL)
write.csv(diagnostics, file.path(output_dir, "07a_MCMC_diagnostics.csv"),
          row.names = FALSE)
write.csv(G_RIL, file.path(output_dir, "07a_G_among_RILs.csv"))
write.csv(R_traits, file.path(output_dir, "07a_G_genetic_correlation.csv"))

## Credible intervals on correlations, not only the posterior-mean matrix.
## A strongly positive posterior genetic variance can still be uncertain.
corr_summary <- do.call(rbind, lapply(seq_len(n_t), function(i)
  do.call(rbind, lapply(i:n_t, function(j) {
    r <- vapply(seq_len(dim(G_posterior)[3]), function(k) {
      a <- G_posterior[, , k]
      a[i, j] / sqrt(a[i, i] * a[j, j])
    }, numeric(1))
    data.frame(trait1 = traits[i], trait2 = traits[j],
      correlation_mean = mean(r),
      lower95 = as.numeric(quantile(r, .025)),
      upper95 = as.numeric(quantile(r, .975)))
  })))
))
write.csv(corr_summary,
          file.path(output_dir, "07a_G_correlation_credible_intervals.csv"),
          row.names = FALSE)

cat("G eigenvalues:", paste(signif(eigen(G_RIL, symmetric = TRUE)$values, 3),
                              collapse = ", "), "\n")
cat("Maximum Rhat:", max(diagnostics$Rhat, na.rm = TRUE),
    "; minimum ESS:", min(diagnostics$ESS, na.rm = TRUE), "\n")
if (any(diagnostics$Rhat > 1.05, na.rm = TRUE) ||
    any(diagnostics$ESS < 200, na.rm = TRUE))
  warning("MCMC diagnostics indicate inadequate convergence/mixing. Do not interpret G yet.")

pdf(file.path(output_dir, "07a_G_diagnostics.pdf"), width = 8, height = 6)
image(seq_along(traits), seq_along(traits), R_traits,
      xaxt = "n", yaxt = "n", xlab = "", ylab = "",
      main = "Among-RIL trait correlations", zlim = c(-1, 1))
axis(1, seq_along(traits), labels = traits, las = 2)
axis(2, seq_along(traits), labels = traits, las = 2)
plot(diagnostics$Rhat, pch = 16, xlab = "Covariance parameter",
     ylab = "Potential scale reduction (Rhat)", las = 1, bty = "n")
abline(h = 1.05, lty = 2, col = "red")
dev.off()

save(G_RIL, G_RIL_scaled, R_traits, G_posterior, R_posterior,
     traits, diagnostics, trait_mean, trait_sd, pair_n, settings,
     run_mode, residual_structure,
     file = file.path(output_dir, paste0("07a_multivariate_G_", run_mode, ".RData")),
     compress = "gzip")
capture.output(sessionInfo(), file = file.path(output_dir, "07a_sessionInfo.txt"))
cat("Saved multivariate RIL G in", output_dir, "\n")
