## 07b_Refit_pleiotropy_with_multivariate_G.R
## Run after 07a FINAL and original 07_Characterize_pleiotropy.R.
## Replaces BLUE-derived Neff with multivariate among-RIL covariance Neff.
## Writes distinct outputs; does not overwrite script 07 or script 05.

threshold_clust <- 0.5
label <- paste0("R", round(threshold_clust * 10))
out_dir <- file.path("pleiotropy_results", label, "multivariate_G")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_env <- new.env(parent = globalenv())
load("RData_G_NaCl/07a_multivariate_G_final.RData", envir = g_env)
if (any(g_env$diagnostics$Rhat > 1.05, na.rm = TRUE) ||
    any(g_env$diagnostics$ESS < 200, na.rm = TRUE))
  stop("G MCMC convergence insufficient; inspect 07a diagnostics first")

q_env <- new.env(parent = globalenv())
load(file.path("pleiotropy_results", label,
  paste0("GWAS_clusters_with_pleiotropy_", label, ".RData")), envir = q_env)
qtl <- q_env$all_GWAS_clusters
traits <- g_env$traits
sig_cols <- paste0("is_sign_", traits, "_NaCl")
S <- as.matrix(qtl[, sig_cols, drop = FALSE])
S[is.na(S)] <- FALSE
storage.mode(S) <- "logical"

## IMPORTANT: N_eff uses the correlation of the seven-trait RIL random-effect
## covariance, NOT covariance of BLUEs and NOT a diagonal error subtraction.
Neff_fun <- function(affected, R) {
  n <- length(affected)
  if (n == 0) return(NA_real_)
  if (n == 1) return(1)
  eig <- eigen(R[affected, affected, drop = FALSE],
               symmetric = TRUE, only.values = TRUE)$values
  n - mean((eig - 1)^2)
}

qtl$Neff_BLUE <- qtl$Neff
qtl$Neff <- vapply(seq_len(nrow(qtl)), function(i)
  Neff_fun(traits[S[i, ]], g_env$R_traits), numeric(1))

## Draw-wise Neff uncertainty from posterior G (conditional on the selected
## affected-trait sets; does NOT integrate GWAS discovery uncertainty).
set.seed(20261028L)
Gdraw <- g_env$G_posterior
n_draw <- dim(Gdraw)[3]
draw_ids <- seq_len(n_draw)
if (n_draw > 500L) draw_ids <- sample(draw_ids, 500L)
patterns <- apply(S, 1, function(z) paste(which(z), collapse = ","))
unique_patterns <- unique(patterns)
Neff_by_pattern <- matrix(NA_real_, nrow = length(unique_patterns),
                          ncol = length(draw_ids),
                          dimnames = list(unique_patterns, NULL))
for (k in seq_along(draw_ids)) {
  R <- cov2cor(Gdraw[, , draw_ids[k]])
  for (pat in unique_patterns) {
    idx <- if (nzchar(pat)) as.integer(strsplit(pat, ",", fixed = TRUE)[[1]])
           else integer()
    Neff_by_pattern[pat, k] <- Neff_fun(traits[idx], R)
  }
}
qtl$Neff_post_lower <- apply(Neff_by_pattern, 1, quantile, .025,
                              na.rm = TRUE)[patterns]
qtl$Neff_post_upper <- apply(Neff_by_pattern, 1, quantile, .975,
                              na.rm = TRUE)[patterns]

## ------------------------------------------------------------
## Original model, with a common dataset for both distances
## ------------------------------------------------------------
fit_scaling <- function(dat, response, c0) {
  d <- dat[dat$chr_type2 %in% c("arm", "center") &
             dat$Ntraits > 0 & is.finite(dat$Neff) &
             dat$Neff > 0 & is.finite(dat[[response]]) &
             dat[[response]] > 0, , drop = FALSE]
  d$region <- factor(d$chr_type2, levels = c("arm", "center"))
  d$logN <- log(d$Neff)
  d$z <- log(d[[response]]) - c0 * d$logN
  if (any(table(d$region) < 5) ||
      any(tapply(d$Neff, d$region, function(z) length(unique(z))) < 2))
    stop("Insufficient Neff variation for ", response)
  null <- lm(z ~ region, data = d)
  common <- lm(z ~ region + logN, data = d)
  full <- lm(z ~ region * logN, data = d)
  b <- coef(full); V <- vcov(full); df <- df.residual(full)
  contrast <- function(L) {
    est <- sum(L * b)
    se <- sqrt(drop(t(L) %*% V %*% L))
    c(exponent = c0 + est, SE = se,
      lower95 = c0 + est - qt(.95, df) * se,
      p_one_sided = pt(est / se, df, lower.tail = FALSE))
  }
  A <- setNames(rep(0, length(b)), names(b)); A["logN"] <- 1
  C <- A; C["regioncenter:logN"] <- 1
  slopes <- rbind(arm = contrast(A), center = contrast(C))
  list(data = d, null = null, common = common, full = full,
       slopes = slopes,
       p_interaction = anova(common, full)[["Pr(>F)"]][2],
       center_minus_arm = unname(b["regioncenter:logN"]))
}

M <- fit_scaling(qtl, "Manhattan", 1)
E <- fit_scaling(qtl, "Euclidean", .5)
slopes <- rbind(data.frame(metric = "Manhattan", region = rownames(M$slopes),
                           M$slopes, row.names = NULL),
                data.frame(metric = "Euclidean", region = rownames(E$slopes),
                           E$slopes, row.names = NULL))
slopes$p_Holm <- p.adjust(slopes$p_one_sided, "holm")
interactions <- data.frame(metric = c("Manhattan", "Euclidean"),
  difference = c(M$center_minus_arm, E$center_minus_arm),
  p = c(M$p_interaction, E$p_interaction))
interactions$p_Holm <- p.adjust(interactions$p, "holm")

## Propagate only uncertainty in G into fitted exponents (a sensitivity check).
## These are NOT full posterior intervals for the slope, because GWAS and
## regression uncertainty are treated separately.
ng <- min(200L, length(draw_ids))
G_sensitivity <- do.call(rbind, lapply(seq_len(ng), function(k) {
  d <- qtl
  d$Neff <- Neff_by_pattern[patterns, k]
  m <- try(fit_scaling(d, "Manhattan", 1), silent = TRUE)
  e <- try(fit_scaling(d, "Euclidean", .5), silent = TRUE)
  if (inherits(m, "try-error") || inherits(e, "try-error")) return(NULL)
  data.frame(draw = k, metric = rep(c("Manhattan", "Euclidean"), each = 2),
    region = rep(c("arm", "center"), 2),
    exponent = c(m$slopes[, "exponent"], e$slopes[, "exponent"]))
}))

write.csv(qtl, file.path(out_dir, "QTL_pleiotropy_multivariateG.csv"),
          row.names = FALSE)
write.csv(slopes, file.path(out_dir, "Scaling_exponents_multivariateG.csv"),
          row.names = FALSE)
write.csv(interactions, file.path(out_dir, "Scaling_interactions_multivariateG.csv"),
          row.names = FALSE)
write.csv(G_sensitivity, file.path(out_dir, "Scaling_G_uncertainty_sensitivity.csv"),
          row.names = FALSE)

pdf(file.path(out_dir, "07b_pleiotropy_multivariateG.pdf"), width = 12, height = 4)
oldpar <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(qtl$Neff_BLUE, qtl$Neff, pch = 16, cex = .5,
     xlab = "Neff from BLUE correlations", ylab = "Neff from RIL G",
     bty = "n", main = "Covariance-model sensitivity")
abline(0, 1, col = "red", lty = 2)
for (response in c("Manhattan", "Euclidean")) {
  fit <- if (response == "Manhattan") M else E
  c0 <- if (response == "Manhattan") 1 else .5
  d <- fit$data
  colors <- c(arm = "forestgreen", center = "deepskyblue3")
  plot(d$Neff, d[[response]], pch = 21,
       bg = adjustcolor(colors[as.character(d$region)], alpha.f = .5),
       xlab = "Neff (RIL multivariate G)", ylab = paste(response, "effect"),
       main = paste0(response, " (null = ", c0, ")"), bty = "n")
  x <- seq(min(d$Neff), max(d$Neff), length.out = 200)
  for (r in levels(d$region)) {
    new <- data.frame(logN = log(x),
                      region = factor(rep(r, length(x)), levels = levels(d$region)))
    lines(x, exp(predict(fit$full, newdata = new) + c0 * log(x)),
          col = colors[r], lwd = 3)
    lines(x, exp(predict(fit$null, newdata = new) + c0 * log(x)),
          col = colors[r], lwd = 2, lty = 2)
  }
}
par(oldpar); dev.off()
save(qtl, M, E, slopes, interactions, G_sensitivity,
     file = file.path(out_dir, "07b_pleiotropy_multivariateG.RData"))
capture.output(sessionInfo(), file = file.path(out_dir, "07b_sessionInfo.txt"))
cat("Multivariate-G pleiotropy analysis saved to", out_dir, "\n")
