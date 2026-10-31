## 14_Figure_multivariate_selection.R
## Manuscript figure: phenotypic selection gradients, GWAS--evolution
## associations, leave-one-chromosome-out prediction, and robustness.
##
## RUN AFTER script 13. This script only reads CSVs, makes no model fits,
## and does not modify the original results.
##
## Uses R5 LD clusters, correctly aligned Pool-seq samples, the one-G0
## drift-aware allele-frequency slope, and GWAS effects oriented to the
## cluster alleles.  All predictors are compared on the SAME 675 clusters.
##
## Figure inference:
## - Panel B: 95% chromosome:2-Mb block-bootstrap intervals (not i.i.d.).
## - Panels C/D: leave-one-chromosome-out sensitivity, not validation on
##   independent experimental evolution populations.
## - Regression/correlation is associative and does not prove causal selection.

## ------------------------------------------------------------
## 0. Settings
## ------------------------------------------------------------
threshold_clust <- 0.5
label <- paste0("R", round(threshold_clust * 10))
base_dir <- file.path("selection_results", "drift_oneG0_nullfit", label,
                      "multivariate_prediction")
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "14_Multivariate_selection_main_figure"
make_pdf <- TRUE
make_png <- TRUE
make_svg <- TRUE                 # editable vector artwork
png_dpi <- 400
width_inches <- 13.2
height_inches <- 9.6
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")
chromosomes <- c("I", "II", "III", "IV", "V", "X")

## Accessible colors on a white background
col_three <- "#648B8B"        # reduced SF + SB + area.F model
col_full <- "#166F88"         # reconstructed full seven-trait model
col_mapped <- "#8E6BA6"       # directly mapped composite beta GWAS
col_SB <- "#B07A42"           # single-trait benchmark
col_negative <- "#B45850"     # negative selection gradient
col_neutral <- "#69737E"      # other estimated gradients
col_light <- "#E8ECEF"
col_grid <- "#E2E6E9"

read_result <- function(name) {
  f <- file.path(base_dir, name)
  if (!file.exists(f))
    stop("Missing output of script 13: ", f,
         "\nRun script 13 from the project root first, or update base_dir.")
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
require_cols <- function(d, cols, file_name) {
  missing <- setdiff(cols, names(d))
  if (length(missing))
    stop(file_name, " missing columns: ", paste(missing, collapse = ", "))
}

## ------------------------------------------------------------
## 1. Load the exact outputs from script 13
## ------------------------------------------------------------
assoc <- read_result("13_Association_summary.csv")
boot <- read_result("13_Genomic_block_bootstrap.csv")
loco <- read_result("13_Leave_chromosome_out_prediction.csv")
fitted <- read_result("13_Leave_chromosome_out_fitted.csv")

require_cols(assoc, c("cohort", "predictor", "n", "r"),
             "13_Association_summary.csv")
require_cols(boot, c("cohort", "metric", "lower95", "upper95", "estimate"),
             "13_Genomic_block_bootstrap.csv")
require_cols(loco, c("model", "n", "R2_LOCO"),
             "13_Leave_chromosome_out_prediction.csv")
require_cols(fitted, c("group", "chromosome", "response",
                       "three_trait", "full_reconstructed", "full_mapped",
                       "SB"), "13_Leave_chromosome_out_fitted.csv")

## Gradient file saved by script 13. The original gradient CSV is an
## alternate input if this particular QC file is unavailable.
grad_qc_file <- file.path(base_dir, "13_Gradient_weight_QC.csv")
if (file.exists(grad_qc_file)) {
  grad <- read.csv(grad_qc_file, stringsAsFactors = FALSE)
  require_cols(grad, c("trait", "beta_BLUE", "beta_SE"),
               "13_Gradient_weight_QC.csv")
  grad_est <- as.numeric(grad$beta_BLUE)
  grad_se <- as.numeric(grad$beta_SE)
} else {
  grad_file <- "RData/selection_gradients_NaCl.csv"
  if (!file.exists(grad_file))
    stop("Neither ", grad_qc_file, " nor ", grad_file,
         " exists. Script 13 normally writes the first file.")
  grad <- read.csv(grad_file, stringsAsFactors = FALSE)
  require_cols(grad, c("trait", "beta", "SE"), grad_file)
  grad_est <- as.numeric(grad$beta)
  grad_se <- as.numeric(grad$SE)
}
if (anyDuplicated(grad$trait) || !all(traits %in% grad$trait))
  stop("Selection-gradient table needs one row per each of the seven traits")
i_grad <- match(traits, grad$trait)
G <- data.frame(trait = traits, estimate = grad_est[i_grad],
                SE = grad_se[i_grad], stringsAsFactors = FALSE)
G$lower95 <- G$estimate - 1.96 * G$SE
G$upper95 <- G$estimate + 1.96 * G$SE
if (any(!is.finite(as.matrix(G[, c("estimate", "SE", "lower95", "upper95")]))))
  stop("Nonfinite selection-gradient estimate or SE; inspect script-00 output")
if (any(G$SE < 0)) stop("Selection-gradient SE cannot be negative")
if (G$estimate[G$trait == "SF"] >= 0 ||
    G$estimate[G$trait == "SB"] <= 0 ||
    G$estimate[G$trait == "area.F"] <= 0)
  warning("Current script-00 gradients do not match the expected SF negative, ",
          "SB positive, area.F positive directions; plot shows actual values.")

## Primary all-common set, same clusters in each comparison
A <- subset(assoc, cohort == "all_common")
B <- subset(boot, cohort == "all_common")
name_by_model <- c(three_trait = "SF + SB + area.F",
                   full_mapped = "Mapped composite beta",
                   full_reconstructed = "Full seven-trait index",
                   SB = "SB alone")
color_by_model <- c(three_trait = col_three,
                    full_mapped = col_mapped,
                    full_reconstructed = col_full,
                    SB = col_SB)
boot_by_model <- c(three_trait = "three", full_mapped = "full",
                   full_reconstructed = "reconstructed")

model_names <- c("three_trait", "full_mapped", "full_reconstructed", "SB")
model_summary <- do.call(rbind, lapply(model_names, function(m) {
  aa <- A[A$predictor == m, , drop = FALSE]
  ll <- loco[loco$model == m, , drop = FALSE]
  if (nrow(aa) != 1L || nrow(ll) != 1L)
    stop("Expected exactly one all_common and one LOCO row for ", m)
  bb <- if (m %in% names(boot_by_model))
    B[B$metric == boot_by_model[[m]], , drop = FALSE] else B[FALSE, ]
  if (m != "SB" && nrow(bb) != 1L)
    stop("Expected one block-bootstrap row for ", m)
  data.frame(model = m, label = unname(name_by_model[[m]]),
             n = aa$n, r = aa$r,
             boot_lower95 = if (nrow(bb)) bb$lower95 else NA_real_,
             boot_upper95 = if (nrow(bb)) bb$upper95 else NA_real_,
             R2_LOCO = ll$R2_LOCO,
             stringsAsFactors = FALSE)
}))

if (length(unique(model_summary$n)) != 1L ||
    any(model_summary$n != nrow(fitted)))
  stop("Association and LOCO files have different cluster counts. ",
       "Re-run script 13 to produce a consistent set of outputs.")
if (anyDuplicated(fitted$group) || any(!is.finite(fitted$response)))
  stop("Held-out fitted file contains duplicate clusters or nonfinite response")
if (!all(chromosomes %in% fitted$chromosome))
  stop("Expected all six C. elegans nuclear chromosomes in held-out file")

## One point per chromosome: decrease in mean squared prediction error
## compared with the reduced three-trait predictor. Positive = full better.
chr_error <- do.call(rbind, lapply(chromosomes, function(chr) {
  d <- fitted[fitted$chromosome == chr, , drop = FALSE]
  keep <- is.finite(d$response) & is.finite(d$three_trait) &
    is.finite(d$full_reconstructed)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("No usable held-out clusters on chromosome ", chr)
  mse_three <- mean((d$response - d$three_trait)^2)
  mse_full <- mean((d$response - d$full_reconstructed)^2)
  data.frame(chromosome = chr, n_clusters = nrow(d),
             MSE_three = mse_three, MSE_full = mse_full,
             delta_MSE = mse_three - mse_full,
             stringsAsFactors = FALSE)
}))
chr_error$delta_MSE_times_1e8 <- chr_error$delta_MSE * 1e8

Delta <- B[B$metric == "delta_reconstructed_minus_three", , drop = FALSE]
if (nrow(Delta) != 1L)
  warning("Block bootstrap contrast full reconstructed - three-trait missing")

## Write plotting-data tables (without replacing any analysis outputs)
write.csv(G, file.path(fig_dir, "14_Figure_selection_gradients.csv"),
          row.names = FALSE)
write.csv(model_summary, file.path(fig_dir, "14_Figure_prediction_comparison.csv"),
          row.names = FALSE)
write.csv(chr_error, file.path(fig_dir, "14_Figure_by_chromosome.csv"),
          row.names = FALSE)

## ------------------------------------------------------------
## 2. Four-panel publication figure (base R; no package dependencies)
## ------------------------------------------------------------
plot_figure <- function() {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(2, 2),
      mar = c(4.8, 6.5, 3.4, 1.6),
      oma = c(3.1, 1.0, 2.8, 1.0),
      family = "sans", las = 1, cex.axis = 0.92, cex.lab = 1.05,
      bg = "white", fg = "#263645", col.axis = "#3E4853",
      col.lab = "#263645")

  panel_title <- function(label, title) {
    mtext(paste0(label, "   ", title), side = 3, line = 1.2,
          adj = 0, font = 2, cex = 1.06, col = "#243646")
  }
  nice_axis <- function(side = 1, ...) {
    axis(side, ..., col = "#C0C7CE", col.axis = "#354451",
         lwd = 0.7, lwd.ticks = 0.7)
  }
  thin_grid <- function(...) abline(..., col = col_grid, lwd = 0.8)

  ## A. Selection gradients. These are phenotypic selection gradients
  ## from the current script-00 regression, not conditional GWAS estimates.
  y <- rev(seq_along(traits))
  grad_rng <- range(c(0, G$lower95, G$upper95))
  pad <- max(0.08 * diff(grad_rng), 0.015)
  plot(NA, xlim = grad_rng + c(-pad, pad), ylim = c(.5, 7.5),
       type = "n", xlab = "Phenotypic selection gradient (beta)",
       ylab = "", yaxt = "n", xaxt = "n", bty = "n")
  thin_grid(v = pretty(grad_rng, n = 5))
  abline(v = 0, col = "#70818B", lty = 2, lwd = 1.3)
  nice_axis(1)
  nice_axis(2, at = y, labels = G$trait, tick = FALSE)
  cols <- ifelse(G$trait == "SF", col_negative,
                 ifelse(G$trait %in% c("SB", "area.F"), col_full, col_neutral))
  segments(G$lower95, y, G$upper95, y, col = cols, lwd = 2.4)
  segments(G$lower95, y - .10, G$lower95, y + .10, col = cols, lwd = 1.2)
  segments(G$upper95, y - .10, G$upper95, y + .10, col = cols, lwd = 1.2)
  points(G$estimate, y, pch = 21, bg = cols, col = "white", cex = 1.4, lwd = .7)
  panel_title("A", "Phenotypic selection gradients")
  mtext("Points: estimates; whiskers: 95% Wald intervals", side = 3,
        line = 0.12, adj = 0, cex = .75, col = col_neutral)

  ## B. Pearson correlations, same clusters, 2-Mb genomic-block bootstrap.
  comparison <- model_summary[match(c("full_reconstructed", "three_trait",
                                        "full_mapped"), model_summary$model), ]
  yy <- rev(seq_len(nrow(comparison)))
  xr <- range(c(0, comparison$r, comparison$boot_lower95,
                comparison$boot_upper95), na.rm = TRUE)
  xr <- c(min(-.01, xr[1] - .025), max(.35, xr[2] + .06))
  plot(NA, xlim = xr, ylim = c(.5, length(yy) + .5),
       xlab = "Pearson correlation with frequency change", ylab = "",
       yaxt = "n", xaxt = "n", bty = "n")
  thin_grid(v = pretty(xr, n = 5))
  abline(v = 0, col = "#9AA5AD", lwd = 1)
  nice_axis(1)
  nice_axis(2, at = yy, labels = comparison$label, tick = FALSE,
            cex.axis = .82)
  for (i in seq_along(yy)) {
    cc <- color_by_model[[comparison$model[i]]]
    segments(comparison$boot_lower95[i], yy[i],
             comparison$boot_upper95[i], yy[i], col = cc, lwd = 3)
    points(comparison$r[i], yy[i], pch = 21, bg = cc,
           col = "white", cex = 1.5, lwd = .6)
  }
  panel_title("B", "Multivariate indices predict evolution")
  mtext("95% genomic-block bootstrap intervals; all common SNP clusters",
        side = 3, line = .12, adj = 0, cex = .75, col = col_neutral)
  if (nrow(Delta) == 1L && all(is.finite(unlist(Delta[1,
       c("estimate", "lower95", "upper95")])))) {
    mtext(sprintf("Full minus reduced: delta r = %.3f [%.3f, %.3f]",
                  Delta$estimate, Delta$lower95, Delta$upper95),
          side = 1, line = 3.35, cex = .79, col = col_full)
  }

  ## C. Held-out predictive R2, with no false claim of independent validation.
  yy <- rev(seq_len(nrow(model_summary)))
  xr <- range(c(0, model_summary$R2_LOCO))
  xmax <- max(.075, xr[2] + .015)
  plot(NA, xlim = c(0, xmax), ylim = c(.5, length(yy) + .5),
       xlab = "Held-out predictive R-squared", ylab = "",
       yaxt = "n", xaxt = "n", bty = "n")
  thin_grid(v = pretty(c(0, xmax), n = 5))
  nice_axis(1)
  nice_axis(2, at = yy, labels = model_summary$label, tick = FALSE,
            cex.axis = .82)
  for (i in seq_along(yy)) {
    cc <- color_by_model[[model_summary$model[i]]]
    segments(0, yy[i], model_summary$R2_LOCO[i], yy[i], col = cc, lwd = 5)
    points(model_summary$R2_LOCO[i], yy[i], pch = 21, bg = cc,
           col = "white", cex = 1.5, lwd = .6)
    text(model_summary$R2_LOCO[i] + .003, yy[i],
         sprintf("%.3f", model_summary$R2_LOCO[i]),
         pos = 4, cex = .79, col = cc)
  }
  panel_title("C", "Prediction on held-out chromosomes")
  mtext("Train on five chromosomes; evaluate on sixth", side = 3,
        line = .12, adj = 0, cex = .75, col = col_neutral)

  ## D. Per-chromosome MSE difference: positive = full model better.
  yy <- rev(seq_len(nrow(chr_error)))
  x <- chr_error$delta_MSE_times_1e8
  lim <- max(abs(x), .7) * 1.20
  plot(NA, xlim = c(-lim, lim), ylim = c(.5, 6.5),
       xlab = "Reduction in MSE with full index (x 1e-8)",
       ylab = "", yaxt = "n", xaxt = "n", bty = "n")
  thin_grid(v = pretty(c(-lim, lim), n = 5))
  abline(v = 0, col = "#6E7C89", lwd = 1.25)
  nice_axis(1)
  nice_axis(2, at = yy,
            labels = paste0(chr_error$chromosome, "  (n=", chr_error$n_clusters, ")"),
            tick = FALSE, cex.axis = .84)
  cc <- ifelse(x >= 0, col_full, col_negative)
  segments(0, yy, x, yy, col = cc, lwd = 3.0)
  points(x, yy, pch = 21, bg = cc, col = "white", cex = 1.4, lwd = .5)
  panel_title("D", "Consistency across chromosomes")
  mtext(sprintf("Full index has lower squared error on %d of %d chromosomes",
                sum(chr_error$delta_MSE > 0), nrow(chr_error)),
        side = 3, line = .12, adj = 0, cex = .75, col = col_neutral)

  ## Full-canvas heading and caveat
  mtext("Pleiotropic genetic effects and evolutionary allele-frequency change",
        side = 3, outer = TRUE, adj = .5, line = 1.0,
        cex = 1.25, font = 2, col = "#213746")
  mtext(sprintf("R5 LD clusters | n = %d | common ancestral-frequency SNPs | drift-aware response | Ne ~ 1,000",
                nrow(fitted)),
        side = 1, outer = TRUE, adj = .5, line = 1.55, cex = .84,
        col = "#405260")
  mtext("Bootstrap and chromosome-held-out results are sensitivity analyses; chromosome folds are not independent evolution experiments.",
        side = 1, outer = TRUE, adj = .5, line = .40, cex = .73,
        col = col_neutral)
}

if (make_pdf) {
  pdf(file.path(fig_dir, paste0(prefix, ".pdf")),
      width = width_inches, height = height_inches, family = "Helvetica",
      paper = "special", useDingbats = FALSE)
  plot_figure()
  dev.off()
}
if (make_png) {
  png(file.path(fig_dir, paste0(prefix, ".png")),
      width = width_inches, height = height_inches, units = "in",
      res = png_dpi, bg = "white")
  plot_figure()
  dev.off()
}
if (make_svg) {
  tryCatch({
    svg(file.path(fig_dir, paste0(prefix, ".svg")),
        width = width_inches, height = height_inches,
        pointsize = 12, bg = "white")
    plot_figure()
    dev.off()
  }, error = function(e) {
    if (dev.cur() > 1L) dev.off()
    warning("SVG output failed (PDF and PNG are unaffected): ",
            conditionMessage(e))
  })
}

cat("\nMain figure saved in: ", normalizePath(fig_dir), "\n", sep = "")
cat("Model comparison (all-common clusters):\n")
print(model_summary[, c("model", "n", "r", "boot_lower95",
                         "boot_upper95", "R2_LOCO")], row.names = FALSE)
cat("\nChromosomes favoring full reconstructed over three-trait index: ",
    sum(chr_error$delta_MSE > 0), "/", nrow(chr_error), "\n", sep = "")
