## 32_Local_QTL_competition_figures.R
## Produce figure-ready PDFs summarising local regional competition.
## Run after 31_Local_QTL_competition_summary.R.

choose_dir <- function(){
  cand <- c('historical_QTL_R5/deep_dive_extended',
            'historical_QTL_R5/deep_dive',
            'deep_dive')
  hit <- cand[file.exists(cand)][1]
  if(is.na(hit)) stop('Cannot find a deep-dive directory')
  hit
}
deep_dir <- choose_dir()
reg_dir <- file.path(deep_dir, 'regional_models')
evo_dir <- file.path(deep_dir, 'evolution')
comp_dir <- file.path(deep_dir, 'competition')
fig_dir  <- file.path(comp_dir, 'figures')
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

eff <- read.csv(file.path(reg_dir, '27_Effects_founder_then_variants_all8.csv'), stringsAsFactors = FALSE)
dec <- read.csv(file.path(reg_dir, '27_Beta_weighted_trait_contributions.csv'), stringsAsFactors = FALSE)
pairs <- read.csv(file.path(comp_dir, '31_Within_region_competition_summary.csv'), stringsAsFactors = FALSE)
rank <- read.csv(file.path(comp_dir, '31_Region_priority_ranking.csv'), stringsAsFactors = FALSE)
traj <- if(file.exists(file.path(evo_dir, '28_ALT_read_counts_all_pools.csv'))) read.csv(file.path(evo_dir, '28_ALT_read_counts_all_pools.csv'), stringsAsFactors = FALSE) else data.frame()

traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
eff$effect_SD <- with(eff, effect / trait_sd)

## Figure 1: heatmaps of conditional effects.
pdf(file.path(fig_dir, '32_Figure_heatmaps_effect_profiles.pdf'), width = 10, height = 7)
for(region in unique(eff$region)){
  z <- eff[eff$region == region &
             eff$stage %in% c('A0_founder_proxy_among_focal_noncarriers_only',
                              'A_founder_proxy_before_variants',
                              'B_variant_conditional_on_founder_proxies_and_other_variants'), , drop = FALSE]
  if(!nrow(z)) next
  terms <- unique(z$term)
  mat <- matrix(NA_real_, nrow = length(terms), ncol = length(traits),
                dimnames = list(terms, traits))
  sig <- mat
  for(i in seq_len(nrow(z))){
    if(z$trait[i] %in% traits && z$term[i] %in% terms){
      mat[z$term[i], z$trait[i]] <- z$effect_SD[i]
      sig[z$term[i], z$trait[i]] <- is.finite(z$q_BH_eight_traits[i]) && z$q_BH_eight_traits[i] < 0.10
    }
  }
  rr <- range(mat, finite = TRUE)
  lim <- max(abs(rr), na.rm = TRUE); if(!is.finite(lim) || lim == 0) lim <- 1
  par(mar = c(6, 8, 4, 2))
  image(x = seq_len(ncol(mat)), y = seq_len(nrow(mat)), z = t(mat[nrow(mat):1, , drop = FALSE]),
        axes = FALSE, xlab = '', ylab = '', col = colorRampPalette(c('#2166ac','white','#b2182b'))(101),
        zlim = c(-lim, lim), main = paste(region, '- conditional effect profiles'))
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2)
  axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2)
  box()
  for(i in seq_len(nrow(mat))) for(j in seq_len(ncol(mat))){
    y0 <- nrow(mat) - i + 1
    if(isTRUE(sig[i, j])) text(j, y0, labels = '*', cex = 1.3)
    if(is.na(mat[i, j])) text(j, y0, labels = '·', col = 'grey50')
  }
  mtext('* q < 0.10 within region/stage/term across eight traits', side = 1, line = 4, cex = 0.8)
}
dev.off()

## Figure 2: beta decompositions.
pdf(file.path(fig_dir, '32_Figure_beta_decomposition.pdf'), width = 10, height = 7)
for(region in unique(dec$region)){
  z <- dec[dec$region == region, , drop = FALSE]
  if(!nrow(z)) next
  terms <- unique(z$term)
  cols <- c('SF'='#1b9e77','SB'='#d95f02','FS'='#7570b3','FB'='#e7298a',
            'BS'='#66a61e','BF'='#e6ab02','area.F'='#a6761d')
  vals <- split(z, z$term)
  par(mar = c(8, 4, 4, 1))
  ymax <- max(abs(tapply(z$weighted_beta_contribution, z$term, sum)), na.rm = TRUE)
  bp <- barplot(rep(0, length(terms)), ylim = c(-1.2 * ymax, 1.2 * ymax),
                names.arg = terms, las = 2, main = paste(region, '- beta decomposition'),
                ylab = 'Weighted contribution to beta')
  abline(h = 0, lty = 2, col = 'grey60')
  for(k in seq_along(terms)){
    zz <- vals[[terms[k]]]
    zz <- zz[match(names(cols), zz$trait, nomatch = 0), , drop = FALSE]
    pos <- neg <- 0
    for(i in seq_len(nrow(zz))){
      h <- zz$weighted_beta_contribution[i]
      if(!is.finite(h)) next
      if(h >= 0){
        rect(bp[k] - 0.35, pos, bp[k] + 0.35, pos + h, col = cols[zz$trait[i]], border = NA)
        pos <- pos + h
      } else {
        rect(bp[k] - 0.35, neg, bp[k] + 0.35, neg + h, col = cols[zz$trait[i]], border = NA)
        neg <- neg + h
      }
    }
  }
  legend('topright', legend = names(cols), fill = cols, bty = 'n', cex = 0.75)
}
dev.off()

## Figure 3: pairwise competition summary.
pdf(file.path(fig_dir, '32_Figure_pairwise_competition_summary.pdf'), width = 9, height = 6)
if(nrow(pairs)){
  x <- pairs$corr_effect_profile
  y <- pairs$fraction_opposite_sign
  labs <- paste(pairs$region, '\n', pairs$term1, ' vs ', pairs$term2, sep = '')
  col <- ifelse(pairs$competition_label == 'redundant_or_same_haplotype', 'grey60',
         ifelse(pairs$competition_label == 'opposite_beta_directions', '#b2182b',
         ifelse(pairs$competition_label == 'different_trait_axes', '#2166ac', '#762a83')))
  cex <- 0.8 + pairs$competition_score / max(pairs$competition_score, na.rm = TRUE)
  plot(x, y, pch = 19, col = col, cex = cex,
       xlab = 'Correlation of local multivariate effect profiles',
       ylab = 'Fraction of traits with opposite effect sign',
       main = 'Local QTL competition summary')
  abline(v = 0, h = 0.5, lty = 2, col = 'grey70')
  text(x, y, labels = labs, pos = 3, cex = 0.6)
  legend('bottomleft', legend = c('redundant/same haplotype','opposite beta directions','different trait axes','mixed'),
         col = c('grey60','#b2182b','#2166ac','#762a83'), pch = 19, bty = 'n', cex = 0.75)
} else {
  plot.new(); text(.5,.5,'No pairwise competition table found')
}
dev.off()

## Figure 4: simple trajectory preview for exact focal or tag SNPs, if present.
if(nrow(traj)){
  pdf(file.path(fig_dir, '32_Figure_selected_trajectories_preview.pdf'), width = 8.5, height = 5.5)
  keep_regions <- head(rank$region, 6)
  for(region in keep_regions){
    z <- traj[traj$region == region & is.finite(traj$ALT_frequency), , drop = FALSE]
    if(!nrow(z)) next
    ids <- unique(z$SNP_ID)
    ids <- ids[seq_len(min(length(ids), 4L))]
    for(sid in ids){
      x <- z[z$SNP_ID == sid, , drop = FALSE]
      if(!nrow(x)) next
      reps <- unique(x$replicate[x$generations > 0])
      cols <- setNames(rainbow(length(reps)), reps)
      par(mar = c(4,4,3,1))
      plot(NA, xlim = c(0, max(x$generations, na.rm = TRUE)), ylim = c(0,1),
           xlab = 'Generation', ylab = 'ALT frequency',
           main = paste(region, '-', sid))
      for(r in reps){
        zz <- x[x$replicate == r & x$generations > 0, , drop = FALSE]
        zz <- zz[order(zz$generations), , drop = FALSE]
        lines(zz$generations, zz$ALT_frequency, type = 'b', pch = 16, col = cols[r])
      }
      g0 <- x[x$generations == 0, , drop = FALSE]
      if(nrow(g0)) points(g0$generations, g0$ALT_frequency, pch = 1, cex = 1.2)
      legend('topleft', legend = reps, col = cols, pch = 16, bty = 'n', cex = 0.8)
    }
  }
  dev.off()
}

cat('Figures written in ', fig_dir, '\n', sep = '')
