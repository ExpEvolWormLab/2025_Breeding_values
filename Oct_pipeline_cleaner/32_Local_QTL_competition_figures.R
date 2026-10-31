## 32_Local_QTL_competition_figures.R
## Produce figure-ready PDFs summarising local regional competition.
## Run after 31_Local_QTL_competition_summary.R.

choose_dir <- function(){
  cand <- c('historical_QTL_R5/deep_dive_extended',
            'historical_QTL_R5/deep_dive',
            'deep_dive')
  required <- c('regional_models/27_Effects_founder_then_variants_all8.csv',
                'regional_models/27_Beta_weighted_trait_contributions.csv',
                'competition/31_Within_region_competition_summary.csv',
                'evolution/28_Exact_variant_GWAS_and_selection_summary.csv',
                'evolution/28_ALT_read_counts_all_pools.csv')
  good <- vapply(cand,function(d) all(file.exists(file.path(d,required))),logical(1))
  if(!any(good)) stop('No complete deep-dive dataset found. Expected these files ',
                     'in historical_QTL_R5/deep_dive_extended/: ',
                     paste(required,collapse=', '))
  cand[which(good)[1]]
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

## Figure 4: selected trajectories, ONE PAGE PER REGION.
## Script 30c exports 28_ALT_read_counts_all_pools.csv without a region
## column in earlier versions. Rejoin the authoritative SNP_ID-to-region
## mapping from 28_Exact_variant_GWAS_and_selection_summary.csv.
## Note that a SNP may belong to more than one nominated historical region.
## This is a presentation join only; no evolutionary data are re-estimated.

if(!nrow(traj)){
  stop('No Pool-seq trajectory rows found at ',
       file.path(evo_dir, '28_ALT_read_counts_all_pools.csv'),
       '. Run 30c or check deep_dir before trying to draw trajectories.')
}
needed <- c('SNP_ID','sample_id','generations','replicate','ALT_frequency')
missing <- setdiff(needed,names(traj))
if(length(missing))stop('Trajectory CSV missing required columns: ',
                        paste(missing,collapse=', '))

site_map_file <- file.path(evo_dir, '28_Exact_variant_GWAS_and_selection_summary.csv')
site_map <- NULL
if(file.exists(site_map_file)){
  site_map <- read.csv(site_map_file, stringsAsFactors = FALSE)
  if(!all(c('SNP_ID','region') %in% names(site_map)))
    stop('Exact-variant summary lacks SNP_ID/region: ',site_map_file)
  keep <- intersect(c('SNP_ID','region','source'),names(site_map))
  site_map <- unique(site_map[,keep,drop=FALSE])
  site_map <- site_map[!is.na(site_map$region) & nzchar(site_map$region),,drop=FALSE]
  ## If one SNP is listed as several sources within one region, prefer a
  ## nominated focal variant to a tag, then historical exact variant.
  if(!'source' %in% names(site_map))site_map$source <- 'unknown'
  site_map$source_priority <- match(site_map$source,
                                    c('focal','R5_tag','historical_exact'))
  site_map$source_priority[is.na(site_map$source_priority)] <- 9L
  site_map <- site_map[order(site_map$source_priority),,drop=FALSE]
  site_map <- site_map[!duplicated(site_map[,c('SNP_ID','region')]),,drop=FALSE]
  old_rows <- nrow(traj)
  if('region'%in%names(traj))traj$region <- NULL
  if('source'%in%names(traj))traj$source <- NULL
  traj <- merge(traj,site_map[,c('SNP_ID','region','source','source_priority'),drop=FALSE],
                by='SNP_ID',all.x=FALSE,sort=FALSE)
  cat('Trajectory region mapping: ',old_rows,' input pool rows -> ',nrow(traj),
      ' region-matched rows. Multiple region assignments may duplicate a SNP.\n',sep='')
} else if('region'%in%names(traj)){
  warning('Missing exact-site region map; using the region column already ',
          'present in the trajectory CSV')
  traj$source <- 'unknown'
  traj$source_priority <- 9L
} else {
  stop('Cannot map trajectories to QTL regions. Neither traj$region nor ',
       site_map_file,' is available. Re-run 30c or restore its exact-site summary.')
}

## Recover exact focal-site trajectories excluded by the corrected 08/09
## SNP universe (especially IV:8,892,663; IV:10,499,135; and the three
## chromosome V positions). These counts were extracted by script 28b.
## The extra rows are descriptive frequencies, not drift-calibrated tests.
raw_candidates <- c(file.path(evo_dir,'raw_sync','28b_Exact_raw_SYNC_ALT_trajectories.csv'),
  file.path('historical_QTL_R5','deep_dive','evolution','raw_sync',
            '28b_Exact_raw_SYNC_ALT_trajectories.csv'),
  file.path('deep_dive','evolution','raw_sync',
            '28b_Exact_raw_SYNC_ALT_trajectories.csv'))
raw_found <- raw_candidates[file.exists(raw_candidates)]
if(length(raw_found) && !is.null(site_map)){
  raw <- read.csv(raw_found[1],stringsAsFactors=FALSE)
  if(all(c('SNP_ID','sample_id','generation','replicate','ALT_frequency') %in% names(raw))){
    if('REF_agrees'%in%names(raw)){
      bad <- !is.na(raw$REF_agrees) & !as.logical(raw$REF_agrees)
      if(any(bad)){
        warning('Omitting ',sum(bad),' raw-sync observations where VCF and SYNC REF disagree')
        raw <- raw[!bad,,drop=FALSE]
      }
    }
    raw$generations <- raw$generation
    raw <- merge(raw,site_map[,c('SNP_ID','region','source','source_priority'),drop=FALSE],
                 by='SNP_ID',all.x=FALSE,sort=FALSE)
    if(nrow(raw)){
      status_lookup <- unique(traj[,intersect(c('sample_id','GA_status'),names(traj)),drop=FALSE])
      if('GA_status'%in%names(status_lookup)){
        status_lookup <- status_lookup[!duplicated(status_lookup$sample_id),,drop=FALSE]
        raw$GA_status <- status_lookup$GA_status[match(raw$sample_id,status_lookup$sample_id)]
      } else raw$GA_status <- NA_character_
      ## Equalize columns and deduplicate by SNP+region+pool, preferring raw
      ## recovery for exact sites when both inputs provide the same observation.
      all_cols <- union(names(raw),names(traj))
      for(nm in setdiff(all_cols,names(raw)))raw[[nm]] <- NA
      for(nm in setdiff(all_cols,names(traj)))traj[[nm]] <- NA
      combined <- rbind(raw[,all_cols,drop=FALSE],traj[,all_cols,drop=FALSE])
      key_pool <- paste(combined$SNP_ID,combined$region,combined$sample_id,sep='::')
      traj <- combined[!duplicated(key_pool),,drop=FALSE]
      cat('Included ',nrow(raw),' region-assigned raw-SYNC pool observations from ',
          raw_found[1],'\n',sep='')
    }
  } else warning('Raw-SYNC trajectory file has unexpected columns: ',raw_found[1])
} else {
  warning('No script-28b raw-SYNC output found. Rare focal variants may be ',
          'absent from the trajectories; plot will show available R5 tags.')
}

traj$generations <- suppressWarnings(as.numeric(traj$generations))
traj$ALT_frequency <- suppressWarnings(as.numeric(traj$ALT_frequency))
traj <- traj[!is.na(traj$region) & nzchar(traj$region) &
               is.finite(traj$generations) & is.finite(traj$ALT_frequency),,drop=FALSE]
## Use regions from the current phenotype models, not just the ranked top six.
## This ensures CHR IV, CHR I, CHR V and other modeled loci all appear.
model_regions <- unique(as.character(eff$region))
plot_regions <- model_regions[model_regions %in% unique(traj$region)]
missing_regions <- setdiff(model_regions,plot_regions)
if(length(missing_regions))
  warning('No mapped, finite ALT-frequency trajectories for regions: ',
          paste(missing_regions,collapse=', '))
if(!length(plot_regions))
  stop('No trajectory records match model regions. Compare eff$region with ',
       'region values in 28_Exact_variant_GWAS_and_selection_summary.csv.')

## Select at most four SNPs per region, prioritizing explicit focal variants,
## then R5 tags. Explicitly report which SNPs were used.
max_sites_per_region <- 4L
selected_sites <- list()
for(region in plot_regions){
  z <- traj[traj$region == region,,drop=FALSE]
  sites <- unique(z[,c('SNP_ID','source','source_priority'),drop=FALSE])
  sites <- sites[order(sites$source_priority,sites$SNP_ID),,drop=FALSE]
  sites <- sites[!duplicated(sites$SNP_ID),,drop=FALSE]
  if(!nrow(sites))next
  sites <- head(sites,max_sites_per_region)
  selected_sites[[region]] <- data.frame(region=region,SNP_ID=sites$SNP_ID,
                                        source=sites$source,
                                        stringsAsFactors=FALSE)
}
selected_sites <- do.call(rbind,selected_sites)
write.csv(selected_sites,file.path(comp_dir,'32_Trajectory_plotted_SNPs.csv'),row.names=FALSE)

output_pdf <- file.path(fig_dir,'32_Figure_selected_trajectories_preview.pdf')
if(file.exists(output_pdf))file.remove(output_pdf)
pdf(output_pdf,width=11,height=8.5,onefile=TRUE)
n_pages <- 0L
for(region in plot_regions){
  these <- selected_sites$SNP_ID[selected_sites$region == region]
  if(!length(these))next
  n_panels <- length(these)
  nr <- if(n_panels <= 2)1L else 2L
  nc <- if(n_panels == 1)1L else 2L
  par(mfrow=c(nr,nc),mar=c(4.1,4.3,3.5,1.0),oma=c(1,0,3,0))
  for(sid in these){
    x <- traj[traj$region == region & traj$SNP_ID == sid,,drop=FALSE]
    ## One record per SNP/pool even if two mapping rows had identical labels.
    x <- x[!duplicated(x$sample_id),,drop=FALSE]
    evolved <- x[x$generations > 0,,drop=FALSE]
    reps <- unique(as.character(evolved$replicate))
    reps <- reps[!is.na(reps) & nzchar(reps)]
    lim <- max(x$ALT_frequency,na.rm=TRUE)
    y_lim <- c(0,min(1,max(0.08,1.12*lim)))
    max_gen <- max(x$generations,na.rm=TRUE)
    plot(NA,xlim=c(0,max(2,max_gen)),ylim=y_lim,
         xlab='Salt generations',ylab='GWAS ALT frequency',
         main=sid)
    abline(h=0,lty=3,col='grey85')
    cols <- if(length(reps)) setNames(grDevices::hcl.colors(length(reps),palette='Dark 3'),reps) else character()
    for(rep in reps){
      zz <- evolved[as.character(evolved$replicate)==rep,,drop=FALSE]
      zz <- zz[order(zz$generations),,drop=FALSE]
      if(!nrow(zz))next
      stat <- if('GA_status'%in%names(zz))tolower(as.character(zz$GA_status[1])) else ''
      lty <- if(stat%in%c('extinct','unknown'))2 else 1
      lines(zz$generations,zz$ALT_frequency,type='b',pch=16,
            col=cols[[rep]],lwd=1.8,lty=lty)
    }
    g0 <- x[x$generations==0,,drop=FALSE]
    if(nrow(g0)){
      for(k in seq_len(nrow(g0))){
        is_A00 <- as.character(g0$sample_id[k])=='A00'
        points(g0$generations[k],g0$ALT_frequency[k],
               pch=if(is_A00)19 else 1,cex=1.3,col='black')
      }
      A00 <- g0$ALT_frequency[as.character(g0$sample_id)=='A00']
      if(length(A00)==1L&&is.finite(A00))abline(h=A00,lty=3,col='grey55')
    }
    if(length(reps)){
      ## replicate is an internal index (3..6), not GA1..GA4. Label using
      ## actual sample IDs: GA110, GA1G36, GA150 -> GA1.
      lbl <- vapply(reps,function(rep){
        sam <- as.character(evolved$sample_id[as.character(evolved$replicate)==rep])
        ga <- sub('^GA([0-9]+).*$', 'GA\\1', sam)
        ga <- ga[grepl('^GA[0-9]+$',ga)]
        if(length(ga))ga[1] else paste('replicate',rep)
      },character(1))
      legend('topright',legend=lbl,col=unname(cols),pch=16,lty=1,
             cex=.72,bty='n',ncol=min(2,length(reps)))
    }
  }
  mtext(paste('Region',region,'— exact-site ALT trajectories (descriptive)'),
        outer=TRUE,side=3,line=1,cex=1.1,font=2)
  n_pages <- n_pages+1L
}
dev.off()
if(n_pages<1L)stop('Internal error: no trajectory PDF pages generated')
cat('Created ',output_pdf,' with ',n_pages,' region pages and ',
    nrow(selected_sites),' site plots.\n',sep='')
cat('Site manifest: ',file.path(comp_dir,'32_Trajectory_plotted_SNPs.csv'),'\n',sep='')

cat('Figures written in ', fig_dir, '\n', sep = '')
