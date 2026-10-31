## 31_Local_QTL_competition_summary.R
## Summarise evidence for local competition / trait trade-offs among terms
## within each nominated region. Uses outputs of scripts 27 and 28.
## Works on the current deep-dive or on an extended rerun.

choose_dir <- function(){
  cand <- c('historical_QTL_R5/deep_dive_extended',
            'historical_QTL_R5/deep_dive',
            'deep_dive')
  hit <- cand[file.exists(cand)][1]
  if(is.na(hit)) stop('Cannot find a deep-dive directory')
  hit
}
deep_dir <- choose_dir()
reg_dir  <- file.path(deep_dir, 'regional_models')
evo_dir  <- file.path(deep_dir, 'evolution')
out_dir  <- file.path(deep_dir, 'competition')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

effects_file <- file.path(reg_dir, '27_Effects_founder_then_variants_all8.csv')
decomp_file  <- file.path(reg_dir, '27_Beta_weighted_trait_contributions.csv')
qc_file      <- file.path(reg_dir, '27_Nested_model_tests_and_identifiability.csv')
epi_file     <- file.path(reg_dir, '27_Pairwise_epistasis_all8.csv')
exact_file   <- file.path(evo_dir, '28_Exact_variant_GWAS_and_selection_summary.csv')
region_sel_file <- file.path(evo_dir, '28_R5_region_full_vs_retained_GA_overlap.csv')

req <- c(effects_file, decomp_file, qc_file, exact_file, region_sel_file)
if(any(!file.exists(req)))stop('Missing files: ', paste(req[!file.exists(req)], collapse = ', '))

eff <- read.csv(effects_file, stringsAsFactors = FALSE)
dec <- read.csv(decomp_file,  stringsAsFactors = FALSE)
qc  <- read.csv(qc_file,      stringsAsFactors = FALSE)
ex  <- read.csv(exact_file,   stringsAsFactors = FALSE)
rsel<- read.csv(region_sel_file, stringsAsFactors = FALSE)
epi <- if(file.exists(epi_file)) read.csv(epi_file, stringsAsFactors = FALSE) else data.frame()

eff$effect_SD <- with(eff, effect / trait_sd)
eff$term_type <- ifelse(grepl('^H_', eff$term), 'founder_proxy',
                  ifelse(grepl('^F', eff$term), 'focal_variant',
                  ifelse(grepl('^R5_', eff$term), 'R5_tag', 'other')))

trim_term <- function(x) sub('^(H_|R5_|F)', '', x)

## Per-term summary
term_key <- unique(eff[, c('region','stage','term','term_type')])
term_out <- vector('list', nrow(term_key))
for(i in seq_len(nrow(term_key))){
  z <- eff[eff$region == term_key$region[i] & eff$stage == term_key$stage[i] & eff$term == term_key$term[i], , drop = FALSE]
  z <- z[is.finite(z$effect_SD), , drop = FALSE]
  if(!nrow(z)) next
  ord <- order(z$p, abs(z$effect_SD), decreasing = FALSE, na.last = TRUE)
  best <- z[ord[1], , drop = FALSE]
  zz <- z[order(abs(z$effect_SD), decreasing = TRUE), , drop = FALSE]
  beta_row <- z[z$trait == 'beta', , drop = FALSE]
  if(nrow(beta_row) == 0) beta_row <- z[0, ]
  d <- dec[dec$region == term_key$region[i] & dec$term == term_key$term[i], , drop = FALSE]
  top_pos_trait <- top_neg_trait <- top_pos_contrib <- top_neg_contrib <- NA
  if(nrow(d)){
    dp <- d[order(d$weighted_beta_contribution, decreasing = TRUE), , drop = FALSE]
    dn <- d[order(d$weighted_beta_contribution, decreasing = FALSE), , drop = FALSE]
    top_pos_trait <- dp$trait[1]
    top_neg_trait <- dn$trait[1]
    top_pos_contrib <- dp$weighted_beta_contribution[1]
    top_neg_contrib <- dn$weighted_beta_contribution[1]
  }
  term_out[[i]] <- data.frame(
    region = term_key$region[i],
    stage = term_key$stage[i],
    term = term_key$term[i],
    term_type = term_key$term_type[i],
    n_traits = nrow(z),
    best_trait_by_p = best$trait[1],
    best_p = best$p[1],
    best_q8 = best$q_BH_eight_traits[1],
    largest_abs_trait = zz$trait[1],
    largest_abs_effect_SD = zz$effect_SD[1],
    beta_effect = if(nrow(beta_row)) beta_row$effect_SD[1] else NA_real_,
    beta_p = if(nrow(beta_row)) beta_row$p[1] else NA_real_,
    top_positive_beta_component_trait = top_pos_trait,
    top_positive_beta_component = top_pos_contrib,
    top_negative_beta_component_trait = top_neg_trait,
    top_negative_beta_component = top_neg_contrib,
    stringsAsFactors = FALSE
  )
}
term_out <- do.call(rbind, term_out)

## Map selection evidence back to each term.
term_out$selection_label <- NA_character_
term_out$phase_aligned_full_slope <- NA_real_
term_out$phase_aligned_retained_slope <- NA_real_
term_out$n_selected_common_q01 <- NA_real_
term_out$n_selected_survivor <- NA_real_
term_out$selected_full <- NA
term_out$selected_retained <- NA
for(i in seq_len(nrow(term_out))){
  tt <- term_out$term[i]
  if(grepl('^F', tt)){
    pos <- suppressWarnings(as.integer(sub('^F', '', tt)))
    rx <- ex[ex$region == term_out$region[i] & grepl(paste0('_', pos, '$'), ex$SNP_ID), , drop = FALSE]
    if(nrow(rx)){
      term_out$selection_label[i] <- rx$selection_status[1]
      if('phase_aligned_full_slope' %in% names(rx)) term_out$phase_aligned_full_slope[i] <- rx$phase_aligned_full_slope[1]
      if('phase_aligned_retained_slope' %in% names(rx)) term_out$phase_aligned_retained_slope[i] <- rx$phase_aligned_retained_slope[1]
      if('all_GA_p_drift_BH_common' %in% names(rx) && is.finite(rx$all_GA_p_drift_BH_common[1]))
        term_out$selected_full[i] <- rx$all_GA_p_drift_BH_common[1] < 0.01
      if('retained_GA_p_drift_BH_common' %in% names(rx) && is.finite(rx$retained_GA_p_drift_BH_common[1]))
        term_out$selected_retained[i] <- rx$retained_GA_p_drift_BH_common[1] < 0.01
    }
  } else {
    g <- trim_term(tt)
    rg <- rsel[rsel$region == term_out$region[i] & as.character(rsel$R5_group) == g, , drop = FALSE]
    if(nrow(rg)){
      if('full_slope_median_common' %in% names(rg)) term_out$phase_aligned_full_slope[i] <- rg$full_slope_median_common[1]
      if('slope_survivor_median_common' %in% names(rg)) term_out$phase_aligned_retained_slope[i] <- rg$slope_survivor_median_common[1]
      if('n_selected_common_q01' %in% names(rg)) term_out$n_selected_common_q01[i] <- rg$n_selected_common_q01[1]
      if('n_selected_survivor' %in% names(rg)) term_out$n_selected_survivor[i] <- rg$n_selected_survivor[1]
      if('selected_full' %in% names(rg)) term_out$selected_full[i] <- as.logical(rg$selected_full[1])
      if('selected_retained' %in% names(rg)) term_out$selected_retained[i] <- as.logical(rg$selected_retained[1])
      term_out$selection_label[i] <- 'R5_group_summary'
    }
  }
}

## Pairwise local competition summaries.
get_vec <- function(region, term){
  z <- eff[eff$region == region & eff$term == term, c('trait','effect_SD','p','q_BH_eight_traits'), drop = FALSE]
  z <- z[!duplicated(z$trait), , drop = FALSE]
  rownames(z) <- z$trait
  z
}
pairs_out <- list(); counter <- 0L
for(region in unique(term_out$region)){
  here <- term_out[term_out$region == region, , drop = FALSE]
  ## Compare conditional marker terms plus founder terms among focal noncarriers if available.
  keep <- here$stage %in% c('B_variant_conditional_on_founder_proxies_and_other_variants',
                            'A0_founder_proxy_among_focal_noncarriers_only',
                            'A_founder_proxy_before_variants')
  here <- here[keep, , drop = FALSE]
  here <- here[order(here$stage, here$term), , drop = FALSE]
  if(nrow(here) < 2) next
  idx <- combn(seq_len(nrow(here)), 2, simplify = FALSE)
  for(pair in idx){
    a <- here[pair[1], , drop = FALSE]
    b <- here[pair[2], , drop = FALSE]
    va <- get_vec(region, a$term)
    vb <- get_vec(region, b$term)
    tr <- intersect(rownames(va), rownames(vb))
    if(length(tr) < 4) next
    xa <- va[tr, 'effect_SD']
    xb <- vb[tr, 'effect_SD']
    ok <- is.finite(xa) & is.finite(xb)
    if(sum(ok) < 4) next
    xa <- xa[ok]; xb <- xb[ok]; tr2 <- tr[ok]
    corr <- suppressWarnings(cor(xa, xb))
    opp <- mean(sign(xa) != sign(xb))
    mag <- mean(abs(xa - xb))
    topa <- a$largest_abs_trait; topb <- b$largest_abs_trait
    bopp <- is.finite(a$beta_effect) && is.finite(b$beta_effect) && sign(a$beta_effect) != sign(b$beta_effect)
    topdiff <- !is.na(topa) && !is.na(topb) && topa != topb
    score <- opp + 0.5 * as.numeric(topdiff) + 0.5 * as.numeric(isTRUE(bopp)) + 0.5 * as.numeric(is.finite(corr) && corr < 0)
    label <- if(is.finite(corr) && corr > 0.75) 'redundant_or_same_haplotype'
      else if(isTRUE(bopp)) 'opposite_beta_directions'
      else if(topdiff && is.finite(corr) && corr < 0.25) 'different_trait_axes'
      else if(is.finite(corr) && corr < 0) 'antagonistic_multitrait_profile'
      else 'mixed_or_unclear'
    counter <- counter + 1L
    pairs_out[[counter]] <- data.frame(
      region = region,
      term1 = a$term,
      stage1 = a$stage,
      best_trait1 = a$largest_abs_trait,
      beta_effect1 = a$beta_effect,
      term2 = b$term,
      stage2 = b$stage,
      best_trait2 = b$largest_abs_trait,
      beta_effect2 = b$beta_effect,
      n_traits_compared = length(tr2),
      corr_effect_profile = corr,
      fraction_opposite_sign = opp,
      mean_abs_profile_difference = mag,
      top_trait_differs = topdiff,
      beta_directions_opposite = isTRUE(bopp),
      competition_score = score,
      competition_label = label,
      stringsAsFactors = FALSE
    )
  }
}
pairs_out <- if(length(pairs_out)) do.call(rbind, pairs_out) else data.frame()

## Region-level ranking.
region_rank <- lapply(split(pairs_out, pairs_out$region), function(z){
  z <- z[order(-z$competition_score, abs(z$corr_effect_profile), na.last = TRUE), , drop = FALSE]
  data.frame(region = z$region[1],
             top_pair = paste(z$term1[1], z$term2[1], sep = ' vs '),
             top_label = z$competition_label[1],
             top_competition_score = z$competition_score[1],
             top_fraction_opposite = z$fraction_opposite_sign[1],
             top_corr = z$corr_effect_profile[1],
             n_tested_pairs = nrow(z),
             stringsAsFactors = FALSE)
})
region_rank <- if(length(region_rank)) do.call(rbind, region_rank) else data.frame()

## Add a compact epistasis flag.
epi_sum <- NULL
if(nrow(epi)){
  tmp <- epi[epi$identifiable %in% TRUE, , drop = FALSE]
  if(nrow(tmp)){
    epi_sum <- aggregate(q_BH_all_testable_interactions ~ region, data = tmp, FUN = function(x) min(x, na.rm = TRUE))
    names(epi_sum)[2] <- 'best_interaction_q'
    region_rank <- merge(region_rank, epi_sum, by = 'region', all.x = TRUE)
  }
}

write.csv(term_out,  file.path(out_dir, '31_Term_trait_profile_summary.csv'), row.names = FALSE)
write.csv(pairs_out, file.path(out_dir, '31_Within_region_competition_summary.csv'), row.names = FALSE)
write.csv(region_rank, file.path(out_dir, '31_Region_priority_ranking.csv'), row.names = FALSE)
cat('Competition summaries written in ', out_dir, '\n', sep = '')
