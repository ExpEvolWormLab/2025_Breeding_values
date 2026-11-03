## 58 — Robustness of beta-vs-PC1 to omitting surviving GA replicates
## RECONSTRUCTS original frequency slopes using the corrected 09f GLS model,
## original A00 baseline and the FROZEN script-43 coverage-eligible SNP set.
## It then aggregates the median phase-aligned SNP slope within each R5 cluster,
## exactly as script 43. This is a sensitivity analysis, not independent
## prospective cross-population forecasting or a new selection significance test.
## Run from the project root after scripts 43, 45, and 53b.
options(stringsAsFactors=FALSE)
set.seed(20261102L)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
effect_names <- paste0('beta_',traits,'_NaCl')
out <- 'selection_results/chrV_deep_dive/58_replicate_sensitivity'
dir.create(out,showWarnings=FALSE,recursive=TRUE)
W <- function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na='')
fpath <- function(paths){z<-paths[file.exists(paths)];if(length(z))z[[1L]]else NA_character_}
need <- function(paths,why){p<-fpath(paths);if(is.na(p))stop('Missing ',why,': ',paste(paths,collapse=' ; '),call.=FALSE);message(why,': ',p);p}
read_tab <- function(path)read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
source(need(c('Oct_pipeline_cleaner/43_Unified_config_and_helpers.R',
              '43_Unified_config_and_helpers.R'),'43 settings'),local=TRUE)
cfg <- uf_config
source(need(c('Oct_pipeline_cleaner/09f_drift_functions.R',
              '09f_drift_functions.R'),'09f corrected drift GLS'),local=TRUE)
source(need(c('Oct_pipeline_cleaner/43x_Integrate_rawSYNC_R5_SNPs.R',
              '43x_Integrate_rawSYNC_R5_SNPs.R'),'raw-SYNC integration'),local=TRUE)
phase <- read_tab(need(c(file.path(cfg$prefix,'43_All_SNP_phase_slope_and_inclusion_QC.csv')),
                       'frozen R5 phase membership and SNP eligibility'))
E <- read_tab(need(c(file.path(cfg$prefix,'45_All_R5_clusters_GWAS_evolution_predictions.csv'),
                     '45_All_R5_clusters_GWAS_evolution_predictions.csv'),'original evolutionary predictions'))
Q <- read_tab(need(c(cfg$r5_gwas_file,'11_Clusters_selection_and_GWAS_09f.csv'),
                       'complete R5 GWAS effect catalogue'))
weights <- read_tab(need(c(cfg$weights_file,'13_Gradient_weight_QC.csv'),'fitness weights'))
E <- E[E$scenario=='primary_survivors'&E$reason=='eligible',,drop=FALSE]
if(anyDuplicated(E$group)||anyDuplicated(Q$group))stop('Duplicated R5 group IDs')
if(!all(c('rs','group','phase','passes_coverage','slope_phase_main')%in%names(phase)))
  stop('Script-43 phase table schema mismatch')
if(anyDuplicated(phase$rs))stop('Duplicated SNPs in frozen R5 phase table')
if(!all(c('group','chr','ps','block','slope','pred_full','pred_SB',
          'LOCO_full','LOCO_SB')%in%names(E)))stop('45 prediction schema mismatch')
idx <- match(as.character(E$group),as.character(Q$group))
if(anyNA(idx))stop('An eligible cluster lacks R5 GWAS effects')
Q <- Q[idx,,drop=FALSE]
if(any(as.character(Q$chr)!=as.character(E$chr)))stop('Chromosome mismatch')
if(!all(effect_names%in%names(Q)))stop('Seven-trait GWAS-effect matrix unavailable')
X <- as.matrix(Q[,effect_names,drop=FALSE]);storage.mode(X)<-'double';colnames(X)<-traits
if(any(!is.finite(X)))stop('Nonfinite GWAS effects')
weights <- weights[match(traits,as.character(weights$trait)),,drop=FALSE]
if(anyNA(weights$trait))stop('Missing selection-gradient weight')
w <- as.numeric(weights$weight_full_GWAS_scale)
if(any(!is.finite(w)))stop('Nonfinite gradient weights')
if(max(abs(drop(X%*%w)-E$pred_full))>1e-9 ||
   max(abs(X[,'SB']-E$pred_SB))>1e-9)
  stop('GWAS effect/weight orientation does not reproduce frozen script 45')
if(nrow(E)!=802L)warning('Expected 802 eligible R5 groups; found ',nrow(E))

## Use ONE shared ancestral A00; never treat it as replicated per GA.
p <- uf_pool_with_R5_raw_recovery(cfg)
if(!identical(as.character(p$pop$sample_id),colnames(p$k_alt)) ||
   !identical(dim(p$k_alt),dim(p$depth)))stop('Pool and manifest alignment failed')
status <- read_tab(cfg$status_file)
if(!all(c('replicate','population','status')%in%names(status)))stop('GA status schema mismatch')
status$population <- as.character(status$population)
status$status <- tolower(as.character(status$status))
expected <- c(GA1='survived',GA2='survived',GA3='extinct',GA4='survived')
if(!identical(unname(status$status[match(names(expected),status$population)]),unname(expected)))
  stop('Unexpected GA survival status; verify corrected sample mapping')
pop <- p$pop
pop$population <- status$population[match(as.character(pop$replicate),as.character(status$replicate))]
pop$population[pop$generations==0] <- 'A00'
if(anyNA(pop$population))stop('Unmapped sequencing population')
if(sum(pop$sample_id=='A00')!=1 || sum(pop$generations==0)!=1 ||
   any(pop$sample_id=='A6G140'))stop('Incorrect A00 / A6G140 baseline configuration')
if(!identical(as.character(pop$population[pop$generations>0])%in%c('GA1','GA2','GA3','GA4'),
              rep(TRUE,sum(pop$generations>0))))stop('Unexpected population labels')
W(pop,'58_Corrected_sample_manifest.csv')

phase$group <- as.character(phase$group)
phase$phase <- as.numeric(phase$phase)
selected <- which(phase$passes_coverage %in% TRUE &
                  phase$group%in%as.character(E$group))
if(!length(selected))stop('No frozen eligible SNPs found')
if(any(!phase$phase[selected]%in%c(-1,1)))stop('Non- +/-1 phase in frozen SNP catalogue')
kk <- match(as.character(phase$rs[selected]),as.character(p$sel$ID))
if(anyNA(kk))stop('Previously eligible SNPs missing from integrated Pool-seq: ',
                   paste(head(phase$rs[selected][is.na(kk)]),collapse=', '))
K <- p$k_alt[kk,,drop=FALSE];D <- p$depth[kk,,drop=FALSE]
if(any(!is.finite(phase$slope_phase_main[selected])))
  stop('Frozen eligible SNP has missing original full-scenario slope')
## Each two-replicate analysis keeps the 43 GWAS-defined R5 member universe
## and original coverage eligibility; failures to fit are explicit missing data.
subsets <- list(
  'GA1_GA2_GA4_original'=c('GA1','GA2','GA4'),
  'GA1_GA2_no_GA4'=c('GA1','GA2'),
  'GA1_GA4_no_GA2'=c('GA1','GA4'),
  'GA2_GA4_no_GA1'=c('GA2','GA4'))

## Verify original 43 slope orientation against the exact recovered counts.
set.seed(20261102L)
chk <- unique(round(seq(1,length(selected),length.out=min(120L,length(selected)))))
orig_checks <- vapply(chk,function(j){
  keep <- pop$population%in%c('A00','GA1','GA2','GA4')
  z <- fit_drift_site_null(as.numeric(K[j,keep]),as.numeric(D[j,keep]),
                           pop[keep,,drop=FALSE],Ne=cfg$Ne,pool_individuals=NA_real_)
  if(!is.finite(z$slope))return(NA_real_)
  z$slope*phase$phase[selected[j]] - phase$slope_phase_main[selected[j]]
},numeric(1))
if(any(!is.finite(orig_checks))||max(abs(orig_checks))>1e-8)
  stop('New 09f GLS does not reproduce original 43 slopes; check scripts/inputs. Max diff: ',
       suppressWarnings(max(abs(orig_checks),na.rm=TRUE)))

## Dynamic GLS for two-replicate scenarios. Do NOT use endpoint frequency
## differences or the frozen three-replicate estimates as subset outcomes.
SNP <- data.frame(rs=as.character(phase$rs[selected]),
                  group=phase$group[selected],
                  phase=phase$phase[selected],
                  full_original=phase$slope_phase_main[selected],
                  stringsAsFactors=FALSE)
fit_qc <- list()
for(nm in names(subsets)[-1]){
  keep <- pop$population%in%c('A00',subsets[[nm]])
  pp <- pop[keep,,drop=FALSE]
  sl <- rep(NA_real_,nrow(K))
  for(i in seq_len(nrow(K))){
    ff <- fit_drift_site_null(as.numeric(K[i,keep]),as.numeric(D[i,keep]),pp,
                              Ne=cfg$Ne,pool_individuals=NA_real_)
    if(is.finite(ff$slope))sl[i] <- ff$slope*SNP$phase[i]
    if(i %% 4000L==0L)cat(nm,': fitted ',i,' / ',nrow(K),' SNPs\n')
  }
  SNP[[nm]] <- sl
  fit_qc[[length(fit_qc)+1L]] <- data.frame(scenario=nm,
    n_original_eligible_SNPs=nrow(SNP),n_subset_fit_SNPs=sum(is.finite(sl)),
    n_subset_missing_SNPs=sum(!is.finite(sl)))
}
W(SNP,'58_SNP_gls_phase_aligned_replicate_subset_slopes.csv')
W(do.call(rbind,fit_qc),'58_SNP_fitting_QC.csv')

## Cluster aggregation: median of ORIGINAL coverage-qualified member SNPs.
## No new trait significance filtering and no choosing sites by response.
bygroup <- split(seq_len(nrow(SNP)),SNP$group)
scenarios <- names(subsets)
Y <- data.frame(group=as.character(E$group),chr=as.character(E$chr),
                ps=as.numeric(E$ps),block=as.character(E$block),
                n_original_eligible=as.integer(E$n_eligible),stringsAsFactors=FALSE)
for(nm in scenarios){
  cn <- if(nm=='GA1_GA2_GA4_original')'full_original' else nm
  out_y <- rep(NA_real_,nrow(E)); num <- integer(nrow(E))
  for(i in seq_len(nrow(E))){
    pos <- bygroup[[as.character(E$group[i])]]
    if(is.null(pos))next
    v <- SNP[[cn]][pos];v <- v[is.finite(v)]
    num[i] <- length(v)
    if(length(v))out_y[i]<-median(v)
  }
  Y[[nm]]<-out_y
  Y[[paste0('n_snps_',nm)]]<-num
}
if(max(abs(Y$GA1_GA2_GA4_original-E$slope))>1e-8)
  stop('Frozen 43 median SNP slope does not reproduce script 45 evolutionary response')
W(Y,'58_R5_cluster_outcomes_fixed_membership.csv')

## Refit y ~ score calibration on five training chromosomes for every subset.
## PC1 itself comes only from TRAIN-chromosome GWAS-effect matrix.
fit_cal <- function(xtr,xte,ytr){
  ok<-is.finite(xtr)&is.finite(ytr)
  if(sum(ok)<30L || sd(xtr[ok])<1e-12)return(rep(NA_real_,length(xte)))
  fit <- lm.fit(cbind(1,xtr[ok]),ytr[ok])
  if(any(!is.finite(fit$coefficients)))return(rep(NA_real_,length(xte)))
  drop(cbind(1,xte)%*%fit$coefficients)
}
B <- drop(X%*%w)
models <- c('beta','SB','PC1_cov','PC1_corr')
run_scenario <- function(outcome,name){
  pred <- data.frame(group=E$group,chr=E$chr,ps=E$ps,block=E$block,
                     scenario=name,observed_slope=outcome,
                     beta=NA_real_,SB=NA_real_,PC1_cov=NA_real_,PC1_corr=NA_real_)
  for(chr_holdout in unique(as.character(E$chr))){
    tr <- as.character(E$chr)!=chr_holdout
    te <- !tr
    pred$beta[te]<-fit_cal(B[tr],B[te],outcome[tr])
    pred$SB[te]<-fit_cal(X[tr,'SB'],X[te,'SB'],outcome[tr])
    for(kind in c('cov','corr')){
      pc <- prcomp(X[tr,,drop=FALSE],center=TRUE,scale.=(kind=='corr'))
      ## The PCA calibration absorbs PC sign; no orientation is chosen from y.
      xs <- predict(pc,newdata=as.data.frame(X))[,1L]
      pred[te,paste0('PC1_',kind)] <- fit_cal(xs[tr],xs[te],outcome[tr])
    }
  }
  perf <- do.call(rbind,lapply(models,function(mod){
    ok <- is.finite(pred$observed_slope)&is.finite(pred[[mod]])
    y <- pred$observed_slope[ok];h <- pred[[mod]][ok]
    denom <- sum((y-mean(y))^2)
    data.frame(scenario=name,model=mod,n_eligible_catalogue=nrow(E),
      n_evaluated=sum(ok),n_chromosomes=length(unique(pred$chr[ok])),
      R2=if(denom>0)1-sum((y-h)^2)/denom else NA_real_,
      MSE=if(length(y))mean((y-h)^2)else NA_real_,
      response_sd=if(length(y)>1)sd(y)else NA_real_,
      correlation=if(length(y)>2)cor(y,h)else NA_real_,
      sign_accuracy=if(length(y))mean(sign(y)==sign(h))else NA_real_)
  }))
  list(pred=pred,metrics=perf)
}
ans <- lapply(scenarios,function(nm)run_scenario(Y[[nm]],nm));names(ans)<-scenarios
FULL <- ans[['GA1_GA2_GA4_original']]$pred
if(max(abs(FULL$beta-E$LOCO_full))>1e-9 ||
   max(abs(FULL$SB-E$LOCO_SB))>1e-9)
  stop('Rebuilt full beta/SB LOCO predictions disagree with script 45')
Pred <- do.call(rbind,lapply(ans,`[[`,'pred')); Metrics<-do.call(rbind,lapply(ans,`[[`,'metrics'))
W(Pred,'58_LOCO_predictions_by_population_subset.csv')
W(Metrics,'58_LOCO_metrics_by_population_subset.csv')

## Fair paired block comparison within each scenario; fitted models remain
## fixed within each bootstrap draw, so intervals are DESCRIPTIVE.
boot <- list();set.seed(20261102L)
for(nm in scenarios){
  d <- ans[[nm]]$pred
  for(comp in c('PC1_cov','SB')){
    ok <- is.finite(d$observed_slope)&is.finite(d$beta)&is.finite(d[[comp]])
    z <- d[ok,,drop=FALSE]
    delta <- (z$observed_slope-z[[comp]])^2-(z$observed_slope-z$beta)^2
    ix <- split(seq_len(nrow(z)),as.character(z$block))
    if(length(ix)<8L)next
    draws <- replicate(3000L,{
      jj<-sample(names(ix),length(ix),replace=TRUE)
      mean(delta[unlist(ix[jj],use.names=FALSE)])
    })
    boot[[length(boot)+1L]]<-data.frame(scenario=nm,comparator=comp,
      n_blocks=length(ix),n_clusters=nrow(z),delta_MSE_beta_advantage=mean(delta),
      CI_low=as.numeric(quantile(draws,.025)),CI_high=as.numeric(quantile(draws,.975)),
      bootstrap_fraction_positive=mean(draws>0))
  }
}
if(length(boot))W(do.call(rbind,boot),'58_Exploratory_paired_block_bootstrap.csv')

## Region-specific direction diagnostics; do not count clustered SNPs as
## independent evidence or evaluate significance using n=4.
focal_groups <- c('848','851','853','854')
focal <- Pred[as.character(Pred$group)%in%focal_groups,,drop=FALSE]
focal$beta_direction_correct <- with(focal,sign(observed_slope)==sign(beta))
focal$PC1_direction_correct <- with(focal,sign(observed_slope)==sign(PC1_cov))
focal$delta_SE_beta_advantage <- with(focal,(observed_slope-PC1_cov)^2-(observed_slope-beta)^2)
W(focal,'58_ChrV_four_focal_groups_population_omission.csv')
## Include original exact-tag trajectories as a separate raw-count visualization.
exact_path<-fpath(c('selection_results/chrV_deep_dive/56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv'))
if(!is.na(exact_path)){
  TT<-read_tab(exact_path)
  a<-TT[as.character(TT$group)%in%focal_groups & TT$sample_id=='A00',]
  a<-a[!duplicated(a$SNP_ID),]
  if(nrow(a))W(a[,intersect(c('group','SNP_ID','depth','ALT_reads','ALT_frequency','phase'),names(a)),drop=FALSE],
               '58_Focal_exact_tag_A00_QC.csv')
}
W(data.frame(n_catalogue=nrow(E),
             n_frozen_SNP=nrow(SNP),n_A00=sum(pop$sample_id=='A00'),
             n_primary_populations=3L,
             max_full_GLSlopes_check=max(abs(orig_checks)),
             max_full_cluster_median_diff=max(abs(Y$GA1_GA2_GA4_original-E$slope)),
             max_full_beta_LOCO_diff=max(abs(FULL$beta-E$LOCO_full)),
             max_full_SB_LOCO_diff=max(abs(FULL$SB-E$LOCO_SB))),
  '58_Immutable_input_and_reproduction_QC.csv')
writeLines(c('Sensitivity uses EXACT original script-43 eligible SNPs, phase, and R5 groups.',
  'Each two-population outcome is a new 09f neutral-covariance GLS slope.',
  'All population subsets reuse one shared A00 baseline, not independent baselines.',
  'PCA and effect-outcome regression are trained on the other five chromosomes.',
  'The target outcome and fitted calibration change with population omission.',
  'This is a robustness analysis of a genome-wide signal; TWO evolved replicates do not',
  'constitute robust independent prospective validation.',
  'Neither the 802 clusters nor the four chromosome-V groups are independent selected mutations.',
  'Do not use replicate omission to recalibrate candidate selection p-values.',
  'GA3 remains excluded because it is confirmed extinct.'),
  file.path(out,'58_Interpretation_cautions.txt'))
cat('Script 58 completed: ',out,'\n')
print(Metrics[,c('scenario','model','n_evaluated','R2','sign_accuracy')],row.names=FALSE)
