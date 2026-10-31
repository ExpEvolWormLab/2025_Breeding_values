## 20_Survivor_only_drift_and_GWAS.R
## Matched sensitivity: exclude non-surviving/unknown GA lines from Pool-seq
## RESPONSE estimation, not from the independent ancestral RIL GWAS.
## Uses the corrected script-08 matrix, script-09 ALT orientation and the
## exact same SNP universe as the full corrected script 09f.
## Explicitly compares the full and survivor-only responses on IDENTICAL QTLs.
## Does not infer extinction from a missing last sequencing time point.

suppressPackageStartupMessages(library(parallel))
helper <- if(file.exists('09f_drift_functions.R')) '09f_drift_functions.R' else
  'drift_nullfit_update/09f_drift_functions.R'
source(helper)

status_file <- 'Oct_pipeline_cleaner/20_GA_population_status.csv'
pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
orientation_file <- 'selection_results/09_Selection_SNP_results.csv'
full_selection_file <- 'selection_results/drift_oneG0_nullfit/09f_SNP_drift_Ne1000.csv'
full_cluster_file <- 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv'
full_predictors_file <- paste0('selection_results/drift_oneG0_nullfit/R5/',
             'multivariate_prediction/13_Cluster_predictors_and_response.csv')
phase_file <- 'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv'
out_dir <- 'selection_results/survivors_only_Ne1000'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

Ne <- 1000L
exclude_G0 <- 'A6G140' # Same baseline policy as corrected full analysis.
chunk_size <- 500L
mc_cores <- 6L
common_range <- c(.10,.90)
n_boot <- 1000L
n_neutral_loci <- 40L
n_sim_per_locus <- 250L
set.seed(20261029)
if(.Platform$OS.type=='windows') mc_cores<-1L

status <- read.csv(status_file,stringsAsFactors=FALSE)
req <- c('replicate','population','status')
if(!all(req %in% names(status))) stop('Status file needs: ',paste(req,collapse=', '))
if(anyDuplicated(as.character(status$replicate)))stop('Duplicate replicate status')
status$status <- tolower(trimws(status$status))
if(!all(status$status %in% c('survived','extinct','unknown')))
  stop('Status must be survived, extinct or unknown')
if(any(status$status=='unknown'))
  warning('Unknown-status population(s) are excluded from this CONSERVATIVE ',
          'survivor-only sensitivity, but MUST NOT be described as extinct: ',
          paste(status$population[status$status=='unknown'],collapse=', '))

## Input reader checks the original 13-column alignment and removes A6G140.
p <- read_drift_inputs(pool_file,orientation_file,exclude_G0=exclude_G0)
pop_full <- p$pop
reps <- unique(as.character(pop_full$replicate[pop_full$generations>0]))
if(!setequal(reps,as.character(status$replicate)))
 stop('Status table must list exactly the four evolving replicate IDs: ',
      paste(reps,collapse=', '))
keep_reps <- as.character(status$replicate[status$status=='survived'])
if(length(keep_reps)<2L) stop('At least two surviving evolved populations required')
keep <- pop_full$generations==0 | as.character(pop_full$replicate)%in%keep_reps
pop <- pop_full[keep,,drop=FALSE]
K <- p$k_alt[,keep,drop=FALSE]
D <- p$depth[,keep,drop=FALSE]
sel <- p$sel
write.csv(status,file.path(out_dir,'20_GA_population_status_used.csv'),row.names=FALSE)
write.csv(pop,file.path(out_dir,'20_Retained_sample_manifest.csv'),row.names=FALSE)
cat('Included evolved populations: ',paste(status$population[status$status=='survived'],collapse=', '),'\n')
cat('Sample counts: old ',nrow(pop_full),'; survivor-only ',nrow(pop),'\n')

## Genome-wide GLS screening with unchanged allele orientation, Ne and
## ancestral frequency definition. The original 09f file remains untouched.
chunks <- split(seq_len(nrow(K)), ceiling(seq_len(nrow(K))/chunk_size))
ans <- parallel::mclapply(chunks,function(ii){
 t(vapply(ii,function(i){
  z <- fit_drift_site_null(K[i,],D[i,],pop,Ne,
                            pool_individuals=NA_real_)
  c(p0=z$p0, p0_null=z$p0_null, slope=z$slope,
    SE=z$SE,z=z$z,p=z$p,n_pools=z$n_pools,
    n_G0=z$n_G0,n_evolved_reps=z$n_evolved_reps)
 },FUN.VALUE=c(p0=0,p0_null=0,slope=0,SE=0,z=0,p=0,
               n_pools=0,n_G0=0,n_evolved_reps=0)))
},mc.cores=mc_cores,mc.preschedule=TRUE)
ans <- do.call(rbind,ans)
if(nrow(ans)!=nrow(sel))stop('SNP order mismatch in survivor refit')
res <- data.frame(rs=sel$rs,ID=sel$ID,CHROM=sel$CHROM,POS=sel$POS,
 orientation=sel$orientation,p0_hat=ans[,'p0'],
 p0_null=ans[,'p0_null'],slope_freq=ans[,'slope'],
 SE_freq=ans[,'SE'],Z_drift=ans[,'z'],p_drift=ans[,'p'],
 n_pools=as.integer(ans[,'n_pools']),n_G0=as.integer(ans[,'n_G0']),
 n_evolved_reps=as.integer(ans[,'n_evolved_reps']),
 stringsAsFactors=FALSE)
res$common_G0 <- is.finite(res$p0_hat) & res$p0_hat>=common_range[1] &
 res$p0_hat<=common_range[2]
res$p_drift_BH_all <- p.adjust(res$p_drift,'BH')
res$p_drift_BY_all <- p.adjust(res$p_drift,'BY')
res$p_drift_BH_common <- NA_real_
ix <- which(res$common_G0 & is.finite(res$p_drift))
res$p_drift_BH_common[ix] <- p.adjust(res$p_drift[ix],'BH')
res$Ne <- Ne
res$excluded_G0 <- exclude_G0
res$retained_G0 <- pop$sample_id[pop$generations==0]
res$survival_policy <- 'confirmed_survivors_only;unknown_excluded'
write.csv(res,file.path(out_dir,'20_SNP_drift_Ne1000_survivors.csv'),row.names=FALSE)

full <- read.csv(full_selection_file,stringsAsFactors=FALSE)
if(!identical(as.character(res$ID),as.character(full$ID)))
 stop('Full-vs-survivor SNP identities/order are inconsistent')
snp_compare <- data.frame(ID=res$ID,CHROM=res$CHROM,POS=res$POS,
 p0_hat=res$p0_hat,common_G0=res$common_G0,
 full_slope_freq=full$slope_freq,survivor_slope_freq=res$slope_freq,
 full_p=full$p_drift,survivor_p=res$p_drift,
 full_BH_common=full$p_drift_BH_common,
 survivor_BH_common=res$p_drift_BH_common,
 stringsAsFactors=FALSE)
write.csv(snp_compare,file.path(out_dir,'20_Compared_SNP_slopes.csv'),row.names=FALSE)

## Phase-align the survivor-only SNP trajectories to the EXISTING R5 cluster
## dosages and summarize the same SNPs and same LD groups as before.
phase <- read.csv(phase_file,stringsAsFactors=FALSE)
if(!all(c('rs','group','phase')%in%names(phase)))
 stop('Phase file needs rs/group/phase')
if(anyDuplicated(phase$rs)||anyNA(phase$phase)||any(!phase$phase%in%c(-1,1)))
 stop('Invalid phase data: expected unique SNPs with phases +/-1')
i <- match(as.character(phase$rs),as.character(res$ID))
phase$survivor_pool_covered <- !is.na(i)
phase$survivor_common <- !is.na(res$common_G0[i]) & res$common_G0[i]
phase$survivor_slope <- phase$phase * res$slope_freq[i]
phase$survivor_slope_common <- ifelse(phase$survivor_common,
                                     phase$survivor_slope,NA_real_)
phase$survivor_q <- res$p_drift_BH_common[i]
phase$survivor_selected <- phase$survivor_common &
  is.finite(phase$survivor_q) & phase$survivor_q<.01
write.csv(phase,file.path(out_dir,'20_SNP_phase_oriented_survivor_slopes.csv'),
          row.names=FALSE)

safe_median <- function(z) {z<-z[is.finite(z)];if(length(z))median(z)else NA_real_}
S <- do.call(rbind,lapply(split(phase,as.character(phase$group)),function(x){
  z<-x$survivor_slope_common
  data.frame(group=as.character(x$group[1]),
   n_survivor_common=sum(is.finite(z)),
   slope_survivor_median_common=safe_median(z),
   n_selected_survivor=sum(x$survivor_selected,na.rm=TRUE),
   any_selected_survivor=any(x$survivor_selected,na.rm=TRUE),
   stringsAsFactors=FALSE)
}))
row.names(S)<-NULL
write.csv(S,file.path(out_dir,'20_Survivor_cluster_slope_summary.csv'),row.names=FALSE)
write.csv(data.frame(
  n_tested_genomewide=sum(is.finite(res$p_drift)),
  n_common_genomewide=sum(res$common_G0 & is.finite(res$p_drift)),
  n_selected_common_BH01_genomewide=sum(res$p_drift_BH_common<.01,na.rm=TRUE),
  n_selected_common_BH05_genomewide=sum(res$p_drift_BH_common<.05,na.rm=TRUE),
  n_GWAS_clusters_with_common_slope=sum(S$n_survivor_common>0),
  n_GWAS_clusters_with_BH01_SNP=sum(S$any_selected_survivor),
  included_replicates=paste(keep_reps,collapse=';'),
  status_unknown_excluded=paste(status$population[status$status=='unknown'],collapse=';'),
  stringsAsFactors=FALSE),
  file.path(out_dir,'20_Survivor_overlap_summary.csv'),row.names=FALSE)

## Previously selected full-data cluster sets and independently mapped
## predictors are fixed. This is essential to a FAIR survivorship comparison.
base_qtl <- read.csv(full_cluster_file,stringsAsFactors=FALSE)
preds <- read.csv(full_predictors_file,stringsAsFactors=FALSE)
if(anyDuplicated(S$group)||anyDuplicated(base_qtl$group)||anyDuplicated(preds$group))
 stop('Duplicated QTL group identifiers')
use <- match(as.character(preds$group),as.character(base_qtl$group))
if(anyNA(use))stop('Original predictor groups missing from GWAS cluster file')
old <- base_qtl[use,,drop=FALSE]
j <- match(as.character(preds$group),as.character(S$group))
if(anyNA(j))stop('Survivor cluster summary does not include every predictor cluster')
dat <- preds
## This is the same response as used to create the predictors in script 13.
dat$response_full <- preds$response
dat$response_survivors <- S$slope_survivor_median_common[j]
dat$n_survivor_common <- S$n_survivor_common[j]
dat$n_full_common <- old$n_slope_common
if('chromosome'%in%names(dat)) dat$chromosome <- as.character(dat$chromosome)
if(!'block'%in%names(dat))
 dat$block <- paste(dat$chromosome,floor(dat$position/2e6),sep=':')
write.csv(dat,file.path(out_dir,'20_Matched_cluster_predictors_and_responses.csv'),
          row.names=FALSE)

model_cols <- c(full_mapped='pred_full_GWAS',
 full_reconstructed='pred_full_reconstructed',
 three_trait='pred_three_trait',
 SF='beta_SF_NaCl',SB='beta_SB_NaCl',FS='beta_FS_NaCl',
 FB='beta_FB_NaCl',BS='beta_BS_NaCl',BF='beta_BF_NaCl',
 area.F='beta_area.F_NaCl')
if(!all(model_cols %in% names(dat)))
 stop('Script 13 predictor CSV missing expected columns: ',
      paste(setdiff(model_cols,names(dat)),collapse=', '))
valid <- is.finite(dat$response_full) & is.finite(dat$response_survivors) &
 dat$n_survivor_common>0 & dat$n_full_common>0 &
 rowSums(!is.finite(as.matrix(dat[,model_cols,drop=FALSE])))==0
z <- dat[valid,,drop=FALSE]
if(nrow(z)<30)stop('Fewer than 30 matched clusters in survival sensitivity')
safe_cor <- function(a,b)if(length(a)>=5&&sd(a)>0&&sd(b)>0)
  cor(a,b) else NA_real_
assoc <- do.call(rbind,lapply(names(model_cols),function(m){
 x <- z[[model_cols[[m]]]]
 data.frame(model=m,n_clusters=nrow(z),
 full_r=safe_cor(x,z$response_full),
 survivor_r=safe_cor(x,z$response_survivors),
 delta_r_survivor_minus_full=safe_cor(x,z$response_survivors)-
                         safe_cor(x,z$response_full),
 stringsAsFactors=FALSE)
}))
write.csv(assoc,file.path(out_dir,'20_Matched_cluster_associations.csv'),row.names=FALSE)

## Leave-one-CHROMOSOME-out regression calibration: two independent response
## definitions, same held-out groups. Do not reuse full-response fitted slopes.
cv_predict <- function(y,x,chr) {
 pred<-rep(NA_real_,length(y))
 for(c in unique(chr)) {
  tr<-chr!=c;te<-!tr
  if(sum(tr)<20||sum(te)<2)next
  if(sd(x[tr])<1e-12){pred[te]<-mean(y[tr]);next}
  fit<-lm(y[tr]~x[tr])
  pred[te]<-coef(fit)[1]+coef(fit)[2]*x[te]
 }
 pred
}
cv <- list()
for (m in names(model_cols)) {
 x<-z[[model_cols[[m]]]]
 for (response_type in c('full','survivors')) {
  y<-z[[paste0('response_',if(response_type=='survivors') 'survivors' else 'full')]]
  pred<-cv_predict(y,x,z$chromosome)
  ok<-is.finite(pred)&is.finite(y)
  cv[[paste(m,response_type,sep=':')]]<-data.frame(model=m,
   response_set=response_type,n_clusters=sum(ok),
   R2_LOCO=1-sum((y[ok]-pred[ok])^2)/sum((y[ok]-mean(y[ok]))^2),
   r_LOCO=safe_cor(y[ok],pred[ok]),
   MSE=mean((y[ok]-pred[ok])^2),stringsAsFactors=FALSE)
 }
}
cv<-do.call(rbind,cv)
write.csv(cv,file.path(out_dir,'20_Matched_LOCO_prediction.csv'),row.names=FALSE)

pdf(file.path(out_dir,'20_Survivor_sensitivity_figure.pdf'),width=10,height=4.5)
op<-par(mfrow=c(1,2),mar=c(4.5,4.7,3,1))
plot(z$response_full,z$response_survivors,pch=16,cex=.65,
     col=adjustcolor('steelblue4',alpha.f=.45),bty='n',
     xlab='Allele-frequency slope: all GA populations',
     ylab='Allele-frequency slope: retained GA populations',
     main='Matched R5 cluster responses')
abline(0,1,col='grey50',lty=2)
plot(z$pred_full_reconstructed,z$response_survivors,pch=16,cex=.65,
     col=adjustcolor('steelblue4',alpha.f=.4),bty='n',
     xlab='Reconstructed full fitness-index GWAS effect',
     ylab='Survivor-only frequency slope',
     main='Fitness prediction: retained populations')
abline(lm(response_survivors~pred_full_reconstructed,data=z),col='firebrick3',lwd=2)
legend('topleft',bty='n',cex=.8,
       legend=paste0('r = ',signif(safe_cor(z$pred_full_reconstructed,
                                         z$response_survivors),3)))
par(op);dev.off()

## Resample the same 2-Mb blocks with both responses AND all predictors
## kept together. Bootstrap intervals are genomic sensitivity, not iid-SNP CI.
blocks <- split(seq_len(nrow(z)),z$block)
set.seed(20261029)
boot <- vector('list',n_boot)
for(b in seq_len(n_boot)) {
 ii<-unlist(blocks[sample(names(blocks),length(blocks),replace=TRUE)],
            use.names=FALSE)
 if(length(ii)<10L)next
 y1<-z$response_full[ii]; y2<-z$response_survivors[ii]
 x<-z$pred_full_reconstructed[ii]
 x3<-z$pred_three_trait[ii]
 boot[[b]]<-data.frame(
  delta_r_full_index=safe_cor(x,y2)-safe_cor(x,y1),
  delta_r_three_index=safe_cor(x3,y2)-safe_cor(x3,y1),
  full_minus_three_survivor=safe_cor(x,y2)-safe_cor(x3,y2),
  full_minus_three_all=safe_cor(x,y1)-safe_cor(x3,y1))
}
boot<-do.call(rbind,Filter(Negate(is.null),boot))
write.csv(boot,file.path(out_dir,'20_Paired_2Mb_block_bootstrap_draws.csv'),row.names=FALSE)
ci <- data.frame(metric=names(boot),estimate=NA_real_,low95=NA_real_,high95=NA_real_)
ci$estimate <- c(
  safe_cor(z$pred_full_reconstructed,z$response_survivors)-
    safe_cor(z$pred_full_reconstructed,z$response_full),
  safe_cor(z$pred_three_trait,z$response_survivors)-
    safe_cor(z$pred_three_trait,z$response_full),
  safe_cor(z$pred_full_reconstructed,z$response_survivors)-
    safe_cor(z$pred_three_trait,z$response_survivors),
  safe_cor(z$pred_full_reconstructed,z$response_full)-
    safe_cor(z$pred_three_trait,z$response_full))
for(k in seq_along(ci$metric)) {
 q<-boot[[k]]
 ci$low95[k]<-quantile(q,.025,na.rm=TRUE)
 ci$high95[k]<-quantile(q,.975,na.rm=TRUE)
}
write.csv(ci,file.path(out_dir,'20_Paired_2Mb_block_bootstrap_summary.csv'),row.names=FALSE)

## Optional conditional neutral calibration of the SURVIVOR sampling design.
## 40 loci x 250 sims is a diagnostic of gross miscalibration only, NOT
## validation of genome-wide FDR tails. Initial SNP frequencies are refit per
## simulation using the same 09f null GLS statistic.
valid_idx<-which(res$common_G0 & is.finite(res$p_drift) &
                 res$p0_hat>.05 & res$p0_hat<.95)
if(n_neutral_loci>0 && n_sim_per_locus>0 && length(valid_idx)) {
 fbins<-cut(res$p0_hat[valid_idx],breaks=c(.05,.1,.25,.5,.75,.9,.95),
            include.lowest=TRUE)
 bybin<-split(valid_idx,fbins,drop=TRUE)
 chosen<-unlist(lapply(bybin,function(ii)sample(ii,min(length(ii),
          ceiling(n_neutral_loci/length(bybin))))),use.names=FALSE)
 chosen<-head(chosen,n_neutral_loci)
 neutral<-list()
 for(i in chosen) {
  fit<-fit_drift_site_null(K[i,],D[i,],pop,Ne)
  if(!is.finite(fit$z))next
  sim<-simulate_wf_site(fit$p0,D[i,],pop,Ne,n_sim=n_sim_per_locus)
  ss<-drop(sim%*%fit$weights)
  pnull<-drop(sim%*%fit$null_weights)
  good<-is.finite(K[i,])&is.finite(D[i,])&D[i,]>0
  lower<-.5/(sum(D[i,good])+1)
  pnull<-pmax(lower,pmin(1-lower,pnull))
  zs<-ss/(fit$slope_scale*sqrt(pnull*(1-pnull)))
  pp<-2*pnorm(abs(zs),lower.tail=FALSE)
  neutral[[length(neutral)+1L]]<-data.frame(ID=sel$ID[i],
   p0_hat=fit$p0,Ne=Ne,n_sim=length(pp),
   frac_p_lt_005=mean(pp<.05),frac_p_lt_001=mean(pp<.01),
   frac_p_lt_0001=mean(pp<.001))
 }
 if(length(neutral)) {
  nullout<-do.call(rbind,neutral)
  write.csv(nullout,file.path(out_dir,'20_Survivor_neutral_pilot.csv'),row.names=FALSE)
  print(colMeans(nullout[,c('frac_p_lt_005','frac_p_lt_001',
                            'frac_p_lt_0001')],na.rm=TRUE))
 }
}

cat('\nCompleted survivor-only sensitivity in: ',out_dir,'\n',sep='')
cat('Matched QTL clusters:',nrow(z),'\n')
print(assoc,row.names=FALSE)
cat('IMPORTANT: missing G50 alone is not proof of extinction. This analysis ',
    'tests survivors/provisionally complete GA trajectories. It conditions ',
    'on post-treatment survival, so it is a SENSITIVITY, not a replacement ',
    'for the all-population response.\n')
cat('Genome-wide SNP discovery p-values remain provisional; neutral pilot ',
    'does not calibrate extreme BH tails.\n')
