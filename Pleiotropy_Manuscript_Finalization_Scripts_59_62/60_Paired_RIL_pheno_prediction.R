# 60_Paired_RIL_pheno_prediction.R
# Compare held-out (5-fold random-RIL) GBLUP vs fixed 75-marker models on
# EXACTLY matched individuals/folds. Legacy in-sample score is deliberately
# excluded from inferential tests.
# Uses script 51 RIL-level predictions; optionally joins script 40 grouping.
# The 75-marker K was fixed in the earlier analysis, not optimized here.
# Only base R required. Run from your project root.

CONFIG <- list(
 prediction='phenotype_prediction/three_panel_comparison/51_Complete_three_method_RIL_predictions.csv',
 GBLUP='phenotype_prediction/GBLUP_rrBLUP/40_Out_of_fold_GBLUP_predictions.csv',
 out='phenotype_prediction/manuscript_finalization/60_paired_CV',
 n_boot=4000L,seed=20261103L
)
resolve <- function(x,optional=FALSE){
 for(f in c(x,basename(x)))if(file.exists(f))return(f)
 if(optional)return(NA_character_)
 stop('Missing ',x,'. Run script 51 or supply the CSV.')
}
f<-resolve(CONFIG$prediction)
z<-read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
need<-c('pop_label','fold','trait','observed','predicted_top75','predicted_GBLUP')
if(length(setdiff(need,names(z))))stop('Script-51 file missing: ',paste(setdiff(need,names(z)),collapse=','))
if('scheme'%in%names(z)) z<-z[z$scheme=='random_RIL',,drop=FALSE]
if(anyDuplicated(paste(z$pop_label,z$trait,sep='::')))stop('Repeated RIL/trait')
if(!all(z$fold%in%1:5))stop('Unexpected fold labels')
if(!'source_population'%in%names(z))z$source_population<-NA_character_
# Validate source-population grouping against individual-level script 40.
f40 <- resolve(CONFIG$GBLUP,optional=TRUE)
if(!is.na(f40)){
  s<-read.csv(f40,stringsAsFactors=FALSE,check.names=FALSE)
  if(all(c('pop_label','trait','fold','scheme','source_population','observed','predicted_GBLUP')%in%names(s))){
    s<-s[s$scheme=='random_RIL',,drop=FALSE]
    if(anyDuplicated(paste(s$pop_label,s$trait,sep='::')))stop('Duplicate script-40 keys')
    key<-paste(z$pop_label,z$trait,sep='::'); j<-match(key,paste(s$pop_label,s$trait,sep='::'))
    if(anyNA(j))stop('Some script-51 individuals absent from script 40')
    if(any(z$fold!=s$fold[j]))stop('Fold mismatch with script 40')
    ok <- is.finite(z$observed)&is.finite(s$observed[j])
    if(any(ok)&max(abs(z$observed[ok]-s$observed[j[ok]]))>1e-6)stop('Phenotype mismatch')
    ok <- is.finite(z$predicted_GBLUP)&is.finite(s$predicted_GBLUP[j])
    if(any(ok)&max(abs(z$predicted_GBLUP[ok]-s$predicted_GBLUP[j[ok]]))>1e-6)stop('GBLUP mismatch')
    z$source_population<-as.character(s$source_population[j])
  }
}
set.seed(CONFIG$seed)
dir.create(CONFIG$out,recursive=TRUE,showWarnings=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
results<-list();draws_out<-list();preds<-list()
for(trait in traits){
 d<-z[z$trait==trait,,drop=FALSE]
 d<-d[is.finite(d$observed)&is.finite(d$predicted_top75)&is.finite(d$predicted_GBLUP),,drop=FALSE]
 if(nrow(d)<50L)stop('Insufficient matched predictions for ',trait)
 # All resampling within each trait uses the exact SAME rows for both models.
 y<-d$observed; A<-d$predicted_GBLUP; B<-d$predicted_top75
 if(anyDuplicated(d$pop_label))stop('Duplicate RILs for ',trait)
 groups<-as.character(d$source_population)
 use_group <- all(!is.na(groups)&nzchar(groups)&groups!='unknown')&&
    length(unique(groups))>=8L
 units<-if(use_group)split(seq_len(nrow(d)),groups)else as.list(seq_len(nrow(d)))
 if(length(units)<20L)warning('Few resampling units for ',trait)
 calc <- function(ii){
   yy<-y[ii];a<-A[ii];b<-B[ii]
   denom<-sum((yy-mean(yy))^2)
   c(R2_GBLUP=if(denom>0)1-sum((yy-a)^2)/denom else NA_real_,
     R2_K75=if(denom>0)1-sum((yy-b)^2)/denom else NA_real_,
     MSE_GBLUP=mean((yy-a)^2),MSE_K75=mean((yy-b)^2),
     delta_MSE_K75_minus_GBLUP=mean((yy-b)^2-(yy-a)^2))
 }
 obs<-calc(seq_len(nrow(d)))
 boot<-replicate(CONFIG$n_boot,{
   choice<-sample.int(length(units),length(units),replace=TRUE)
   ix<-unlist(units[choice],use.names=FALSE)
   calc(ix)
 })
 ci<-t(apply(boot,1,function(v)unname(quantile(v,c(.025,.975),na.rm=TRUE))))
 results[[trait]]<-data.frame(trait=trait,n_RIL=nrow(d),
   n_resampling_units=length(units),
   resampling_unit=if(use_group)'source_population' else 'RIL_only',
   metric=names(obs),estimate=as.numeric(obs),lo95=ci[,1],hi95=ci[,2],
   interpretation='paired bootstrap of held-out predictions; conditional on original model fits',
   stringsAsFactors=FALSE)
 draws_out[[trait]]<-data.frame(trait=trait,iteration=seq_len(CONFIG$n_boot),
   delta_MSE_K75_minus_GBLUP=boot['delta_MSE_K75_minus_GBLUP',])
 preds[[trait]]<-d[,c('pop_label','trait','fold','observed','predicted_GBLUP','predicted_top75','source_population')]
}
res<-do.call(rbind,results)
write.csv(res,file.path(CONFIG$out,'60_Paired_RIL_accuracy_CIs.csv'),row.names=FALSE)
write.csv(do.call(rbind,draws_out),file.path(CONFIG$out,'60_Paired_bootstrap_draws.csv'),row.names=FALSE)
write.csv(do.call(rbind,preds),file.path(CONFIG$out,'60_Matched_prediction_QC.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(CONFIG$out,'60_sessionInfo.txt'))
pdf(file.path(CONFIG$out,'60_Paired_heldout_R2.pdf'),width=9,height=5)
sub<-res[res$metric%in%c('R2_GBLUP','R2_K75'),,drop=FALSE]
ylim<-range(c(sub$lo95,sub$hi95,0),finite=TRUE)
plot(NA,xlim=c(.5,length(traits)+.5),ylim=ylim,xaxt='n',xlab='',
     ylab='Out-of-fold predictive R-squared (95% paired bootstrap CI)',
     main='RIL prediction: genome-wide GBLUP versus 75-marker GWAS')
axis(1,at=seq_along(traits),labels=traits);abline(h=0,lty=2,col='gray65')
for(k in 1:2){
 mod<-c('R2_GBLUP','R2_K75')[k]
 p<-sub[sub$metric==mod,,drop=FALSE];p<-p[match(traits,p$trait),]
 x<-seq_along(traits)+(k-1.5)*.18
 cc<-c('steelblue4','darkorange3')[k]
 arrows(x,p$lo95,x,p$hi95,code=3,angle=90,length=.025,col=cc)
 points(x,p$estimate,pch=19,col=cc)
}
legend('topright',legend=c('GBLUP','75-marker GWAS'),
  col=c('steelblue4','darkorange3'),pch=19,bty='n')
dev.off()
cat('60 COMPLETE: ',CONFIG$out,'\n')
if(any(res$resampling_unit=='RIL_only'))
 warning('At least one trait uses RIL-level resampling; related RILs can reduce effective sample size.')
