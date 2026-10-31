## 41_GBLUP_comparisons_inference_figures.R
## Manuscript-oriented inference using OUT-OF-FOLD predictions from script 40.
## Paired source-population block bootstrap where source-pop IDs are complete;
## otherwise RIL resampling (less independent if relatives are clustered).
## Comparison with script-16 top-75-SNP prediction is fixed/pre-specified:
## K=75 is NOT optimized after viewing validation results.

pred_dir <- 'phenotype_prediction/GBLUP_rrBLUP'
old_dir <- 'phenotype_prediction/CV_RILs'
out_dir <- pred_dir
n_boot <- 2500L
seed <- 20261029L
set.seed(seed)

pred_file <- file.path(pred_dir,'40_Out_of_fold_GBLUP_predictions.csv')
if(!file.exists(pred_file))stop('Run script 40 first')
zz <- read.csv(pred_file,stringsAsFactors=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
need <- c('pop_label','scheme','fold','trait','observed','predicted_GBLUP',
          'baseline_training_mean','source_population')
if(!all(need%in%names(zz)))stop('Script-40 prediction file malformed')
zz$source_population <- as.character(zz$source_population)
zz$source_population[zz$source_population%in%c('','NA','unknown')]<-NA_character_

## Use the SAME script-16 RIL folds and test individuals for a paired analysis.
old_file <- file.path(old_dir,'16_Heldout_predictions_K75.csv')
if(file.exists(old_file)){
  old <- read.csv(old_file,stringsAsFactors=FALSE,check.names=FALSE)
  if(!all(c('pop_label','fold')%in%names(old))||anyDuplicated(old$pop_label))
     stop('Script-16 predictions have malformed RIL IDs')
  ix <- match(zz$pop_label,as.character(old$pop_label))
  zz$predicted_topK75 <- NA_real_
  for(t in intersect(traits,names(old))){
    ii <- which(zz$trait==t&zz$scheme=='random_RIL'&!is.na(ix))
    if(length(ii)){
      if(any(zz$fold[ii]!=old$fold[ix[ii]]))stop('Script-16 folds do not match script-40 random folds')
      zz$predicted_topK75[ii] <- as.numeric(old[[t]][ix[ii]])
    }
  }
} else {
  zz$predicted_topK75 <- NA_real_
  warning('Script-16 K75 predictions missing: top-QTL comparison skipped')
}

bootstrap_performance <- function(d,model_col){
  valid <- is.finite(d$observed)&is.finite(d[[model_col]])
  d <- d[valid,,drop=FALSE]
  if(nrow(d)<20)return(NULL)
  y <- d$observed; pr <- d[[model_col]]
  denom <- sum((y-mean(y))^2)
  R2 <- if(denom>0)1-sum((y-pr)^2)/denom else NA_real_
  MSE <- mean((y-pr)^2)
  ## Group bootstrap when all groups known and >=8 groups represented.
  group_ok <- all(!is.na(d$source_population))&&length(unique(d$source_population))>=8
  units <- if(group_ok)split(seq_len(nrow(d)),d$source_population) else
    as.list(seq_len(nrow(d)))
  B <- n_boot
  draws <- matrix(NA_real_,nrow=B,ncol=2)
  for(b in seq_len(B)){
    ii <- unlist(units[sample.int(length(units),length(units),replace=TRUE)],use.names=FALSE)
    yy <- y[ii];pp <- pr[ii];de <- sum((yy-mean(yy))^2)
    draws[b,1] <- if(de>0)1-sum((yy-pp)^2)/de else NA_real_
    draws[b,2] <- mean((yy-pp)^2)
  }
  data.frame(n=length(y),n_independent_resampling_units=length(units),
             bootstrap_unit=if(group_ok)'source_population' else 'RIL',
             R2=R2,R2_lo=unname(quantile(draws[,1],.025,na.rm=TRUE)),
             R2_hi=unname(quantile(draws[,1],.975,na.rm=TRUE)),
             MSE=MSE,MSE_lo=unname(quantile(draws[,2],.025,na.rm=TRUE)),
             MSE_hi=unname(quantile(draws[,2],.975,na.rm=TRUE)))
}

## Paired MSE improvement; positive means GBLUP is more accurate.
paired_model_comparison <- function(d,model_col){
  d <- d[is.finite(d$observed)&is.finite(d$predicted_GBLUP)&
         is.finite(d[[model_col]]),,drop=FALSE]
  if(nrow(d)<20)return(NULL)
  y <- d$observed
  diff <- (y-d[[model_col]])^2-(y-d$predicted_GBLUP)^2
  group_ok <- all(!is.na(d$source_population))&&length(unique(d$source_population))>=8
  units <- if(group_ok)split(seq_len(nrow(d)),d$source_population) else
    as.list(seq_len(nrow(d)))
  boots <- replicate(n_boot,{
    ix <- unlist(units[sample.int(length(units),length(units),replace=TRUE)],use.names=FALSE)
    mean(diff[ix])
  })
  data.frame(n=nrow(d),boot_unit=if(group_ok)'source_population'else'RIL',
    mean_MSE_improvement_GBLUP=mean(diff),
    CI_lo=unname(quantile(boots,.025)),CI_hi=unname(quantile(boots,.975)),
    two_sided_bootstrap_p=min(1,2*min(mean(boots<=0),mean(boots>=0))))
}

perf <- list();compare <- list()
for(scheme in unique(zz$scheme))for(t in traits){
  d <- zz[zz$scheme==scheme&zz$trait==t,,drop=FALSE]
  if(!nrow(d))next
  for(m in c('predicted_GBLUP','baseline_training_mean','predicted_topK75')){
    pp <- bootstrap_performance(d,m)
    if(is.null(pp))next
    pp$scheme <- scheme;pp$trait <- t;pp$model <- m
    perf[[length(perf)+1L]] <- pp
  }
  for(m in c('baseline_training_mean','predicted_topK75')){
    pp <- paired_model_comparison(d,m)
    if(is.null(pp))next
    pp$scheme <- scheme;pp$trait <- t;pp$comparator <- m
    compare[[length(compare)+1L]] <- pp
  }
}
if(!length(perf))stop('No valid prediction performance calculations')
perf <- do.call(rbind,perf)
write.csv(perf,file.path(out_dir,'41_Bootstrap_GBLUP_accuracy_CIs.csv'),row.names=FALSE)
if(length(compare)){
  compare <- do.call(rbind,compare)
  write.csv(compare,file.path(out_dir,'41_Paired_GBLUP_vs_baselines.csv'),row.names=FALSE)
}

## Predicted-value covariance is NOT genetic covariance G.
corr <- list()
for(scheme in unique(zz$scheme)){
  d <- zz[zz$scheme==scheme,,drop=FALSE]
  wide <- list();ob <- list()
  for(t in traits){
    s <- d[d$trait==t,c('pop_label','predicted_GBLUP','observed'),drop=FALSE]
    if(anyDuplicated(s$pop_label))stop('More than one prediction per RIL/trait')
    wide[[t]] <- s
  }
  id <- Reduce(intersect,lapply(wide,function(s)as.character(s$pop_label)))
  if(length(id)<20)next
  P <- O <- matrix(NA_real_,nrow=length(id),ncol=length(traits),
                  dimnames=list(id,traits))
  for(t in traits){
    j <- match(id,wide[[t]]$pop_label)
    P[,t] <- wide[[t]]$predicted_GBLUP[j]
    O[,t] <- wide[[t]]$observed[j]
  }
  for(i in seq_along(traits))for(j in seq_along(traits)){
    if(j<=i)next
    c1 <- suppressWarnings(cor(O[,i],O[,j],use='pairwise.complete.obs'))
    c2 <- suppressWarnings(cor(P[,i],P[,j],use='pairwise.complete.obs'))
    corr[[length(corr)+1L]] <- data.frame(scheme=scheme,trait1=traits[i],
      trait2=traits[j],r_observed_BLUE=c1,r_out_of_fold_GBLUP_predictions=c2,
      difference=c2-c1)
  }
}
if(length(corr)){
  corr <- do.call(rbind,corr)
  write.csv(corr,file.path(out_dir,'41_Observed_vs_predicted_trait_correlations.csv'),row.names=FALSE)
}

## Main figure: fixed K75 comparator, not optimally selected K.
pdf(file.path(out_dir,'41_GBLUP_vs_top75_marker_prediction.pdf'),width=11,height=6)
for(scheme in unique(perf$scheme)){
  z <- perf[perf$scheme==scheme,,drop=FALSE]
  lab <- traits[traits%in%z$trait]
  mm <- c('predicted_GBLUP','predicted_topK75','baseline_training_mean')
  xx <- seq_along(lab)
  plot(NA,xlim=c(.6,length(lab)+.4),ylim=range(c(z$R2_lo,z$R2_hi,0),finite=TRUE),
     xaxt='n',xlab='',ylab='Held-out R2 (bootstrap 95% CI)',
     main=paste('Genomic prediction:',scheme))
  axis(1,at=xx,labels=lab)
  abline(h=0,lty=2,col='grey70')
  col <- c('steelblue4','darkorange3','grey60')
  for(k in seq_along(mm)){
    s <- z[z$model==mm[k],,drop=FALSE]
    jj <- match(lab,s$trait)
    at <- xx+(k-2)*.18
    keep <- !is.na(jj)
    j <- jj[keep];x <- at[keep]
    arrows(x,s$R2_lo[j],x,s$R2_hi[j],code=3,angle=90,length=.025,col=col[k])
    points(x,s$R2[j],pch=19,col=col[k])
  }
  legend('topright',legend=c('GBLUP','Top 75 train-ranked SNPs','Fold training mean'),
     col=col,pch=19,bty='n',cex=.8)
}
dev.off()

if(exists('corr')&&is.data.frame(corr)&&nrow(corr)){
  pdf(file.path(out_dir,'41_Multitrait_correlation_recovery.pdf'),width=7,height=6)
  for(scheme in unique(corr$scheme)){
    z <- corr[corr$scheme==scheme,,drop=FALSE]
    plot(z$r_observed_BLUE,z$r_out_of_fold_GBLUP_predictions,
       xlim=c(-1,1),ylim=c(-1,1),pch=19,col='steelblue4',
       xlab='Observed RIL BLUE pairwise correlation',
       ylab='GBLUP predicted-value pairwise correlation',
       main=paste('Trait covariance pattern:',scheme))
    abline(0,1,lty=2,col='grey50')
  }
  dev.off()
}
cat('Script 41 complete. Bootstrap and model comparisons are exploratory, paired when possible.\n')
