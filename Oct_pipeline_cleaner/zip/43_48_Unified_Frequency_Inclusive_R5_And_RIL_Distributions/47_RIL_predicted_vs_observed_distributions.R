## 47_RIL_predicted_vs_observed_distributions.R
## ANSWERS: do independently GWAS-selected QTLs reproduce RIL phenotypes
## and their distributions as reliably as GBLUP?
## Only matched RANDOM-RIL out-of-fold validation. Never use source-held-out
## predictions in this manuscript comparison. No new GWAS/GBLUP runs.

gblup_file <- 'phenotype_prediction/GBLUP_rrBLUP/40_Out_of_fold_GBLUP_predictions.csv'
top75_file <- 'phenotype_prediction/CV_RILs/16_Heldout_predictions_K75.csv'
output_dir <- 'phenotype_prediction/GBLUP_rrBLUP'
if(!file.exists(gblup_file)||!file.exists(top75_file))
 stop('Missing the original GBLUP out-of-fold predictions or script-16 K75 predictions')
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
a <- read.csv(gblup_file,stringsAsFactors=FALSE)
b <- read.csv(top75_file,stringsAsFactors=FALSE)
a <- a[a$scheme=='random_RIL',,drop=FALSE]
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
if(anyDuplicated(paste(a$pop_label,a$trait))||anyDuplicated(b$pop_label))
  stop('Non-unique individual/trait predictions')
need <- c('pop_label','trait','fold','observed','predicted_GBLUP')
if(!all(need%in%names(a))||!all(c('pop_label','fold')%in%names(b)))
  stop('Required prediction columns missing')
matchid <- match(a$pop_label,b$pop_label)
if(anyNA(matchid))stop('Some RILs missing from original top75 file')
if(any(a$fold!=b$fold[matchid]))stop('GBLUP and top75 fold assignments differ')
a$predicted_top75 <- NA_real_
for(tr in setdiff(traits,'beta')){
 ix <- which(a$trait==tr)
 if(!tr%in%names(b))stop('Top75 predictions missing trait ',tr)
 a$predicted_top75[ix] <- as.numeric(b[[tr]][matchid[ix]])
}
write.csv(a,file.path(output_dir,'47_RIL_observed_and_two_CV_predictions.csv'),row.names=FALSE)

## Summary of accuracy and variance shrinkage of OUT-OF-FOLD predictions.
performance <- function(y,p){ok <- is.finite(y)&is.finite(p)
 if(sum(ok)<20)return(c(n=sum(ok),r=NA_real_,R2=NA_real_,predicted_SD=NA_real_,
    observed_SD=NA_real_,ratio_SD=NA_real_,RMSE=NA_real_))
 y<-y[ok];p<-p[ok]
 c(n=length(y),r=cor(y,p),
   R2=if(var(y)>0)1-sum((y-p)^2)/sum((y-mean(y))^2)else NA_real_,
   predicted_SD=sd(p),observed_SD=sd(y),ratio_SD=sd(p)/sd(y),
   RMSE=sqrt(mean((y-p)^2)))}
res <- list()
for(tr in traits)for(m in c('predicted_GBLUP','predicted_top75')){
 z <- a[a$trait==tr,,drop=FALSE]
 if(tr=='beta'&&m=='predicted_top75')next
 v <- performance(z$observed,z[[m]])
 res[[length(res)+1L]] <- data.frame(trait=tr,method=m,t(v),row.names=NULL)
}
write.csv(do.call(rbind,res),file.path(output_dir,'47_Trait_distribution_and_prediction_accuracy.csv'),row.names=FALSE)

pdf(file.path(output_dir,'47_Observed_vs_predicted_RIL_distributions.pdf'),
    width=12,height=8,onefile=TRUE)
for(tr in traits){
 z <- a[a$trait==tr & is.finite(a$observed)&is.finite(a$predicted_GBLUP),,drop=FALSE]
 if(!nrow(z))next
 par(mfrow=c(1,2),mar=c(5,4.5,4,1.5))
 cols <- c(observed='black',GBLUP='#2b70ac',top75='#c16a42')
 rg <- range(c(z$observed,z$predicted_GBLUP,z$predicted_top75),na.rm=TRUE)
 plot(z$observed,z$predicted_GBLUP,pch=16,col=adjustcolor(cols['GBLUP'],alpha.f=.35),
      xlab='Observed RIL phenotype (adjusted BLUE)',
      ylab='Out-of-fold genomic prediction',
      xlim=rg,ylim=rg,main=paste(tr,'— held-out individual predictions'))
 abline(a=0,b=1,col='grey40',lty=2)
 if(tr!='beta'){
  ok <- is.finite(z$predicted_top75)
  points(z$observed[ok],z$predicted_top75[ok],pch=1,
         col=adjustcolor(cols['top75'],alpha.f=.45),cex=.65)
 }
 legend('topleft',legend=c('GBLUP',if(tr!='beta')'75 GWAS markers' else NULL,
        'Perfect prediction'),
        col=c(cols['GBLUP'],if(tr!='beta')cols['top75'] else NULL,'grey40'),
        pch=c(16,if(tr!='beta')1 else NULL,NA),
        lty=c(NA,if(tr!='beta')NA else NULL,2),cex=.75,bty='n')
 dens <- density(z$observed,na.rm=TRUE)
 plot(dens,lwd=2,col=cols['observed'],main=paste(tr,'— phenotypic distributions'),
      xlab='Adjusted phenotype',ylab='Density',xlim=rg)
 lines(density(z$predicted_GBLUP,na.rm=TRUE),lwd=2,col=cols['GBLUP'])
 if(tr!='beta'&&sum(is.finite(z$predicted_top75))>20){
    lines(density(z$predicted_top75,na.rm=TRUE),lwd=2,col=cols['top75'])
 }
 legend('topright',legend=c('Observed','GBLUP out-of-fold',
        if(tr!='beta')'75-marker out-of-fold' else NULL),
        col=c(cols['observed'],cols['GBLUP'],
              if(tr!='beta')cols['top75'] else NULL),lwd=2,bty='n',cex=.76)
 mtext('Random-RIL CV only; beta 75-marker prediction unavailable',side=1,outer=FALSE,
       line=4.1,cex=.65)
}
dev.off()
cat('47 DONE: paired observed/predicted phenotypic distributions (random folds only).\n')
