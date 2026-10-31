# 51_Phenotype_three_panel_manuscript_figure.R
# ------------------------------------------------------------
# Three columns x multiple traits for phenotype prediction ONLY.
#   Column 1 = historical in-sample significant-QTL score (legacy R2=0.2);
#   Column 2 = training-fold selected top-75 GWAS markers (5-fold CV);
#   Column 3 = genome-wide GBLUP predictions (same 5 folds).
#
# Beta's missing 75-marker predictions are supplied by script 50.
# Main text: SF, SB, area.F, beta. Supplement: all 8 traits.
# Horizontal and vertical values are standardized relative to the
# observed phenotype *within each dataset*, with in-sample scores
# converted to phenotype units by a FULL-DATA regression for
# DESCRIPTION ONLY. This is not held-out calibration.
# No source-population-held-out results are used.
# ------------------------------------------------------------

CONFIG <- list(
 old_49 = 'phenotype_prediction/three_panel_comparison/49_Individual_RIL_scores_and_predictions.csv',
 beta_50 = 'phenotype_prediction/three_panel_comparison/50_Beta_top75_out_of_fold.csv',
 out = 'phenotype_prediction/three_panel_comparison',
 main_traits = c('SF','SB','area.F','beta'),
 all_traits = c('SF','SB','FS','FB','BS','BF','area.F','beta'),
 n_boot = 1000L,
 seed = 20261030L
)
dir.create(CONFIG$out,recursive=TRUE,showWarnings=FALSE)
for(f in c(CONFIG$old_49,CONFIG$beta_50)) if(!file.exists(f))
 stop('Required input missing: ',f,'\nRun script 49 (for its individual-level CSV), then script 50.')

z <- read.csv(CONFIG$old_49,check.names=FALSE,stringsAsFactors=FALSE)
need <- c('pop_label','trait','observed','original_observed',
          'in_sample_GWAS_score','predicted_top75','predicted_GBLUP','fold')
if(!all(need%in%names(z)))
 stop('Script 49 CSV is missing: ',paste(setdiff(need,names(z)),collapse=', '))
if('scheme'%in%names(z)) z<-z[z$scheme=='random_RIL',,drop=FALSE]
z <- z[z$trait%in%CONFIG$all_traits,,drop=FALSE]
if(anyDuplicated(paste(z$pop_label,z$trait,sep='::')))
 stop('Duplicate RIL/trait rows in script 49 predictions')
if(!setequal(unique(z$trait),CONFIG$all_traits))
 stop('Not all eight phenotypes present in script 49 predictions')
for(nm in c('observed','original_observed','in_sample_GWAS_score',
            'predicted_top75','predicted_GBLUP'))z[[nm]]<-as.numeric(z[[nm]])

be <- read.csv(CONFIG$beta_50,stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c('pop_label','fold','observed','predicted_top75')%in%names(be)))
 stop('Script 50 beta CSV has unexpected columns')
if(anyDuplicated(be$pop_label)) stop('Duplicate beta RIL IDs')
ix <- which(z$trait=='beta')
j <- match(z$pop_label[ix],be$pop_label)
if(anyNA(j)) stop('Some beta RILs not covered by script 50')
if(!identical(as.integer(z$fold[ix]),as.integer(be$fold[j])))
 stop('Script 49 and script 50 use different folds')
good <- is.finite(z$observed[ix]) & is.finite(be$observed[j])
if(any(good) && max(abs(z$observed[ix][good]-be$observed[j][good]))>1e-6)
 stop('Beta targets differ between script 49 and script 50')
if(any(xor(is.finite(z$observed[ix]),is.finite(be$observed[j]))))
 stop('Beta missingness differs between script 49 and script 50')
if(any(is.finite(z$predicted_top75[ix])))
 warning('Overwriting existing beta top-75 values with independently fitted script 50 results')
z$predicted_top75[ix] <- as.numeric(be$predicted_top75[j])
write.csv(z,file.path(CONFIG$out,'51_Complete_three_method_RIL_predictions.csv'),row.names=FALSE)

# Compare legacy and current observed phenotype definitions on shared RILs.
qc <- do.call(rbind,lapply(CONFIG$all_traits,function(tr){
 d <- z[z$trait==tr,,drop=FALSE]
 ok <- is.finite(d$original_observed)&is.finite(d$observed)
 yy <- d$original_observed[ok]; vv <- d$observed[ok]
 data.frame(trait=tr,n_paired=sum(ok),
            correlation_original_vs_current=if(length(yy)>3&&sd(yy)>0&&sd(vv)>0)
              cor(yy,vv)else NA_real_,
            max_absolute_original_vs_current=if(length(yy))max(abs(yy-vv))else NA_real_,
            mean_legacy=if(length(yy))mean(yy)else NA_real_,
            sd_legacy=if(length(yy)>1)sd(yy)else NA_real_,
            mean_current=if(length(vv))mean(vv)else NA_real_,
            sd_current=if(length(vv)>1)sd(vv)else NA_real_)
}))
write.csv(qc,file.path(CONFIG$out,'51_Legacy_vs_current_phenotype_QC.csv'),row.names=FALSE)
if(any(is.finite(qc$correlation_original_vs_current) &
       qc$correlation_original_vs_current < 0.95))
 warning('Some original and current phenotype BLUEs correlate below 0.95. ',
         'See 51_Legacy_vs_current_phenotype_QC.csv; plots standardize separately.')

# In-sample score calibrated on the SAME complete-case phenotype data.
# This makes all scatterplots comparable visually (observed-phenotype SD units),
# but must never be interpreted as cross-validated accuracy.
z$legacy_calibrated_fit <- NA_real_
legacy_fit <- list()
for(tr in CONFIG$all_traits){
 ii <- which(z$trait==tr)
 good <- is.finite(z$original_observed[ii])&is.finite(z$in_sample_GWAS_score[ii])
 if(sum(good)<10 || sd(z$in_sample_GWAS_score[ii][good])<1e-12)
  stop('Insufficient legacy score variation for ',tr)
 fit <- lm(z$original_observed[ii][good] ~ z$in_sample_GWAS_score[ii][good])
 cf <- coef(fit)
 z$legacy_calibrated_fit[ii] <- unname(cf[1])+unname(cf[2])*z$in_sample_GWAS_score[ii]
 legacy_fit[[length(legacy_fit)+1L]] <- data.frame(trait=tr,
   n_legacy=sum(good),full_data_intercept=unname(cf[1]),
   full_data_slope=unname(cf[2]),in_sample_r2=summary(fit)$r.squared)
}
write.csv(do.call(rbind,legacy_fit),
          file.path(CONFIG$out,'51_Legacy_full_data_calibration_descriptive.csv'),row.names=FALSE)

model_names <- c('Legacy significant-QTL score','75 GWAS markers','GBLUP genome-wide')
model_keys <- c('legacy_calibrated_fit','predicted_top75','predicted_GBLUP')
model_types <- c('in_sample_full_data_fit','random_RIL_out_of_fold','random_RIL_out_of_fold')

score <- function(y,p,insample){
 ok <- is.finite(y)&is.finite(p)
 y<-y[ok];p<-p[ok]
 if(length(y)<10||sd(y)<=0 || sd(p)<=0)
  return(c(n=length(y),Pearson_r=NA_real_,in_sample_r2=NA_real_,
           CV_R2=NA_real_,RMSE=NA_real_,predicted_SD_over_observed_SD=NA_real_))
 c(n=length(y),Pearson_r=cor(y,p),
   in_sample_r2=if(insample)cor(y,p)^2 else NA_real_,
   CV_R2=if(!insample)1-sum((y-p)^2)/sum((y-mean(y))^2)else NA_real_,
   RMSE=sqrt(mean((y-p)^2)),predicted_SD_over_observed_SD=sd(p)/sd(y))
}

summary_rows <- list()
for(tr in CONFIG$all_traits){
 d <- z[z$trait==tr,,drop=FALSE]
 for(m in seq_along(model_keys)){
  y <- if(m==1L)d$original_observed else d$observed
  ss <- score(y,d[[model_keys[m]]],m==1L)
  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    trait=tr,model=model_names[m],validation=model_types[m],
    t(as.matrix(ss)),row.names=NULL,check.names=FALSE)
 }
}
metrics <- do.call(rbind,summary_rows)
write.csv(metrics,file.path(CONFIG$out,'51_Three_panel_manuscript_metrics.csv'),row.names=FALSE)

# Paired difference in held-out predictive R2 (GWAS - GBLUP), resampled
# by random CV fold. These CIs condition on the fitted models and do not
# capture uncertainty in phenotype construction or RIL family dependencies.
set.seed(CONFIG$seed)
paired <- list()
for(tr in CONFIG$all_traits){
 d <- z[z$trait==tr & is.finite(z$observed) &
         is.finite(z$predicted_top75)&is.finite(z$predicted_GBLUP),,drop=FALSE]
 if(nrow(d)<30) next
 calc_R2 <- function(o,p){
  den <- sum((o-mean(o))^2)
  if(den<1e-10) return(NA_real_)
  1-sum((o-p)^2)/den
 }
 B <- rep(NA_real_,CONFIG$n_boot)
 folds <- split(seq_len(nrow(d)),as.character(d$fold))
 for(b in seq_len(CONFIG$n_boot)){
  indices <- unlist(lapply(folds,function(ii)sample(ii,length(ii),replace=TRUE)),use.names=FALSE)
  o <- d$observed[indices]
  B[b] <- calc_R2(o,d$predicted_top75[indices])-calc_R2(o,d$predicted_GBLUP[indices])
 }
 B <- B[is.finite(B)]
 paired[[length(paired)+1L]] <- data.frame(trait=tr,n_RIL=nrow(d),
   R2_GWAS=calc_R2(d$observed,d$predicted_top75),
   R2_GBLUP=calc_R2(d$observed,d$predicted_GBLUP),
   difference_GWAS_minus_GBLUP=calc_R2(d$observed,d$predicted_top75)-
     calc_R2(d$observed,d$predicted_GBLUP),
   CI_low=if(length(B)>20)unname(quantile(B,.025))else NA_real_,
   CI_high=if(length(B)>20)unname(quantile(B,.975))else NA_real_)
}
write.csv(do.call(rbind,paired),
          file.path(CONFIG$out,'51_Paired_GWAS_vs_GBLUP_fold_stratified_bootstrap.csv'),row.names=FALSE)

# Return a standardized scatter panel: all x and y values are in units of
# the observed standard deviation, not standardized separately to unit SD.
# The legacy score is first CALIBRATED IN-SAMPLE; held-out predictions are
# not recalibrated, so regression-to-mean remains visible.
cols <- c('#344054','#CB7830','#1C7D8D')
plot_panel <- function(tr,m,limits){
 d <- z[z$trait==tr,,drop=FALSE]
 observed <- if(m==1L)d$original_observed else d$observed
 pred <- d[[model_keys[m]]]
 valid <- is.finite(observed)&is.finite(pred)
 if(sum(valid)<10){plot.new();title(main=paste(tr,'— insufficient data'));return(invisible(NULL))}
 mu <- mean(observed[valid]); sig <- sd(observed[valid])
 if(!is.finite(sig)||sig<=0)stop('Zero phenotype SD for ',tr)
 obs <- (observed[valid]-mu)/sig
 pp <- (pred[valid]-mu)/sig
 lab <- score(observed,pred,m==1L)
 label <- if(m==1L)sprintf('In-sample r² = %.3f',lab[['in_sample_r2']])
          else sprintf('Held-out R² = %.3f',lab[['CV_R2']])
 plot(obs,pp,xlim=limits,ylim=limits,pch=16,cex=.45,
      col=grDevices::adjustcolor(cols[m],alpha.f=.28),
      xlab='',ylab='',main=paste(tr,'  |  ',model_names[m]),
      cex.main=.88,font.main=2,axes=TRUE,bty='l')
 abline(0,1,col='grey62',lty=2,lwd=1.2)
 abline(h=0,v=0,col='grey90',lwd=.6)
 legend('topleft',legend=label,bty='n',cex=.81,text.col=cols[m])
}

get_limits <- function(tr){
 d <- z[z$trait==tr,,drop=FALSE]
 zz <- list()
 for(m in seq_along(model_keys)){
  y <- if(m==1L)d$original_observed else d$observed
  p <- d[[model_keys[m]]]
  ok <- is.finite(y)&is.finite(p)
  if(sum(ok)<10)next
  mu <- mean(y[ok]);sig <- sd(y[ok])
  zz[[length(zz)+1L]] <- c((y[ok]-mu)/sig,(p[ok]-mu)/sig)
 }
 values <- unlist(zz,use.names=FALSE)
 if(!length(values)) return(c(-3,3))
 lo <- min(-2.5, floor(min(values)*2)/2)
 hi <- max(2.5, ceiling(max(values)*2)/2)
 c(lo,hi)
}

draw_page <- function(traits_page,full_title){
 nr <- length(traits_page)
 op <- par(mfrow=c(nr,3),
           mar=c(3.1,3.55,2.25,0.65),
           oma=c(2.0,2.15,2.8,.2),mgp=c(1.9,.65,0),tcl=-.3)
 on.exit(par(op),add=TRUE)
 for(tr in traits_page){
  lim <- get_limits(tr)
  for(m in 1:3) plot_panel(tr,m,lim)
 }
 mtext('Observed adjusted RIL phenotype (within-trait SD units)',
       side=1,outer=TRUE,line=.55,cex=.95)
 mtext('Genetic score / prediction (within-trait SD units)',
       side=2,outer=TRUE,line=.55,cex=.95)
 mtext(full_title,side=3,outer=TRUE,line=.8,cex=1.1,font=2)
}

main_pdf <- file.path(CONFIG$out,'51_MAIN_FIGURE_four_traits_three_models.pdf')
pdf(main_pdf,width=15.5,height=14.3,onefile=TRUE,useDingbats=FALSE)
draw_page(CONFIG$main_traits,'RIL phenotype reconstruction versus held-out genomic prediction')
dev.off()

# Render an optional publication-sized raster from the same plotting function.
# Failure to open a PNG device must not prevent the vector PDF from being made.
tryCatch({
  png(file.path(CONFIG$out,'51_MAIN_FIGURE_four_traits_three_models.png'),
      width=4650,height=4290,res=300)
  draw_page(CONFIG$main_traits,'RIL phenotype reconstruction versus held-out genomic prediction')
  dev.off()
},error=function(e){
  warning('PNG rendering failed; use the vector PDF: ',conditionMessage(e))
  if(dev.cur()>1L)try(dev.off(),silent=TRUE)
})

supp_pdf <- file.path(CONFIG$out,'51_SUPPLEMENT_all_eight_traits.pdf')
pdf(supp_pdf,width=15.5,height=8.3,onefile=TRUE,useDingbats=FALSE)
for(ind in split(seq_along(CONFIG$all_traits),ceiling(seq_along(CONFIG$all_traits)/2))){
 draw_page(CONFIG$all_traits[ind],
           'RIL phenotypes: discovery-panel QTL score versus held-out predictors')
}
dev.off()

cat('\nSaved publication figure: ',main_pdf,
    '\nSupplement: ',supp_pdf,
    '\nMetrics: ',file.path(CONFIG$out,'51_Three_panel_manuscript_metrics.csv'),'\n',sep='')
cat('IMPORTANT: column 1 is full-data, legacy R2=0.2, in-sample.\n',
    'Columns 2/3 are random-RIL held-out.\n',
    'The standardized visual scale does not change the underlying CV R2.\n')
