## 47_RIL_predicted_vs_observed_distributions.R
## Revised 2026-10-29. No model refitting is performed.
## Compares the ORIGINAL random-RIL out-of-fold predictions, not in-sample fits.
## Intentionally EXCLUDES source-population-held-out validation.
## Requires only base R and the existing script-16 / script-40 CSV files.

## Run from project root: source('47_RIL_predicted_vs_observed_distributions.R')

candidate_paths <- function(paths, description) {
  found <- paths[file.exists(paths)]
  if (!length(found)) stop(
    'Cannot find ', description, '. Checked:\n - ',
    paste(paths, collapse='\n - '),
    '\nRun from your project root or edit the paths at the top of script 47.',
    call.=FALSE)
  found[1L]
}

GBLUP_FILE <- candidate_paths(c(
  'phenotype_prediction/GBLUP_rrBLUP/40_Out_of_fold_GBLUP_predictions.csv',
  'GBLUP_rrBLUP/40_Out_of_fold_GBLUP_predictions.csv',
  '40_Out_of_fold_GBLUP_predictions.csv'), 'rrBLUP random-RIL predictions')
TOP75_FILE <- candidate_paths(c(
  'phenotype_prediction/CV_RILs/16_Heldout_predictions_K75.csv',
  'phenotype_prediction/16_Heldout_predictions_K75.csv',
  '16_Heldout_predictions_K75.csv'), 'script-16 top-75-marker predictions')
OUT_DIR <- 'phenotype_prediction/GBLUP_rrBLUP'
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
a <- read.csv(GBLUP_FILE, stringsAsFactors=FALSE, check.names=FALSE)
b <- read.csv(TOP75_FILE, stringsAsFactors=FALSE, check.names=FALSE)
required_a <- c('pop_label','scheme','fold','trait','observed','predicted_GBLUP')
required_b <- c('pop_label','fold',setdiff(traits,'beta'))
if (!all(required_a %in% names(a)))
  stop('GBLUP CSV lacks: ',paste(setdiff(required_a,names(a)),collapse=', '))
if (!all(required_b %in% names(b)))
  stop('Top-75 CSV lacks: ',paste(setdiff(required_b,names(b)),collapse=', '))

## Use only matched random-RIL folds; NEVER plot source-population validation.
a <- a[a$scheme == 'random_RIL' & a$trait %in% traits, , drop=FALSE]
if (!nrow(a)) stop('No random_RIL predictions in script-40 input.')
a$pop_label <- as.character(a$pop_label)
b$pop_label <- as.character(b$pop_label)
if (anyDuplicated(paste(a$pop_label,a$trait,sep='|||')))
  stop('Duplicate RIL/trait rows in random_RIL GBLUP file')
if (anyDuplicated(b$pop_label)) stop('Duplicate RILs in top-75 predictions')
if (!all(traits %in% a$trait))
  stop('GBLUP predictions missing traits: ',paste(setdiff(traits,a$trait),collapse=', '))
j <- match(a$pop_label,b$pop_label)
if (anyNA(j))
  stop(sum(is.na(j)), ' GBLUP RIL IDs absent from script-16 prediction CSV')
if (any(as.character(a$fold)!=as.character(b$fold[j])))
  stop('Cross-validation fold assignments disagree between scripts 16 and 40')

## Values in both input CSVs are out-of-fold individual predictions.
a$predicted_top75 <- NA_real_
for (tr in setdiff(traits,'beta')) {
  ix <- which(a$trait==tr)
  a$predicted_top75[ix] <- as.numeric(b[[tr]][j[ix]])
}
a$observed <- as.numeric(a$observed)
a$predicted_GBLUP <- as.numeric(a$predicted_GBLUP)
write.csv(a,file.path(OUT_DIR,'47_Matched_random_RIL_out_of_fold_predictions.csv'),row.names=FALSE)

metrics <- function(y,p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]
  if (length(y)<10L || !is.finite(sd(y)) || sd(y)==0)
    return(data.frame(n=length(y),R2=NA_real_,correlation=NA_real_,
      RMSE=NA_real_,MAE=NA_real_,observed_SD=NA_real_,predicted_SD=NA_real_,
      predicted_over_observed_SD=NA_real_,calibration_intercept=NA_real_,
      calibration_slope=NA_real_))
  cal <- if(sd(p)>0)coef(lm(y~p)) else c(NA_real_,NA_real_)
  data.frame(n=length(y),R2=1-sum((y-p)^2)/sum((y-mean(y))^2),
    correlation=if(sd(p)>0)cor(y,p) else NA_real_,
    RMSE=sqrt(mean((y-p)^2)),MAE=mean(abs(y-p)),
    observed_SD=sd(y),predicted_SD=sd(p),
    predicted_over_observed_SD=sd(p)/sd(y),
    calibration_intercept=unname(cal[1]),calibration_slope=unname(cal[2]))
}
summary <- list()
for (tr in traits) {
  z <- a[a$trait==tr,,drop=FALSE]
  for (model in c('GBLUP',if(tr!='beta') 'Top75_GWAS_markers')) {
    p <- if(model=='GBLUP') z$predicted_GBLUP else z$predicted_top75
    summary[[length(summary)+1L]] <- cbind(
      data.frame(trait=tr,model=model,stringsAsFactors=FALSE),
      metrics(z$observed,p))
  }
}
summary <- do.call(rbind,summary)
rownames(summary) <- NULL
write.csv(summary,file.path(OUT_DIR,'47_Random_RIL_distribution_prediction_accuracy.csv'),row.names=FALSE)

## The observed distribution includes all phenotyped RILs. Prediction
## distributions use the identical RIL subset for each corresponding model.
## The blue/orange curves are predictions of genetic values, so a smaller
## SD than the observed phenotype is expected: this is NOT a formal
## goodness-of-fit test of full phenotypic distribution equality.
col_obs <- '#222222'; col_g <- '#2778B5'; col_t <- '#D57739'
fig <- file.path(OUT_DIR,'47_Observed_vs_predicted_RIL_distributions.pdf')
pdf(fig,width=12,height=7.6,onefile=TRUE)
for (tr in traits) {
  z <- a[a$trait==tr & is.finite(a$observed) &
           is.finite(a$predicted_GBLUP),,drop=FALSE]
  if(nrow(z)<10)next
  y <- z$observed; g <- z$predicted_GBLUP
  has_t <- tr!='beta' && sum(is.finite(z$predicted_top75))>=10L
  t <- if(has_t)z$predicted_top75 else numeric(0)
  range_all <- range(c(y,g,t),finite=TRUE)
  expand <- diff(range_all)*.035
  if(!is.finite(expand)||expand==0) expand <- 0.01
  rg <- range_all + c(-expand,expand)
  s_g <- summary[summary$trait==tr & summary$model=='GBLUP',,drop=FALSE]
  s_t <- summary[summary$trait==tr & summary$model=='Top75_GWAS_markers',,drop=FALSE]

  par(mfrow=c(1,2),mar=c(5,4.8,4.6,1.2),oma=c(1,0,0,0))
  plot(y,g,xlim=rg,ylim=rg,pch=16,cex=.55,
       col=adjustcolor(col_g,alpha.f=.35),
       xlab='Observed adjusted RIL phenotype',
       ylab='Out-of-fold predicted value',
       main=sprintf('%s: held-out RIL predictions (n=%d)',tr,nrow(z)))
  abline(0,1,lty=2,col='grey40',lwd=1.4)
  if(has_t) {
    good <- is.finite(t)
    points(y[good],t[good],pch=1,cex=.56,
           col=adjustcolor(col_t,alpha.f=.47))
  }
  leg <- c(sprintf('GBLUP (R2 = %.3f)',s_g$R2[1]))
  if(has_t) leg <- c(leg,sprintf('Top 75 GWAS (R2 = %.3f)',s_t$R2[1]))
  leg <- c(leg,'1:1 line')
  legend('topleft',legend=leg,bty='n',cex=.73,
    col=c(col_g,if(has_t)col_t,'grey40'),
    pch=c(16,if(has_t)1,NA_integer_),
    lty=c(NA_integer_,if(has_t)NA_integer_,2))

  if(length(unique(y))>2L && sd(y)>0) {
    d_y <- density(y,na.rm=TRUE)
    d_g <- if(length(unique(g))>2L && sd(g)>0)density(g,na.rm=TRUE) else NULL
    d_t <- if(has_t && length(unique(t[is.finite(t)]))>2L &&
              sd(t[is.finite(t)])>0)density(t,na.rm=TRUE) else NULL
    max_d <- max(c(d_y$y,if(!is.null(d_g))d_g$y,
                 if(!is.null(d_t))d_t$y),finite=TRUE)
    plot(d_y,col=col_obs,lwd=2.3,xlim=rg,ylim=c(0,1.06*max_d),
         xlab='Adjusted phenotype',ylab='Kernel density',
         main=sprintf('%s: observed and predicted distributions',tr))
    if(!is.null(d_g))lines(d_g,col=col_g,lwd=2.1)
    if(!is.null(d_t))lines(d_t,col=col_t,lwd=2.1)
    legend('topright',legend=c('Observed','GBLUP',if(has_t)'Top 75 GWAS'),
           col=c(col_obs,col_g,if(has_t)col_t),lwd=2.1,bty='n',cex=.75)
  }else plot.new()
  mtext('Same random-RIL cross-validation folds; no source-population-held-out results',
        outer=TRUE,side=1,line=-.5,cex=.8)
}
dev.off()
cat('Script 47 finished. Figure: ',fig,'\n',
    'Prediction/phenotype table: ',file.path(OUT_DIR,'47_Random_RIL_distribution_prediction_accuracy.csv'),
    '\n',sep='')
