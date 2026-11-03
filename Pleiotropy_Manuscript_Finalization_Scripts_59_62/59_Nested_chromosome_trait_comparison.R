# 59_Nested_chromosome_trait_comparison.R
# Aim: a fair, chromosome-held-out comparison of the fixed fitness-weighted
# seven-trait score and the best single-trait score SELECTED ONLY USING
# evolution on the OTHER chromosomes (nested chromosome CV).
# The selected single trait is NOT a prospective, evolution-free predictor.
# Uses the previously fixed 802-region catalogue; never selects by q-value.
# R >= 4.1; base R only. Run from the project root.

CONFIG <- list(
  qtl = 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
  evolution = 'selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
  weights = 'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv',
  out = 'selection_results/manuscript_finalization/59_nested_comparison',
  n_boot = 4000L, seed = 20261103L
)
resolve <- function(x) {
  if (file.exists(x)) return(x)
  y <- basename(x)
  if (file.exists(y)) return(y)
  stop('Missing required input: ', x, '\nProvide it at the original path or in the working directory.')
}
dir.create(CONFIG$out, recursive=TRUE, showWarnings=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
cols <- paste0('beta_',traits,'_NaCl')
Q <- read.csv(resolve(CONFIG$qtl), check.names=FALSE, stringsAsFactors=FALSE)
E <- read.csv(resolve(CONFIG$evolution), check.names=FALSE, stringsAsFactors=FALSE)
W <- read.csv(resolve(CONFIG$weights), check.names=FALSE, stringsAsFactors=FALSE)
check <- function(x,need,what) {
  z <- setdiff(need,names(x))
  if(length(z)) stop(what,' lacks: ', paste(z,collapse=', '))
}
check(Q,c('group','chr',cols),'GWAS table')
check(E,c('group','chr','ps','block','scenario','reason','slope','pred_full','pred_SB','LOCO_full','LOCO_SB'),'Evolution table')
check(W,c('trait','weight_full_GWAS_scale'),'Weight table')
if(anyDuplicated(Q$group)||anyDuplicated(W$trait)) stop('Duplicated GWAS group or trait weight')
E <- E[E$scenario=='primary_survivors' & E$reason=='eligible' & is.finite(E$slope),,drop=FALSE]
if(anyDuplicated(E$group)) stop('Duplicated primary group')
if(nrow(E)!=802L) warning('Expected 802 clusters; matched ',nrow(E))
j <- match(as.character(E$group),as.character(Q$group))
if(anyNA(j)) stop('No GWAS row for some R5 groups')
Q <- Q[j,,drop=FALSE]
if(any(as.character(E$chr)!=as.character(Q$chr))) stop('Chromosome mismatch')
X <- as.matrix(Q[,cols,drop=FALSE]);storage.mode(X)<-'double';colnames(X)<-traits
if(any(!is.finite(X)))stop('Nonfinite seven-trait effects')
W <- W[match(traits,W$trait),,drop=FALSE]
if(anyNA(W$trait)||any(!is.finite(W$weight_full_GWAS_scale)))stop('Missing/invalid gradient weights')
w <- W$weight_full_GWAS_scale
if(max(abs(drop(X%*%w)-E$pred_full))>1e-9)stop('Fitness score mismatch versus script 45')
if(max(abs(X[,'SB']-E$pred_SB))>1e-9)stop('SB score mismatch versus script 45')
if(any(!is.finite(E$slope)))stop('Nonfinite response')

pred_lm <- function(xtr,ytr,xte){
  if(length(ytr)<20L||sd(xtr)<=1e-12)return(rep(mean(ytr),length(xte)))
  m <- lm.fit(cbind(1,xtr),ytr)
  if(any(!is.finite(m$coefficients)))stop('Singular calibration')
  drop(cbind(1,xte)%*%m$coefficients)
}
# Inner LOCO within the five outer-training chromosomes: choose the trait
# minimizing their pooled squared prediction errors (never outer-test values).
run_nested <- function(d,X,w,scenario){
  ch <- as.character(d$chr); y <- as.numeric(d$slope)
  if(length(unique(ch))<5L)stop('Need >=5 chromosomes for nested LOCO')
  beta <- drop(X%*%w)
  P <- data.frame(group=d$group,chr=ch,ps=d$ps,block=d$block,slope=y,
                  raw_beta=beta,raw_SB=X[,'SB'],
                  beta=NA_real_,SB=NA_real_,nested_single=NA_real_,
                  chosen_trait=NA_character_,stringsAsFactors=FALSE)
  inner_rows <- list()
  for(cc in unique(ch)){
    train <- which(ch!=cc); test <- which(ch==cc)
    if(length(test)<1 || length(train)<40)stop('Empty chromosome fold')
    P$beta[test] <- pred_lm(beta[train],y[train],beta[test])
    P$SB[test] <- pred_lm(X[train,'SB'],y[train],X[test,'SB'])
    # For every candidate trait, build out-of-chromosome predictions inside
    # the outer training set. No test-chromosome outcomes influence selection.
    cv_mse <- vapply(traits,function(tr){
      pi <- rep(NA_real_,length(train))
      for(c_inner in unique(ch[train])){
        a <- which(ch[train]!=c_inner)
        b <- which(ch[train]==c_inner)
        pi[b] <- pred_lm(X[train[a],tr],y[train[a]],X[train[b],tr])
      }
      if(any(!is.finite(pi)))stop('Nonfinite inner CV scores')
      mean((y[train]-pi)^2)
    },numeric(1))
    winner <- traits[which.min(cv_mse)]
    P$nested_single[test] <- pred_lm(X[train,winner],y[train],X[test,winner])
    P$chosen_trait[test] <- winner
    inner_rows[[length(inner_rows)+1L]] <- data.frame(
      scenario=scenario,heldout_chr=cc,trait=traits,inner_MSE=as.numeric(cv_mse),
      chosen=traits==winner,n_inner_chromosomes=length(unique(ch[train])))
  }
  # Foldwise null prediction uses the training-response mean, not the test mean.
  P$training_null <- NA_real_
  for(cc in unique(ch))P$training_null[ch==cc] <- mean(y[ch!=cc])
  models <- c('beta','SB','nested_single','training_null')
  metric <- do.call(rbind,lapply(models,function(model){
    pred <- P[[model]]
    data.frame(scenario=scenario,model=model,n_clusters=length(y),
      n_chromosomes=length(unique(ch)),n_blocks=length(unique(d$block)),
      MSE=mean((y-pred)^2),
      predictive_R2=1-sum((y-pred)^2)/sum((y-mean(y))^2),
      calibrated_direction_accuracy=mean(sign(y)==sign(pred)),
      stringsAsFactors=FALSE)
  }))
  list(pred=P,metric=metric,inner=do.call(rbind,inner_rows))
}

base <- run_nested(E,X,w,'all_eligible')
# Pre-specified primary evidence is all 802; these EXPLORATORY omissions were
# chosen after seeing regional outcomes and cannot serve as independent tests.
sensitivity <- list(base)
if(any(E$block=='V_3')){
  k <- E$block!='V_3'
  sensitivity[[length(sensitivity)+1L]] <- run_nested(E[k,,drop=FALSE],X[k,,drop=FALSE],w,'excluding_V_3_refit')
}
if(any(as.character(E$chr)=='V')){
  k <- as.character(E$chr)!='V'
  sensitivity[[length(sensitivity)+1L]] <- run_nested(E[k,,drop=FALSE],X[k,,drop=FALSE],w,'excluding_V_refit')
}
metric <- do.call(rbind,lapply(sensitivity,`[[`,'metric'))
inner <- do.call(rbind,lapply(sensitivity,`[[`,'inner'))
P <- base$pred
# Basic reproduction check of beta/SB coefficients from script 45.
repro <- data.frame(n_groups=nrow(E),n_blocks=length(unique(E$block)),
 max_diff_LOCO_beta=max(abs(P$beta-E$LOCO_full)),
 max_diff_LOCO_SB=max(abs(P$SB-E$LOCO_SB)))
if(repro$max_diff_LOCO_beta>1e-9||repro$max_diff_LOCO_SB>1e-9)
 stop('LOCO scores disagree with original script 45')

set.seed(CONFIG$seed)
blocks <- split(seq_len(nrow(P)),as.character(P$block))
if(length(blocks)<10L)stop('Too few genomic blocks to bootstrap')
comp <- c('SB','nested_single','training_null')
boot <- do.call(rbind,lapply(comp,function(m){
  yy<-P$slope
  delta <- (yy-P[[m]])^2-(yy-P$beta)^2
  dd <- replicate(CONFIG$n_boot,{
    ix <- unlist(blocks[sample.int(length(blocks),length(blocks),replace=TRUE)],use.names=FALSE)
    mean(delta[ix])
  })
  data.frame(comparator=m,positive_means_beta_better=mean(delta),
    CI_low=unname(quantile(dd,.025)),CI_high=unname(quantile(dd,.975)),
    fraction_bootstrap_positive=mean(dd>0),n_boot=length(dd),
    interpretation='conditional block bootstrap; nested selection not refitted in bootstrap')
}))
# Calibration-free directional accuracy: uses phenotypic selection weights,
# and NOT the evolution-fitted calibration coefficient of LOCO.
P$raw_beta_correct <- sign(P$raw_beta)==sign(P$slope)
P$raw_SB_correct <- sign(P$raw_SB)==sign(P$slope)
raw_dir <- data.frame(n=nrow(P),accuracy_raw_beta=mean(P$raw_beta_correct),
 accuracy_raw_SB=mean(P$raw_SB_correct),n_discordant=sum(sign(P$raw_beta)!=sign(P$raw_SB)))
write.csv(P,file.path(CONFIG$out,'59_Matched_predictions.csv'),row.names=FALSE)
write.csv(metric,file.path(CONFIG$out,'59_Accuracy_and_refitted_sensitivity.csv'),row.names=FALSE)
write.csv(inner,file.path(CONFIG$out,'59_Inner_trait_selection.csv'),row.names=FALSE)
write.csv(boot,file.path(CONFIG$out,'59_Paired_block_bootstrap.csv'),row.names=FALSE)
write.csv(raw_dir,file.path(CONFIG$out,'59_Calibration_free_direction_QC.csv'),row.names=FALSE)
write.csv(repro,file.path(CONFIG$out,'59_Reproduction_QC.csv'),row.names=FALSE)
pdf(file.path(CONFIG$out,'59_Nested_trait_selection_figures.pdf'),width=10,height=4.7)
par(mfrow=c(1,2),mar=c(5,4.3,3,1))
prim<-metric[metric$scenario=='all_eligible',,drop=FALSE]
barplot(prim$predictive_R2,names.arg=prim$model,las=2,ylab='LOCO predictive R-squared',
 main='Nested trait choice versus fixed beta',col='gray75')
abline(h=0,lty=2,col='gray50')
plot(P$beta,P$slope,pch=16,cex=.6,col=rgb(.2,.2,.2,.35),
 xlab='Full beta: held-out predictions',ylab='Observed cluster slope',
 main='Fitness-informed evolutionary prediction')
abline(0,1,lty=2,col='gray40')
dev.off()
writeLines(capture.output(sessionInfo()),file.path(CONFIG$out,'59_sessionInfo.txt'))
print(metric)
cat('59 COMPLETE. Nested single-trait selection used only training chromosomes.\n')
