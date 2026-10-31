## 40_rrBLUP_GBLUP_CV.R
## GBLUP of seven adjusted RIL BLUEs + fixed composite beta index.
## rrBLUP::mixed.solve fits variance components ONLY on the training RILs.
## Held-out genomic predictions use K_test,train (NOT rrBLUP's in-sample u).
## K built with training-fold allele frequencies; no outcome leakage.
## Reuses script-16 folds for a paired top-QTL comparison.
## Optional leave-entire-source-population-out sensitivity.

suppressPackageStartupMessages(library(rrBLUP))
preparation_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
fold_file <- 'phenotype_prediction/CV_RILs/16_RIL_folds.csv'
origin_file <- '19_RIL_source_populations.csv'  # optional
out_dir <- 'phenotype_prediction/GBLUP_rrBLUP'
marker_set <- 'pruned'   # 'pruned' recommended; 'full' sensitivity (slower)
min_train_MAF <- 0.01
max_train_missing <- 0.20
block_size <- 1500L
n_group_folds <- 5L
run_source_population_CV <- TRUE
seed <- 20261029L
set.seed(seed)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

if(!file.exists(preparation_file))stop('Missing: ',preparation_file)
e <- new.env(parent=baseenv());load(preparation_file,envir=e)
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
if(!all(c('genotype_matrix','genotype_pruned_matrix','phen_NaCl')%in%ls(e)))
  stop('RData lacks phenotype or SNP genotypes')
ph <- e$phen_NaCl
if(!all(c('pop_label',traits[1:7])%in%names(ph)))stop('Phenotype BLUE table malformed')
if(anyDuplicated(ph$pop_label))stop('Expected one BLUE per RIL')
if(!'beta'%in%names(ph)){
  need <- c('beta_NaCl','BLUE_mean','BLUE_sd')
  if(!all(need%in%ls(e)))stop('Cannot reconstruct beta from saved selection gradients')
  ## Index formed from the ORIGINAL BLUE units and pre-estimated fitness weights.
  ## This index is not a newly cross-fitted selection-gradient estimate.
  mu <- e$BLUE_mean[traits[1:7]]
  ss <- e$BLUE_sd[traits[1:7]]
  ww <- e$beta_NaCl[traits[1:7]]
  if(any(!is.finite(mu))||any(!is.finite(ss))||any(ss<=0)||any(!is.finite(ww)))
    stop('Invalid fixed beta index coefficients or BLUE scaling')
  xx <- sweep(as.matrix(ph[,traits[1:7],drop=FALSE]),2,mu,'-')
  xx <- sweep(xx,2,ss,'/')
  ph$beta <- drop(xx %*% ww)
}
G <- if(marker_set=='pruned')e$genotype_pruned_matrix else if(marker_set=='full')e$genotype_matrix else stop('marker_set must be pruned or full')
if(!is.matrix(G))G<-as.matrix(G)
storage.mode(G)<-'double'
if(anyDuplicated(rownames(G))||anyDuplicated(colnames(G)))stop('Duplicate genotype IDs')
ids <- intersect(rownames(G),as.character(ph$pop_label))
if(length(ids)<100L)stop('Too few matching genotyped RILs')
G <- G[ids,,drop=FALSE]
Y <- as.matrix(ph[match(ids,ph$pop_label),traits,drop=FALSE]);storage.mode(Y)<-'double'
if(any(G[is.finite(G)]<0|G[is.finite(G)]>2))stop('Expected 0/1/2 ALT dosages')
if(!length(colnames(G)))stop('Genotype matrix needs SNP IDs')
rm(e,ph);invisible(gc())

origin <- rep(NA_character_,length(ids));names(origin)<-ids
if(file.exists(origin_file)){
  oo <- read.csv(origin_file,stringsAsFactors=FALSE)
  if(!all(c('pop_label','source_population')%in%names(oo))||anyDuplicated(oo$pop_label))
    stop('Origin file needs unique pop_label and source_population columns')
  origin <- as.character(oo$source_population[match(ids,oo$pop_label)])
  origin[is.na(origin)|!nzchar(origin)|origin=='unknown']<-NA_character_
}

if(file.exists(fold_file)){
  fo <- read.csv(fold_file,stringsAsFactors=FALSE)
  if(!all(c('pop_label','fold')%in%names(fo))||anyDuplicated(fo$pop_label))stop('Malformed script-16 folds')
  random_folds <- as.integer(fo$fold[match(ids,fo$pop_label)])
  if(anyNA(random_folds))stop('Script 16 folds do not cover all RILs; rerun 16 or supply the exact folds')
} else {
  warning('No script-16 fold file: generating random folds, so top-QTL paired comparison is unavailable')
  random_folds <- sample(rep(seq_len(5L),length.out=length(ids)))
}

## Group-held-out CV is possible only for documented populations. Do not put
## all unknown-origin RILs into an artificial huge source-population group.
group_folds <- rep(NA_integer_,length(ids))
if(run_source_population_CV){
  if(anyNA(origin)){
    warning('RIL source-population CV skipped: unannotated origins (file: ',origin_file,')')
  } else {
    tab <- sort(table(origin),decreasing=TRUE)
    if(length(tab)<n_group_folds){
      warning('RIL source-population CV skipped: too few distinct origin groups')
    } else {
      sizes <- integer(n_group_folds)
      assign <- setNames(integer(length(tab)),names(tab))
      for(g in names(tab)){
        f <- which.min(sizes)
        assign[g] <- f
        sizes[f] <- sizes[f]+as.integer(tab[g])
      }
      group_folds <- unname(assign[origin])
    }
  }
}
fold_table <- data.frame(pop_label=ids,random_fold=random_folds,
                         source_pop_fold=group_folds,source_population=origin)
write.csv(fold_table,file.path(out_dir,'40_RIL_folds_and_origins.csv'),row.names=FALSE)

## This is a cross-kinship implementation of VanRaden's K=WW'/sum[2p(1-p)].
## Every p and missing-data imputation mean comes exclusively from training RILs.
## Use all marker genotypes (not phenotypes) to compute K_test,train.
make_fold_kinship <- function(G,tr,te,min_MAF,max_miss,block_size){
  Ktr <- matrix(0,length(tr),length(tr));Kte <- matrix(0,length(te),length(tr))
  used <- 0L;denom <- 0;total <- ncol(G)
  for(first in seq.int(1L,total,by=block_size)){
    jj <- first:min(first+block_size-1L,total)
    Xtr <- G[tr,jj,drop=FALSE]
    Xte <- G[te,jj,drop=FALSE]
    missing <- colMeans(!is.finite(Xtr))
    means <- colMeans(Xtr,na.rm=TRUE)
    p <- means/2
    good <- is.finite(p)&is.finite(means)&missing<=max_miss&
           pmin(p,1-p)>=min_MAF & p>0 & p<1
    if(!any(good))next
    Xtr <- Xtr[,good,drop=FALSE]
    Xte <- Xte[,good,drop=FALSE]
    mu <- means[good]; pp <- p[good]
    bad <- which(!is.finite(Xtr),arr.ind=TRUE)
    if(nrow(bad))Xtr[bad] <- mu[bad[,2]]
    bad <- which(!is.finite(Xte),arr.ind=TRUE)
    if(nrow(bad))Xte[bad] <- mu[bad[,2]]
    Xtr <- sweep(Xtr,2,mu,'-');Xte <- sweep(Xte,2,mu,'-')
    Ktr <- Ktr + tcrossprod(Xtr)
    Kte <- Kte + Xte%*%t(Xtr)
    denom <- denom + sum(2*pp*(1-pp))
    used <- used+ncol(Xtr)
    if(first%%(block_size*30L)==1L)cat('  GRM markers processed:',first,'/',total,'\n')
  }
  if(used<100L||denom<=0)stop('Too few valid markers for genomic kinship')
  list(train=Ktr/denom,test_train=Kte/denom,n_markers=used,denom=denom)
}

pred_rows <- list();fit_rows <- list();qc_rows <- list()
for(scheme in c('random_RIL','source_population')){
  ff <- if(scheme=='random_RIL')random_folds else group_folds
  if(all(is.na(ff)))next
  if(anyNA(ff))stop('Partial fold assignment in ',scheme)
  for(f in sort(unique(ff))){
    tr <- which(ff!=f);te <- which(ff==f)
    if(length(tr)<100L||length(te)<5L)stop('Training/testing fold size too small')
    cat('\n',scheme,' fold ',f,': ',length(tr),' training / ',length(te),' test lines\n',sep='')
    KK <- make_fold_kinship(G,tr,te,min_train_MAF,max_train_missing,block_size)
    qc_rows[[length(qc_rows)+1L]] <- data.frame(scheme=scheme,fold=f,
      n_train=length(tr),n_test=length(te),n_genome_markers=KK$n_markers,
      kinship_normalization=KK$denom,marker_set=marker_set)
    for(trait in traits){
      ytr <- Y[tr,trait];tte <- Y[te,trait]
      o <- which(is.finite(ytr))
      if(length(o)<60L){warning('Insufficient training phenotype: ',scheme,' ',trait);next}
      Ktr <- KK$train[o,o,drop=FALSE]
      Kcross <- KK$test_train[,o,drop=FALSE]
      fit <- tryCatch(rrBLUP::mixed.solve(y=ytr[o],K=Ktr,method='REML'),error=function(err){
        warning('GBLUP failed ',scheme,' fold ',f,' trait ',trait,': ',conditionMessage(err));NULL
      })
      if(is.null(fit))next
      mu <- as.numeric(fit$beta[1])
      lambda <- as.numeric(fit$Ve/fit$Vu)
      if(!is.finite(lambda)||lambda<=0)stop('Invalid rrBLUP Ve/Vu')
      A <- Ktr + diag(lambda+1e-9,nrow(Ktr))
      alpha <- tryCatch(solve(A,ytr[o]-mu),error=function(err){
        warning('Cross-kinship solve failed: ',conditionMessage(err));rep(NA_real_,length(o))
      })
      predicted <- as.vector(mu+Kcross%*%alpha)
      baseline <- rep(mean(ytr[o]),length(te))
      pred_rows[[length(pred_rows)+1L]] <- data.frame(
        pop_label=ids[te],scheme=scheme,fold=f,trait=trait,
        observed=tte,predicted_GBLUP=predicted,
        baseline_training_mean=baseline,
        source_population=origin[te],stringsAsFactors=FALSE)
      fit_rows[[length(fit_rows)+1L]] <- data.frame(scheme=scheme,fold=f,trait=trait,
        n_train_complete=length(o),Vu=as.numeric(fit$Vu),Ve=as.numeric(fit$Ve),
        variance_ratio=lambda,pseudo_PVE_BLUE=as.numeric(fit$Vu/(fit$Vu+fit$Ve)))
    }
    rm(KK);invisible(gc())
  }
}
if(!length(pred_rows))stop('No genomic predictions generated')
out <- do.call(rbind,pred_rows)
write.csv(out,file.path(out_dir,'40_Out_of_fold_GBLUP_predictions.csv'),row.names=FALSE)
write.csv(do.call(rbind,fit_rows),file.path(out_dir,'40_Foldwise_variance_components.csv'),row.names=FALSE)
write.csv(do.call(rbind,qc_rows),file.path(out_dir,'40_Kinship_marker_QC.csv'),row.names=FALSE)

performance <- list()
for(k in unique(paste(out$scheme,out$trait,sep='::'))){
  zz <- out[paste(out$scheme,out$trait,sep='::')==k,,drop=FALSE]
  zz <- zz[is.finite(zz$observed)&is.finite(zz$predicted_GBLUP),,drop=FALSE]
  if(nrow(zz)<20L)next
  ob <- zz$observed;pr <- zz$predicted_GBLUP;null <- zz$baseline_training_mean
  denom <- sum((ob-mean(ob))^2)
  performance[[length(performance)+1L]] <- data.frame(
    scheme=zz$scheme[1],trait=zz$trait[1],n=nrow(zz),
    R2=1-sum((ob-pr)^2)/denom,
    R2_training_mean_baseline=1-sum((ob-null)^2)/denom,
    MSE_GBLUP=mean((ob-pr)^2),MSE_training_mean=mean((ob-null)^2),
    MSE_improvement_vs_null=mean((ob-null)^2-(ob-pr)^2),
    Pearson_r=if(sd(pr)>0&&sd(ob)>0)cor(ob,pr)else NA_real_)
}
write.csv(do.call(rbind,performance),file.path(out_dir,'40_GBLUP_predictive_performance.csv'),row.names=FALSE)
cat('Script 40 complete. Out-of-fold predictions are for RIL BLUEs and fixed beta index.\n')
