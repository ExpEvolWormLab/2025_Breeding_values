# 50_Beta_top75_CV_and_K_sensitivity.R
# ------------------------------------------------------------
# Extends the design of script 16 to the fixed fitness index beta.
# Uses the SAME pre-existing 5 random-RIL folds as scripts 16 and 40.
# Marker rankings, LD clumping and lambda estimation are derived
# from OUTER TRAINING RILs only. Never select K based on test R2.
# Outputs a fixed K=75 cross-validated beta predictor, plus an
# exploratory K-sensitivity curve (K=1...200).
#
# Requires no packages beyond base R. Run from project root.
# ------------------------------------------------------------

CONFIG <- list(
  prep = 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData',
  folds = 'phenotype_prediction/CV_RILs/16_RIL_folds.csv',
  gblup = 'phenotype_prediction/GBLUP_rrBLUP/40_Out_of_fold_GBLUP_predictions.csv',
  out = 'phenotype_prediction/three_panel_comparison',
  K_grid = c(1L, 2L, 5L, 10L, 20L, 40L, 75L, 125L, 200L),
  main_K = 75L,
  lambda_grid = c(1,10,100,1000),
  ld_r2 = .80,
  max_ranked_to_prune = 4000L,
  rank_block = 3000L,
  seed = 20261029L,
  min_train = 40L,
  min_test = 5L
)
stopifnot(CONFIG$main_K %in% CONFIG$K_grid)
dir.create(CONFIG$out, recursive=TRUE, showWarnings=FALSE)
for (f in c(CONFIG$prep, CONFIG$folds, CONFIG$gblup))
  if (!file.exists(f)) stop('Missing required input: ',f)

e <- new.env(parent=baseenv())
load(CONFIG$prep, envir=e)
if (!all(c('phen_NaCl','genotype_matrix') %in% ls(e)))
  stop('Preparatory RData lacks phen_NaCl or genotype_matrix')
ph <- e$phen_NaCl
G <- e$genotype_matrix
if(!is.matrix(G)) G <- as.matrix(G)
storage.mode(G) <- 'double'
if (is.null(rownames(G)) || is.null(colnames(G)) ||
    anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)))
  stop('Genotype matrix must have unique RIL and SNP IDs')
if(!'pop_label' %in% names(ph) || anyDuplicated(ph$pop_label))
  stop('Phenotypes need unique pop_label IDs')
ph$pop_label <- as.character(ph$pop_label)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')

# Exactly match script 40's fixed composite index and standardization.
if (!'beta' %in% names(ph)) {
  if (!all(c('beta_NaCl','BLUE_mean','BLUE_sd') %in% ls(e)))
    stop('beta missing and beta_NaCl/BLUE_mean/BLUE_sd unavailable')
  mu <- e$BLUE_mean[traits]
  ss <- e$BLUE_sd[traits]
  ww <- e$beta_NaCl[traits]
  if(any(!is.finite(mu)) || any(!is.finite(ss)) || any(ss<=0) ||
     any(!is.finite(ww))) stop('Invalid fixed beta coefficients or standardization')
  xx <- sweep(as.matrix(ph[,traits,drop=FALSE]),2,mu,'-')
  xx <- sweep(xx,2,ss,'/')
  ph$beta <- as.vector(xx %*% ww)
}
ph$beta <- as.numeric(ph$beta)

folds <- read.csv(CONFIG$folds,stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c('pop_label','fold') %in% names(folds)) ||
   anyDuplicated(folds$pop_label)) stop('Malformed script 16 fold table')
folds$pop_label <- as.character(folds$pop_label)
ids <- intersect(rownames(G),ph$pop_label)
if(length(ids)<100) stop('Too few matching RIL genotypes/phenotypes')
if(!all(ids %in% folds$pop_label)) stop('Some RILs lack script 16 fold labels')
G <- G[ids,,drop=FALSE]
y <- ph$beta[match(ids,ph$pop_label)]
fold <- as.integer(folds$fold[match(ids,folds$pop_label)])
if(anyNA(fold) || length(unique(fold))!=5L ||
   !setequal(sort(unique(fold)),1:5)) stop('Expected five folds numbered 1...5')
if(any(G[is.finite(G)]<0 | G[is.finite(G)]>2))
  stop('Expect 0/1/2 genotype dosage as in script 16')
rm(e,ph);invisible(gc())
cat('Beta model: ',length(ids),' RILs; ',ncol(G),' genotype SNPs; ',
    sum(is.finite(y)),' measured beta values\n',sep='')

# Verify that the response and fold assignment match the completed GBLUP.
g <- read.csv(CONFIG$gblup, stringsAsFactors=FALSE)
if(!all(c('pop_label','scheme','fold','trait','observed') %in% names(g)))
  stop('Unexpected format in script 40 predictions')
g <- g[g$scheme=='random_RIL' & g$trait=='beta',,drop=FALSE]
if(anyDuplicated(g$pop_label) || !length(g$pop_label))
  stop('Missing/duplicate random-RIL beta GBLUP predictions')
ix <- match(ids,g$pop_label)
if(anyNA(ix)) stop('Script 40 GBLUP does not contain every genotype RIL')
if(!identical(as.integer(g$fold[ix]),fold))
  stop('Script 16 versus script 40 fold mismatch; cannot do paired CV')
ok <- is.finite(y) & is.finite(g$observed[ix])
if(any(ok) && max(abs(y[ok]-g$observed[ix][ok]))>1e-6)
  stop('The beta target differs from GBLUP. Check index construction.')
if(any(xor(is.finite(y),is.finite(g$observed[ix]))))
  stop('Beta missingness differs from GBLUP beta observations')

rank_train <- function(G,tr,y,block_size) {
  yy <- y-mean(y); yss <- sum(yy^2)
  if(!is.finite(yss) || yss<1e-10) return(integer())
  z <- rep(NA_real_,ncol(G))
  for(first in seq.int(1L,ncol(G),by=block_size)){
    jj <- first:min(first+block_size-1L,ncol(G))
    xx <- G[tr,jj,drop=FALSE]
    m <- colMeans(xx,na.rm=TRUE);m[!is.finite(m)]<-0
    bad <- which(!is.finite(xx),arr.ind=TRUE)
    if(nrow(bad)) xx[bad] <- m[bad[,2]]
    xx <- sweep(xx,2,m,'-')
    sxx <- colSums(xx^2)
    r <- as.vector(crossprod(yy,xx))/sqrt(pmax(sxx*yss,1e-30))
    r[sxx<1e-10] <- NA_real_
    z[jj] <- r
  }
  order(abs(z),decreasing=TRUE,na.last=NA)
}

prune_ranked <- function(G,tr,ord,top_n,r2_threshold,cap) {
  sel <- integer(0)
  for(j in head(ord,cap)){
    x <- G[tr,j]
    if(sum(is.finite(x))<10L || sd(x,na.rm=TRUE)==0) next
    if(!length(sel)){
      sel <- j
    } else {
      mat <- G[tr,sel,drop=FALSE]
      complete <- is.finite(x) & rowSums(is.finite(mat))==length(sel)
      if(sum(complete)<10L) next
      r <- suppressWarnings(cor(x[complete],mat[complete,,drop=FALSE]))
      if(!any(is.finite(r) & r^2 >= r2_threshold)) sel <- c(sel,j)
    }
    if(length(sel)>=top_n) break
  }
  sel
}

standardize <- function(A,B){
  m <- colMeans(A,na.rm=TRUE);m[!is.finite(m)]<-0
  bad <- which(!is.finite(A),arr.ind=TRUE)
  if(nrow(bad)) A[bad] <- m[bad[,2]]
  bad <- which(!is.finite(B),arr.ind=TRUE)
  if(nrow(bad)) B[bad] <- m[bad[,2]]
  s <- apply(A,2,sd);s[!is.finite(s)|s<1e-10]<-1
  list(train=sweep(sweep(A,2,m,'-'),2,s,'/'),
       test=sweep(sweep(B,2,m,'-'),2,s,'/'))
}

ridge_predict <- function(X,y,B,lambda){
  yy <- mean(y)
  if(ncol(X)<1) return(rep(yy,nrow(B)))
  b <- tryCatch(solve(crossprod(X)+diag(lambda,ncol(X)),
                      crossprod(X,y-yy)),error=function(e)rep(0,ncol(X)))
  as.vector(yy+B%*%b)
}

choose_lambda <- function(X,y,lambda_grid,seed_i){
  # Mirrors script 16's inner validation, using outer-training RILs only.
  set.seed(seed_i)
  in_fold <- sample(rep(1:3,length.out=length(y)))
  errors <- rep(NA_real_,length(lambda_grid))
  for(h in seq_along(lambda_grid)){
    pp <- rep(NA_real_,length(y))
    for(j in 1:3){
      tr <- which(in_fold!=j);te <- which(in_fold==j)
      st <- standardize(X[tr,,drop=FALSE],X[te,,drop=FALSE])
      pp[te] <- ridge_predict(st$train,y[tr],st$test,lambda_grid[h])
    }
    errors[h] <- mean((y-pp)^2)
  }
  lambda_grid[which.min(errors)]
}

Kvals <- CONFIG$K_grid
P <- matrix(NA_real_,length(ids),length(Kvals),
            dimnames=list(ids,paste0('K',Kvals)))
log_rows <- list();mark_rows <- list()
for(f in 1:5){
  tr <- which(fold!=f & is.finite(y))
  te <- which(fold==f & is.finite(y))
  if(length(tr)<CONFIG$min_train || length(te)<CONFIG$min_test)
    stop('Insufficient measured beta for fold ',f)
  cat('Fitting held-out beta: fold ',f,' / 5; ',length(tr),' training / ',
      length(te),' test RILs\n',sep='')
  ord <- rank_train(G,tr,y[tr],CONFIG$rank_block)
  leads <- prune_ranked(G,tr,ord,max(Kvals),CONFIG$ld_r2,
                        CONFIG$max_ranked_to_prune)
  if(!length(leads)) stop('No usable beta markers in fold ',f)
  mark_rows[[length(mark_rows)+1L]] <- data.frame(
    fold=f,rank=seq_along(leads),SNP=colnames(G)[leads])
  for(h in seq_along(Kvals)){
    k <- min(Kvals[h],length(leads))
    jj <- leads[seq_len(k)]
    st <- standardize(G[tr,jj,drop=FALSE],G[te,jj,drop=FALSE])
    lambda <- choose_lambda(st$train,y[tr],CONFIG$lambda_grid,
                            CONFIG$seed+10000L*8L+100L*f+h)
    P[te,h] <- ridge_predict(st$train,y[tr],st$test,lambda)
    log_rows[[length(log_rows)+1L]] <- data.frame(
      fold=f,K_requested=Kvals[h],K_used=k,lambda=lambda,
      n_train=length(tr),n_test=length(te))
  }
}

pred <- data.frame(pop_label=ids,fold=fold,trait='beta',observed=y,
                   predicted_top75=P[,which(Kvals==CONFIG$main_K)],
                   check.names=FALSE)
write.csv(pred,file.path(CONFIG$out,'50_Beta_top75_out_of_fold.csv'),row.names=FALSE)
write.csv(data.frame(pop_label=ids,fold=fold,observed=y,P,check.names=FALSE),
          file.path(CONFIG$out,'50_Beta_all_K_out_of_fold.csv'),row.names=FALSE)
write.csv(do.call(rbind,log_rows),file.path(CONFIG$out,'50_Beta_fold_tuning.csv'),row.names=FALSE)
write.csv(do.call(rbind,mark_rows),file.path(CONFIG$out,'50_Beta_ranked_markers.csv'),row.names=FALSE)

metrics <- do.call(rbind,lapply(seq_along(Kvals),function(h){
  ok <- is.finite(y)&is.finite(P[,h])
  yy <- y[ok];pp <- P[ok,h]
  data.frame(trait='beta',K=Kvals[h],n=length(yy),
             r=if(length(yy)>3&&sd(pp)>0)cor(yy,pp)else NA_real_,
             R2=if(length(yy)>3)1-sum((yy-pp)^2)/sum((yy-mean(yy))^2)else NA_real_)
}))
write.csv(metrics,file.path(CONFIG$out,'50_Beta_prediction_vs_K.csv'),row.names=FALSE)

old_k_file <- 'phenotype_prediction/CV_RILs/16_Prediction_vs_number_of_SNPs.csv'
old_k <- if(file.exists(old_k_file))read.csv(old_k_file)else NULL
all_k <- if(!is.null(old_k))rbind(old_k,metrics)else metrics
write.csv(all_k,file.path(CONFIG$out,'50_All_traits_K_sensitivity_descriptive.csv'),row.names=FALSE)

pdf(file.path(CONFIG$out,'50_Marker_count_sensitivity_supplement.pdf'),width=8.5,height=5.5)
tr_names <- unique(all_k$trait)
cols <- setNames(grDevices::hcl.colors(length(tr_names),'Dark 3'),tr_names)
plot(NA,xlim=range(Kvals),ylim=range(c(all_k$R2,0),na.rm=TRUE),
     xlab='Number of training-ranked, LD-pruned GWAS markers (K)',
     ylab='Random-RIL held-out predictive R-squared',
     main='Sensitivity to marker count (exploratory)',bty='n')
abline(h=0,col='grey60',lty=3)
abline(v=CONFIG$main_K,col='grey50',lty=2)
for(t in tr_names){
 z <- all_k[all_k$trait==t,,drop=FALSE]
 lines(z$K,z$R2,type='b',lwd=1.5,pch=16,cex=.6,col=cols[t])
}
legend('topright',legend=tr_names,col=cols[tr_names],lwd=1.5,pch=16,
       cex=.77,ncol=2,bty='n')
mtext('Do not select K by the outer-fold test R-squared shown here.',
      side=1,line=3.6,cex=.72,col='grey40')
dev.off()
cat('\nFinished beta top-75 predictions and K sensitivity: ',CONFIG$out,'\n',sep='')
cat('WARNING: K=75 is a fixed historical comparator, not an optimized K.\n')
