## 16_Crossvalidated_phenotype_prediction.R
## Nested marker ranking within 5 held-out RIL folds, no use of full-data GWAS p-values.
## Primary markers: the full SNP matrix, not clusters discovered on all phenotypes.
## Standardized ridge effects on top K LD-pruned, TRAIN-ranked SNPs.
## Correlations among out-of-fold genetic predictions are compared with BLUE correlations.
## Run from project root after 00. Requires the preparatory RData only.

traits <- c('SF','SB','FS','FB','BS','BF','area.F')
K_values <- c(1L, 2L, 5L, 10L, 20L, 40L, 75L, 125L, 200L)
n_outer <- 5L
seed <- 20261029L
max_leads <- 4000L             # SNP candidates inspected for train-set LD clumping
ld_prune_r2 <- .80
lambda_grid <- c(1, 10, 100, 1000)  # penalty on standardized X'X scale
output_dir <- 'phenotype_prediction/CV_RILs'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

p <- new.env(parent = globalenv())
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData', envir = p)
ph <- p$phen_NaCl
ph$pop_label <- as.character(ph$pop_label)
## genotype_matrix: RILs x SNPs. Markers are selected within each outer training set.
G <- p$genotype_matrix
if (!is.matrix(G)) G <- as.matrix(G)
storage.mode(G) <- 'numeric'
if (anyDuplicated(rownames(G)) || anyDuplicated(colnames(G)) ||
    anyDuplicated(ph$pop_label)) stop('Duplicate genotype or phenotype IDs')
ids <- intersect(rownames(G), ph$pop_label)
if (length(ids) < 30L) stop('Too few RILs matching genotypes and BLUEs')
G <- G[ids, , drop = FALSE]
Y <- as.matrix(ph[match(ids, ph$pop_label), traits, drop = FALSE])
storage.mode(Y) <- 'numeric'
rm(p, ph)
n <- nrow(G)
cat('RILs:',n,' SNPs:',ncol(G),'\n')

## Same folds for all traits, so genetic predictions are paired across traits.
fold <- sample(rep(seq_len(n_outer), length.out = n))
write.csv(data.frame(pop_label=ids,fold=fold), file.path(output_dir,'16_RIL_folds.csv'),row.names=FALSE)

## Data-frame extraction for a small selected SNP set. Missing genotypes are
## imputed by the TRAIN mean and training SDs are used for the held-out set.
standardize_markers <- function(Gtr, Gte) {
  mu <- colMeans(Gtr,na.rm=TRUE)
  mu[!is.finite(mu)] <- 0
  for (j in seq_along(mu)) {
    Gtr[!is.finite(Gtr[,j]),j] <- mu[j]
    Gte[!is.finite(Gte[,j]),j] <- mu[j]
  }
  sig <- apply(Gtr,2,sd)
  sig[!is.finite(sig)|sig < 1e-10] <- 1
  list(train=sweep(sweep(Gtr,2,mu,'-'),2,sig,'/'),
       test=sweep(sweep(Gte,2,mu,'-'),2,sig,'/'))
}

## Training-only marginal correlation rankings. Chunking avoids creating an
## RIL x genome-wide-SNP standardized matrix in memory.
rank_train <- function(G, tr, y, block_size=3000L) {
  z <- rep(NA_real_,ncol(G))
  yy <- y - mean(y)
  yss <- sum(yy^2)
  if (!is.finite(yss) || yss < 1e-10) return(integer())
  for (first in seq.int(1L,ncol(G),by=block_size)) {
    ii <- first:min(first+block_size-1L,ncol(G))
    x <- G[tr,ii,drop=FALSE]
    ## matrix can contain NAs; fill by train means, no test-data leakage
    mx <- colMeans(x,na.rm=TRUE)
    mx[!is.finite(mx)] <- 0
    bad <- which(!is.finite(x),arr.ind=TRUE)
    if (nrow(bad)) x[bad] <- mx[bad[,2]]
    x <- sweep(x,2,mx,'-')
    sxx <- colSums(x^2)
    v <- as.vector(crossprod(yy,x))/sqrt(pmax(sxx*yss,1e-30))
    v[sxx < 1e-10] <- NA_real_
    z[ii] <- v
  }
  order(abs(z),decreasing=TRUE,na.last=NA)
}

## Greedy clumping from ranked candidate SNPs using training RIL genotypes.
## If too few candidates survive, inspect additional ranked SNPs.
prune_ranked <- function(G,tr,ord,n_target,threshold,cap) {
  if (!length(ord)) return(integer())
  selected <- integer(0)
  for (j in head(ord,cap)) {
    x <- G[tr,j]
    if (sum(is.finite(x)) < 10) next
    if (!length(selected)) { selected <- j; if(length(selected)>=n_target)break; next }
    mat <- G[tr,selected,drop=FALSE]
    keep <- is.finite(x) & rowSums(is.finite(mat)) == ncol(mat)
    if (sum(keep) < 10L) next
    r <- suppressWarnings(cor(x[keep],mat[keep,,drop=FALSE]))
    if (!any(is.finite(r) & r^2 >= threshold)) selected <- c(selected,j)
    if (length(selected) >= n_target) break
  }
  selected
}

## Small ridge regression; penalty λ||b||² with a free intercept.
fit_predict <- function(X,y,Xtest,lam) {
  ym <- mean(y)
  if (!ncol(X)) return(rep(ym,nrow(Xtest)))
  b <- tryCatch(solve(crossprod(X)+diag(lam,ncol(X)),crossprod(X,y-ym)),
                error=function(e) rep(0,ncol(X)))
  as.vector(ym + Xtest %*% b)
}

## Inner validation tunes shrinkage, using only outer-training data.
choose_lambda <- function(X,y,lambda_grid,seed_i) {
  if(ncol(X)==0L) return(lambda_grid[1])
  set.seed(seed_i)
  inner <- sample(rep(1:3,length.out=length(y)))
  mse <- numeric(length(lambda_grid))
  for(k in seq_along(lambda_grid)){
    pe <- rep(NA_real_,length(y))
    for(f in 1:3) {
      a <- which(inner != f); b <- which(inner==f)
      ## Training fold means/SDs are re-used within each inner fold.
      xc <- X[a,,drop=FALSE]; xv <- X[b,,drop=FALSE]
      mu <- colMeans(xc); sig <- apply(xc,2,sd)
      sig[!is.finite(sig)|sig<1e-10] <- 1
      xc <- sweep(sweep(xc,2,mu,'-'),2,sig,'/')
      xv <- sweep(sweep(xv,2,mu,'-'),2,sig,'/')
      pe[b] <- fit_predict(xc,y[a],xv,lambda_grid[k])
    }
    mse[k] <- mean((y-pe)^2)
  }
  lambda_grid[which.min(mse)]
}

pred <- array(NA_real_,dim=c(n,length(traits),length(K_values)),
              dimnames=list(ids,traits,paste0('K',K_values)))
rank_log <- NULL
for(t in seq_along(traits)) {
  trait <- traits[t]
  cat('\nTrait',trait,'\n')
  ok_y <- which(is.finite(Y[,t]))
  for(f in seq_len(n_outer)){
    tr <- intersect(which(fold!=f),ok_y)
    te <- intersect(which(fold==f),ok_y)
    if(length(tr)<30 ||length(te)<3)next
    ord <- rank_train(G,tr,Y[tr,t])
    leads <- prune_ranked(G,tr,ord,max(K_values),ld_prune_r2,max_leads)
    cat(' fold',f,':',length(leads),'LD-pruned markers\n')
    rank_log <- rbind(rank_log,data.frame(trait=trait,fold=f,
                                         rank=seq_along(leads),SNP=colnames(G)[leads]))
    for(h in seq_along(K_values)){
      k <- min(K_values[h],length(leads))
      if(k<1)next
      cols <- leads[seq_len(k)]
      ## Original train-only matrix for tuning, with train-only imputation.
      XX <- G[tr,cols,drop=FALSE]; Xt <- G[te,cols,drop=FALSE]
      st <- standardize_markers(XX,Xt)
      lambda <- choose_lambda(st$train,Y[tr,t],lambda_grid,
                              seed+10000*t+100*f+h)
      pred[te,t,h] <- fit_predict(st$train,Y[tr,t],st$test,lambda)
    }
  }
}
write.csv(rank_log,file.path(output_dir,'16_Train_fold_SNP_ranks.csv'),row.names=FALSE)

## Out-of-fold R2, not in-sample reconstruction R2. Denominator uses a
## constant (complete-case observed mean) to permit comparisons across K.
performance <- NULL
for(t in seq_along(traits))for(h in seq_along(K_values)){
  yy<-Y[,t]; pp<-pred[,t,h]; ok<-is.finite(yy)&is.finite(pp)
  yy<-yy[ok];pp<-pp[ok]
  r<-if(length(yy)>5 && sd(yy)>0 && sd(pp)>0) cor(yy,pp) else NA_real_
  performance<-rbind(performance,data.frame(trait=traits[t],K=K_values[h],
   n=length(yy),r=r,R2=if(length(yy)>5)1-sum((yy-pp)^2)/sum((yy-mean(yy))^2) else NA_real_))
}
write.csv(performance,file.path(output_dir,'16_Prediction_vs_number_of_SNPs.csv'),row.names=FALSE)

## Correlations of predicted genetic values between phenotypes, compared with
## observed adjusted RIL BLUE correlations. Estimates are NOT G correlations.
chosen_K <- 75L
h <- which.min(abs(K_values-chosen_K))
P <- pred[,,h]
phen_corr <- cor(Y,use='pairwise.complete.obs')
pred_corr <- cor(P,use='pairwise.complete.obs')
cor_table <- do.call(rbind,lapply(seq_along(traits),function(i)
 do.call(rbind,lapply(seq_along(traits),function(j){
  if(j<=i)return(NULL)
  data.frame(trait1=traits[i],trait2=traits[j],K=K_values[h],
             r_observed_BLUE=phen_corr[i,j],r_predicted_genetic=pred_corr[i,j],
             difference=pred_corr[i,j]-phen_corr[i,j])
 }))))
write.csv(cor_table,file.path(output_dir,'16_Observed_vs_predicted_trait_correlations.csv'),row.names=FALSE)
write.csv(data.frame(pop_label=ids,fold=fold,Y,check.names=FALSE),
          file.path(output_dir,'16_Observed_BLUEs.csv'),row.names=FALSE)
write.csv(data.frame(pop_label=ids,fold=fold,P,check.names=FALSE),
          file.path(output_dir,'16_Heldout_predictions_K75.csv'),row.names=FALSE)

pdf(file.path(output_dir,'16_Phenotypic_prediction_and_correlations.pdf'),
    width=11,height=5)
par(mfrow=c(1,2),mar=c(4.2,4.4,2.5,1))
cols<-grDevices::hcl.colors(length(traits),'Dark 3')
plot(NA,xlim=range(K_values),ylim=range(c(-.1,performance$R2),na.rm=TRUE),
     xlab='Top LD-pruned SNPs (ranked within training RILs)',
     ylab='Held-out R2',main='Out-of-fold prediction',bty='n')
for(t in seq_along(traits)){
 z<-subset(performance,trait==traits[t]); lines(z$K,z$R2,type='b',pch=16,
 col=cols[t],lwd=1.5,cex=.6)
}
abline(h=0,lty=3)
legend('bottomright',legend=traits,col=cols,lty=1,pch=16,cex=.7,bty='n',ncol=2)
plot(cor_table$r_observed_BLUE,cor_table$r_predicted_genetic,
     pch=16,col='#506f87',xlab='Observed BLUE correlation',
     ylab='Out-of-fold predicted-value correlation',
     main=paste0('Trait correlations (K=',K_values[h],')'),bty='n')
abline(0,1,lty=2,col='gray50')
par(mfrow=c(1,1));dev.off()
cat('\nSaved cross-validated prediction results to',output_dir,'\n')
cat('Note: random RIL folds are not necessarily family-independent;\n',
    'the polygenic predictions target adjusted RIL BLUEs, not individual phenotypes.\n')
