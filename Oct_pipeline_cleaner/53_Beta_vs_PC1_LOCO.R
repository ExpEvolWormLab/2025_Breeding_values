## 53_Beta_vs_PC1_LOCO.R
## Compare the selection-gradient beta predictor with a fitness-blind PC1
## of the seven trait-specific GWAS-effect vectors. All calibration and
## PCA estimation is confined to training chromosomes.
##
## Run from the project root. Uses script 45's primary GA1/GA2/GA4 response.
## Prefer full R5 GWAS table (802 testable groups); if absent, falls back to
## script 13's complete seven-trait effect table (~674 matched testable groups).
## No packages are required.

set.seed(20261102L)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
cols <- paste0('beta_', traits, '_NaCl')

first_present <- function(paths) {
  good <- paths[file.exists(paths)]
  if (length(good)) good[[1L]] else NA_character_
}

evo_path <- first_present(c(
  'selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
  '45_All_R5_clusters_GWAS_evolution_predictions.csv',
  '/mnt/data/pleiotropy_new_results/45_All_R5_clusters_GWAS_evolution_predictions.csv'))
weights_path <- first_present(c(
  'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv',
  '13_Gradient_weight_QC.csv', '/mnt/data/13_Gradient_weight_QC.csv'))
full_path <- first_present(c(
  'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv'))
older_path <- first_present(c(
  'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Cluster_predictors_and_response.csv',
  '13_Cluster_predictors_and_response.csv', '/mnt/data/13_Cluster_predictors_and_response.csv'))
if (is.na(evo_path) || is.na(weights_path))
  stop('Missing evolutionary predictions from script 45 or script-13 weight QC table')
if (is.na(full_path) && is.na(older_path))
  stop('Missing both the R5 GWAS cluster table and the older script-13 effect table')

read_input <- function(f) read.csv(f, stringsAsFactors=FALSE,check.names=FALSE)
E <- read_input(evo_path)
E <- E[E$scenario == 'primary_survivors' & E$reason == 'eligible' &
         is.finite(E$slope),,drop=FALSE]
if (anyDuplicated(as.character(E$group))) stop('Duplicate evolutionary group IDs')
W <- read_input(weights_path)
W <- W[match(traits,W$trait),,drop=FALSE]
if (anyNA(W$trait)) stop('Missing a trait in weight table')
w <- as.numeric(W$weight_full_GWAS_scale)
if (any(!is.finite(w))) stop('Non-finite beta weights')

if (!is.na(full_path)) {
  Q <- read_input(full_path)
  input_kind <- 'full_R5'
  chr_col <- 'chr'
} else {
  Q <- read_input(older_path)
  input_kind <- 'older_13'
  chr_col <- 'chromosome'
}
if (!all(c('group', chr_col, cols) %in% names(Q)))
  stop('QTL table does not contain group, chromosome and all seven GWAS effect columns')
if (anyDuplicated(as.character(Q$group))) stop('Duplicate GWAS group IDs')
j <- match(as.character(E$group), as.character(Q$group))
E <- E[!is.na(j),,drop=FALSE]
j <- j[!is.na(j)]
Q <- Q[j,,drop=FALSE]
if (any(as.character(Q[[chr_col]]) != as.character(E$chr)))
  stop('Mismatch in GWAS/evolution chromosome labels')
X <- as.matrix(Q[,cols,drop=FALSE]); storage.mode(X) <- 'double'
keep <- rowSums(is.finite(X)) == length(traits) & is.finite(E$slope) &
  !is.na(E$chr) & !is.na(E$block)
E <- E[keep,,drop=FALSE]; X <- X[keep,,drop=FALSE]
colnames(X) <- traits
if (length(unique(E$chr)) != 6L) stop('Expected six chromosomes')
if (nrow(E) < 30L) stop('Too few common groups')
y <- as.numeric(E$slope)
chrom <- as.character(E$chr)

## Verify exact consistency with selection-weighted 45 predictor.
beta_effect <- drop(X %*% w)
if (any(abs(beta_effect-as.numeric(E$pred_full)) > 1e-8))
  stop('Selection-gradient projection is inconsistent with script 45')

## Out-of-chromosome models: calibration regression y ~ effect.
## Exactly the same form as script 45. Note that training uses evolutionary
## slopes from OTHER chromosomes; this is not prediction without calibration.
fit_calibration <- function(train_score, test_score, ytrain) {
  fit <- lm.fit(cbind(1,train_score), ytrain)
  if (any(!is.finite(fit$coefficients))) stop('Nonestimable calibration')
  drop(cbind(1,test_score) %*% fit$coefficients)
}

pred <- data.frame(group=E$group, chromosome=chrom, block=E$block,
  slope=y, raw_beta=beta_effect, raw_SB=X[,'SB'],
  beta=NA_real_, SB=NA_real_, PC1_cov=NA_real_, PC1_corr=NA_real_)
loadings <- list()
for (cc in unique(chrom)) {
  tr <- chrom != cc; te <- !tr
  pred$beta[te] <- fit_calibration(beta_effect[tr],beta_effect[te],y[tr])
  pred$SB[te] <- fit_calibration(X[tr,'SB'],X[te,'SB'],y[tr])
  for (kind in c('cov','corr')) {
    use_corr <- kind == 'corr'
    pca <- prcomp(X[tr,,drop=FALSE],center=TRUE,scale.=use_corr)
    v <- pca$rotation[,1L]
    mu <- pca$center
    scl <- if (use_corr) pca$scale else rep(1,length(mu))
    scores_tr <- drop(sweep(sweep(X[tr,,drop=FALSE],2,mu,'-'),2,scl,'/') %*% v)
    scores_te <- drop(sweep(sweep(X[te,,drop=FALSE],2,mu,'-'),2,scl,'/') %*% v)
    pred[te,paste0('PC1_',kind)] <- fit_calibration(scores_tr,scores_te,y[tr])
    frac <- pca$sdev[1L]^2/sum(pca$sdev^2)
    loadings[[length(loadings)+1L]] <- data.frame(
      holdout_chr=cc, kind=kind, trait=traits,
      loading=as.numeric(v),explained_fraction=frac)
  }
}

r2_cv <- function(z) 1 - sum((y-z)^2)/sum((y-mean(y))^2)
metrics <- do.call(rbind,lapply(c('beta','SB','PC1_cov','PC1_corr'),function(nm) {
  z <- pred[[nm]]
  data.frame(input_kind=input_kind,n_clusters=length(y),n_blocks=length(unique(E$block)),
    model=nm, LOCO_R2=r2_cv(z), LOCO_MSE=mean((y-z)^2),
    prediction_vs_response_r=cor(z,y))
}))

## Paired block-bootstrap uncertainty conditional on the fitted LOCO models.
## This is an exploratory uncertainty interval, not a formal permutation p.
set.seed(20261102L)
block_idx <- split(seq_along(y),as.character(E$block))
block_keys <- names(block_idx)
B <- 4000L
boot <- do.call(rbind,lapply(c('SB','PC1_cov','PC1_corr'), function(other) {
  delta <- (y-pred[[other]])^2-(y-pred$beta)^2
  est <- mean(delta)
  draws <- replicate(B, {
    idx <- unlist(block_idx[sample(block_keys,length(block_keys),replace=TRUE)],
                  use.names=FALSE)
    mean(delta[idx])
  })
  data.frame(input_kind=input_kind,comparator=other,n_clusters=length(y),
    n_blocks=length(block_keys),delta_MSE_beta_advantage=est,
    CI_low=unname(quantile(draws,.025)),CI_high=unname(quantile(draws,.975)),
    bootstrap_sign_p_exploratory=min(1,2*min(mean(draws<=0),mean(draws>=0))))
}))
chr_rows <- do.call(rbind,lapply(unique(chrom),function(cc) {
  z <- chrom==cc
  data.frame(holdout_chr=cc,n_clusters=sum(z),
    R2_beta_on_chr=1-sum((y[z]-pred$beta[z])^2)/sum((y[z]-mean(y[z]))^2),
    R2_PC1_cov_on_chr=1-sum((y[z]-pred$PC1_cov[z])^2)/sum((y[z]-mean(y[z]))^2),
    mean_delta_MSE_beta_vs_PC1_cov=mean((y[z]-pred$PC1_cov[z])^2-(y[z]-pred$beta[z])^2),
    mean_delta_MSE_beta_vs_PC1_corr=mean((y[z]-pred$PC1_corr[z])^2-(y[z]-pred$beta[z])^2))
}))

out <- 'selection_results/beta_vs_PC1'
dir.create(out,showWarnings=FALSE,recursive=TRUE)
write.csv(metrics,file.path(out,'53_LOCO_metrics.csv'),row.names=FALSE)
write.csv(boot,file.path(out,'53_Paired_block_bootstrap.csv'),row.names=FALSE)
write.csv(do.call(rbind,loadings),file.path(out,'53_PC1_loadings.csv'),row.names=FALSE)
write.csv(chr_rows,file.path(out,'53_Chromosome_sensitivity.csv'),row.names=FALSE)
write.csv(pred,file.path(out,'53_Matched_predictions.csv'),row.names=FALSE)
cat('\nInput source:',input_kind,' | matched clusters:',nrow(E),'\n')
print(metrics,row.names=FALSE)
print(boot,row.names=FALSE)
if (input_kind == 'older_13') message('The 802-cluster analysis needs the full R5 GWAS table at: ',
  'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv')
