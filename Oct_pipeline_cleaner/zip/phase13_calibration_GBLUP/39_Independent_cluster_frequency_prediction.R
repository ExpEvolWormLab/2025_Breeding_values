## 39_Independent_cluster_frequency_prediction.R
## Assign EACH R5 cluster to ONE ancestral-MAF class using the median
## MAF of its available GWAS SNPs (no selection statistics involved).
## Re-estimate its phase-aligned slope ONCE, from all eligible SNPs, then
## compare fixed RIL predictors with held-out-chromosome evaluation.
## Main population set = GA1,GA2,GA4. GA3-containing set = supplementary.
## Paired chromosome-aware genomic block bootstrap for full-vs-reduced/SB.

phase_file <- 'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv'
cluster_file <- 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv'
weights_file <- 'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv'
base_dir <- 'selection_results/frequency_strata'
output_dir <- file.path(base_dir,'unique_cluster_classes')
min_snp_per_cluster <- 2L
min_clusters_per_model <- 15L
min_chromosomes <- 3L
block_bp <- 2000000L
n_boot <- 2000L
seed <- 20261029L
set.seed(seed)
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

required <- c(phase_file,cluster_file,weights_file,
 file.path(base_dir,'33_All_SNPs_retained_GA_Ne1000.csv'))
if(!all(file.exists(required)))stop('Missing prerequisite: ',paste(required[!file.exists(required)],collapse=', '))
phase <- read.csv(phase_file,stringsAsFactors=FALSE)
qtl <- read.csv(cluster_file,stringsAsFactors=FALSE,check.names=FALSE)
W <- read.csv(weights_file,stringsAsFactors=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
trait_effect_cols <- paste0('beta_',traits,'_NaCl')
stopifnot(all(c('rs','group','phase')%in%names(phase)))
stopifnot(all(c('group','chr','beta_beta_NaCl','beta_SB_NaCl',trait_effect_cols)%in%names(qtl)))
stopifnot(all(c('trait','weight_full_GWAS_scale','weight_three_GWAS_scale')%in%names(W)))
if(anyDuplicated(phase$rs)||anyDuplicated(qtl$group))stop('Duplicate SNP/cluster IDs')
if(any(!phase$phase%in%c(-1,1)))stop('Phase must be +/-1')
W <- W[match(traits,W$trait),,drop=FALSE]
if(anyNA(W$trait))stop('Missing index weights')
B <- as.matrix(qtl[,trait_effect_cols,drop=FALSE]);storage.mode(B)<-'double'
qtl$pred_full <- drop(B %*% W$weight_full_GWAS_scale)
qtl$pred_reduced <- drop(B %*% W$weight_three_GWAS_scale)
qtl$pred_mapped_beta <- as.numeric(qtl$beta_beta_NaCl)
qtl$pred_SB <- as.numeric(qtl$beta_SB_NaCl)
qtl$group <- as.character(qtl$group)
phase$group <- as.character(phase$group)
models <- c(full='pred_full',reduced='pred_reduced',SB='pred_SB',mapped_beta='pred_mapped_beta')
classes <- c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50')

safe_cor <- function(x,y){
 if(length(x)<5 || sd(x)<=0 || sd(y)<=0)return(NA_real_)
 cor(x,y)
}
LOCO <- function(dat,pred_col){
 pred <- rep(NA_real_,nrow(dat))
 for(chr in unique(dat$CHROM)){
   te <- which(dat$CHROM==chr);tr <- which(dat$CHROM!=chr)
   if(length(tr)<10 || length(te)<1 || sd(dat[[pred_col]][tr])<1e-10)next
   fit <- lm(y~x,data=data.frame(y=dat$slope[tr],x=dat[[pred_col]][tr]))
   pred[te] <- unname(coef(fit)[1]+coef(fit)[2]*dat[[pred_col]][te])
 }
 pred
}

bootstrap_one <- function(d, comp, B=n_boot){
 # Resample genomic 2-Mb blocks, preserving within-block LD.
 valid <- is.finite(d$slope)&is.finite(d$LOCO_full)&is.finite(d[[paste0('LOCO_',comp)]])
 x <- d[valid,,drop=FALSE]
 blocks <- split(seq_len(nrow(x)),x$block)
 if(length(blocks)<6L)return(data.frame(comparison=comp,n_blocks=length(blocks),
   delta_MSE=NA_real_,CI_lo=NA_real_,CI_hi=NA_real_,p_two_sided=NA_real_))
 mse_full <- mean((x$slope-x$LOCO_full)^2)
 mse_other <- mean((x$slope-x[[paste0('LOCO_',comp)]])^2)
 delta <- mse_other-mse_full  # positive = full better
 boots <- replicate(B,{
   take <- unlist(blocks[sample.int(length(blocks),length(blocks),replace=TRUE)],use.names=FALSE)
   mean((x$slope[take]-x[[paste0('LOCO_',comp)]][take])^2-
        (x$slope[take]-x$LOCO_full[take])^2)
 })
 data.frame(comparison=comp,n_blocks=length(blocks),delta_MSE=delta,
   CI_lo=unname(quantile(boots,.025)),CI_hi=unname(quantile(boots,.975)),
   p_two_sided=min(1,2*min(mean(boots<=0),mean(boots>=0))))
}

out_rows <- list();metrics <- list();boot_rows <- list();qc <- list()
for(scenario in c('retained_GA','all_GA')){
  fn <- file.path(base_dir,paste0('33_All_SNPs_',scenario,'_Ne1000.csv'))
  if(!file.exists(fn)){
    if(scenario=='retained_GA')stop('Missing primary survivor-only script33 output')
    next
  }
  snp <- read.csv(fn,stringsAsFactors=FALSE)
  if(anyDuplicated(snp$ID))stop('Duplicate IDs in ',fn)
  i <- match(phase$rs,as.character(snp$ID))
  good <- !is.na(i)
  x <- phase[good,c('rs','group','phase'),drop=FALSE]
  sm <- snp[i[good],,drop=FALSE]
  x$MAF <- as.numeric(sm$MAF_G0)
  x$slope_phase <- x$phase*as.numeric(sm$slope_freq_per_generation)
  x$SNP_POS <- suppressWarnings(as.numeric(sub('^.*_','',x$rs)))
  x <- x[is.finite(x$MAF)&is.finite(x$slope_phase),,drop=FALSE]
  all <- split(x,x$group)
  res <- lapply(all,function(z){
    if(nrow(z)<min_snp_per_cluster)return(NULL)
    data.frame(group=as.character(z$group[1]),
      n_snp=nrow(z),
      median_MAF=median(z$MAF),
      min_MAF=min(z$MAF),max_MAF=max(z$MAF),
      fraction_snp_MAF_below_0.10=mean(z$MAF<.10),
      slope=median(z$slope_phase),
      median_pos=median(z$SNP_POS,na.rm=TRUE),
      stringsAsFactors=FALSE)
  })
  res <- do.call(rbind,res)
  if(is.null(res)||!nrow(res))stop('No eligible R5 clusters in ',scenario)
  j <- match(res$group,qtl$group)
  res$CHROM <- as.character(qtl$chr[j])
  for(k in names(models))res[[models[[k]]]] <- as.numeric(qtl[[models[[k]]]][j])
  res <- res[!is.na(res$CHROM)&is.finite(res$median_MAF),,drop=FALSE]
  res$MAF_class <- as.character(cut(res$median_MAF,c(-Inf,.01,.05,.10,.25,Inf),
     right=FALSE,labels=classes))
  res$scenario <- scenario
  res$block <- paste0(res$CHROM,'_',floor(res$median_pos/block_bp))
  qc[[length(qc)+1L]] <- data.frame(scenario=scenario,
    n_unique_clusters=nrow(res),n_unique_group_ids=length(unique(res$group)),
    all_bins_disjoint=!anyDuplicated(res$group),n_blocks=length(unique(res$block)))
  if(anyDuplicated(res$group))stop('R5 group counted more than once')
  for(bin in c(classes,'ALL')){
    z <- if(bin=='ALL')res else res[res$MAF_class==bin,,drop=FALSE]
    if(!nrow(z))next
    nchr <- length(unique(z$CHROM))
    for(model in names(models)){
      col <- models[[model]]
      val <- is.finite(z[[col]])&is.finite(z$slope)
      zz <- z[val,,drop=FALSE]
      if(nrow(zz)<min_clusters_per_model||nchr<min_chromosomes)next
      pred <- LOCO(zz,col)
      key <- paste0('LOCO_',model)
      z[[key]] <- NA_real_
      z[[key]][val] <- pred
      valid <- is.finite(pred)
      r <- safe_cor(zz[[col]],zz$slope)
      denom <- sum((zz$slope[valid]-mean(zz$slope[valid]))^2)
      R2 <- if(sum(valid)>10&&denom>0)1-sum((zz$slope[valid]-pred[valid])^2)/denom else NA_real_
      metrics[[length(metrics)+1L]] <- data.frame(scenario=scenario,MAF_class=bin,
        model=model,n_clusters=nrow(zz),n_chromosomes=nchr,
        n_blocks=length(unique(zz$block)),Pearson_r=r,LOCO_R2=R2,
        LOCO_MSE=if(any(valid))mean((zz$slope[valid]-pred[valid])^2)else NA_real_)
    }
    z$scenario <- scenario;z$reported_class <- bin
    out_rows[[length(out_rows)+1L]] <- z
    if(all(paste0('LOCO_',c('full','reduced','SB'))%in%names(z))){
      for(comp in c('reduced','SB')){
        b <- bootstrap_one(z,comp)
        b$scenario <- scenario;b$MAF_class <- bin
        boot_rows[[length(boot_rows)+1L]] <- b
      }
    }
  }
}
write.csv(do.call(rbind,qc),file.path(output_dir,'39_Disjoint_cluster_membership_QC.csv'),row.names=FALSE)
if(!length(metrics))stop('No class had enough independently assigned clusters')
write.csv(do.call(rbind,metrics),file.path(output_dir,'39_Prediction_by_unique_MAF_class.csv'),row.names=FALSE)
write.csv(do.call(rbind,out_rows),file.path(output_dir,'39_Cluster_predictions_and_residuals.csv'),row.names=FALSE)
if(length(boot_rows))write.csv(do.call(rbind,boot_rows),
  file.path(output_dir,'39_Paired_block_bootstrap_full_vs_comparators.csv'),row.names=FALSE)

## Matched genomic-block bootstrap tests whether the FULL-index correlation
## differs across *disjoint* frequency classes. Resample the SAME genomic
## 2-Mb block identifiers across all classes to retain cross-class linkage.
heterogeneity <- list()
raw <- do.call(rbind,out_rows)
raw <- raw[raw$reported_class != 'ALL' & is.finite(raw$pred_full) &
           is.finite(raw$slope), , drop=FALSE]
reference_class <- '0.10-0.25'
for(scenario in unique(raw$scenario)){
  ds <- raw[raw$scenario==scenario,,drop=FALSE]
  bk <- split(seq_len(nrow(ds)),ds$block)
  if(length(bk)<8L)next
  for(alt in setdiff(classes,reference_class)){
    rr <- ds[ds$MAF_class==reference_class,,drop=FALSE]
    aa <- ds[ds$MAF_class==alt,,drop=FALSE]
    if(nrow(rr)<min_clusters_per_model||nrow(aa)<min_clusters_per_model)next
    obs <- safe_cor(rr$pred_full,rr$slope)-safe_cor(aa$pred_full,aa$slope)
    boots <- replicate(n_boot,{
      take <- unlist(bk[sample.int(length(bk),length(bk),replace=TRUE)],use.names=FALSE)
      sampled <- ds[take,,drop=FALSE]
      r1 <- sampled[sampled$MAF_class==reference_class,,drop=FALSE]
      r2 <- sampled[sampled$MAF_class==alt,,drop=FALSE]
      if(nrow(r1)<10L||nrow(r2)<10L)NA_real_ else
        safe_cor(r1$pred_full,r1$slope)-safe_cor(r2$pred_full,r2$slope)
    })
    ok <- boots[is.finite(boots)]
    if(length(ok)<n_boot/2)next
    heterogeneity[[length(heterogeneity)+1L]] <- data.frame(
      scenario=scenario,reference_class=reference_class,comparison_class=alt,
      n_ref=nrow(rr),n_comparison=nrow(aa),n_blocks=length(bk),
      correlation_difference=obs,CI_lo=unname(quantile(ok,.025)),
      CI_hi=unname(quantile(ok,.975)),
      two_sided_bootstrap_p=min(1,2*min(mean(ok<=0),mean(ok>=0))),
      usable_bootstraps=length(ok))
  }
}
if(length(heterogeneity))write.csv(do.call(rbind,heterogeneity),
  file.path(output_dir,'39_Frequency_class_heterogeneity_block_bootstrap.csv'),row.names=FALSE)

mm <- do.call(rbind,metrics)
pdf(file.path(output_dir,'39_Unique_frequency_class_prediction.pdf'),width=10,height=6)
for(scenario in unique(mm$scenario)){
  z <- mm[mm$scenario==scenario&mm$MAF_class!='ALL',,drop=FALSE]
  mat <- xtabs(LOCO_R2~MAF_class+model,z)
  vals <- as.matrix(mat[intersect(classes,rownames(mat)),,drop=FALSE])
  barplot(t(vals),beside=TRUE,ylab='Leave-one-chromosome-out R2',
     main=paste('Disjoint R5 frequency classes:',scenario),
     las=2,col=c('steelblue4','grey60','tan3','darkseagreen4'))
  abline(h=0,lty=2,col='grey40')
  legend('topright',legend=colnames(vals),fill=c('steelblue4','grey60','tan3','darkseagreen4'),bty='n',cex=.75)
}
dev.off()
cat('Script 39 complete: each R5 group belongs to exactly ONE MAF class.\n')
