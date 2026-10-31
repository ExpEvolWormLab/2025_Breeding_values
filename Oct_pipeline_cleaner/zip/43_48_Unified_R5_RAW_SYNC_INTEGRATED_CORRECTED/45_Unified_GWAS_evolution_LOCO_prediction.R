## 45_Unified_GWAS_evolution_LOCO_prediction.R
## Main-text genome-wide prediction, ALL eligible ancestral frequencies.
## No SNP/cluster is filtered by selection p-value or observed response.
## GA1, GA2, GA4 primary; extinct GA3 supplementary only.
source('43_Unified_config_and_helpers.R')
cfg <- uf_config
out <- cfg$prefix
need <- c(cfg$r5_gwas_file,cfg$weights_file,
          file.path(out,'43_Unified_R5_cluster_registry.csv'),
          file.path(out,'44_primary_survivors_WF_cluster_selection.csv'))
uf_check(need)
registry <- uf_read(file.path(out,'43_Unified_R5_cluster_registry.csv'))
qtl <- uf_read(cfg$r5_gwas_file)
w <- uf_read(cfg$weights_file)
sel <- uf_read(file.path(out,'44_primary_survivors_WF_cluster_selection.csv'))
if(anyDuplicated(registry$group)||anyDuplicated(qtl$group)||anyDuplicated(sel$group))
  stop('Duplicate cluster ID')
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
bcols <- paste0('beta_',traits,'_NaCl')
req <- c('group','chr','beta_beta_NaCl','beta_SB_NaCl',bcols)
if(!all(req%in%names(qtl)))stop('Script-11 R5 GWAS columns not found')
if(!all(c('trait','weight_full_GWAS_scale','weight_three_GWAS_scale')%in%names(w)))
  stop('Script-13 selection-gradient weight columns not found')
w <- w[match(traits,w$trait),,drop=FALSE]
if(anyNA(w$trait))stop('Missing one or more weights for the seven traits')
w1 <- as.numeric(w$weight_full_GWAS_scale)
w3 <- as.numeric(w$weight_three_GWAS_scale)
if(any(!is.finite(w1))||any(!is.finite(w3)))stop('Invalid selection gradient weights')
B <- as.matrix(qtl[,bcols,drop=FALSE]);storage.mode(B) <- 'double'
qtl$pred_full <- as.numeric(B%*%w1)
qtl$pred_reduced <- as.numeric(B%*%w3)
qtl$pred_SB <- as.numeric(qtl$beta_SB_NaCl)
qtl$pred_direct_beta <- as.numeric(qtl$beta_beta_NaCl)
models <- c(full='pred_full', reduced='pred_reduced', SB='pred_SB',
            direct_beta='pred_direct_beta')
x <- registry
jj <- match(as.character(x$group),as.character(qtl$group))
if(anyNA(jj))stop('Registry group missing from original R5 GWAS')
for(m in unname(models))x[[m]] <- qtl[[m]][jj]
x$chr <- as.character(x$chr)
x$block <- paste0(x$chr,'_',floor(as.numeric(x$ps)/cfg$block_bp))
k <- match(as.character(x$group),as.character(sel$group))
for(n in c('p_MC_upper_simultaneous','q_BY_conservative','q_BH_conservative',
           'significant_BY_05','tag_SNP_ID')){
 col <- paste0('selection_',n)
 x[[col]] <- sel[[n]][k]
}
if(anyNA(x$group) || anyDuplicated(x$group))stop('Bad registry membership')

safe_cor <- function(a,b){ok <- is.finite(a)&is.finite(b)
 if(sum(ok)<10L||sd(a[ok])<=0||sd(b[ok])<=0)return(NA_real_)
 unname(cor(a[ok],b[ok]))}
loco <- function(d, response, predictor){
 observed <- d[[response]];effect <- d[[predictor]];chr <- d$chr
 output <- rep(NA_real_,nrow(d))
 for(cc in unique(chr[!is.na(chr)])){
   te <- which(chr==cc & is.finite(observed)&is.finite(effect))
   tr <- which(chr!=cc & is.finite(observed)&is.finite(effect))
   if(!length(te)||length(tr)<25L||length(unique(chr[tr]))<3L||
      sd(effect[tr])<1e-12)next
   f <- lm(observed[tr]~effect[tr])
   cf <- coef(f)
   if(all(is.finite(cf)))output[te] <- unname(cf[1]+cf[2]*effect[te])
 }
 output
}
allpred <- list();metrics <- list();boots <- list()
scenarios <- c(primary_survivors='primary_slope_median_phase',
               supplementary_all_GA='supp_slope_median_phase')
for(scenario in names(scenarios)){
 scol <- scenarios[[scenario]]
 if(!scol%in%names(x))stop('Missing scenario slope ',scol)
 z <- x
 z$scenario <- scenario
 z$slope <- as.numeric(z[[scol]])
 for(model in names(models)){
   cn <- models[[model]]
   pred <- loco(z,'slope',cn)
   z[[paste0('LOCO_',model)]] <- pred
   good <- is.finite(z$slope)&is.finite(z[[cn]])&is.finite(pred)
   valid <- z[good,,drop=FALSE]
   if(nrow(valid)<30L||length(unique(valid$chr))<4L){
     metrics[[length(metrics)+1L]] <- data.frame(scenario=scenario,model=model,
       n_clusters=nrow(valid),n_chromosomes=length(unique(valid$chr)),
       n_blocks=length(unique(valid$block)),Pearson_r=NA_real_,LOCO_R2=NA_real_,
       LOCO_MSE=NA_real_,LOCO_null_MSE=NA_real_)
   } else {
     err <- valid$slope-valid[[paste0('LOCO_',model)]]
     denom <- sum((valid$slope-mean(valid$slope))^2)
     metrics[[length(metrics)+1L]] <- data.frame(scenario=scenario,model=model,
       n_clusters=nrow(valid),n_chromosomes=length(unique(valid$chr)),
       n_blocks=length(unique(valid$block)),
       Pearson_r=safe_cor(valid$slope,valid[[cn]]),
       LOCO_R2=if(denom>0)1-sum(err^2)/denom else NA_real_,
       LOCO_MSE=mean(err^2),LOCO_null_MSE=mean((valid$slope-mean(valid$slope))^2))
   }
 }
 for(m in c('reduced','SB','direct_beta')){
   z1 <- z[,c('slope','block',paste0('LOCO_',c('full',m))),drop=FALSE]
   bb <- uf_block_boot(z1,paste0('LOCO_full'),paste0('LOCO_',m),
                       B=cfg$bootstrap_B,seed=cfg$seed+match(m,names(models))+
                         if(scenario=='primary_survivors')0L else 100L)
   bb$scenario <- scenario;bb$model_A <- 'full';bb$model_B <- m
   boots[[length(boots)+1L]] <- bb
 }
 allpred[[scenario]] <- z
}
y <- do.call(rbind,allpred);rownames(y) <- NULL
if(anyDuplicated(paste(y$scenario,y$group)))stop('More than one prediction per cluster/scenario')
uf_write(y,file.path(out,'45_All_R5_clusters_GWAS_evolution_predictions.csv'))
uf_write(do.call(rbind,metrics),file.path(out,'45_LOCO_performance_main_and_supp.csv'))
uf_write(do.call(rbind,boots),file.path(out,'45_Paired_genomic_block_bootstrap.csv'))

## Genomic-block bootstrap CIs for the original GWAS-effect/evolution
## correlation, with each chromosome spatially segmented into 2Mb blocks.
set.seed(cfg$seed+45L)
ci <- list()
for(s in names(scenarios))for(m in names(models)){
 zz <- y[y$scenario==s & is.finite(y$slope) &
           is.finite(y[[models[[m]]]]) & !is.na(y$block),,drop=FALSE]
 blk <- split(seq_len(nrow(zz)),zz$block)
 if(length(blk)<10L)next
 obs <- safe_cor(zz$slope,zz[[models[[m]]]])
 boot <- replicate(cfg$bootstrap_B,{
   idx <- unlist(blk[sample.int(length(blk),length(blk),replace=TRUE)],use.names=FALSE)
   safe_cor(zz$slope[idx],zz[[models[[m]]]][idx])
 })
 boot <- boot[is.finite(boot)]
 if(length(boot)<100L)next
 ci[[length(ci)+1L]] <- data.frame(scenario=s,model=m,n_clusters=nrow(zz),
    n_blocks=length(blk),r=obs,CI_lower=unname(quantile(boot,.025)),
    CI_upper=unname(quantile(boot,.975)),n_bootstrap=length(boot))
}
if(length(ci))uf_write(do.call(rbind,ci),file.path(out,'45_Block_bootstrap_correlations.csv'))

## Explicitly prohibit selection-enriched subsets as PRIMARY prediction:
## a supplementary descriptive table only; NOT separate CV fits.
qc <- data.frame(n_catalogue=nrow(registry),
 n_primary_with_phase_slope=sum(is.finite(registry$primary_slope_median_phase)),
 n_primary_with_testable_tag=sum(!is.na(registry$tag_SNP_ID)),
 n_primary_selected_BY05=sum(sel$significant_BY_05,na.rm=TRUE),
 n_primary_with_full_GWAS_effect=sum(is.finite(x$pred_full)),
 n_all_trait_effects_finite=sum(rowSums(is.finite(B))==7L),
 primary_includes_non_significant=TRUE,
 main_population='GA1+GA2+GA4',supp_population='GA1+GA2+GA3+GA4')
uf_write(qc,file.path(out,'45_Unified_prediction_QC.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'45_sessionInfo.txt'))
cat('45 DONE: ',qc$n_primary_with_phase_slope,' eligible R5 responses (no MAF cutoff).\n',sep='')
