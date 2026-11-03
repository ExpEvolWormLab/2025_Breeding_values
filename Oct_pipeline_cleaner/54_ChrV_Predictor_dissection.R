## 54_ChrV_Predictor_dissection.R
## Purpose: audit the fitness-beta versus PC1 difference in chr V, 6--8 Mb.
## No genome-wide parameter is refitted using a selected positive result.
## Requires summary results of 11, 13, 53; optionally 17 and 45.
## Run from the original project ROOT, not from the Oct_pipeline_cleaner subfolder.
## Base R only. No RIL genotype or PoolSeq input required.

options(stringsAsFactors=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
effect_cols <- paste0('beta_', traits, '_NaCl')
out_dir <- 'selection_results/chrV_deep_dive'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
first_existing <- function(x) { x <- x[file.exists(x)]; if(length(x)) x[1] else NA_character_ }
read_required <- function(paths,description) {
  x <- first_existing(paths)
  if(is.na(x))stop('Cannot find ',description,':\n',paste(paths,collapse='\n'),call.=FALSE)
  message(description,': ',x)
  read.csv(x,stringsAsFactors=FALSE,check.names=FALSE)
}
write_out <- function(d,name) write.csv(d,file.path(out_dir,name),row.names=FALSE,na='')
q <- read_required(c('selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
                     '11_Clusters_selection_and_GWAS_09f.csv'),'R5 cluster GWAS')
p <- read_required(c('selection_results/beta_vs_PC1/53_Matched_predictions.csv',
                     '53_Matched_predictions.csv',
                     'selection_results/beta_vs_PC1/53b_Full802_matched_predictions.csv',
                     '53b_Full802_matched_predictions.csv',
                     '53_Full802_Matched_predictions.csv'), 'Matched LOCO predictions')
wdf <- read_required(c('selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv',
                       '13_Gradient_weight_QC.csv'), 'Fitness weights')
required_q <- c('group','chr','ps','n_SNPs',effect_cols)
required_p <- c('group','block','slope','raw_beta','raw_SB','beta','SB','PC1_cov')
if(!all(required_q %in% names(q))||!all(required_p %in% names(p)))stop('Unrecognized 11 or 53 column names')
if(!'chromosome'%in%names(p)) {
  if('chr'%in%names(p))p$chromosome<-p$chr else stop('No chromosome in 53 file')
}
if(anyDuplicated(q$group)||anyDuplicated(p$group))stop('Duplicated R5 group in input')
if(!all(traits%in%wdf$trait))stop('Missing trait weights')
w <- setNames(as.numeric(wdf$weight_full_GWAS_scale[match(traits,wdf$trait)]),traits)
if(any(!is.finite(w)))stop('Missing or invalid weights')
ix <- match(as.character(p$group),as.character(q$group))
if(anyNA(ix))stop('53 contains groups missing from 11')
q <- q[ix,,drop=FALSE]
if(any(as.character(q$chr)!=as.character(p$chromosome)))stop('Chromosome mismatch between 11 and 53')
X <- as.matrix(q[,effect_cols,drop=FALSE]); storage.mode(X)<-'double'; colnames(X)<-traits
if(any(!is.finite(X)) || any(!is.finite(p$slope)))stop('Nonfinite GWAS effects/responses')
raw_beta <- drop(X %*% w)
if(max(abs(raw_beta-p$raw_beta))>1e-8 || max(abs(X[,'SB']-p$raw_SB))>1e-8)
  stop('Fitness-beta or SB orientation/scale mismatch: do NOT interpret decomposition')
if(nrow(p)!=802L)warning('Expected 802 matched clusters, found ',nrow(p))

## Full V chromosome is held out, never train PCA or calibration on V data.
tr <- p$chromosome!='V'; te <- !tr
if(sum(te)==0 || sum(tr)==0)stop('No chromosome V or no training chromosomes')
fit_cal <- function(a,y) {
  cc <- coef(lm(y ~ a)); if(any(!is.finite(cc)))stop('Nonestimable calibration');cc
}
bfit <- fit_cal(raw_beta[tr],p$slope[tr])
pcfit <- prcomp(X[tr,,drop=FALSE],center=TRUE,scale.=FALSE)
loadings <- pcfit$rotation[,1]
center <- pcfit$center
Z <- sweep(X,2,center,'-')
train_pc <- drop(Z[tr,,drop=FALSE]%*%loadings)
if(cor(train_pc,raw_beta[tr])<0) { loadings <- -loadings;train_pc <- -train_pc }
pcc <- fit_cal(train_pc,p$slope[tr])
pred_beta <- unname(bfit[1]+bfit[2]*raw_beta[te])
pred_pc <- unname(pcc[1]+pcc[2]*drop(Z[te,,drop=FALSE]%*%loadings))
qc <- data.frame(check=c('n_matched','n_chrV','n_V_3','beta_raw_max_abs_difference',
  'beta_LOCO_max_abs_difference_chrV','PC1_cov_LOCO_max_abs_difference_chrV',
  'PCA_fraction_of_GWAS_variance','train_beta_calibration_slope','train_PC1_calibration_slope'),
  value=c(nrow(p),sum(te),sum(te & p$block=='V_3'), max(abs(raw_beta-p$raw_beta)),
          max(abs(pred_beta-p$beta[te])),max(abs(pred_pc-p$PC1_cov[te])),
          pcfit$sdev[1]^2/sum(pcfit$sdev^2),unname(bfit[2]),unname(pcc[2])))
write_out(qc,'54_Validation_and_model_QC.csv')
if(max(abs(pred_beta-p$beta[te]))>1e-8||max(abs(pred_pc-p$PC1_cov[te]))>1e-8)
  stop('Could not reproduce chr V held-out predictions: inspect script 53 before proceeding')

## Baseline-relative, per-trait contributions in the V-held-out models.
## All effects remain on the SAME input GWAS scales; do not standardize them again.
Cbeta <- sweep(X[te,,drop=FALSE],2,center,'-')
Cbeta <- sweep(Cbeta,2,w,'*') * unname(bfit[2])
Cpc <- sweep(Z[te,,drop=FALSE],2,loadings,'*') * unname(pcc[2])
base_beta <- unname(bfit[1]+bfit[2]*sum(center*w))
base_pc <- unname(pcc[1])
if(max(abs(rowSums(Cbeta)+base_beta-p$beta[te]))>1e-8 ||
   max(abs(rowSums(Cpc)+base_pc-p$PC1_cov[te]))>1e-8)
  stop('Per-trait contribution reconstruction failed')

P <- p[te,,drop=FALSE]; Q <- q[te,,drop=FALSE]
P$ps <- Q$ps
P$n_SNPs <- Q$n_SNPs
P$n_poolseq_covered <- if('n_poolseq_covered'%in%names(Q))Q$n_poolseq_covered else NA_integer_
P$observed_sign <- sign(P$slope)
P$beta_sign_correct <- sign(P$beta)==P$observed_sign
P$PC1_sign_correct <- sign(P$PC1_cov)==P$observed_sign
P$delta_SE_beta_advantage <- (P$slope-P$PC1_cov)^2-(P$slope-P$beta)^2
P$effect_distance_beta_PC1 <- abs(P$beta-P$PC1_cov)
P$in_V3 <- P$block=='V_3'
P$in_6_to_8Mb <- P$ps>=6000000 & P$ps<8000000
P$beta_predicted_base <- base_beta
P$PC1_predicted_base <- base_pc
P$neutral_test_tag <- NA_character_
P$test_q_BY <- NA_real_
P$test_selected_BY05 <- NA
P$median_A00_MAF <- NA_real_

## Optional unified catalogue: do not use selection-tag evidence to choose a predictor.
epath <- first_existing(c('selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
                          '45_All_R5_clusters_GWAS_evolution_predictions.csv'))
if(!is.na(epath)) {
  e <- read.csv(epath,check.names=FALSE,stringsAsFactors=FALSE)
  e <- e[e$scenario=='primary_survivors' & e$reason=='eligible',,drop=FALSE]
  if(anyDuplicated(e$group))stop('Duplicated eligible records in script 45')
  m <- match(P$group,e$group)
  if(anyNA(m))stop('Missing script-45 records for chr V groups')
  if(max(abs(e$slope[m]-P$slope))>1e-8)stop('Mismatch in script-45 evolutionary response')
  P$neutral_test_tag <- e$selection_tag_SNP_ID[m]
  P$test_q_BY <- e$selection_q_BY_conservative[m]
  P$test_selected_BY05 <- e$selection_significant_BY_05[m]
  P$median_A00_MAF <- e$median_A00_MAF[m]
}
lf <- first_existing(c('pleiotropy_results/R5/linkage_diagnostics/17_Within_cluster_lead_comparisons.csv',
                       '17_Within_cluster_lead_comparisons.csv'))
if(!is.na(lf)) {
  ld <- read.csv(lf,stringsAsFactors=FALSE,check.names=FALSE)
  ld <- ld[ld$group%in%P$group,,drop=FALSE]
  write_out(ld,'54_Existing_trait_pair_linkage_diagnostics_chrV.csv')
  P$n_candidate_distinct_trait_pairs <- vapply(P$group,function(g)
    sum(ld$group==g & ld$classification=='candidate_distinct_linked_signals'),integer(1))
}
write_out(P[order(P$ps),,drop=FALSE],'54_All_tested_chrV_cluster_prediction_audit.csv')
F <- P[P$in_V3,,drop=FALSE];F <- F[order(F$ps),,drop=FALSE]
write_out(F,'54_V3_12_cluster_prediction_audit.csv')

D <- do.call(rbind,lapply(seq_len(nrow(P)),function(i) {
  data.frame(group=P$group[i],chr='V',block=P$block[i],ps=P$ps[i],
     trait=traits,GWAS_effect=as.numeric(X[which(te)[i],]),
     fitness_gradient_weight=as.numeric(w),
     beta_raw_trait_product=as.numeric(X[which(te)[i],]*w),
     beta_calibrated_centered_contribution=as.numeric(Cbeta[i,]),
     PC1_loading=as.numeric(loadings),
     PC1_calibrated_centered_contribution=as.numeric(Cpc[i,]),
     stringsAsFactors=FALSE)
}))
write_out(D,'54_chrV_Per_trait_beta_and_PC1_contributions.csv')
write_out(D[D$group%in%F$group,], '54_V3_Per_trait_beta_and_PC1_contributions.csv')
write_out(data.frame(trait=traits,fitness_weight=as.numeric(w),
                     trained_chrV_holdout_PC1_loading=as.numeric(loadings),
                     GWAS_effect_training_mean=as.numeric(center),
                     PC1_GWAS_effect_variance_fraction=pcfit$sdev[1]^2/sum(pcfit$sdev^2)),
          '54_ChrV_heldout_PCA_loadings_and_weights.csv')

## Squared-error contributions sum across CLUSTERS, but these clusters are correlated.
all_delta <- (p$slope-p$PC1_cov)^2-(p$slope-p$beta)^2
V_delta <- sum(P$delta_SE_beta_advantage)
V3_delta <- sum(F$delta_SE_beta_advantage)
summary <- data.frame(region=c('genome_all','chrV_all','V_3','chrV_excluding_V_3'),
  n_clusters=c(nrow(p),nrow(P),nrow(F),sum(!P$in_V3)),
  sum_SE_adv_beta=c(sum(all_delta),V_delta,V3_delta,V_delta-V3_delta),
  mean_SE_adv_beta=c(mean(all_delta),mean(P$delta_SE_beta_advantage),
                     mean(F$delta_SE_beta_advantage),mean(P$delta_SE_beta_advantage[!P$in_V3])),
  beta_sign_accuracy=c(mean(sign(p$beta)==sign(p$slope)),mean(P$beta_sign_correct),
                       mean(F$beta_sign_correct),mean(P$beta_sign_correct[!P$in_V3])),
  PC1_sign_accuracy=c(mean(sign(p$PC1_cov)==sign(p$slope)),mean(P$PC1_sign_correct),
                       mean(F$PC1_sign_correct),mean(P$PC1_sign_correct[!P$in_V3])))
summary$V3_share_of_genomic_advantage <- NA_real_
if(abs(sum(all_delta))>1e-20)summary$V3_share_of_genomic_advantage[summary$region=='V_3'] <- V3_delta/sum(all_delta)
write_out(summary,'54_Regional_descriptive_summary_NO_independent_cluster_inference.csv')

## FIGURES: no false visual suggestion that 12 clusters are independent loci.
pdf(file.path(out_dir,'54_V3_Observed_vs_beta_and_PC1_predictions.pdf'),width=11,height=7)
old <- par(mar=c(5,5,4,1))
plot(F$ps/1e6,F$slope, type='b',pch=16,lwd=2,col='black',
     ylim=range(c(F$slope,F$beta,F$PC1_cov),finite=TRUE),
     xlab='R5 cluster representative position on chromosome V (Mb)',
     ylab='Allele-frequency slope / generation',
     main='Chr V, 6–8 Mb (V_3): observed and chromosome-held-out predictions')
lines(F$ps/1e6,F$beta,type='b',pch=17,col='#2166AC')
lines(F$ps/1e6,F$PC1_cov,type='b',pch=15,col='#B2182B')
abline(h=0,lty=3)
text(F$ps/1e6,F$slope,labels=F$group,pos=3,cex=.72)
legend('topright',legend=c('Observed phase-aligned','Fitness beta LOCO','Genetic PC1 LOCO'),
       col=c('black','#2166AC','#B2182B'),pch=c(16,17,15),lty=1,cex=.8,bty='n')
mtext('Cluster representative positions; physical span and RIL LD require script 55',side=1,line=4,cex=.7)
par(old);dev.off()

pdf(file.path(out_dir,'54_V3_Seven_trait_contributions_and_effects.pdf'),width=12,height=9)
show_heat <- function(M,title,legend) {
  if(!is.matrix(M))M<-as.matrix(M)
  vals<-max(abs(M),na.rm=TRUE);if(!is.finite(vals)||vals==0)vals<-1
  palette <- colorRampPalette(c('#2166AC','white','#B2182B'))(101)
  par(mar=c(5,7,3,4))
  image(seq_len(nrow(M)),seq_len(ncol(M)),M,zlim=c(-vals,vals),col=palette,
        axes=FALSE,xlab='R5 cluster group',ylab='Trait',main=title)
  axis(1,at=seq_len(nrow(M)),labels=rownames(M),las=2,cex.axis=.78)
  axis(2,at=seq_len(ncol(M)),labels=colnames(M),las=2,cex.axis=.8)
  box();mtext(legend,3,line=.25,cex=.68)
}
ord <- match(F$group,P$group)
M1 <- Cbeta[ord,,drop=FALSE];rownames(M1)<-as.character(F$group);colnames(M1)<-traits
M2 <- Cpc[ord,,drop=FALSE];rownames(M2)<-as.character(F$group);colnames(M2)<-traits
show_heat(M1,'Fitness beta: per-trait contributions to held-out predictions',
          'Centered on training chromosome means; contribution sums + intercept = prediction')
show_heat(M2,'PC1: per-trait contributions to held-out predictions',
          'Centered GWAS effects multiplied by PC1 loading and held-out calibration slope')
par(mar=c(6,7,3,1))
barplot(t(M1),beside=FALSE,legend.text=traits,args.legend=list(x='topright',cex=.75),
        main='Regional beta prediction, centered trait components',
        ylab='Contribution to allele-frequency slope',xlab='R5 group',las=2)
abline(h=0,lty=3);dev.off()

cat('\nSCRIPT 54 COMPLETED. Output directory: ',out_dir,'\n',sep='')
print(summary,row.names=FALSE)
cat('\nCAUTION: post hoc region of interest; these are correlated R5 LD clusters,\n',
    'not 12 independent selected mutations. Selection p-values are conditional\n',
    'on the Ne and null-model calibration assumptions.\n',sep='')
