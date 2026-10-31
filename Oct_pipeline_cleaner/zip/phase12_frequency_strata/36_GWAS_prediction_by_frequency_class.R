## 36_GWAS_prediction_by_frequency_class.R
## Does a pre-existing, independently built RIL GWAS predictor describe
## phase-aligned SNP trajectories differently across G0 frequency strata?
## The same R5 SNP->group->phase assignments are used for all comparisons.
## Exploratory: bins of SNPs are disjoint; bins of CLUSTERS need not be.

phase_file <- 'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv'
cluster_file <- 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv'
weight_file <- 'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv'
predictor_QC_file <- 'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Cluster_predictors_and_response.csv'
base_dir <- 'selection_results/frequency_strata'
min_snp_per_cluster <- 2L
min_clusters_to_report <- 20L

req<-c(phase_file,cluster_file,weight_file,
 file.path(base_dir,'33_All_SNPs_all_GA_Ne1000.csv'),
 file.path(base_dir,'33_All_SNPs_retained_GA_Ne1000.csv'))
if(any(!file.exists(req)))stop('Missing data; run 33 and check phase/predictor paths: ',
 paste(req[!file.exists(req)],collapse=', '))
phase<-read.csv(phase_file,stringsAsFactors=FALSE)
qtl<-read.csv(cluster_file,stringsAsFactors=FALSE,check.names=FALSE)
weights<-read.csv(weight_file,stringsAsFactors=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
effect_cols <- paste0('beta_',traits,'_NaCl')
if(!all(c('rs','group','phase')%in%names(phase)))stop('Invalid phase file')
if(anyDuplicated(phase$rs) || any(!phase$phase%in%c(-1,1)))
 stop('Phase rows must be unique; phase = +/-1')
if(!all(c('group','chr','beta_beta_NaCl',effect_cols)%in%names(qtl)))
 stop('Script-11 cluster file missing GWAS effects')
if(!all(c('trait','weight_full_GWAS_scale','weight_three_GWAS_scale')%in%names(weights)))
 stop('Script-13 gradient QC lacks original reconstructed-index weights')
weights <- weights[match(traits,weights$trait),,drop=FALSE]
if(anyNA(weights$trait))stop('Missing one of the seven trait weights')
if(anyDuplicated(qtl$group))stop('Duplicate GWAS cluster groups')
B<-as.matrix(qtl[,effect_cols,drop=FALSE]);storage.mode(B)<-'double'
qtl$pred_full_reconstructed <- drop(B%*%weights$weight_full_GWAS_scale)
qtl$pred_full_GWAS <- as.numeric(qtl$beta_beta_NaCl)
qtl$pred_three_trait <- drop(B%*%weights$weight_three_GWAS_scale)
qtl$chromosome <- as.character(qtl$chr)
preds <- qtl
phase$group<-as.character(phase$group)
preds$group<-as.character(preds$group)
preds$chromosome<-as.character(preds$chromosome)
## Validation: the 675 old common-frequency cluster predictors must match.
if(file.exists(predictor_QC_file)){
 old<-read.csv(predictor_QC_file,stringsAsFactors=FALSE)
 jj<-match(as.character(old$group),preds$group)
 if(anyNA(jj))stop('Old 13 predictors do not occur in full GWAS cluster file')
 for(col in c('pred_full_reconstructed','pred_full_GWAS','pred_three_trait')){
   diff<-abs(preds[[col]][jj]-old[[col]])
   if(any(is.finite(diff))&&max(diff[is.finite(diff)])>1e-7)
      stop('Reconstructed predictor differs from script 13: ',col)
 }
}
classes <- c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50')
results<-list();cluster_rows<-list()
models <- c(full_index='pred_full_reconstructed',direct_beta='pred_full_GWAS',
            three_trait='pred_three_trait',SB_only='beta_SB_NaCl')
safe_cor<-function(x,y)if(length(x)>=5&&is.finite(sd(x))&&is.finite(sd(y))&&sd(x)>0&&sd(y)>0)
 cor(x,y)else NA_real_
cv_LOCO<-function(x,y,chr){
 pred<-rep(NA_real_,length(y))
 for(cc in unique(chr)){
  tr<-chr!=cc;te<-!tr
  if(sum(tr)<15||sum(te)<1||!is.finite(sd(x[tr]))||sd(x[tr])<1e-12)next
  fit<-lm(y[tr]~x[tr])
  pred[te]<-unname(coef(fit)[1]+coef(fit)[2]*x[te])
 }
 good<-is.finite(pred)&is.finite(y)
 if(sum(good)<10||sd(y[good])<=0)return(NA_real_)
 1-sum((y[good]-pred[good])^2)/sum((y[good]-mean(y[good]))^2)
}
for(scenario in c('all_GA','retained_GA')){
 file<-file.path(base_dir,paste0('33_All_SNPs_',scenario,'_Ne1000.csv'))
 snp<-read.csv(file,stringsAsFactors=FALSE)
 if(anyDuplicated(snp$ID))stop('Unexpected duplicated SNP ID in 33 output')
 ii<-match(as.character(phase$rs),as.character(snp$ID))
 a<-phase[!is.na(ii),c('rs','group','phase'),drop=FALSE]
 ss<-snp[ii[!is.na(ii)],c('MAF_class','slope_freq_per_generation'),drop=FALSE]
 a$MAF_class<-ss$MAF_class
 a$slope_phase_aligned<-a$phase*ss$slope_freq_per_generation
 a<-a[is.finite(a$slope_phase_aligned)&!is.na(a$MAF_class),,drop=FALSE]
 for(bin in classes){
  zz<-a[a$MAF_class==bin,,drop=FALSE]
  if(!nrow(zz))next
  gg<-split(zz$slope_phase_aligned,zz$group)
  sub<-data.frame(group=names(gg),
       n_sites=vapply(gg,length,integer(1)),
       median_phase_slope=vapply(gg,median,numeric(1)),stringsAsFactors=FALSE)
  sub<-sub[sub$n_sites>=min_snp_per_cluster,,drop=FALSE]
  if(!nrow(sub))next
  j<-match(sub$group,preds$group)
  sub<-cbind(sub,preds[j,c('chromosome',unname(models)),drop=FALSE])
  sub$scenario<-scenario;sub$MAF_class<-bin
  cluster_rows[[length(cluster_rows)+1L]]<-sub
  for(m in names(models)){
    pred<-sub[[models[[m]]]];y<-sub$median_phase_slope;chr<-sub$chromosome
    valid<-is.finite(pred)&is.finite(y)&!is.na(chr)
    r<-if(sum(valid)>=min_clusters_to_report)safe_cor(pred[valid],y[valid])else NA_real_
    R2<-if(sum(valid)>=min_clusters_to_report&&length(unique(chr[valid]))>=3)
      cv_LOCO(pred[valid],y[valid],chr[valid]) else NA_real_
    results[[length(results)+1L]]<-data.frame(
      scenario=scenario,MAF_class=bin,model=m,n_clusters=sum(valid),
      median_n_SNP_per_cluster=median(sub$n_sites),
      n_chromosomes=length(unique(chr[valid])),
      correlation_r=r,R2_LOCO=R2,stringsAsFactors=FALSE)
  }
 }
}
if(!length(results))stop('No eligible clusters; inspect ID/phase mapping and thresholds')
out<-do.call(rbind,results)
write.csv(out,file.path(base_dir,'36_Prediction_by_MAF_class.csv'),row.names=FALSE)
write.csv(do.call(rbind,cluster_rows),
 file.path(base_dir,'36_Cluster_frequency_stratified_responses.csv'),row.names=FALSE)

cols<-c(full_index='#31688e',direct_beta='#35b779',three_trait='#de8f05',SB_only='#bc4b51')
pdf(file.path(base_dir,'36_Quantitative_prediction_across_MAF_classes.pdf'),width=10,height=7)
par(mfrow=c(2,2),mar=c(7,4.7,3,1))
for(s in c('all_GA','retained_GA'))for(metric in c('correlation_r','R2_LOCO')){
  plot(NA,xlim=c(.8,length(classes)+.2),ylim=c(-.5,.8),xaxt='n',
       xlab='',ylab=metric,main=paste(s,metric))
  axis(1,at=seq_along(classes),labels=classes,las=2,cex.axis=.8)
  abline(h=0,lty=2,col='grey60')
  for(m in names(cols)){
    zz<-out[out$scenario==s&out$model==m,,drop=FALSE]
    j<-match(classes,zz$MAF_class)
    yy<-zz[[metric]][j]
    lines(seq_along(classes),yy,type='b',col=cols[[m]],pch=16,lwd=1.7)
  }
  if(s=='all_GA'&&metric=='correlation_r')
    legend('topright',legend=names(cols),col=cols,pch=16,lty=1,bty='n',cex=.72)
}
dev.off()
cat('Script 36 complete. Interpretation is exploratory; no bin-specific q-values\n',
    'and no claim of independence among LD clusters counted in multiple bins.\n',sep='')
