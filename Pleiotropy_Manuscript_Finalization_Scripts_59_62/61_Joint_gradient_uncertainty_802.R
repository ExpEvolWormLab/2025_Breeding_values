# 61_Joint_gradient_uncertainty_802.R
# Aim: propagate the FULL JOINT COVARIANCE of the phenotypic fitness-gradient
# estimates into 802-cluster GWAS-effect predictions and LOCO performance.
# The covariance of lm(model_beta) comes from the ORIGINAL independent RIL
# fitness-phenotype fit (script 00), not from genomic evolutionary responses.
# NO rerunning MCMCglmm or GWAS. No GWAS or sampling-error propagation.
# R >= 4.1, base R only. Run from project root.

CONFIG <- list(
 prep='RData/Preparatory_GLMM_NaCl_only_withBeta.RData',
 gradient='RData/selection_gradients_NaCl.csv',
 weight='selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv',
 qtl='selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
 evolution='selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
 out='selection_results/manuscript_finalization/61_gradient_uncertainty',
 n_gradient_draws=2000L,seed=20261103L
)
resolve <- function(f){for(x in c(f,basename(f)))if(file.exists(x))return(x)
 stop('Missing required input ',f)}
dir.create(CONFIG$out,recursive=TRUE,showWarnings=FALSE)
traits<-c('SF','SB','FS','FB','BS','BF','area.F')
eff_cols<-paste0('beta_',traits,'_NaCl')
G<-read.csv(resolve(CONFIG$gradient),stringsAsFactors=FALSE,check.names=FALSE)
W<-read.csv(resolve(CONFIG$weight),stringsAsFactors=FALSE,check.names=FALSE)
Q<-read.csv(resolve(CONFIG$qtl),stringsAsFactors=FALSE,check.names=FALSE)
E<-read.csv(resolve(CONFIG$evolution),stringsAsFactors=FALSE,check.names=FALSE)
req <- function(z,n,name){if(length(setdiff(n,names(z))))stop(name,' missing ',paste(setdiff(n,names(z)),collapse=', '))}
req(G,c('trait','beta','BLUE_sd','replicate_sd'),'gradient')
req(W,c('trait','GWAS_pheno_SD','weight_full_GWAS_scale'),'weight')
req(Q,c('group','chr',eff_cols),'GWAS')
req(E,c('group','chr','ps','scenario','reason','slope','block','pred_full','LOCO_full'),'evolution')
if(anyDuplicated(G$trait)||anyDuplicated(W$trait)||anyDuplicated(Q$group))stop('Duplicate IDs')
G<-G[match(traits,G$trait),,drop=FALSE]
W<-W[match(traits,W$trait),,drop=FALSE]
if(anyNA(G$trait)||anyNA(W$trait))stop('Missing trait weights')
E<-E[E$scenario=='primary_survivors' & E$reason=='eligible' & is.finite(E$slope),,drop=FALSE]
if(anyDuplicated(E$group))stop('Duplicate eligible group')
if(nrow(E)!=802L)warning('Only ',nrow(E),' eligible groups, not 802')
i<-match(as.character(E$group),as.character(Q$group))
if(anyNA(i))stop('Cannot match R5 GWAS groups')
Q<-Q[i,,drop=FALSE]
if(any(as.character(E$chr)!=as.character(Q$chr)))stop('Chromosome mismatch')
X<-as.matrix(Q[,eff_cols,drop=FALSE]);storage.mode(X)<-'double';colnames(X)<-traits
if(any(!is.finite(X)))stop('Nonfinite GWAS effects')

prep<-new.env(parent=globalenv());load(resolve(CONFIG$prep),envir=prep)
if(!exists('model_beta',envir=prep,inherits=FALSE))stop('Preparation RData lacks original model_beta')
mb<-prep$model_beta
V<-stats::vcov(mb)[traits,traits,drop=FALSE]
V<-(V+t(V))/2
if(any(!is.finite(V))||inherits(try(chol(V),silent=TRUE),'try-error'))
 stop('The full joint selection-gradient covariance is not positive definite')
mu<-as.numeric(G$beta)
if(any(!is.finite(mu)))stop('Nonfinite reported selection gradients')
if(max(abs(unname(stats::coef(mb)[traits])-mu))>1e-6)
 stop('Original model_beta coefficients disagree with gradient CSV')
# The component-GWAS effects were mapped on standardized phenotype scales.
# Convert fitness gradients using the previously documented script-13 scaling.
conv<-as.numeric(G$replicate_sd/G$BLUE_sd)*as.numeric(W$GWAS_pheno_SD)
if(any(!is.finite(conv))||any(conv<=0))stop('Invalid GWAS scale factors')
w_ref<-as.numeric(W$weight_full_GWAS_scale)
w_unscaled<-mu*conv
scale_factor<-sum(w_unscaled*w_ref)/sum(w_unscaled^2)
if(!is.finite(scale_factor)||scale_factor<=0 ||
   max(abs(scale_factor*w_unscaled-w_ref))>1e-6)
 stop('Cannot reproduce frozen beta weights from original fitness model')
if(max(abs(drop(X%*%w_ref)-E$pred_full))>1e-9)
 stop('Frozen beta score does not match script 45')

set.seed(CONFIG$seed)
Z<-matrix(rnorm(CONFIG$n_gradient_draws*length(traits)),
          nrow=CONFIG$n_gradient_draws,ncol=length(traits))
BETA<-sweep(Z%*%chol(V),2,mu,'+')
colnames(BETA)<-traits
WG<-sweep(BETA,2,conv*scale_factor,'*')
# Each column is a candidate fitness-derived raw score. No evolutionary
# information used in the weights, only their subsequent calibration.
SCORE<-X%*%t(WG)
ch<-as.character(E$chr);y<-as.numeric(E$slope)
pp<-matrix(NA_real_,nrow=nrow(E),ncol=CONFIG$n_gradient_draws)
for(cc in unique(ch)){
 tr<-which(ch!=cc);te<-which(ch==cc)
 xtr<-SCORE[tr,,drop=FALSE];xte<-SCORE[te,,drop=FALSE]
 yy<-y[tr]; xm<-colMeans(xtr);ym<-mean(yy)
 xx<-sweep(xtr,2,xm,'-'); yc<-yy-ym
 v<-colSums(xx^2); b<-colSums(xx*yc)/v
 if(any(!is.finite(b)))stop('Nonfinite LOCO coefficients for some gradient draws')
 a<-ym-b*xm
 pp[te,]<-sweep(sweep(xte,2,b,'*'),2,a,'+')
}
if(any(!is.finite(pp)))stop('Incomplete chromosome-held-out predictions')
sst<-sum((y-mean(y))^2)
cor_vector <- function(M,y){
 my<-mean(y);dy<-y-my;mc<-colMeans(M)
 dx<-sweep(M,2,mc,'-');den<-sqrt(colSums(dx^2)*sum(dy^2))
 colSums(dx*dy)/den
}
correct_sign <- colMeans(sign(SCORE)==matrix(sign(y),nrow=length(y),
 ncol=CONFIG$n_gradient_draws))
pred_r2<-1-colSums((pp-y)^2)/sst
raw_r <- cor_vector(SCORE,y)
## Important: these across-gradient distributions are conditional on the
## observed genome, fixed allele orientation, GWAS effects and evolved data.
summary_draws<-data.frame(draw=seq_len(CONFIG$n_gradient_draws),
  GWAS_evolution_r=raw_r,LOCO_R2=pred_r2,
  raw_direction_accuracy=correct_sign,
  V3_raw_direction_accuracy=if(any(E$block=='V_3'))
   colMeans(sign(SCORE[E$block=='V_3',,drop=FALSE])==
    matrix(sign(y[E$block=='V_3']),sum(E$block=='V_3'),CONFIG$n_gradient_draws))
   else NA_real_)
stat_one <- function(a,name) data.frame(metric=name,
 point_frozen=NA_real_,mean=mean(a),median=median(a),
 q025=unname(quantile(a,.025)),q975=unname(quantile(a,.975)),
 fraction_positive=mean(a>0))
summary_stats <- rbind(stat_one(raw_r,'GWAS_evolution_r'),
 stat_one(pred_r2,'LOCO_R2'),stat_one(correct_sign,'raw_direction_accuracy'))
summary_stats$point_frozen<-c(cor(E$pred_full,y),
  1-sum((E$LOCO_full-y)^2)/sst,mean(sign(E$pred_full)==sign(y)))
# Pointwise sign stability among the focal clusters and for all 802 groups.
ss<-sign(SCORE)==matrix(sign(y),length(y),CONFIG$n_gradient_draws)
cluster_summary<-data.frame(group=E$group,chr=E$chr,ps=E$ps,
 block=E$block,mean_raw_score=rowMeans(SCORE),
 raw_score_q025=apply(SCORE,1,quantile,.025),
 raw_score_q975=apply(SCORE,1,quantile,.975),
 frac_correct_direction=rowMeans(ss),
 fraction_raw_score_positive=rowMeans(SCORE>0),
 observed_slope=y,frozen_raw_score=E$pred_full)
weight_summary<-data.frame(trait=traits,frozen_weight=w_ref,
  mean=colMeans(WG),lo95=apply(WG,2,quantile,.025),
  hi95=apply(WG,2,quantile,.975),
  fraction_positive=colMeans(WG>0))
write.csv(summary_draws,file.path(CONFIG$out,'61_Gradient_draw_metrics.csv'),row.names=FALSE)
write.csv(summary_stats,file.path(CONFIG$out,'61_Gradient_uncertainty_summary.csv'),row.names=FALSE)
write.csv(cluster_summary,file.path(CONFIG$out,'61_Cluster_direction_stability.csv'),row.names=FALSE)
write.csv(weight_summary,file.path(CONFIG$out,'61_Gradient_weight_posteriors.csv'),row.names=FALSE)
write.csv(V,file.path(CONFIG$out,'61_Selection_gradient_covariance.csv'))
pdf(file.path(CONFIG$out,'61_Gradient_uncertainty_figures.pdf'),width=10,height=4.8)
par(mfrow=c(1,2),mar=c(5,4.2,3,1))
lim<-range(c(weight_summary$lo95,weight_summary$hi95,0),finite=TRUE)
plot(seq_along(traits),weight_summary$frozen_weight,pch=19,xaxt='n',
 xlab='Trait',ylab='Fitness-gradient weight (GWAS scale)',ylim=lim,
 main='Joint gradient uncertainty')
axis(1,at=seq_along(traits),labels=traits)
arrows(seq_along(traits),weight_summary$lo95,seq_along(traits),weight_summary$hi95,
  code=3,angle=90,length=.04,col='gray50')
abline(h=0,lty=2,col='gray70')
hist(pred_r2,breaks=35,xlab='Chromosome-held-out predictive R-squared',
 main='Sensitivity to fitness-gradient uncertainty',col='gray80',border='white')
abline(v=summary_stats$point_frozen[summary_stats$metric=='LOCO_R2'],
 col='steelblue4',lwd=2)
legend('topright',legend='Frozen beta',lwd=2,col='steelblue4',bty='n')
dev.off()
writeLines(capture.output(sessionInfo()),file.path(CONFIG$out,'61_sessionInfo.txt'))
cat('61 COMPLETE. This propagates joint selection-gradient uncertainty ONLY.\n')
print(summary_stats)
