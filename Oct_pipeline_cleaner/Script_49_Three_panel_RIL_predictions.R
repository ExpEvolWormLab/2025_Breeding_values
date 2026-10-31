## 49 — Three-column, eight-trait manuscript comparison
## Column 1: original legacy R2=0.2 significant-marker score, IN SAMPLE
## Column 2: top-75 GWAS marker prediction, OUT OF FOLD
## Column 3: rrBLUP GBLUP prediction, OUT OF FOLD
## The first column is NOT a cross-validated predictor and is NOT R5-based.
## Run from project root: source('49_Three_panel_GWAS_score_vs_CV_predictions.R')

options(stringsAsFactors=FALSE)
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
legacy_threshold <- 0.2
output_dir <- 'phenotype_prediction/three_panel_comparison'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

## --- 1. Load the original full-data GWAS marker effects and RIL genotypes ----
cluster_file <- sprintf('RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R%d_with_residuals.RData',
                        round(10*legacy_threshold))
geno_file <- sprintf('data_BLUE_GridLMM_NaCl/genotypes_matrix_clusters_R%d.txt',
                     round(10*legacy_threshold))
phenotype_file <- 'RData_old/Gmatrices_P.RData'
for(f in c(cluster_file,geno_file,phenotype_file))
  if(!file.exists(f))stop('Required original analysis input is missing: ',f)

gwas_env <- new.env(parent=globalenv()); load(cluster_file,envir=gwas_env)
if(!exists('all_GWAS_clusters',envir=gwas_env,inherits=FALSE))
  stop('The legacy GWAS RData must contain all_GWAS_clusters')
clust <- get('all_GWAS_clusters',envir=gwas_env)
if(!'rs'%in%names(clust))stop('all_GWAS_clusters has no rs (marker ID) column')
clust$rs <- as.character(clust$rs)
if(anyNA(clust$rs)||anyDuplicated(clust$rs))
  stop('Missing/duplicate marker IDs in all_GWAS_clusters$rs; resolve before scoring')

## Original file layout: rows = SNPs; columns = RILs; rownames are SNP IDs.
G <- as.matrix(read.table(geno_file,check.names=FALSE,
                          stringsAsFactors=FALSE))
if(is.null(rownames(G))||is.null(colnames(G)))stop('Legacy genotype matrix lacks SNP/RIL IDs')
if(anyDuplicated(rownames(G))||anyDuplicated(colnames(G)))
  stop('Duplicate SNP or RIL names in legacy genotype matrix')
if(!all(clust$rs%in%rownames(G)))
  stop('Some legacy GWAS SNP IDs do not occur in genotype matrix (check file header and marker IDs): ',
       paste(head(setdiff(clust$rs,rownames(G)),10),collapse=', '))
G <- G[clust$rs,,drop=FALSE]  # explicit SNP-ID match, never rely on sorting
suppressWarnings(storage.mode(G) <- 'double')
if(any(!is.finite(G)))stop('Legacy genotype matrix has missing/non-numeric values: inspect coding')

phen_env <- new.env(parent=globalenv());load(phenotype_file,envir=phen_env)
if(!exists('all_traits_corrected',envir=phen_env,inherits=FALSE))
  stop('Phenotype RData must contain all_traits_corrected')
original_pheno <- get('all_traits_corrected',envir=phen_env)
if(!'pop_label'%in%names(original_pheno))stop('Missing pop_label in all_traits_corrected')
original_pheno$pop_label <- as.character(original_pheno$pop_label)
if(anyDuplicated(original_pheno$pop_label))stop('Duplicate RIL IDs in original phenotypes')

## Construct original genetic scores (no phenotype intercept); missing phenotypes
## do NOT affect scoring. The score is purposely NOT refitted or rescaled.
score <- matrix(NA_real_,nrow=ncol(G),ncol=length(traits),
                dimnames=list(colnames(G),traits))
selection_count <- setNames(integer(length(traits)),traits)
for(tr in traits){
  sig_col <- paste0('is_sign_',tr,'_NaCl')
  effect_col <- paste0('beta_',tr,'_NaCl')
  if(!all(c(sig_col,effect_col)%in%names(clust)))
    stop('Missing full-data GWAS column(s) for ',tr)
  keep <- as.logical(clust[[sig_col]])
  if(anyNA(keep))stop('NA significance flags for ',tr)
  effects <- suppressWarnings(as.numeric(as.character(clust[[effect_col]])))
  if(any(keep & !is.finite(effects)))stop('Nonfinite selected effect estimates: ',tr)
  selection_count[tr] <- sum(keep)
  score[,tr] <- if(sum(keep)==0) rep(0,ncol(G)) else
    as.vector(crossprod(G[keep,,drop=FALSE],effects[keep]))
}
write.csv(data.frame(trait=traits,n_selected_marker_rows=as.integer(selection_count),
                     cluster_definition='legacy_R2_0.2_full_data'),
          file.path(output_dir,'49_Legacy_marker_selection_counts.csv'),row.names=FALSE)

## --- 2. Load already matched, genuinely held-out predictions ----------------
## Prefer the table produced by script 47; can also read a standalone copy.
candidates <- c('phenotype_prediction/GBLUP_rrBLUP/47_Matched_random_RIL_out_of_fold_predictions.csv',
 'phenotype_prediction/GBLUP_rrBLUP/47_RIL_observed_and_two_CV_predictions.csv',
 'phenotype_prediction/GBLUP_rrBLUP/40_Out_of_fold_GBLUP_predictions.csv')
cv_file <- candidates[file.exists(candidates)][1]
if(is.na(cv_file))stop('No matched script-47 or GBLUP held-out prediction CSV found. See README.')
cv <- read.csv(cv_file,check.names=FALSE,stringsAsFactors=FALSE)
if('scheme'%in%names(cv))cv <- cv[cv$scheme=='random_RIL',,drop=FALSE]
if(!all(c('pop_label','trait','observed','predicted_GBLUP')%in%names(cv)))
  stop('Held-out CSV lacks required RIL, trait, observed, or GBLUP columns')
if(!'predicted_top75'%in%names(cv)){
  stop('CSV has no predicted_top75 column. Re-run script 47 or supply its matched output.')
}
cv$pop_label <- as.character(cv$pop_label)
cv$trait <- as.character(cv$trait)
cv <- cv[cv$trait%in%traits,,drop=FALSE]
if(!nrow(cv)||anyDuplicated(paste(cv$pop_label,cv$trait,sep=':::')))
  stop('Empty or duplicated held-out RIL/trait rows')
cv$observed <- as.numeric(cv$observed)
cv$predicted_GBLUP <- as.numeric(cv$predicted_GBLUP)
cv$predicted_top75 <- as.numeric(cv$predicted_top75)
if(!all(unique(cv$pop_label)%in%rownames(score)))
  stop('Held-out RIL IDs absent from original legacy genotypes. Examples: ',
       paste(head(setdiff(unique(cv$pop_label),rownames(score))),collapse=', '))
if(!all(unique(cv$pop_label)%in%original_pheno$pop_label))
  stop('Held-out RIL IDs absent from original phenotypes')

## Add in-sample score by explicit RIL IDs, and independently verify observed
## phenotypes refer to the same data as the original legacy scatterplot.
cv$in_sample_GWAS_score <- mapply(function(id,tr)score[id,tr],cv$pop_label,cv$trait)
cv$original_observed <- NA_real_
for(tr in traits){
  ids <- which(cv$trait==tr)
  ph_col <- paste0(tr,'_NaCl')
  if(!ph_col%in%names(original_pheno))stop('Missing original phenotype column ',ph_col)
  cv$original_observed[ids] <- as.numeric(original_pheno[[ph_col]][
                          match(cv$pop_label[ids],original_pheno$pop_label)])
}
## Check and report, not stop: adjusted phenotype versions may differ.
paired_ok <- is.finite(cv$observed)&is.finite(cv$original_observed)
delta <- if(any(paired_ok))max(abs(cv$observed[paired_ok]-cv$original_observed[paired_ok])) else NA_real_
cat('Maximum difference between legacy and CV observed phenotypes: ',delta,'\n')
if(is.finite(delta)&&delta>0.01)
 warning('The old and new phenotype values differ; panels will use their respective observed values. ',
         'Interpret comparisons as partly different phenotype definitions.')
write.csv(cv,file.path(output_dir,'49_Individual_RIL_scores_and_predictions.csv'),row.names=FALSE)

## --- 3. Statistics: distinguish fitted correlation from held-out R2 ----------
measure <- function(y,p,heldout){
 ok <- is.finite(y)&is.finite(p);y<-y[ok];p<-p[ok]
 if(length(y)<5||sd(y)==0||sd(p)==0)
  return(c(n=length(y),cor=NA_real_,cor_squared=NA_real_,R2_predictive=NA_real_,
           sd_ratio=NA_real_,RMSE=NA_real_))
 c(n=length(y),cor=cor(y,p),cor_squared=cor(y,p)^2,
   R2_predictive=if(heldout)1-sum((y-p)^2)/sum((y-mean(y))^2)else NA_real_,
   sd_ratio=sd(p)/sd(y),RMSE=sqrt(mean((p-y)^2)))
}
res <- list()
for(tr in traits){
 z <- cv[cv$trait==tr,,drop=FALSE]
 for(model in c('in_sample_GWAS_score','predicted_top75','predicted_GBLUP')){
  if(tr=='beta'&&model=='predicted_top75')next
  y <- if(model=='in_sample_GWAS_score')z$original_observed else z$observed
  res[[length(res)+1L]] <- data.frame(trait=tr,method=model,
     t(measure(y,z[[model]],heldout=model!='in_sample_GWAS_score')),
     validation=if(model=='in_sample_GWAS_score')'in_sample' else 'random_RIL_out_of_fold',
     row.names=NULL)
 }
}
metrics <- do.call(rbind,res)
write.csv(metrics,file.path(output_dir,'49_Three_panel_accuracy_metrics.csv'),row.names=FALSE)

## --- 4. Plot eight traits x three columns, over several legible PDF pages ----
## Important: first panel has an independent y scale because the original
## score has no intercept/calibration. CV panels use identity axes.
color <- c(in_sample_GWAS_score='#454545',predicted_top75='#d17a36',predicted_GBLUP='#2878a8')
model_titles <- c(in_sample_GWAS_score='Original GWAS score (in-sample)',
                  predicted_top75='75 GWAS markers (held-out)',
                  predicted_GBLUP='GBLUP (held-out)')
plot_one <- function(z,tr,model){
 y <- if(model=='in_sample_GWAS_score')z$original_observed else z$observed
 p <- z[[model]]
 good <- is.finite(y)&is.finite(p);y<-y[good];p<-p[good]
 if(length(y)<5){plot.new();title(main=paste(tr,' — unavailable'));return(invisible(NULL))}
 if(model=='in_sample_GWAS_score'){
  plot(y,p,pch=16,cex=.4,col=adjustcolor(color[model],alpha.f=.27),
       xlab='Observed RIL phenotype',ylab='Genetic score (uncalibrated)',
       main=paste(tr,'—',model_titles[model]),cex.main=.8)
  abline(lm(p~y),col=color[model],lwd=1.7)
  lab <- sprintf('In-sample r² = %.3f',cor(y,p)^2)
 }else{
  axes <- range(c(y,p),finite=TRUE)
  plot(y,p,pch=16,cex=.4,col=adjustcolor(color[model],alpha.f=.30),
       xlim=axes,ylim=axes,xlab='Observed RIL phenotype',
       ylab='Held-out prediction',main=paste(tr,'—',model_titles[model]),cex.main=.8)
  abline(0,1,lty=2,col='grey45')
  r2 <- 1-sum((y-p)^2)/sum((y-mean(y))^2)
  lab <- sprintf('Held-out R² = %.3f',r2)
 }
 legend('topleft',legend=lab,bty='n',cex=.78,text.col=color[model])
}

pdf(file.path(output_dir,'49_Three_panel_eight_traits.pdf'),
    width=14.4,height=9.2,onefile=TRUE)
for(batch in split(traits,ceiling(seq_along(traits)/2))){
 op <- par(mfrow=c(2,3),mar=c(4.25,4.6,3.2,1.1),oma=c(1.1,.6,1.3,.3))
 for(tr in batch){
  z <- cv[cv$trait==tr,,drop=FALSE]
  for(model in names(model_titles)){
   if(tr=='beta'&&model=='predicted_top75'){
    plot.new();text(.5,.55,'No matched out-of-fold\n75-marker beta model',cex=1.1)
    title(main=paste(tr,'—',model_titles[model]),cex.main=.8)
   }else plot_one(z,tr,model)
  }
 }
 mtext('Left: full-data legacy R²=0.2 GWAS score (in-sample)  |  Middle/right: random-RIL held-out prediction',
       side=3,outer=TRUE,line=.2,cex=.9,font=2)
 par(op)
}
dev.off()

## Additionally create a single-page four-trait selection for main-text layout.
focus <- c('SF','SB','area.F','beta')
pdf(file.path(output_dir,'49_Three_panel_selected_traits.pdf'),width=14.4,height=14.3)
op <- par(mfrow=c(4,3),mar=c(4,4.4,3,1),oma=c(.4,.5,1.3,.2))
for(tr in focus){
 z<-cv[cv$trait==tr,,drop=FALSE]
 for(model in names(model_titles)){
  if(tr=='beta'&&model=='predicted_top75'){
   plot.new();text(.5,.5,'No held-out top-75\nbeta predictions');title(main=paste(tr,'—',model_titles[model]),cex.main=.8)
  }else plot_one(z,tr,model)
 }
}
mtext('GWAS phenotype reconstruction versus independent genomic prediction',side=3,outer=TRUE,line=.2,font=2)
par(op);dev.off()
cat('Completed three-panel comparisons. Outputs in: ',output_dir,'\n',sep='')
