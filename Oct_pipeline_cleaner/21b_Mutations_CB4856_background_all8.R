## 21_Mutations_CB4856_background_all8.R
## Two chrIV variants, CB4856 local ancestry, and their interactions across
## SF SB FS FB BS BF area.F and the derived fitness index beta.
## This is a targeted, kinship-adjusted BLUE-level analysis, not proof of causality.
## Run from your project root. All outputs are new; nothing upstream is overwritten.

suppressPackageStartupMessages(library(AGHmatrix))

## ------------------------- USER SETTINGS -----------------------------
preparation_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
ancestry_file <- '21_CB4856_local_ancestry.csv'  # Optional. See README.
origin_file <- '19_RIL_source_populations.csv' # Optional. TRUE provenance only.
gradient_file <- 'RData/selection_gradients_NaCl.csv'
output_dir <- 'mutation_effects_chrIV/CB4856_all8'
positions <- c(8892663L,10499135L)
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
min_lines <- 40L
min_background_lines <- 6L  # QC flag; not a magic criterion for statistical power
alpha <- 0.05
set.seed(20261029)
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

if (!file.exists(preparation_file) &&
    file.exists(paste0(preparation_file,'.RData')))
  preparation_file <- paste0(preparation_file,'.RData')
if (!file.exists(preparation_file)) stop('Cannot find preparation RData')

e <- new.env(parent=emptyenv())
load(preparation_file,envir=e)
## Prefer ONE adjusted RIL BLUE row per line. The mapping table may contain
## repeated phenotype observations and must not be mistaken for line-level BLUEs.
ph_name <- intersect(c('phen_NaCl','all_BLUEs_genotyped_NaCl',
                       'data_for_GWAS_NaCl_genotyped_NaCl'),ls(e))
if (!length(ph_name)) stop('No NaCl phenotype table in preparation RData')
ph <- e[[ph_name[1]]]
seven <- setdiff(traits,'beta')
if (!all(c('pop_label',seven) %in% names(ph)))
  stop('Phenotype object lacks the seven traits or pop_label')
if(anyDuplicated(as.character(ph$pop_label))) {
  warning('One-row-per-RIL BLUEs unavailable: collapsing replicate phenotype ',
          'observations by mean. These are NOT adjusted line BLUEs.')
  ph <- aggregate(ph[,intersect(traits,names(ph)),drop=FALSE],
                  by=list(pop_label=as.character(ph$pop_label)),
                  FUN=function(x)if(any(is.finite(x)))mean(x,na.rm=TRUE)else NA_real_)
  if('pop_label'%in%names(ph))names(ph)[names(ph)=='pop_label'] <- 'pop_label'
}
## Reconstruct beta for line-level BLUEs using the EXACT weights and
## GWAS-preparation standardization when it was saved in script 00.
if(!'beta'%in%names(ph)) {
  needed<-c('beta_for_mapping','replicate_mean','replicate_sd')
  if(!all(needed%in%ls(e)))
    stop('One-row-per-RIL BLUEs lack beta; need beta_for_mapping, ',
         'replicate_mean and replicate_sd from recomputed script 00')
  weights<-e$beta_for_mapping[seven]
  means<-e$replicate_mean[seven]
  scales<-e$replicate_sd[seven]
  if(any(!is.finite(c(weights,means,scales)))||any(scales<=0))
    stop('Invalid stored beta reconstruction weights')
  zz<-sweep(as.matrix(ph[,seven,drop=FALSE]),2,means,'-')
  zz<-sweep(zz,2,scales,'/')
  ph$beta<-drop(zz%*%weights)
  cat('Reconstructed beta for adjusted RIL BLUEs from saved GWAS weights\n')
}
if (!('genotype_matrix' %in% ls(e)) || !('genotype_pruned_matrix' %in% ls(e)))
  stop('Preparation requires genotype_matrix and genotype_pruned_matrix')

G <- e$genotype_matrix
P <- e$genotype_pruned_matrix
ph$pop_label <- as.character(ph$pop_label)
if (anyDuplicated(ph$pop_label)) stop('Expected one BLUE row per RIL: duplicate pop_label')
ids <- intersect(rownames(G),intersect(rownames(P),ph$pop_label))
if (length(ids)<min_lines) stop('Insufficient RILs shared by phenotype/genotype')
ph <- ph[match(ids,ph$pop_label),,drop=FALSE]
P <- P[ids,,drop=FALSE]
find_marker <- function(pos) {
  target <- paste0('IV_',pos)
  if (target %in% colnames(G)) return(target)
  options <- grep(paste0('^IV[_:.]',pos,'($|[_:.])'),colnames(G),value=TRUE)
  if(length(options)!=1L)stop('Cannot uniquely locate IV:',pos,' in full RIL matrix')
  options
}
markers <- vapply(positions,find_marker,character(1))
raw <- G[ids,markers,drop=FALSE]
rm(G);e$genotype_matrix<-NULL;gc(FALSE)
colnames(raw) <- c('mut1','mut2')
if(any(is.finite(raw)&(raw<0|raw>2)))stop('Mutation genotypes outside 0..2')
if(any(is.finite(raw)&abs(raw-round(raw))>1e-5))
  warning('Mutation dosages are noninteger; binary genotype-class QC is approximate')
## Dosage 0/2 in inbred lines -> homozygote-ALT contrast 0/1.
mut <- raw/2

## TRUE RIL origin labels, if supplied. Never infer origin from genotype.
origin <- rep('unknown',length(ids))
if(file.exists(origin_file)) {
  d <- read.csv(origin_file,stringsAsFactors=FALSE)
  if (!all(c('pop_label','source_population') %in% names(d)) || anyDuplicated(d$pop_label))
    stop('Invalid RIL origin CSV: need unique pop_label,source_population')
  origin <- as.character(d$source_population[match(ids,d$pop_label)])
  origin[is.na(origin)|!nzchar(origin)] <- 'unknown'
} else warning('RIL origin not supplied: population-origin confounding cannot be checked')

## Local CB4856 ancestry. Two separate measurements are important because
## a recombinant segment may switch founder identity between the two positions.
H <- matrix(NA_real_,nrow=length(ids),ncol=2,
            dimnames=list(ids,c('CB4856_8892663','CB4856_10499135')))
H_available <- file.exists(ancestry_file)
if(H_available) {
  a <- read.csv(ancestry_file,stringsAsFactors=FALSE,check.names=FALSE)
  req <- c('pop_label',colnames(H))
  if (!all(req %in% names(a)) || anyDuplicated(a$pop_label))
    stop('Ancestry CSV must have unique pop_label and two CB4856_ columns')
  for(j in seq_len(2))H[,j]<-as.numeric(a[[colnames(H)[j]]][match(ids,a$pop_label)])
  if(any(is.finite(H)&(H<0|H>1)))
    stop('CB4856 local ancestry must be 0..1 (indicator or validated probability)')
} else {
  warning('No CB4856 ancestry file. Mutation-only tests run; founder-background ',
          'tests are marked not_estimable. A per-RIL template will be created.')
  write.csv(data.frame(pop_label=ids,CB4856_8892663=NA_real_,
    CB4856_10499135=NA_real_),
    file.path(output_dir,'21_CB4856_local_ancestry_TEMPLATE.csv'),row.names=FALSE)
}

## Save genotype/origin/ancestry counts, including missing classes.
q <- data.frame(pop_label=ids,source_population=origin,
                mut1_dosage=raw[,1],mut2_dosage=raw[,2],
                CB4856_8892663=H[,1],CB4856_10499135=H[,2])
write.csv(q,file.path(output_dir,'21_Genotype_and_ancestry_QC.csv'),row.names=FALSE)
gtab <- as.data.frame(table(factor(round(raw[,1]),levels=0:2),
                            factor(round(raw[,2]),levels=0:2)))
names(gtab) <- c('mut1_dosage','mut2_dosage','n_RILs')
write.csv(gtab,file.path(output_dir,'21_Mutation_two_locus_counts.csv'),row.names=FALSE)
if(H_available) {
  cells <- data.frame(mut1=factor(round(mut[,1])),mut2=factor(round(mut[,2])),
                      H1=factor(ifelse(is.finite(H[,1]),round(H[,1]),NA)),
                      H2=factor(ifelse(is.finite(H[,2]),round(H[,2]),NA)))
  tab <- as.data.frame(table(cells,useNA='ifany'))
  write.csv(tab,file.path(output_dir,'21_Mutation_by_CB4856_cells.csv'),row.names=FALSE)
}

## Kinship excludes every chromosome-IV marker; local ancestry is included as
## an explicit fixed effect rather than absorbed into a chrIV-containing GRM.
keep_P <- !grepl('^IV([_:.]|$)',colnames(P))
if(sum(keep_P)<100L)stop('Fewer than 100 chromosome-IV-excluded pruned markers')
K <- as.matrix(AGHmatrix::Gmatrix(P[,keep_P,drop=FALSE],method='VanRaden'))
if(is.null(rownames(K)) || !all(ids%in%rownames(K)))stop('Missing kinship IDs')
K <- K[ids,ids,drop=FALSE]
K <- (K+t(K))/2
rm(P);e$genotype_pruned_matrix<-NULL;gc(FALSE)

## Fixed variance ratio under a baseline model, then nested GLS tests using
## the same whitening for all compared models (EMMAX-style approximation).
make_whitener <- function(y,K,X0) {
  eig <- eigen(K,symmetric=TRUE)
  U <- eig$vectors;lambda <- pmax(eig$values,0)
  yy <- as.vector(crossprod(U,y))
  xx <- crossprod(U,X0)
  n <- length(y);p <- ncol(X0)
  reml <- function(log_delta){
    v<-lambda+exp(log_delta)
    z<-1/v
    XX<-crossprod(xx,xx*z)
    if(qr(XX)$rank<ncol(XX))return(Inf)
    b<-tryCatch(solve(XX,crossprod(xx,yy*z)),error=function(e) NULL)
    if(is.null(b))return(Inf)
    rss<-sum((yy-drop(xx%*%b))^2*z)
    if(!is.finite(rss)||rss<=0)return(Inf)
    ld<-determinant(XX,logarithm=TRUE)
    if(ld$sign<=0)return(Inf)
    (n-p)*log(rss/(n-p))+sum(log(v))+as.numeric(ld$modulus)
  }
  delta <- exp(optimize(reml,c(-9,9))$minimum)
  scale <- 1/sqrt(lambda+delta)
  list(y=yy*scale,
       transform=function(X) sweep(crossprod(U,X),1,scale,'*'),
       delta=delta)
}
fit_w <- function(y,X) {
  n<-length(y);p<-ncol(X)
  if(n<=p+3||qr(X,tol=1e-8)$rank<p)return(NULL)
  Q<-qr(X,tol=1e-8)
  b<-as.vector(qr.coef(Q,y))
  resid<-y-drop(X%*%b)
  rss<-sum(resid^2)
  if(!is.finite(rss)||rss<=0)return(NULL)
  xx<-crossprod(X)
  inv<-tryCatch(solve(xx),error=function(e)NULL)
  if(is.null(inv))return(NULL)
  se<-sqrt(pmax(0,diag(inv)*rss/(n-p)))
  list(b=b,se=se,rss=rss,df=n-p,p=p,condition=kappa(X))
}
add_test <- function(small,big) {
  if(is.null(small)||is.null(big)||big$p<=small$p)return(NA_real_)
  df1<-big$p-small$p
  f<-max(0,(small$rss-big$rss)/df1)/(big$rss/big$df)
  pf(f,df1,big$df,lower.tail=FALSE)
}
term_est <- function(f,X,name) {
  if(is.null(f)||!name%in%colnames(X))
    return(c(beta=NA_real_,SE=NA_real_,p=NA_real_,min80_SD=NA_real_))
  j<-match(name,colnames(X))
  b<-f$b[j]; se<-f$se[j]
  c(beta=b,SE=se,
    p=if(is.finite(se)&&se>0)2*pt(abs(b/se),f$df,lower.tail=FALSE)else NA_real_,
    min80_SD=NA_real_)
}

## p-values are approximate and require targeted validation. Models with
## unidentifiable interactions report NA, NEVER arbitrarily drop collinear terms.
records <- list();power <- list();decomposition <- list();modelqc<-list()
for(trait in traits) {
  yy<-suppressWarnings(as.numeric(as.character(ph[[trait]])))
  complete_mut <- is.finite(yy)&complete.cases(mut)
  if(sum(complete_mut)<min_lines){warning('Insufficient RILs for ',trait);next}
  ## Mutation-only analysis uses every eligible RIL even if haplotype ancestry
  ## is missing. The background analysis is a separately labeled subset.
  for(analysis in c('mutation_only','with_CB4856')) {
    if(analysis=='with_CB4856'&&!H_available)next
    keep <- complete_mut
    if(analysis=='with_CB4856')keep <- keep&complete.cases(H)
    n <- sum(keep)
    if(n<min_lines){
      warning('Insufficient complete data for ',trait,' ',analysis,' (n=',n,')')
      next
    }
    y<-yy[keep]
    cov<-data.frame(M1=as.vector(mut[keep,1]),M2=as.vector(mut[keep,2]),
                    H1=H[keep,1],H2=H[keep,2])
    origin_ <- factor(origin[keep]); origin_df<-data.frame(origin=origin_)
    base <- if(nlevels(origin_)>1)model.matrix(~origin,data=origin_df) else
      matrix(1,n,1,dimnames=list(NULL,'(Intercept)'))
    ## Ensure comparisons use the exact SAME RILs and kinship for a given
    ## phenotype and analysis family.
    bnames<-colnames(base)
    make_design <- function(terms){
      x<-base
      for(nm in terms){
        x<-cbind(x,switch(nm,
          H1=cov$H1,H2=cov$H2,M1=cov$M1,M2=cov$M2,
          M1xM2=cov$M1*cov$M2,
          M1xH1=cov$M1*cov$H1,M2xH2=cov$M2*cov$H2))
        colnames(x)[ncol(x)]<-nm
      }
      x
    }
    X_base <- make_design(character(0))
    ## Two nearby ancestry columns can be identical or non-variable. Include
    ## the second only if it adds a linearly independent contrast.
    terms0 <- character(0)
    if(analysis=='with_CB4856')for(h in c('H1','H2')) {
      trial <- make_design(c(terms0,h))
      if(qr(trial,tol=1e-8)$rank==ncol(trial)) terms0<-c(terms0,h)
    }
    if(analysis=='with_CB4856'&&!length(terms0))
      warning('No estimable CB4856 ancestry variation for ',trait)

    X_H <- make_design(terms0)
    X_M1 <- make_design(c(terms0,'M1'))
    X_M2 <- make_design(c(terms0,'M2'))
    X_A <- make_design(c(terms0,'M1','M2'))
    X_MM <- make_design(c(terms0,'M1','M2','M1xM2'))
    X_H1M1 <- if(analysis=='with_CB4856')
      make_design(c(terms0,'M1','M2','M1xH1')) else NULL
    X_H2M2 <- if(analysis=='with_CB4856')
      make_design(c(terms0,'M1','M2','M2xH2')) else NULL
    X_HM <- if(analysis=='with_CB4856')
      make_design(c(terms0,'M1','M2',
                    'M1xH1','M2xH2')) else NULL
    X_all <- if(analysis=='with_CB4856')
      make_design(c(terms0,'M1','M2','M1xM2',
                    'M1xH1','M2xH2')) else NULL
    if(qr(X_base)$rank<ncol(X_base))stop('Origin covariates are rank-deficient')
    w <- make_whitener(y,K[keep,keep,drop=FALSE],X_base)
    F <- function(X) if(!is.null(X))fit_w(w$y,w$transform(X))else NULL
    f0<-F(X_base);fH<-F(X_H);f1<-F(X_M1);f2<-F(X_M2)
    fA<-F(X_A);fMM<-F(X_MM);fHM<-F(X_HM);fAll<-F(X_all)
    fH1M1<-F(X_H1M1);fH2M2<-F(X_H2M2)
    ## For each SNP, the conditional estimate requires both SNPs, plus
    ## the CB4856 local backgrounds where available.
    for(term in c('M1_marginal','M2_marginal','M1','M2','H1','H2','M1xM2','M1xH1','M2xH2')) {
      if(term%in%c('H1','H2')&&analysis!='with_CB4856')next
      if(term%in%c('M1xH1','M2xH2')&&analysis!='with_CB4856')next
      if(term=='M1_marginal') {f<-f1;X<-X_M1}
      if(term=='M2_marginal') {f<-f2;X<-X_M2}
      if(term=='M1') { f<-fA;X<-X_A }
      if(term=='M2') { f<-fA;X<-X_A }
      if(term%in%c('H1','H2')){f<-fH;X<-X_H}
      if(term=='M1xM2'){f<-fMM;X<-X_MM}
      if(term=='M1xH1'){f<-fH1M1;X<-X_H1M1}
      if(term=='M2xH2'){f<-fH2M2;X<-X_H2M2}
      actual_term <- sub('_marginal$','',term)
      if(is.null(X))X<-X_A
      zz<-term_est(f,X,actual_term)
      mde<-if(!is.null(f)&&is.finite(zz['SE'])&&zz['SE']>0)
        (qt(1-alpha/2,f$df)+qt(.8,f$df))*zz['SE']/sd(y)else NA_real_
      records[[length(records)+1L]] <- data.frame(
        trait=trait,analysis=analysis,term=term,n_RIL=n,
        beta_as_homozygote_contrast=unname(zz['beta']),
        SE=unname(zz['SE']),p_wald=unname(zz['p']),
        minimum_effect_80pct_power_trait_SD=mde,
        model_estimable=!is.null(f),
        model_condition=if(!is.null(f))f$condition else NA_real_,
        source_pop_adjusted=nlevels(origin_)>1,
        two_locus_r2=if(sd(cov$M1)>0&&sd(cov$M2)>0)
          cor(cov$M1,cov$M2)^2 else NA_real_,
        stringsAsFactors=FALSE)
    }
    modelqc[[length(modelqc)+1L]] <- data.frame(trait=trait,
      analysis=analysis,n_RIL=n,
      n_M1_ALT=sum(cov$M1>0),n_M2_ALT=sum(cov$M2>0),
      n_M1_only=sum(cov$M1>0&cov$M2==0),
      n_M2_only=sum(cov$M2>0&cov$M1==0),
      n_both=sum(cov$M1>0&cov$M2>0),
      n_source_groups=nlevels(origin_),
      CB4856_carriers_site1=if(analysis=='with_CB4856')sum(cov$H1>.5)else NA_integer_,
      CB4856_carriers_site2=if(analysis=='with_CB4856')sum(cov$H2>.5)else NA_integer_,
      correlation_M1_H1=if(analysis=='with_CB4856'&&sd(cov$M1)>0&&sd(cov$H1)>0)
        cor(cov$M1,cov$H1)else NA_real_,
      correlation_M2_H2=if(analysis=='with_CB4856'&&sd(cov$M2)>0&&sd(cov$H2)>0)
        cor(cov$M2,cov$H2)else NA_real_,
      H_model_estimable=!is.null(fH),
      adjusted_two_mutations_estimable=!is.null(fA),
      M1_M2_epistasis_estimable=!is.null(fMM),
      M_by_CB4856_both_estimable=!is.null(fHM),
      full_interaction_model_estimable=!is.null(fAll),
      p_mutations_given_H=add_test(fH,fA),
      p_H_given_origin=if(analysis=='with_CB4856')add_test(f0,fH)else NA_real_,
      p_M1M2_given_additive=add_test(fA,fMM),
      p_MxH_given_additive=add_test(fA,fHM),
      p_all_interactions_given_additive=add_test(fA,fAll),
      estimable_H_covariates=paste(terms0,collapse=';'),
      variance_ratio_delta=w$delta,
      stringsAsFactors=FALSE)
    ## Optional decomposition of beta into its seven trait components.
    ## Use complete cases, identical GLS whitening and design. Since beta
    ## is a linear combination, contributions must sum to the beta coefficient.
    if(trait=='beta') {
      seven<-setdiff(traits,'beta')
      candidates<-list()
      if(all(c('beta_for_mapping','replicate_sd')%in%ls(e)))
        candidates[['saved_beta_for_mapping_BLUE_raw']] <-
          e$beta_for_mapping[seven]/e$replicate_sd[seven]
      if(all(c('beta_NaCl','BLUE_sd')%in%ls(e)))
        candidates[['beta_NaCl_over_BLUE_sd']] <- e$beta_NaCl[seven]/e$BLUE_sd[seven]
      if('beta_for_mapping'%in%ls(e))
        candidates[['saved_beta_for_mapping']]<-e$beta_for_mapping
      if('beta_NaCl'%in%ls(e))
        candidates[['saved_beta_NaCl']]<-e$beta_NaCl
      if('selection_gradients'%in%ls(e) &&
         all(c('trait','beta_for_mapping')%in%names(e$selection_gradients)))
        candidates[['saved_gradient_table']]<-setNames(
          e$selection_gradients$beta_for_mapping,e$selection_gradients$trait)
      if(file.exists(gradient_file)) {
        tab<-tryCatch(read.csv(gradient_file,stringsAsFactors=FALSE),
                      error=function(e)NULL)
        if(!is.null(tab)&&all(c('trait','beta_for_mapping')%in%names(tab)))
          candidates[['gradient_csv']]<-setNames(tab$beta_for_mapping,tab$trait)
      }
      ## Check actual phenotype-level reconstruction rather than assuming the
      ## newest gradient file uses the same scale as a historical beta.
      wt<-NULL
      complete_seven<-all(vapply(seven,function(t)
        all(is.finite(as.numeric(ph[[t]][keep]))),logical(1)))
      mat<-as.matrix(ph[keep,seven,drop=FALSE]);storage.mode(mat)<-'double'
      if(complete_seven)for(nm in names(candidates)) {
        z<-candidates[[nm]]
        if(length(z)!=7)next
        if(!is.null(names(z)) && all(seven%in%names(z)))z<-z[seven]
        else if(is.null(names(z)))next
        z<-as.numeric(z);names(z)<-seven
        resid_index <- drop(mat%*%z)-y
        if(all(is.finite(z))&&all(is.finite(resid_index))&&
           max(abs(resid_index-mean(resid_index)))<1e-5){
          wt<-z
          break
        }
      }
      ok<-!is.null(wt)&&complete_seven
      if(ok){
        mat<-as.matrix(ph[keep,seven,drop=FALSE]); storage.mode(mat)<-'double'
        reconstruction<-drop(mat%*%wt)
        if(max(abs((reconstruction-y)-mean(reconstruction-y)))<1e-5){
          ## Shared whitening is based on beta, not on trait-specific variance
          ## components. This guarantees an exact additive decomposition.
          fx<-if(!is.null(fAll))X_all else if(!is.null(fMM))X_MM else X_A
          fB<-F(fx)
          if(!is.null(fB)) {
            wx<-w$transform(fx)
            for(t in seven){
              fT<-fit_w(w$transform(matrix(mat[,t],ncol=1))[,1],wx)
              if(is.null(fT))next
              for(tt in intersect(c('M1','M2','M1xM2','M1xH1','M2xH2'),colnames(fx))){
                j<-match(tt,colnames(fx))
                decomposition[[length(decomposition)+1L]]<-data.frame(
                  analysis=analysis,term=tt,trait=t,weight=wt[t],
                  component_effect=fT$b[j],weighted_contribution=wt[t]*fT$b[j],
                  beta_total_from_same_GLS=fB$b[j],n_RIL=n)
              }
            }
          }
        } else warning('Cannot reconstruct saved beta on complete rows; skipping beta decomposition')
      } else warning('No compatible beta weights/complete traits for beta decomposition')
    }
  }
  cat('Completed trait: ',trait,'\n',sep='')
}
if(!length(records))stop('No valid trait models')
res<-do.call(rbind,records)
qc<-do.call(rbind,modelqc)
## One BH adjustment over eight traits for EACH prespecified term/family;
## other tests exploratory and not familywise confirmatory.
res$q_BH_eight_traits <- NA_real_
for(fam in unique(paste(res$analysis,res$term,sep='::'))){
  ii<-which(paste(res$analysis,res$term,sep='::')==fam)
  res$q_BH_eight_traits[ii]<-p.adjust(res$p_wald[ii],method='BH')
}
write.csv(res,file.path(output_dir,'21_All_eight_trait_effects.csv'),row.names=FALSE)
write.csv(qc,file.path(output_dir,'21_Model_identifiability_and_power_QC.csv'),row.names=FALSE)
if(length(decomposition)){
  dd<-do.call(rbind,decomposition)
  write.csv(dd,file.path(output_dir,'21_Beta_trait_contributions.csv'),row.names=FALSE)
  agg<-aggregate(weighted_contribution~analysis+term,data=dd,sum)
  totals<-unique(dd[,c('analysis','term','beta_total_from_same_GLS')])
  comparison<-merge(agg,totals,by=c('analysis','term'))
  comparison$abs_difference<-abs(comparison$weighted_contribution-
                                   comparison$beta_total_from_same_GLS)
  write.csv(comparison,file.path(output_dir,'21_Beta_decomposition_check.csv'),row.names=FALSE)
}
## Plot the independently tested eight outcomes: effect estimates with
## approximate confidence intervals, not inferred causal effects.
pdf(file.path(output_dir,'21_Effects_eight_traits.pdf'),width=10,height=7)
par(mfrow=c(2,2),mar=c(4,6,3,1))
for(family in c('mutation_only','with_CB4856')) {
  if(!family%in%res$analysis)next
  for(term in c('M1','M2')){
    z<-subset(res,analysis==family&term==term)
    ok<-is.finite(z$beta_as_homozygote_contrast)&is.finite(z$SE)
    z<-z[ok,]
    if(!nrow(z))next
    xlim<-range(c(0,z$beta_as_homozygote_contrast+1.96*z$SE,
                  z$beta_as_homozygote_contrast-1.96*z$SE))
    plot(z$beta_as_homozygote_contrast,seq_len(nrow(z)),
         yaxt='n',ylim=c(.5,nrow(z)+.5),xlim=xlim,
         pch=16,xlab='Effect (BLUE units; ALT homozygote contrast)',ylab='',
         main=paste(family,term),bty='n')
    arrows(z$beta_as_homozygote_contrast-1.96*z$SE,seq_len(nrow(z)),
           z$beta_as_homozygote_contrast+1.96*z$SE,seq_len(nrow(z)),
           angle=90,code=3,length=.035)
    axis(2,at=seq_len(nrow(z)),labels=z$trait,las=2)
    abline(v=0,lty=2,col='grey50')
  }
}
dev.off()
writeLines(c(paste0('RILs analyzed in genotype table: ',length(ids)),
  paste0('CB4856 background available: ',H_available),
  paste0('Mutation genotype marker IDs: ',paste(markers,collapse=', ')),
  'Mutation dose is ALT allele count divided by two; use homozygote contrast.',
  'Kinship excluded chromosome IV. P-values are approximate BLUE-level GLS tests.',
  'An interaction that is collinear with founder ancestry is not estimable.',
  'SB, area.F and beta are outcomes, NOT proof of mechanistic trait mediation.'),
  file.path(output_dir,'21_Analysis_notes.txt'))
cat('Saved: ',output_dir,'\n',sep='')
