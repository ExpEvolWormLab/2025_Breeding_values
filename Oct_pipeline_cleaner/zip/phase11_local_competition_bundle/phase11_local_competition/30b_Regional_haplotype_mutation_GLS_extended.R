## Auto-generated extended version for phase 11 local-QTL competition analyses.
## 27_Regional_haplotype_mutation_GLS_all8.R
## For EACH historical region: estimate founder-similarity effects first, then
## the independent contribution of nominated variants/R5 tags conditional on
## the same founder terms; screen estimable pairwise genetic interactions.
## Approximate fixed-variance-ratio (EMMAX-style) kinship-adjusted BLUE models.
## Requires output of script 26. These are associations, not causal editing tests.

suppressPackageStartupMessages(library(AGHmatrix))
preparation_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
origin_file <- '19_RIL_source_populations.csv' # optional, documented origins
model_group_file <- '30_Regions_and_R5_groups_extended.csv'
focal_file <- '30_Focal_variants_extended.csv'
script26_dir <- 'historical_QTL_R5/deep_dive_extended'
output_dir <- file.path(script26_dir,'regional_models')
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
minimum_n <- 70L
minimum_cell_n_for_interaction <- 5L
min_ancestry_called_fraction <- 0.55
min_ancestry_class_n <- 15L
maximum_ancestry_terms <- 2L
minimum_nonIV_pruned_SNPs <- 100L
set.seed(20261029)
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

required <- c(preparation_file,model_group_file,focal_file,
 file.path(script26_dir,'26_R5_model_tag_SNPs.csv'),
 file.path(script26_dir,'26_RIL_founder_similarity_PROXIES.csv'))
if(any(!file.exists(required)))stop('Files missing; run 26 first: ',paste(required[!file.exists(required)],collapse=', '))
e<-new.env(parent=baseenv());load(preparation_file,envir=e)
phname<-intersect(c('phen_NaCl','all_BLUEs_genotyped_NaCl','data_for_GWAS_NaCl_genotyped_NaCl'),ls(e))
if(!length(phname))stop('Missing phenotype table in script 00')
ph<-e[[phname[1]]]
if(!all(c('pop_label',traits[1:7])%in%names(ph)))stop('Missing phenotypes')
if(anyDuplicated(ph$pop_label))stop('Use unique adjusted BLUEs, not repeated raw phenotype observations')
if(!'beta'%in%names(ph)){
 if(!all(c('beta_for_mapping','replicate_mean','replicate_sd')%in%ls(e)))stop('Missing saved coefficients to reconstruct beta')
 zz<-sweep(as.matrix(ph[,traits[1:7],drop=FALSE]),2,e$replicate_mean[traits[1:7]],'-')
 zz<-sweep(zz,2,e$replicate_sd[traits[1:7]],'/')
 ph$beta<-drop(zz%*%e$beta_for_mapping[traits[1:7]])
}
if(!all(c('genotype_matrix','genotype_pruned_matrix')%in%ls(e)))stop('Missing RIL genotype matrices')
G<-e$genotype_matrix;P<-e$genotype_pruned_matrix
ids<-intersect(as.character(ph$pop_label),intersect(rownames(G),rownames(P)))
if(length(ids)<minimum_n)stop('Too few genotyped and phenotyped RILs')
ph<-ph[match(ids,ph$pop_label),,drop=FALSE];P<-P[ids,,drop=FALSE]
origin<-setNames(rep('unknown',length(ids)),ids)
if(file.exists(origin_file)){
 oo<-read.csv(origin_file,stringsAsFactors=FALSE)
 if(!all(c('pop_label','source_population')%in%names(oo))||anyDuplicated(oo$pop_label))stop('Invalid origin file')
 origin[]<-as.character(oo$source_population[match(ids,oo$pop_label)])
 origin[is.na(origin)|!nzchar(origin)]<-'unknown'
}
groups<-read.csv(model_group_file,stringsAsFactors=FALSE)
groups<-groups[toupper(as.character(groups$include))=='TRUE',,drop=FALSE]
groups$group<-as.character(groups$group)
focal<-read.csv(focal_file,stringsAsFactors=FALSE)
tags<-read.csv(file.path(script26_dir,'26_R5_model_tag_SNPs.csv'),stringsAsFactors=FALSE)
tags$group<-as.character(tags$group)
proxies<-read.csv(file.path(script26_dir,'26_RIL_founder_similarity_PROXIES.csv'),stringsAsFactors=FALSE)
proxies$group<-as.character(proxies$group)
if(anyDuplicated(paste(proxies$pop_label,proxies$group)))stop('Duplicate RIL-group ancestry proxy')

## Cache chromosome-specific kinship; no focal chromosome in the GRM.
Kcache<-new.env(parent=emptyenv())
getK<-function(chr){
 if(exists(chr,envir=Kcache,inherits=FALSE))return(get(chr,envir=Kcache))
 keep<-!grepl(paste0('^',chr,'([_:.]|$)'),colnames(P))
 if(sum(keep)<minimum_nonIV_pruned_SNPs)stop('Not enough LOCO kinship SNPs for chromosome ',chr)
 kk<-as.matrix(AGHmatrix::Gmatrix(P[,keep,drop=FALSE],method='VanRaden'))
 if(is.null(rownames(kk))||!all(ids%in%rownames(kk)))stop('Kinship RIL labels not aligned')
 kk<-kk[ids,ids,drop=FALSE];kk<-(kk+t(kk))/2
 assign(chr,kk,envir=Kcache);kk
}
fit_w<-function(y,X){
 if(nrow(X)<ncol(X)+5L||qr(X,tol=1e-8)$rank<ncol(X))return(NULL)
 fit<-qr(X,tol=1e-8);b<-as.vector(qr.coef(fit,y))
 rs<-y-drop(X%*%b);rss<-sum(rs^2);df<-length(y)-ncol(X)
 if(!is.finite(rss)||rss<=0||df<5L)return(NULL)
 inv<-tryCatch(solve(crossprod(X)),error=function(err)NULL)
 if(is.null(inv))return(NULL)
 se<-sqrt(pmax(0,diag(inv)*rss/df))
 list(beta=b,se=se,rss=rss,df=df,p=ncol(X),kappa=kappa(X))
}
make_whitener<-function(y,K,X0){
 ee<-eigen(K,symmetric=TRUE);U<-ee$vectors;lambda<-pmax(ee$values,0)
 ey<-drop(crossprod(U,y));xx<-crossprod(U,X0)
 n<-length(y);p<-ncol(X0)
 obj<-function(ld){
   v<-lambda+exp(ld);xxw<-crossprod(xx,xx/v)
   if(qr(xxw)$rank<ncol(xxw))return(Inf)
   bb<-tryCatch(solve(xxw,crossprod(xx,ey/v)),error=function(err)NULL)
   if(is.null(bb))return(Inf)
   rss<-sum((ey-drop(xx%*%bb))^2/v)
   detx<-determinant(xxw,logarithm=TRUE)
   if(!is.finite(rss)||rss<=0||detx$sign<=0)return(Inf)
   (n-p)*log(rss/(n-p))+sum(log(v))+as.numeric(detx$modulus)
 }
 delta<-exp(optimize(obj,c(-9,9))$minimum)
 scale<-1/sqrt(lambda+delta)
 list(y=ey*scale,transform=function(X)sweep(crossprod(U,X),1,scale,'*'),delta=delta)
}
ptest<-function(f0,f1){
 if(is.null(f0)||is.null(f1)||f1$p<=f0$p)return(NA_real_)
 ff<-max(0,(f0$rss-f1$rss)/(f1$p-f0$p))/(f1$rss/f1$df)
 pf(ff,f1$p-f0$p,f1$df,lower.tail=FALSE)
}
getterm<-function(f,X,nm){
 if(is.null(f)||!nm%in%colnames(X))return(c(estimate=NA_real_,SE=NA_real_,p=NA_real_))
 j<-match(nm,colnames(X));bb<-f$beta[j];s<-f$se[j]
 c(estimate=bb,SE=s,p=if(is.finite(s)&&s>0)2*pt(abs(bb/s),f$df,lower.tail=FALSE)else NA_real_)
}

fits<-list(); diagnostics<-list(); ancestryQC<-list(); interaction<-list();decomp<-list()
for(region in unique(groups$region)){
 tab<-groups[groups$region==region,,drop=FALSE]
 chr<-unique(tab$CHROM);if(length(chr)!=1L)stop('One region spans chromosomes')
 kk<-getK(chr)
 tg<-tags[tags$region==region,,drop=FALSE]
 ff<-focal[focal$region==region,,drop=FALSE]
 ## Named variants always have priority; cluster representative tags follow.
 mk<-data.frame(label=character(),SNP_ID=character(),group=character(),source=character())
 if(nrow(ff))for(j in seq_len(nrow(ff))){
   id<-paste0(ff$CHROM[j],'_',ff$POS[j]);gg<-tg$group[match(id,tg$tag_SNP_ID)]
   mk<-rbind(mk,data.frame(label=paste0('F',ff$POS[j]),SNP_ID=id,
                       group=if(length(gg)&&!is.na(gg))gg else NA_character_,source='focal'))
 }
 for(j in seq_len(nrow(tg))){
   id<-as.character(tg$tag_SNP_ID[j])
   if(!id%in%mk$SNP_ID)
      mk<-rbind(mk,data.frame(label=paste0('R5_',tg$group[j]),SNP_ID=id,
                         group=tg$group[j],source='R5_tag'))
 }
 if(any(!mk$SNP_ID%in%colnames(G)))stop('RIL matrix lacks regional markers: ',
                         paste(mk$SNP_ID[!mk$SNP_ID%in%colnames(G)],collapse=','))
 M<-as.matrix(G[ids,mk$SNP_ID,drop=FALSE])/2
 colnames(M)<-mk$label
 if(any(is.finite(M)&(M<0|M>1)))stop('Mutation/tag dosage outside [0,2] at ',region)
 ## Build ancestry proxy candidates using marker windows distinct from focal SNPs.
 ## Only confidently distinguishable founders are used (NAs stay missing).
 aa<-list()
 for(j in seq_len(nrow(tg))){
   g<-tg$group[j];z<-proxies[proxies$group==g,,drop=FALSE]
   if(!nrow(z))next
   pidx<-match(ids,z$pop_label)
   who<-as.character(z$best_founder_proxy[pidx])
   if(region=='beta_I'){
     target<-if(g=='53')c('PX174','RC301')else c('RC301')
   }else if(region=='beta_V')target<-c('CB4507')
   else if(region=='chrIV_CB4856')target<-c('CB4856')
   else {
     tt<-sort(table(who),decreasing=TRUE)
     target<-if(length(tt))names(tt)[1]else character()
   }
   if(!length(target))next
   H<-ifelse(is.na(who),NA_real_,as.numeric(who%in%target))
   na<-sum(is.finite(H));n1<-sum(H==1,na.rm=TRUE);n0<-sum(H==0,na.rm=TRUE)
   qc<-data.frame(region=region,group=g,target_founders=paste(target,collapse=';'),
      n_RILs_called=na,n_target=n1,n_other=n0,call_fraction=na/length(ids),
      enough_information=na/length(ids)>=min_ancestry_called_fraction&&
         n1>=min_ancestry_class_n&&n0>=min_ancestry_class_n)
   ancestryQC[[length(ancestryQC)+1L]]<-qc
   if(qc$enough_information)aa[[paste0('H_',g)]]<-H
 }
 ## Prefer best-covered founder-proxy terms, never force ambiguous ancestry.
 if(length(aa)>maximum_ancestry_terms){
    cov<-vapply(aa,function(h)sum(is.finite(h)),integer(1))
    aa<-aa[order(cov,decreasing=TRUE)[seq_len(maximum_ancestry_terms)]]
 }
 ## To compare H and H+M, use the SAME RILs. Retain no H if a strict
 ## complete-case intersection becomes too small.
 if(length(aa)){
   intersection<-Reduce('&',lapply(aa,is.finite))
   if(sum(intersection)<minimum_n||sum(intersection)<.50*length(ids)){
      aa<-aa[which.max(vapply(aa,function(h)sum(is.finite(h)),integer(1)))]
   }
 }
 H<-if(length(aa))do.call(cbind,aa)else matrix(numeric(0),nrow=length(ids),ncol=0)
 if(ncol(H)==1L)H<-matrix(H,ncol=1L,dimnames=list(NULL,names(aa)))
 if(ncol(H)>1L)colnames(H)<-names(aa)
 ## A staged test that uses all eight outcomes, with the same covariates.
 for(trait in traits){
  y0<-as.numeric(as.character(ph[[trait]]))
  ok<-is.finite(y0)&rowSums(!is.finite(M))==0L
  if(ncol(H))ok<-ok&rowSums(!is.finite(H))==0L
  if(sum(ok)<minimum_n){warning('Insufficient complete data for ',region,':',trait);next}
  y<-y0[ok];Mm<-M[ok,,drop=FALSE];Hh<-H[ok,,drop=FALSE]
  oo<-factor(origin[ids[ok]]);base<-if(nlevels(oo)>1L)
    model.matrix(~ oo) else matrix(1,sum(ok),1L,dimnames=list(NULL,'(Intercept)'))
  if(qr(base)$rank<ncol(base))stop('RIL-origin model rank deficient')
  X0<-base
  chosenH<-character()
  if(ncol(Hh))for(k in colnames(Hh)){
    trial<-cbind(X0,Hh[,k]);colnames(trial)[ncol(trial)]<-k
    if(qr(trial,tol=1e-8)$rank==ncol(trial)){
       X0<-trial;chosenH<-c(chosenH,k)
    }
  }
  ## Prefer the named variants over redundant R5 proxies; save omissions.
  X1<-X0;chosenM<-character();droppedM<-character()
  for(k in colnames(Mm)){
    trial<-cbind(X1,Mm[,k]);colnames(trial)[ncol(trial)]<-k
    if(qr(trial,tol=1e-8)$rank==ncol(trial)){
        X1<-trial;chosenM<-c(chosenM,k)
    }else droppedM<-c(droppedM,k)
  }
  w<-make_whitener(y,kk[ok,ok,drop=FALSE],base)
  f0<-fit_w(w$y,w$transform(base))
  fH<-fit_w(w$y,w$transform(X0))
  fA<-fit_w(w$y,w$transform(X1))
  d<-data.frame(region=region,CHROM=chr,trait=trait,n_RIL=length(y),
     n_source_populations=nlevels(oo),n_H_terms=length(chosenH),
     n_selected_SNP_terms=length(chosenM),
     H_terms=paste(chosenH,collapse=';'),variant_terms=paste(chosenM,collapse=';'),
     dropped_collinear_SNP_terms=paste(droppedM,collapse=';'),
     p_H_given_origin=ptest(f0,fH),
     p_markers_given_H=ptest(fH,fA),
     variance_ratio_delta=w$delta,
     joint_condition=if(!is.null(fA))fA$kappa else NA_real_)
  diagnostics[[length(diagnostics)+1L]]<-d
  ## Estimate the background association among NONCARRIERS of every
  ## explicitly nominated focal mutation in that region. This is distinct
  ## from stage A, which includes focal-allele carriers and can absorb them.
  focal_terms<-intersect(paste0('F',ff$POS),colnames(Mm))
  if(length(chosenH)&&length(focal_terms)){
    no_foc<-rowSums(Mm[,focal_terms,drop=FALSE]>0,na.rm=FALSE)==0
    if(sum(no_foc)>=minimum_n &&
       qr(X0[no_foc,,drop=FALSE],tol=1e-8)$rank==ncol(X0)){
      wN<-make_whitener(y[no_foc],kk[which(ok)[no_foc],which(ok)[no_foc],drop=FALSE],
                        base[no_foc,,drop=FALSE])
      fN<-fit_w(wN$y,wN$transform(X0[no_foc,,drop=FALSE]))
      for(k in chosenH){
        eeN<-getterm(fN,X0,k)
        fits[[length(fits)+1L]]<-data.frame(region=region,trait=trait,
          stage='A0_founder_proxy_among_focal_noncarriers_only',term=k,
          effect=eeN['estimate'],SE=eeN['SE'],p=eeN['p'],
          n_RIL=sum(no_foc),trait_sd=sd(y[no_foc]),
          minimum_detectable_80pct_SD=(qnorm(.975)+qnorm(.8))*eeN['SE']/sd(y[no_foc]),
          stringsAsFactors=FALSE)
      }
    }
  }
  for(k in chosenH){
    est<-getterm(fH,X0,k)
    fits[[length(fits)+1L]]<-data.frame(region=region,trait=trait,
      stage='A_founder_proxy_before_variants',term=k,
      effect=est['estimate'],SE=est['SE'],p=est['p'],n_RIL=length(y),
      trait_sd=sd(y),minimum_detectable_80pct_SD=
       (qnorm(.975)+qnorm(.8))*est['SE']/sd(y),stringsAsFactors=FALSE)
  }
  for(k in chosenM){
    est<-getterm(fA,X1,k)
    fits[[length(fits)+1L]]<-data.frame(region=region,trait=trait,
       stage='B_variant_conditional_on_founder_proxies_and_other_variants',term=k,
       effect=est['estimate'],SE=est['SE'],p=est['p'],n_RIL=length(y),
       trait_sd=sd(y),minimum_detectable_80pct_SD=
       (qnorm(.975)+qnorm(.8))*est['SE']/sd(y),stringsAsFactors=FALSE)
  }
  ## Pairwise epistasis is a SCREEN among distinguishable loci. Require
  ## genotype classes (using homozygous inbred 0/1 coding) to be populated.
  if(length(chosenM)>=2L)for(pair in combn(chosenM,2,simplify=FALSE)){
     m1<-Mm[,pair[1]];m2<-Mm[,pair[2]]
     tab<-table(factor(round(m1),levels=0:1),factor(round(m2),levels=0:1))
     mincell<-min(tab)
     ix<-paste(pair,collapse='x')
     trial<-cbind(X1,m1*m2);colnames(trial)[ncol(trial)]<-ix
     identifiable<-mincell>=minimum_cell_n_for_interaction&&
        qr(trial,tol=1e-8)$rank==ncol(trial)
     fI<-if(identifiable)fit_w(w$y,w$transform(trial))else NULL
     est<-getterm(fI,trial,ix)
     interaction[[length(interaction)+1L]]<-data.frame(region=region,trait=trait,
       locus1=pair[1],locus2=pair[2],n_RIL=length(y),
       n_00=tab[1,1],n_01=tab[1,2],n_10=tab[2,1],n_11=tab[2,2],
       identifiable=identifiable,interaction_beta=est['estimate'],
       interaction_SE=est['SE'],p=est['p'],
       p_added_interaction=ptest(fA,fI),
       minimum_detectable_80pct_SD=(qnorm(.975)+qnorm(.8))*est['SE']/sd(y))
  }
  ## Genetic-marker x founder-proxy interaction is meaningful only if
  ## the focal allele occurs in BOTH ancestry classes with enough RILs.
  ## Otherwise y~M+H+M:H is rank-deficient (as at IV:10,499,135).
  if(length(chosenH)&&length(chosenM))for(hterm in chosenH){
    hvec<-Hh[,hterm]
    for(k in chosenM){
      mv<-Mm[,k];ct<-table(factor(round(mv),levels=0:1),
                           factor(round(hvec),levels=0:1))
      nm<-paste0(k,'x',hterm)
      xx<-cbind(X1,mv*hvec);colnames(xx)[ncol(xx)]<-nm
      good<-min(ct)>=minimum_cell_n_for_interaction &&
            qr(xx,tol=1e-8)$rank==ncol(xx)
      fX<-if(good)fit_w(w$y,w$transform(xx))else NULL
      est<-getterm(fX,xx,nm)
      interaction[[length(interaction)+1L]]<-data.frame(region=region,
        trait=trait,locus1=k,locus2=hterm,n_RIL=length(y),
        n_00=ct[1,1],n_01=ct[1,2],n_10=ct[2,1],n_11=ct[2,2],
        identifiable=good,interaction_beta=est['estimate'],
        interaction_SE=est['SE'],p=est['p'],p_added_interaction=ptest(fA,fX),
        minimum_detectable_80pct_SD=(qnorm(.975)+qnorm(.8))*est['SE']/sd(y))
    }
  }
  ## Exact beta decomposition on the SAME RILs, same model, same whitening.
  if(trait=='beta'&&!is.null(fA)&&all(c('beta_for_mapping','replicate_sd')%in%ls(e))){
    seven<-traits[1:7]
    wt<-as.numeric(e$beta_for_mapping[seven]/e$replicate_sd[seven]);names(wt)<-seven
    Y<-as.matrix(ph[ok,seven,drop=FALSE]);storage.mode(Y)<-'double'
    if(all(is.finite(Y))&&all(is.finite(wt))){
      recon<-drop(Y%*%wt);rr<-recon-y
      if(max(abs(rr-mean(rr)))<1e-5){
       Xt<-w$transform(X1)
       for(t in seven){
         ft<-fit_w(w$transform(matrix(Y[,t],ncol=1))[,1],Xt)
         for(k in chosenM){
           jj<-match(k,colnames(X1))
           decomp[[length(decomp)+1L]]<-data.frame(region=region,term=k,
             trait=t,weight=wt[t],trait_effect=ft$beta[jj],
             weighted_beta_contribution=wt[t]*ft$beta[jj],
             beta_effect_same_GLS=fA$beta[jj])
         }
       }
      } else warning('beta reconstruction failed for ',region,': skipping decomposition')
    }
  }
  cat('Completed ',region,':',trait,' (n=',length(y),')\n',sep='')
 }
}
if(!length(fits))stop('No estimable regional effects: inspect founder-proxy QC')
out<-do.call(rbind,fits)
## BH across traits for each region+stage+term; also global across all tests.
out$q_BH_eight_traits<-NA_real_
for(k in unique(paste(out$region,out$stage,out$term,sep='::'))){
 ix<-which(paste(out$region,out$stage,out$term,sep='::')==k)
 out$q_BH_eight_traits[ix]<-p.adjust(out$p[ix],method='BH')
}
out$q_BH_all_regional_terms<-p.adjust(out$p,method='BH')
write.csv(out,file.path(output_dir,'27_Effects_founder_then_variants_all8.csv'),row.names=FALSE)
diag<-do.call(rbind,diagnostics)
diag$q_H_eight_traits<-diag$q_markers_eight_traits<-NA_real_
for(r in unique(diag$region)){
 ix<-which(diag$region==r)
 diag$q_H_eight_traits[ix]<-p.adjust(diag$p_H_given_origin[ix],method='BH')
 diag$q_markers_eight_traits[ix]<-p.adjust(diag$p_markers_given_H[ix],method='BH')
}
write.csv(diag,file.path(output_dir,'27_Nested_model_tests_and_identifiability.csv'),row.names=FALSE)
if(length(ancestryQC))write.csv(do.call(rbind,ancestryQC),
 file.path(output_dir,'27_Ancestry_information_by_region.csv'),row.names=FALSE)
if(length(interaction)){
 epi<-do.call(rbind,interaction)
 epi$q_BH_all_testable_interactions<-p.adjust(epi$p,method='BH')
 write.csv(epi,file.path(output_dir,'27_Pairwise_epistasis_all8.csv'),row.names=FALSE)
}
if(length(decomp)){
 d<-do.call(rbind,decomp)
 write.csv(d,file.path(output_dir,'27_Beta_weighted_trait_contributions.csv'),row.names=FALSE)
 audit<-aggregate(weighted_beta_contribution~region+term,d,sum)
 tot<-unique(d[,c('region','term','beta_effect_same_GLS')])
 audit<-merge(audit,tot,by=c('region','term'))
 audit$abs_reconstruction_error<-abs(audit$weighted_beta_contribution-audit$beta_effect_same_GLS)
 write.csv(audit,file.path(output_dir,'27_Beta_decomposition_reconstruction_QC.csv'),row.names=FALSE)
}
cat('Regional conditional association analysis complete. Founder assignments are PROXIES.\n')

## Manuscript-oriented exploratory forest plots. All effects are scaled by
## the observed trait SD in the *matched* fitted subset (not effect on fitness).
pdf(file.path(output_dir,'27_Regional_forest_plots.pdf'),width=10,height=7)
for(region in unique(out$region)){
 for(stage in c('A_founder_proxy_before_variants',
                'A0_founder_proxy_among_focal_noncarriers_only',
                'B_variant_conditional_on_founder_proxies_and_other_variants')){
  zz<-out[out$region==region&out$stage==stage&
          is.finite(out$effect)&is.finite(out$SE)&is.finite(out$trait_sd)&out$trait_sd>0,,drop=FALSE]
  if(!nrow(zz))next
  lab<-paste(zz$term,zz$trait,sep=' / ')
  yy<-seq_len(nrow(zz));eff<-zz$effect/zz$trait_sd;se<-zz$SE/zz$trait_sd
  xr<-range(c(0,eff-1.96*se,eff+1.96*se),finite=TRUE)
  par(mar=c(4,10,4,1))
  plot(eff,yy,ylim=c(.5,length(yy)+.5),yaxt='n',xlim=xr,pch=19,
       xlab='Associated effect in phenotype SD units',ylab='',
       main=paste(region,'\n',stage),cex.main=.85)
  arrows(eff-1.96*se,yy,eff+1.96*se,yy,angle=90,code=3,length=.025)
  axis(2,at=yy,labels=lab,las=2,cex.axis=.65)
  abline(v=0,lty=2,col='grey60')
 }
}
dev.off()
