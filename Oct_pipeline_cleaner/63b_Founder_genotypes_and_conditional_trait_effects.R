# 63b_Founder_genotypes_and_conditional_trait_effects.R
# Companion figures to 63. Uses frozen historical script-26 founder genotype
# table and script-27 kinship-adjusted, conditional regional QTL estimates.
# These are genotype backgrounds / ASSOCIATIONS, not proven causal variants.
# Run: source('Oct_pipeline_cleaner/63b_Founder_genotypes_and_conditional_trait_effects.R')
# Base R only.

CFG <- list(
 founder='historical_QTL_R5/deep_dive_extended/26_Founder_allele_signature_all_sites.csv',
 effects='historical_QTL_R5/deep_dive_extended/regional_models/27_Effects_founder_then_variants_all8.csv',
 founder_ids='18_founder_sample_ids.txt',
 out='selection_results/manuscript_finalization/63_three_regions',
 max_founder_markers=35L
)
if(!file.exists(CFG$founder)||!file.exists(CFG$effects)||!file.exists(CFG$founder_ids))
 stop('Requires original scripts 26/27 results plus 18_founder_sample_ids.txt. Missing: ',
      paste(c(CFG$founder,CFG$effects,CFG$founder_ids)[!file.exists(c(CFG$founder,CFG$effects,CFG$founder_ids))],collapse=', '))
S<-read.csv(CFG$founder,stringsAsFactors=FALSE,check.names=FALSE)
E<-read.csv(CFG$effects,stringsAsFactors=FALSE)
founders<-trimws(readLines(CFG$founder_ids,warn=FALSE));founders<-founders[nzchar(founders)]
if(length(founders)!=16L||anyDuplicated(founders))stop('Expected 16 unique VCF founder IDs')
if(!all(founders%in%names(S)))stop('Founder column(s) absent from script 26 output')
if(!all(c('region','CHROM','POS','SNP_ID')%in%names(S)))stop('Unexpected script 26 table')
if(!all(c('region','trait','stage','term','effect','SE','trait_sd')%in%names(E)))stop('Unexpected script 27 table')
dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)
meta<-list(I=list(region='beta_I',focal=c(7167738,7491656,7596359)),
           IV=list(region='chrIV_CB4856',focal=c(8892663,10499135)),
           V=list(region='beta_V',focal=c(6322080,6637695,6826421)))
for(chr in names(meta)){
 r<-meta[[chr]]$region
 s<-S[S$region==r & S$CHROM==chr,,drop=FALSE]
 e<-E[E$region==r & E$stage=='B_variant_conditional_on_founder_proxies_and_other_variants' &
      E$trait%in%c('SF','SB','FS','FB','BS','BF','area.F','beta'),,drop=FALSE]
 if(!nrow(s)||!nrow(e)){
    warning(chr,' lacks founder or conditional-effect data; skipping supplemental plot')
    next
 }
 s<-s[order(s$POS),,drop=FALSE]
 # Prespecify focal exact positions first, then evenly spaced background SNPs
 idx<-which(s$POS%in%meta[[chr]]$focal)
 extra<-unique(round(seq_len(nrow(s))[unique(round(seq(1,nrow(s),length.out=min(nrow(s),CFG$max_founder_markers))))]))
 take<-sort(unique(c(idx,extra)))
 if(length(take)>CFG$max_founder_markers+length(idx))
   take<-sort(unique(c(idx,head(extra,CFG$max_founder_markers))))
 s<-s[take,,drop=FALSE]
 mat<-t(as.matrix(s[,founders,drop=FALSE]));storage.mode(mat)<-'double'
 mat[!is.finite(mat)]<-NA_real_;mat[is.finite(mat)&mat==2]<-1;mat[is.finite(mat)&mat==0]<-0
 mat[is.finite(mat)&!mat%in%c(0,1)]<-NA_real_ # no unphased heterozygotes as binary founder signatures
 # Use the same full rows in effect table and retain named tag terms.
 e$effect_SD <- e$effect / e$trait_sd
 e$SE_SD <- e$SE / e$trait_sd
 e <- e[is.finite(e$effect_SD)&is.finite(e$SE_SD),,drop=FALSE]
 # Restrict forest to explicitly nominated mutations (Fposition) where possible.
 fterms<-paste0('F',meta[[chr]]$focal)
 ee<-e[e$term%in%fterms,,drop=FALSE]
 if(!nrow(ee))ee<-e[order(e$p),,drop=FALSE][seq_len(min(20L,nrow(e))),,drop=FALSE]
 make<-function(){
    op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
    layout(matrix(1:2,nrow=1),widths=c(1.15,1.0))
    par(mar=c(7.7,7,3.7,1))
    image(seq_len(ncol(mat)),seq_len(nrow(mat)),t(mat),axes=FALSE,
          col=c('#F4F1E8','#B03D45'),zlim=c(0,1),
          xlab='',ylab='',main=paste('Chr',chr,'| 16 founder genotype signatures'))
    axis(1,at=seq_len(ncol(mat)),labels=s$POS,las=2,cex.axis=.65)
    axis(2,at=seq_len(nrow(mat)),labels=founders,las=1,cex.axis=.76)
    mtext('1-based WBcel235 position',side=1,line=5.8,cex=.7)
    box()
    par(mar=c(7.7,6,3.7,1.2))
    if(!nrow(ee)){
      plot.new();text(.5,.5,'No estimable conditional effects')
    }else{
      labs<-paste0(ee$term,': ',ee$trait)
      yy<-rev(seq_len(nrow(ee)))
      lo<-ee$effect_SD-1.96*ee$SE_SD;hi<-ee$effect_SD+1.96*ee$SE_SD
      xr<-range(c(lo,hi,0),finite=TRUE);pad<-max(.02,diff(xr)*.12)
      plot(ee$effect_SD,yy,xlim=xr+c(-pad,pad),ylim=c(.5,max(yy)+.7),
           pch=16,col='#24689B',yaxt='n',xlab='Conditional trait effect / trait SD',ylab='',
           main='Kinship-adjusted RIL associations')
      arrows(lo,yy,hi,yy,angle=90,code=3,length=.04,col='#24689B')
      abline(v=0,lty=2,col='grey70')
      axis(2,at=yy,labels=labs,las=1,cex.axis=.58)
    }
    mtext('Founder presence is not causal assignment; conditional tests are exploratory.',outer=TRUE,
          side=1,line=-1,cex=.75,col='grey40')
 }
 pdf(file.path(CFG$out,paste0('63b_Chr',chr,'_founder_and_trait_effects.pdf')),
     width=14,height=8,useDingbats=FALSE);make();dev.off()
 png(file.path(CFG$out,paste0('63b_Chr',chr,'_founder_and_trait_effects.png')),
     width=14,height=8,units='in',res=300);make();dev.off()
 write.csv(ee,file.path(CFG$out,paste0('63b_Chr',chr,'_conditional_effects_plotted.csv')),row.names=FALSE)
 cat('Wrote founder/trait panels for ',chr,'\n')
}
writeLines(c('Supplemental founder backgrounds: allele 0/1 after diploid homozygote coding.',
             'Do not equate VCF allele presence with ancestry if allele orientation differs.',
             '27 conditionals are kinship-adjusted but depend on chosen founder proxy and variant terms.',
             'Unidentifiable interactions should not be inferred.'),
           file.path(CFG$out,'63b_Interpretation.txt'))
