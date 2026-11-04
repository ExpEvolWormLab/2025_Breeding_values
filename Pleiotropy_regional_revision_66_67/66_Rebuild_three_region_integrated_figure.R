# 66_Rebuild_three_region_integrated_figure.R
# Publication-oriented regional figure, descriptive only.
# Run from the CeMEE project root AFTER scripts 63 and 63b.
# Base R, no VEP. Does not make unvalidated selection or causality claims.
# 2026-11-03

CFG <- list(
  source_dir = 'selection_results/manuscript_finalization/63_three_regions',
  out = 'selection_results/manuscript_finalization/66_regional_synthesis',
  historical_reads = 'historical_QTL_R5/deep_dive/regional_models/28_ALT_read_counts_all_pools.csv',
  founder_signatures = 'historical_QTL_R5/deep_dive/regional_models/26_Founder_allele_signature_all_sites.csv',
  nested_tests = 'historical_QTL_R5/deep_dive/regional_models/27_Nested_model_tests_and_identifiability.csv',
  dpi = 300L
)

META <- list(
 I=list(title='Chr I | RC301 / PX174 backgrounds',focus=c(7167738L,7491656L,7596359L),
        term='F7167738',highlight=c(53L,54L)),
 IV=list(title='Chr IV | CB4856 and tax-6 region',focus=c(8892663L,10499135L),
        term='F10499135',highlight=c(727L,731L)),
 V=list(title='Chr V | linked multivariate response',focus=c(6322080L,6637695L,6826421L),
        term='F6637695',highlight=c(848L,851L,853L,854L))
)
TRAITS <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
PALETTE <- c(GA1='#287AB5',GA2='#DD8C23',GA4='#20916F')
FIRST <- function(p) { a<-p[file.exists(p)]; if(length(a)) a[1] else NA_character_ }
if(!dir.exists(CFG$source_dir))stop('Missing script 63/63b folder: ',CFG$source_dir)
dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)

read_input <- function(chr,kind){
 p <- file.path(CFG$source_dir,sprintf('63_%s_%s.csv',chr,kind))
 if(!file.exists(p))stop('Missing script 63 input: ',p)
 read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
}
# 63b names are different from 63 names.
read_effects <- function(chr){
 p<-file.path(CFG$source_dir,sprintf('63b_Chr%s_conditional_effects_plotted.csv',chr))
 if(!file.exists(p))stop('Missing script 63b input: ',p)
 d<-read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
 required<-c('trait','term','effect_SD','SE_SD','q_BH_all_regional_terms')
 if(!all(required%in%names(d)))stop('Unexpected 63b schema in ',p)
 d
}

# For exact frequencies, use the already validated script 63 export.
# Never treat A6G140 as a second G0 ancestor, or GA3 as a survivor.
get_trajectory <- function(chr,focus){
 p<-file.path(CFG$source_dir,paste0('63_',chr,'_focal_rare_allele_trajectories.csv'))
 if(file.exists(p)){
   a<-read.csv(p,stringsAsFactors=FALSE)
   req<-c('SNP_ID','population','generation','rare_frequency')
   if(!all(req%in%names(a)))stop('Bad 63 trajectory schema: ',p)
 } else {
   # Optional fallback for exact read counts. This does not invent unmeasured sites.
   if(!file.exists(CFG$historical_reads))return(data.frame())
   a<-read.csv(CFG$historical_reads,stringsAsFactors=FALSE)
   req<-c('SNP_ID','sample_id','generations','ALT_frequency')
   if(!all(req%in%names(a)))stop('Bad historical read-count schema')
   a<-a[a$SNP_ID%in%paste0(chr,'_',focus),,drop=FALSE]
   if(!nrow(a))return(data.frame())
   a$generation<-a$generations
   a$population<-sub('^GA([124]).*','GA\\1',a$sample_id)
   a$population[a$sample_id=='A00']<-'A00'
   a$rare_frequency<-NA_real_
   for(id in unique(a$SNP_ID)){
     ii<-which(a$SNP_ID==id)
     anc<-a$ALT_frequency[ii[a$population[ii]=='A00' & a$generation[ii]==0]]
     if(length(anc)!=1L||!is.finite(anc))stop('Missing/ambiguous A00 frequency for ',id)
     a$rare_frequency[ii]<-if(anc>0.5)1-a$ALT_frequency[ii]else a$ALT_frequency[ii]
   }
 }
 a<-a[a$SNP_ID%in%paste0(chr,'_',focus) & a$population%in%c('A00','GA1','GA2','GA4'),,drop=FALSE]
 if(anyDuplicated(paste(a$SNP_ID,a$population,a$generation)))
   stop('Duplicate site/population/generation in ',chr,' trajectories')
 a
}

# Validate founder configurations and summarize identification limits.
if(file.exists(CFG$founder_signatures)){
 S<-read.csv(CFG$founder_signatures,stringsAsFactors=FALSE)
 rows<-list()
 for(chr in names(META))for(pos in META[[chr]]$focus){
   d<-S[S$CHROM==chr & S$POS==pos,,drop=FALSE]
   if(!nrow(d))next
   signatures<-unique(paste(d$REF,d$ALT,ifelse(is.na(d$ALT_founder_names),'none called',d$ALT_founder_names)))
   if(length(signatures)!=1L)stop('Inconsistent founder signatures for ',chr,':',pos)
   u<-d[1,,drop=FALSE]
   rows[[length(rows)+1L]]<-data.frame(
     SNP_ID=paste0(chr,'_',pos),REF=u$REF,ALT=u$ALT,
     n_founders_ALT=u$n_founders_ALT,
     founders_ALT=ifelse(is.na(u$ALT_founder_names),'none called',u$ALT_founder_names),
     stringsAsFactors=FALSE)
 }
 if(length(rows))write.csv(do.call(rbind,rows),file.path(CFG$out,'66_focal_founder_signatures.csv'),row.names=FALSE)
}
if(file.exists(CFG$nested_tests)){
 N<-read.csv(CFG$nested_tests,stringsAsFactors=FALSE)
 r<-N[N$CHROM%in%names(META),,drop=FALSE]
 write.csv(r,file.path(CFG$out,'66_original_nested_model_tests.csv'),row.names=FALSE)
}

QC <- list()
draw <- function(){
 op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
 layout(matrix(1:9,nrow=3,byrow=TRUE),widths=c(1.25,1.13,1.12),heights=rep(1,3))
 par(oma=c(4,1,3,1),mgp=c(2.2,.65,0),tcl=-.25)
 for(chr in names(META)){
   reg<-META[[chr]]
   z<-read_input(chr,'R5_prediction_window')
   if(!all(c('ps','group','primary_slope_median_phase','LOCO_full')%in%names(z)))
     stop('Required columns missing from R5 window ',chr)
   z<-z[is.finite(z$ps)&is.finite(z$primary_slope_median_phase)&is.finite(z$LOCO_full),,drop=FALSE]
   e<-read_effects(chr)
   a<-get_trajectory(chr,reg$focus)
   QC[[chr]] <<- data.frame(chr=chr,n_R5_finite=nrow(z),n_exact_sites=length(unique(a$SNP_ID)),
                            n_trajectories=nrow(a),has_exact_counts=nrow(a)>0)
   # A: discrete genomic clusters, NOT a smoothed locus-wise trend.
   par(mar=c(4,4.6,3.5,1))
   xx<-z$ps/1e6
   yy<-1000*z$primary_slope_median_phase
   pr<-1000*z$LOCO_full
   lim<-range(c(yy,pr,0),finite=TRUE)
   lim<-lim+c(-1,1)*max(.2,.12*diff(lim))
   plot(xx,yy,pch=16,col='#287AB5',cex=1.1,ylim=lim,
        xlab='Genomic position (Mb)',ylab='Per-generation slope (x 0.001)',
        main=reg$title)
   points(xx,pr,pch=15,col='#C34D50',cex=.9)
   abline(h=0,col='grey80',lty=2)
   for(k in which(z$group%in%reg$highlight))
     text(xx[k],yy[k],labels=z$group[k],pos=3,cex=.75,font=2)
   if(chr=='I')legend('topleft',c('Observed phase-aligned','LOCO beta prediction'),
                      pch=c(16,15),col=c('#287AB5','#C34D50'),bty='n',cex=.72)
   # B: single *prespecified* focal term. Inference is conditional, exploratory.
   par(mar=c(4.3,4.6,3.5,1))
   f<-e[e$term==reg$term,,drop=FALSE]
   f<-f[match(TRAITS,f$trait),,drop=FALSE]
   f<-f[is.finite(f$effect_SD)&is.finite(f$SE_SD),,drop=FALSE]
   if(!nrow(f)){
     plot.new();title(paste('No focal effects',reg$term))
   }else{
     y<-rev(seq_len(nrow(f)))
     lo<-f$effect_SD-1.96*f$SE_SD;hi<-f$effect_SD+1.96*f$SE_SD
     lim<-range(c(lo,hi,0),finite=TRUE)
     lim<-lim+c(-1,1)*max(.12,.1*diff(lim))
     plot(f$effect_SD,y,ylim=c(.5,length(y)+.5),xlim=lim,
          yaxt='n',ylab='',xlab='Conditional RIL effect / trait SD',
          pch=16,col='#287AB5',main=paste('Focal SNP',substring(reg$term,2)))
     axis(2,at=y,labels=f$trait,las=1,cex.axis=.75)
     segments(lo,y,hi,y,col='#287AB5',lwd=1.25)
     segments(lo,y-.11,lo,y+.11,col='#287AB5')
     segments(hi,y-.11,hi,y+.11,col='#287AB5')
     abline(v=0,lty=2,col='grey80')
     # Asterisk is only for tested focal term, with whole-region correction.
     idx<-which(is.finite(f$q_BH_all_regional_terms) & f$q_BH_all_regional_terms<.05)
     if(length(idx))text(f$effect_SD[idx],y[idx]+.31,'*',cex=1,col='#C34D50')
   }
   # C: exact site trajectories from historical timepoints, no interpolation across sites.
   par(mar=c(4,4.7,3.5,1))
   if(!nrow(a)){
     plot.new();title('Focal exact-site trajectories')
     text(.5,.6,'No exact read counts for focal variants',cex=.95)
     text(.5,.43,'Do not reconstruct from R5 cluster slopes',cex=.7,col='grey45')
   }else{
     maxgen<-max(a$generation,na.rm=TRUE)
     plot(NA,xlim=c(0,maxgen+2),ylim=c(0,1),xlab='Generation',
          ylab='A00-minor allele frequency',main='Historical focal trajectories')
     abline(h=c(0,.5,1),lty=3,col='grey90')
     ids<-paste0(chr,'_',reg$focus)
     for(j in seq_along(ids)){
       b<-a[a$SNP_ID==ids[j],,drop=FALSE]
       if(!nrow(b))next
       b0<-b[b$population=='A00'&b$generation==0,,drop=FALSE]
       if(nrow(b0)!=1L)stop('Expected exactly one A00 at ',ids[j])
       for(pop in names(PALETTE)){
         v<-b[b$population==pop & b$generation>0,,drop=FALSE]
         if(!nrow(v))next
         v<-v[order(v$generation),,drop=FALSE]
         xs<-c(0,v$generation)
         ys<-c(b0$rare_frequency,v$rare_frequency)
         lines(xs,ys,col=PALETTE[[pop]],lty=j,lwd=1.3,type='b',pch=c(16,15,17)[j],cex=.6)
       }
     }
     legend('topleft',names(PALETTE),col=unname(PALETTE),lty=1,
            bty='n',cex=.7,horiz=TRUE)
     legend('topright',as.character(reg$focus),lty=seq_along(reg$focus),
            col='grey35',bty='n',cex=.62,title='SNP position')
   }
 }
 mtext('Three candidate regions: observations, cross-chromosome predictions, and exploratory RIL associations',
       side=3,outer=TRUE,line=1,cex=1.05,font=2)
 mtext('Not causal proof; linked R5 groups are not independent loci. Asterisk: regional-term BH q < 0.05. Historical trajectories end at G50.',
       side=1,outer=TRUE,line=1.1,cex=.72,col='grey30')
}
pdf(file.path(CFG$out,'66_Integrated_three_regions.pdf'),width=17.2,height=12.5,useDingbats=FALSE)
draw();dev.off()
png(file.path(CFG$out,'66_Integrated_three_regions.png'),width=17.2,height=12.5,units='in',res=CFG$dpi)
draw();dev.off()
qc<-do.call(rbind,QC)
write.csv(qc,file.path(CFG$out,'66_trajectory_completeness.csv'),row.names=FALSE)
writeLines(c(
 'Panels A: primary 802-R5 catalogue, phase-aligned median responses and LOCO beta predictions.',
 'Panels B: conditionally estimated *RIL associations* in trait SD units, not variant causation.',
 'Asterisks: q_BH_all_regional_terms<0.05, not significance after all possible follow-up analyses.',
 'Panels C: historical exact-site counts, G0/G10/G36/G50; exclude A6G140 and extinct GA3.',
 'MISSING IV panel C is a data-availability gap, not no evolutionary response.',
 'Do not infer gene consequence from physical proximity or unvalidated online gene predictions.',
 'Selection candidates remain conditional on unresolved extreme-tail neutral-test calibration.'
 ),file.path(CFG$out,'66_figure_interpretation_notes.txt'))
print(qc)
