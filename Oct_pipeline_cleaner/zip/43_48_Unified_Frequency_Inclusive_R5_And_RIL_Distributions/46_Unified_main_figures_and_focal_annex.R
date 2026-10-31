## 46_Unified_main_figures_and_focal_annex.R
## Main: GA1 GA2 GA4, all ancestral frequencies, 1 row per R5 cluster.
## Supplement: GA3-inclusive response, rare focal raw-SYNC variants.
source('43_Unified_config_and_helpers.R')
cfg <- uf_config
out <- cfg$prefix
uf_check(c(file.path(out,'43_Unified_R5_cluster_registry.csv'),
 file.path(out,'44_primary_survivors_WF_cluster_selection.csv'),
 file.path(out,'45_All_R5_clusters_GWAS_evolution_predictions.csv'),
 file.path(out,'45_LOCO_performance_main_and_supp.csv')))
x <- uf_read(file.path(out,'45_All_R5_clusters_GWAS_evolution_predictions.csv'))
met <- uf_read(file.path(out,'45_LOCO_performance_main_and_supp.csv'))
x <- x[x$scenario=='primary_survivors',,drop=FALSE]
if(anyDuplicated(x$group))stop('One row per R5 cluster required')
reg <- uf_read(file.path(out,'43_Unified_R5_cluster_registry.csv'))
test <- uf_read(file.path(out,'44_primary_survivors_WF_cluster_selection.csv'))
if(!identical(sort(as.character(reg$group)),sort(as.character(x$group))))
 stop('Registry/prediction cluster universes are not identical')
red <- reg[order(match(reg$group, x$group)),,drop=FALSE]
if(anyNA(red$group))stop('Registry mismatch')
qc <- data.frame(metric=c('R5_GWAS_total','R5_with_any_Poolseq','R5_response_eligible',
 'R5_with_calibrated_tag_test','R5_selected_BY_q05','R5_zero_A00_ALT_test_tags',
 'R5_non_selected_but_in_prediction'),
 value=c(nrow(reg),sum(reg$n_present>0),sum(reg$n_eligible>0),nrow(test),
 sum(test$significant_BY_05,na.rm=TRUE),
 sum(test$A00_ALT_reads==0,na.rm=TRUE),
 sum(x$n_eligible>0 & !(x$selection_significant_BY_05 %in% TRUE),na.rm=TRUE)))
uf_write(qc,file.path(out,'46_Manuscript_counts.csv'))
if(nrow(test)<50L)warning('Very few tested clusters; inspect 43/44 QC')

pal <- c(full='#3266a8',reduced='#cc8e38',SB='#b65468',direct_beta='#65945b')
fig <- file.path(out,'46_Main_unified_frequency_inclusive_R5.pdf')
pdf(fig,width=11,height=8.3,onefile=TRUE)
par(mfrow=c(2,2),mar=c(4.4,4.5,3,1),oma=c(0,0,2,0))
## panel 1: ALL tested tags, not only positive discoveries
isok <- is.finite(x$selection_p_MC_upper_simultaneous) &
        is.finite(x$median_A00_MAF)
if(sum(isok)>10L){
 z <- x[isok,,drop=FALSE]
 xx <- pmax(z$median_A00_MAF,1e-5)
 yy <- -log10(pmax(as.numeric(z$selection_p_MC_upper_simultaneous),1e-300))
 plot(xx,yy,log='x',pch=16,cex=.53,col=adjustcolor('grey30',alpha.f=.45),
    xlab='Median ancestral MAF of eligible cluster SNPs',ylab='-log10 conservative MC p upper bound',
    main='A. Frequency-inclusive cluster tests')
 hi <- which(z$selection_significant_BY_05 %in% TRUE)
 if(length(hi))points(xx[hi],yy[hi],pch=16,cex=.85,col='#b2182b')
 legend('topright',legend=c('All tested R5 tags','BY q<0.05'),
   col=c('grey30','#b2182b'),pch=16,cex=.75,bty='n')
} else plot.new()
## panel 2: effect response, including all eligible/non-selected clusters
ok <- is.finite(x$pred_full)&is.finite(x$primary_slope_median_phase)
if(sum(ok)>15L){
 z <- x[ok,,drop=FALSE]
 plot(z$pred_full,z$primary_slope_median_phase,pch=16,cex=.64,
  col=adjustcolor(pal['full'],alpha.f=.35),
  xlab='GWAS-reconstructed seven-trait effect on beta',
  ylab='Median phase-aligned allele-frequency slope / generation',
  main='B. Independent GWAS effects vs evolution')
 abline(h=0,v=0,lty=3,col='grey60')
 f <- lm(primary_slope_median_phase~pred_full,data=z)
 abline(f,col=pal['full'],lwd=2)
 legend('topleft',legend=sprintf('n=%d; Pearson r=%.3f',nrow(z),cor(z$pred_full,z$primary_slope_median_phase)),
        bty='n',cex=.77)
} else plot.new()
## panel 3: held-out chromosome performance, no source-pop folds
m <- met[met$scenario=='primary_survivors',,drop=FALSE]
inds <- match(c('full','reduced','SB','direct_beta'),m$model)
m <- m[inds[!is.na(inds)],,drop=FALSE]
if(nrow(m)){
 val <- m$LOCO_R2;names(val) <- m$model
 bp <- barplot(val,col=pal[m$model],border=NA,
   ylim=range(c(-.08,val,0),finite=TRUE)*c(1.25,1.25),
   ylab='Leave-one-chromosome-out predictive R2',
   main='C. Evolution prediction (out-of-chromosome)')
 abline(h=0,col='grey50',lty=2)
 text(bp,val,labels=sprintf('%.3f',val),pos=ifelse(val>=0,3,1),cex=.78)
} else plot.new()
## panel 4: size of full study universe with QC and tests
levels <- c('Original R5 catalogue','Poolseq represented','Usable phase slopes','WF tested tags','BY q<0.05')
n <- c(nrow(reg),sum(reg$n_present>0),sum(reg$n_eligible>0),nrow(test),
       sum(test$significant_BY_05,na.rm=TRUE))
bp <- barplot(n,names.arg=levels,las=2,col=c('grey70','grey55','#4f8baf','#386da6','#b2182b'),
   ylab='Number of distinct LD clusters',main='D. Inclusion and selection accounting',cex.names=.72)
text(bp,n,labels=n,pos=3,cex=.8)
mtext('Primary population set: GA1, GA2, GA4 (GA3 extinct)',outer=TRUE,cex=1.05)
dev.off()

## Separate supplement: GA3-inclusive sensitivity and historical 675 baseline
if('supp_slope_median_phase'%in%names(reg)){
 ss <- x[,c('group','primary_slope_median_phase','supp_slope_median_phase'),drop=FALSE]
 uf_write(ss,file.path(out,'46_Supp_all_GA3_inclusive_cluster_slopes.csv'))
 pdf(file.path(out,'46_Supp_GA3_inclusive_comparison.pdf'),width=8,height=7)
 ok <- is.finite(ss$primary_slope_median_phase)&is.finite(ss$supp_slope_median_phase)
 if(sum(ok)>5L){
  plot(ss$primary_slope_median_phase[ok],ss$supp_slope_median_phase[ok],
    pch=16,col=adjustcolor('steelblue',alpha.f=.4),
    xlab='Primary: GA1 + GA2 + GA4 slope',ylab='Supplement: including extinct GA3',
    main='Effect of including GA3 (not the primary inference)')
  abline(a=0,b=1,lty=2,col='grey50')
 }else plot.new()
 dev.off()
}

## Annex: nominated variants previously omitted by the founder-R5 filter.
## They are NOT added as independent R5 clusters and NEVER enter 45/46
## genome-wide predictions or their BY multiplicity family.
rawf <- cfg$raw_focal_file
tf <- cfg$previous_focal_test_file
if(file.exists(rawf)){
 a <- uf_read(rawf)
 must <- c('SNP_ID','sample_id','generation','ALT_reads','depth','ALT_frequency','REF_agrees')
 if(!all(must%in%names(a)))stop('Malformed focal raw-SYNC file')
 if(any(a$REF_agrees%in%FALSE,na.rm=TRUE))stop('Focal variant reference allele mismatch')
 a <- a[a$sample_id != 'A6G140',,drop=FALSE]
 r5ids <- unique(as.character(reg$tag_SNP_ID[!is.na(reg$tag_SNP_ID)]))
 a$already_R5_primary_selection_tag <- a$SNP_ID%in%r5ids
 if(file.exists(tf)){
   t <- uf_read(tf)
   t <- t[t$scenario=='retained_GA'&as.numeric(t$Ne)==cfg$Ne,,drop=FALSE]
   if(anyDuplicated(t$SNP_ID))stop('Duplicate focal neutral test per SNP/scenario/Ne')
   j <- match(a$SNP_ID,t$SNP_ID)
   a$focal_target_WF_p <- t$p_WF_conditional[j]
   a$focal_target_simulations <- t$n_sim[j]
 }
 uf_write(a,file.path(out,'46_Supp_focal_raw_SYNC_annex_NOT_genomewide_FDR.csv'))
 pdf(file.path(out,'46_Supp_named_focal_trajectories.pdf'),width=11,height=7,onefile=TRUE)
 ids <- unique(a$SNP_ID)
 for(id in ids){
    z <- a[a$SNP_ID==id & a$sample_id!='A6G140',,drop=FALSE]
    z <- z[is.finite(z$ALT_frequency)&is.finite(z$generation),,drop=FALSE]
    if(!nrow(z))next
    par(mar=c(5,5,3,2))
    plot(NA,xlim=range(z$generation),ylim=c(0,max(.05,1.08*max(z$ALT_frequency))),
      xlab='Generation',ylab='ALT-read frequency',main=paste('Targeted focal variant:',id))
    baseline <- z[z$sample_id=='A00',,drop=FALSE]
    if(nrow(baseline))points(0,baseline$ALT_frequency[1],pch=16,col='black',cex=1.25)
    labs <- c('3'='GA1','4'='GA2','5'='GA3 (extinct)','6'='GA4')
    clrs <- c('3'='#2979b8','4'='#d17a25','5'='grey65','6'='#46884d')
    for(r in c('3','4','5','6')){
      h <- z[as.character(z$replicate)==r & z$generation>0,,drop=FALSE]
      if(!nrow(h))next
      h <- h[order(h$generation),,drop=FALSE]
      lines(h$generation,h$ALT_frequency,col=clrs[r],lwd=2,
            lty=if(r=='5')2 else 1,type='b',pch=16)
    }
    legend('topleft',legend=c('A00','GA1','GA2','GA4','GA3 extinct (supplement)'),
      pch=16,lty=c(NA,1,1,1,2),col=c('black',clrs[c('3','4','6','5')]),bty='n',cex=.8)
    mtext('Focal p-values are TARGETED; not R5 genome-wide FDR claims',side=3,line=.1,cex=.72)
 }
 dev.off()
}
writeLines(capture.output(sessionInfo()),file.path(out,'46_sessionInfo.txt'))
cat('46 DONE. Main and supplementary figure PDFs in: ',out,'\n',sep='')
