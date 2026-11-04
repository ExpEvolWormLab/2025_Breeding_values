# 63_Three_regions_manuscript_figures.R
# Manuscript figures for the THREE frozen historical candidate regions.
# No causal inference, linkage independence or selection significance is implied.
# Run from the project root: source('Oct_pipeline_cleaner/63_Three_regions_manuscript_figures.R')
# Dependencies: base R only. Figures: vector PDFs and high-resolution PNGs.

CFG <- list(
  out = 'selection_results/manuscript_finalization/63_three_regions',
  main_candidates = c('selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
                      '45_All_R5_clusters_GWAS_evolution_predictions.csv'),
  gwas_candidates = c('selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
                      'selection_results/drift_oneG0_nullfit/11_Clusters_selection_and_GWAS_09f.csv',
                      '11_Clusters_selection_and_GWAS_09f.csv'),
  trajectory_candidates = c('historical_QTL_R5/deep_dive_extended/evolution/28_ALT_read_counts_all_pools.csv',
                            'historical_QTL_R5/deep_dive_extended/evolution/30c_ALT_read_counts_all_pools.csv',
                            'selection_results/chrV_deep_dive/56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv',
                            'chrV_deep_dive/56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv',
                            '56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv'),
  show_ga3 = FALSE,  # GA3 extinct; use ONLY as separately labelled supplementary analysis
  dpi = 320L
)
regions <- list(
  I  = list(lo=5883574L,hi=7996414L,focus=c(7167738L,7491656L,7596359L),
            tag=c('Shared RC301/PX174','RC301-private','RC301-private'),
            title='Chr I | RC301/PX174-linked alternatives'),
  IV = list(lo=8650000L,hi=10600000L,focus=c(8892663L,10499135L),
            tag=c('Domestication candidate','Domestication mutation (tax-6 region)'),
            title='Chr IV | CB4856-linked background and tax-6 candidate'),
  V  = list(lo=6046042L,hi=7278103L,focus=c(6322080L,6637695L,6826421L),
            tag=c('Post-founder mutation (cest-24)','CB4507-linked','CB4507-linked'),
            title='Chr V | linked multilocus response')
)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
blue <- '#1866A3'; red <- '#B73545'; grey <- '#677585'; green <- '#1F8F79'
first_file <- function(z) { f <- z[file.exists(z)]; if(length(f)) f[1] else NA_character_ }
read_csv <- function(path) read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
main_file <- first_file(CFG$main_candidates)
gwas_file <- first_file(CFG$gwas_candidates)
if(is.na(main_file)||is.na(gwas_file)) stop('Missing script-45 and/or script-11 CSV. Edit CFG file paths.')
X <- read_csv(main_file); G <- read_csv(gwas_file)
X <- X[X$scenario=='primary_survivors',,drop=FALSE]
if(anyDuplicated(X$group)||anyDuplicated(G$group))stop('Duplicated R5 group; cannot make figures')
needed <- c('group','chr','ps','pred_full','primary_slope_median_phase','n_eligible')
if(!all(needed %in% names(X)))stop('Script 45 table lacks expected columns')
cc <- paste0('beta_',traits,'_NaCl')
if(!all(cc %in% names(G)))stop('Script 11 table lacks seven GWAS effect columns')
G <- G[match(as.character(X$group),as.character(G$group)),,drop=FALSE]
if(anyNA(G$group))stop('R5 GWAS/prediction catalogue mismatch')
if(any(abs(X$pred_full - G$beta_beta_NaCl)>1e-7,na.rm=TRUE))
 warning('Original beta score and script-45 beta predictor differ: check orientation before interpreting effects')
# Trait-standardized coloring is descriptive only, never an inference about effect size.
trait_sd <- vapply(cc,function(z) sd(G[[z]],na.rm=TRUE),numeric(1))
if(any(!is.finite(trait_sd)|trait_sd==0)) stop('Invalid trait SD')
dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)

tr_files <- CFG$trajectory_candidates[file.exists(CFG$trajectory_candidates)]
# With heterogeneous CSV schemas, read each source separately and only keep recognized columns.
read_trajectories <- function(chr, ids){
 out <- list()
 for(f in tr_files){
   a <- read_csv(f)
   if(!('SNP_ID'%in%names(a)))next
   a <- a[a$SNP_ID %in% ids,,drop=FALSE]
   if(!nrow(a))next
   if(!'generation'%in%names(a)){
     if('generations'%in%names(a))a$generation<-a$generations else next
   }
   if(!'ALT_frequency'%in%names(a))next
   if(!'population'%in%names(a)){
      if(!'replicate'%in%names(a))next
      a$population<-as.character(a$replicate)
      a$population[a$population=='3']<-'GA1';a$population[a$population=='4']<-'GA2'
      a$population[a$population=='5']<-'GA3';a$population[a$population=='6']<-'GA4'
   }
   a$population[as.character(a$sample_id)=='A00'] <- 'A00'
   a$population <- as.character(a$population)
   if('replicate' %in% names(a)){
      r <- as.character(a$replicate)
      labs <- c('3'='GA1','4'='GA2','5'='GA3','6'='GA4')
      use <- r %in% names(labs)
      a$population[use] <- unname(labs[r[use]])
   }
   if('sample_id' %in% names(a))
      a$population[as.character(a$sample_id)=='A00'] <- 'A00'
   if(!CFG$show_ga3)a<-a[!grepl('GA3',a$population),,drop=FALSE]
   out[[length(out)+1L]]<-a[,c('SNP_ID','population','generation','ALT_frequency'),drop=FALSE]
 }
 if(!length(out))return(data.frame())
 a <- do.call(rbind,out)
 # duplicated exact rows can appear in newer and historical archives: retain one
 a <- a[!duplicated(a),,drop=FALSE]
 # Each exact SNP/pop/time must identify a unique frequency, fail if conflicting.
 key <- paste(a$SNP_ID,a$population,a$generation)
 if(anyDuplicated(key)){
   d <- split(a$ALT_frequency,key)
   if(any(vapply(d,function(t) diff(range(t,na.rm=TRUE))>1e-8,logical(1))))
      stop('Conflicting Pool-seq records across trajectory files for ',chr)
   a <- a[!duplicated(key),,drop=FALSE]
 }
 # Orient toward the allele *initially rare at A00*, not arbitrary ALT orientation.
 a$rare_frequency <- NA_real_
 for(id in unique(a$SNP_ID)){
    idx <- which(a$SNP_ID==id)
    anc <- a$ALT_frequency[idx[a$population[idx]=='A00' & a$generation[idx]==0]]
    if(length(anc)!=1L||!is.finite(anc))next
    a$rare_frequency[idx] <- if(anc>0.5) 1-a$ALT_frequency[idx] else a$ALT_frequency[idx]
 }
 a
}

plot_one <- function(chr, reg) {
  
  lo <- reg$lo
  hi <- reg$hi
  
  # Select evolutionary responses
  z <- X[
    X$chr == chr &
      X$ps >= lo &
      X$ps <= hi,
    ,
    drop = FALSE
  ]
  
  if (!nrow(z)) {
    stop("No R5 clusters within chromosome ", chr)
  }
  
  # Match GWAS effects using chromosome + R5 cluster ID
  key_X <- paste(z$chr, z$group, sep = "_")
  key_G <- paste(G$chr, G$group, sep = "_")
  
  if (anyDuplicated(key_G)) {
    stop("Duplicate R5 cluster IDs in GWAS table G")
  }
  
  idx <- match(key_X, key_G)
  
  if (anyNA(idx)) {
    stop(
      "Missing GWAS information for R5 clusters: ",
      paste(key_X[is.na(idx)], collapse = ", ")
    )
  }
  
  q <- G[idx, , drop = FALSE]
  
  # Verify exact cluster correspondence
  stopifnot(
    identical(
      paste(z$chr, z$group),
      paste(q$chr, q$group)
    )
  )
  
 if(!nrow(z)) stop('No R5 clusters within ',chr)
 pos<-z$ps/1e6
 known<-which(is.finite(z$primary_slope_median_phase))
 n_ok<-length(known)
 focal_ids<-paste0(chr,'_',reg$focus)
 a<-read_trajectories(chr,focal_ids)
 sumrow <- data.frame(chr=chr,region_start=lo,region_end=hi,
     n_R5_groups=nrow(z),n_eligible=sum(z$n_eligible>0,na.rm=TRUE),
     n_response_finite=n_ok,n_exact_focal_sites=if(nrow(a))length(unique(a$SNP_ID))else 0L,
     plot_ready_trajectories=nrow(a)>0,stringsAsFactors=FALSE)
 write.csv(z,file.path(CFG$out,paste0('63_',chr,'_R5_prediction_window.csv')),row.names=FALSE)
 if(nrow(a)) write.csv(a,file.path(CFG$out,paste0('63_',chr,'_focal_rare_allele_trajectories.csv')),row.names=FALSE)
 draw <- function(){
   op<-par(no.readonly=TRUE);on.exit(par(op),add=TRUE)
   par(oma=c(1.7,0,0,0))
   layout(matrix(c(1,1,2,3),nrow=2,byrow=TRUE),heights=c(1.05,1.05))
   par(mar=c(3.4,5,3.6,1.5),mgp=c(2.5,.8,0))
   ylim<-range(c(z$primary_slope_median_phase,0),na.rm=TRUE)
   pad<-max(1e-5,diff(ylim)*.2);ylim<-ylim+c(-pad,pad)
   plot(pos,z$primary_slope_median_phase,pch=ifelse(z$pred_full>0,16,1),cex=1.1,
        col=ifelse(z$pred_full>0,red,blue),xlim=c(lo,hi)/1e6,ylim=ylim,
        xlab='Genomic position (Mb)',ylab='Phase-aligned slope / generation',
        main=paste('A |',reg$title))
   abline(h=0,lty=2,col='grey55')
   abline(v=reg$focus/1e6,lty=3,col='grey70')
   abline(v=pos[z$group%in%c(53,54,727,731,848,851,853,854)],lty=3,col='grey85')
   for(i in seq_len(nrow(z)))if(is.finite(z$primary_slope_median_phase[i]))
      text(pos[i],z$primary_slope_median_phase[i],z$group[i],pos=3,cex=.68,col='#333333')
   legend('topleft',legend=c('Beta effect positive','Beta effect negative'),
     col=c(red,blue),pch=c(16,1),bty='n',cex=.73)
   mtext(sprintf('%d R5 groups; %d finite cluster responses | not independent loci',nrow(z),n_ok),
         side=3,line=.35,cex=.67,col=grey)
   par(mar=c(5.6,5.4,3.2,1))
   # Heatmap rows are seven GWAS effects standardized by genome-wide effect SD.
   zz <- t(as.matrix(q[, cc, drop = FALSE]))
   zz <- sweep(zz, 1, trait_sd, "/")
   
   zz[!is.finite(zz)] <- 0
   
   # Truncate values without dropping matrix dimensions
   zz[zz < -3] <- -3
   zz[zz >  3] <-  3
   
   stopifnot(
     is.matrix(zz),
     nrow(zz) == length(traits),
     ncol(zz) == nrow(q)
   )
   cols <- grDevices::colorRampPalette(c('#28558A','#EAF0F4','#BC3448'))(151)
   image(seq_len(ncol(zz)),seq_len(nrow(zz)),t(zz),col=cols,zlim=c(-3,3),
      axes=FALSE,xlab='',ylab='',main='B | Relative seven-trait GWAS effects')
   axis(1,at=seq_len(ncol(zz)),labels=z$group,las=2,cex.axis=.68)
   axis(2,at=seq_len(nrow(zz)),labels=traits,las=1,cex.axis=.77)
   mtext('R5 cluster ID',side=1,line=3.9,cex=.75)
   box()
   par(mar=c(5.6,4.4,3.2,1.5))
   if(!nrow(a)||all(!is.finite(a$rare_frequency))){
      plot.new(); title('C | Exact Pool-seq trajectories')
      text(.5,.60,'Not available in current inputs',cex=1.05)
      text(.5,.42,'Provide script 28/30c exact-site counts',cex=.8,col=grey)
   }else{
      aa<-a[is.finite(a$rare_frequency),,drop=FALSE]
      plot(NA,xlim=c(0,max(aa$generation,na.rm=TRUE)),ylim=c(0,1),
           xlab='Generation',ylab='Initially rare allele frequency',
           main='C | Focal exact-SNP trajectories')
      abline(h=c(0,.5,1),lty=3,col='grey85')
      palette<-c('GA1'='#2676AA','GA2'='#D28B38','GA4'='#298269','GA3'='#888888')
      ids<-unique(aa$SNP_ID); pchs<-setNames(seq(15,15+length(ids)-1),ids)
      for(id in ids){
         b<-aa[aa$SNP_ID==id,,drop=FALSE]
         b0<-b[b$population=='A00' & b$generation==0,,drop=FALSE]
         if(nrow(b0))points(0,b0$rare_frequency[1],pch=pchs[id],col='black',cex=1.1)
         for(pop in names(palette)){
            v<-b[b$population==pop & b$generation>0,,drop=FALSE]
            if(!nrow(v))next
            v<-v[order(v$generation),,drop=FALSE]
            if(nrow(b0))v<-rbind(b0[1,,drop=FALSE],v)
            lines(v$generation,v$rare_frequency,col=palette[pop],lwd=1.4,
                  lty=match(id,ids),type='b',pch=pchs[id],cex=.55)
         }
      }
      legend('topleft',legend=names(palette)[names(palette)%in%aa$population],
             col=palette[names(palette)%in%aa$population],lwd=2,bty='n',cex=.73)
      legend('bottomright',legend=sub('^[IVX]+_','',ids),pch=unname(pchs),
             bty='n',cex=.68,ncol=1,title='Focal positions')
   }
   mtext('Descriptive regional analysis | R5 clusters can be linked; no causal or independent-selection attribution',
         side=1,outer=TRUE,line=-1,cex=.75,col=grey)
 }
 pdf(file.path(CFG$out,paste0('63_Chr',chr,'_regional_figure.pdf')),
     width=13.2,height=8.8,useDingbats=FALSE);draw();dev.off()
 png(file.path(CFG$out,paste0('63_Chr',chr,'_regional_figure.png')),
     width=13.2,height=8.8,units='in',res=CFG$dpi);draw();dev.off()
 sumrow
}
res<-do.call(rbind,Map(plot_one,names(regions),regions))
write.csv(res,file.path(CFG$out,'63_Figure_input_and_completeness_QC.csv'),row.names=FALSE)
writeLines(c('All three vector figures completed.',
 'Missing trajectory panel is labelled explicitly; no inferred trajectories are generated.',
 'Cluster genomic position in the main panel is the R5 representative position, not the full interval.',
 'Panel A holds the phase-aligned response from the primary GA1/GA2/GA4 analysis.',
 'Panel B standardizes each trait effect by its genome-wide across-cluster SD for visualization ONLY.',
 'Panel C uses counts/ALT allele frequencies and orients to the A00-minor allele; not phase-aligned slopes.',
 'Use complete sample manifest and ascertainment QC in Methods.',
 capture.output(sessionInfo())),file.path(CFG$out,'63_figure_methods_notes.txt'))
print(res)
