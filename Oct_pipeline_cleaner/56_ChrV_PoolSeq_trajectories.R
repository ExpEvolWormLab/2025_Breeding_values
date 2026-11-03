## 56_ChrV_PoolSeq_trajectories.R
## Examine EXACT-SNP Pool-seq trajectories and phase-aligned frequencies in V_3.
## Relies on the corrected 08/09 sample orientation and optional 43x raw-SYNC
## integration. Does NOT perform new selection tests or assume independent loci.
## From the project root, after 54 (and optionally 55).

options(stringsAsFactors=FALSE)
out_dir <- 'selection_results/chrV_deep_dive'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
fpath <- function(paths){x<-paths[file.exists(paths)];if(length(x))x[1L]else NA_character_}
need <- function(paths,what){z<-fpath(paths);if(is.na(z))stop('Missing ',what,': ',paste(paths,collapse='; '));z}
readtab <- function(f)read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
W <- function(x,n)write.csv(x,file.path(out_dir,n),row.names=FALSE,na='')
A <- readtab(need(c(file.path(out_dir,'54_V3_12_cluster_prediction_audit.csv')),'script 54 output'))
phase_path <- need(c('selection_results/unified_frequency_inclusive_R5/43_All_SNP_phase_slope_and_inclusion_QC.csv',
                     'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv'),'R5 phase membership')
phase <- readtab(phase_path)
if(!all(c('rs','group','phase')%in%names(phase))||anyDuplicated(phase$rs))
 stop('Invalid phase table or duplicate sites')
phase$group<-as.character(phase$group)
phase$phase<-as.numeric(phase$phase)
if(any(!phase$phase%in%c(-1,1)))stop('R5 phase must be +/-1')
phase<-phase[phase$group%in%as.character(A$group),,drop=FALSE]
if(!nrow(phase))stop('No phased V_3 SNPs')

source(need(c('Oct_pipeline_cleaner/09f_drift_functions.R','09f_drift_functions.R',
              'drift_nullfit_update/09f_drift_functions.R'),'corrected 09f helper'),local=FALSE)
source(need(c('Oct_pipeline_cleaner/43_Unified_config_and_helpers.R',
              '43_Unified_config_and_helpers.R'),'43 config/helper'),local=FALSE)
source(need(c('Oct_pipeline_cleaner/43x_Integrate_rawSYNC_R5_SNPs.R',
              '43x_Integrate_rawSYNC_R5_SNPs.R'),'43x raw-SYNC integrator'),local=FALSE)
cfg<-uf_config
uf_check(c(cfg$pool_file,cfg$orientation_file,cfg$status_file))

## Recover low-frequency V sites, including 6322080, 6637695 and 6826421
## when original 34a data are available. Original sample columns are checked.
if(file.exists(cfg$raw_focal_file)){
  p<-uf_pool_with_R5_raw_recovery(cfg)
} else {
  warning('34a raw recovered counts missing: rare focal SNPs may be unavailable; ',
          'missing sites will be listed in QC, not silently substituted')
  p<-read_drift_inputs(cfg$pool_file,cfg$orientation_file,exclude_G0='A6G140')
}
if(!identical(as.character(p$pop$sample_id),colnames(p$k_alt)) ||
   !identical(dim(p$k_alt),dim(p$depth)))stop('Pool sample/order validation failed')
if(anyDuplicated(p$sel$ID))stop('Duplicated oriented Poolseq SNP IDs')
if(sum(p$pop$sample_id=='A00')!=1L || sum(p$pop$generations==0)!=1L)
 stop('Expected precisely the A00 ancestral pool after excluding A6G140')
status<-readtab(cfg$status_file)
if(!all(c('replicate','population','status')%in%names(status)))stop('Invalid GA status')
if(!all(c('GA1','GA2','GA3','GA4')%in%status$population))stop('Missing GA statuses')
if(!identical(as.character(status$status[match(c('GA1','GA2','GA3','GA4'),status$population)]),
              c('survived','survived','extinct','survived')))
  stop('Unexpected GA status (main = GA1, GA2, GA4; GA3 is supplementary)')
pop<-p$pop
pop$population<-as.character(status$population[match(pop$replicate,status$replicate)])
pop$status<-as.character(status$status[match(pop$replicate,status$replicate)])
pop$status[pop$generations==0]<-'baseline'
pop$population[pop$generations==0]<-'A00'
if(anyNA(pop$population))stop('Unmapped Poolseq replicate in sample manifest')

## Sample a small, reproducible set of exact SNPs per cluster.
## Selection tags were chosen using GWAS p and coverage, not response sign.
choose<-list(); all_candidates<-list()
extra_markers<-fpath(c(file.path(out_dir,'55_Selected_RIL_markers.csv')))
extra<-if(!is.na(extra_markers))readtab(extra_markers) else data.frame()
for(i in seq_len(nrow(A))){
 g<-as.character(A$group[i]); mem<-phase[phase$group==g,,drop=FALSE]
 if(!nrow(mem))next
 pr<-character()
 if('tag_for_neutral_test'%in%names(mem))
   pr<-c(pr,as.character(mem$rs[mem$tag_for_neutral_test%in%TRUE]))
 if('neutral_test_tag'%in%names(A))pr<-c(pr,as.character(A$neutral_test_tag[i]))
 pr<-c(pr, intersect(c('V_6322080','V_6300718','V_6637695','V_6826421','V_6959758'),mem$rs))
 if(nrow(extra))pr<-c(pr,as.character(extra$rs[as.character(extra$group)==g]))
 pr<-unique(pr[pr%in%mem$rs])
 if('gwas_min_p'%in%names(mem)){
   ranking<-mem$rs[order(mem$gwas_min_p,mem$rs,na.last=TRUE)]
 } else ranking<-mem$rs[order(abs(as.numeric(sub('^.*_','',mem$rs))-A$ps[i]))]
 pr<-unique(c(pr,head(ranking,3L)))
 all_candidates[[length(all_candidates)+1L]]<-data.frame(group=g,rs=pr,stringsAsFactors=FALSE)
 # Limit to 4 sites per group: always retain the prespecified neutral tag
 # and rare variants before other representatives. Export full candidate list.
 choose[[length(choose)+1L]]<-data.frame(group=g,rs=head(pr,4L),stringsAsFactors=FALSE)
}
key<-unique(do.call(rbind,choose));cand<-unique(do.call(rbind,all_candidates))
W(cand,'56_Candidate_exact_SNP_manifest.csv')
key$phase<-phase$phase[match(key$rs,phase$rs)]
key$pos<-as.numeric(sub('^.*_','',key$rs))
key$is_wf_test_tag<-if('tag_for_neutral_test'%in%names(phase))
 phase$tag_for_neutral_test[match(key$rs,phase$rs)] else
 key$rs%in%as.character(A$neutral_test_tag)
key$in_oriented_poolseq<-key$rs%in%as.character(p$sel$ID)
W(key,'56_Exact_SNP_coverage_and_phase_QC.csv')
missing<-key$rs[!key$in_oriented_poolseq]
if(length(missing))warning('Sites without oriented PoolSeq counts: ',paste(missing,collapse=', '))
key<-key[key$in_oriented_poolseq,,drop=FALSE]
if(!nrow(key))stop('None of the nominated exact sites is present in corrected Poolseq')

traj<-vector('list',nrow(key))
perrep<-list()
for(j in seq_len(nrow(key))){
 id<-key$rs[j];idx<-match(id,as.character(p$sel$ID))
 d<-as.numeric(p$depth[idx,]);k<-as.numeric(p$k_alt[idx,])
 ok<-is.finite(d)&is.finite(k)&d>0&k>=0&k<=d
 f<-ifelse(ok,k/d,NA_real_)
 aligned<-if(key$phase[j]==1)f else 1-f
 df<-data.frame(group=key$group[j],SNP_ID=id,position=key$pos[j],
  phase=key$phase[j],is_wf_test_tag=key$is_wf_test_tag[j],
  sample_id=pop$sample_id,replicate=pop$replicate,population=pop$population,
  status=pop$status,generation=pop$generations,depth=d,ALT_reads=k,
  ALT_frequency=f,phase_aligned_frequency=aligned,
  approx_sampling_SE=sqrt(f*(1-f)/d),valid_depth_and_count=ok,
  stringsAsFactors=FALSE)
 traj[[j]]<-df
 base<-df$phase_aligned_frequency[df$sample_id=='A00']
 if(length(base)!=1L)stop('A00 baseline not unique')
 for(ga in c('GA1','GA2','GA3','GA4')){
   z<-df[df$population==ga&df$generation>0&is.finite(df$phase_aligned_frequency),,drop=FALSE]
   if(!nrow(z))next
   last<-z[which.max(z$generation),,drop=FALSE]
   perrep[[length(perrep)+1L]]<-data.frame(group=key$group[j],SNP_ID=id,
     population=ga,status=z$status[1],last_generation=last$generation[1],
     n_evolved_timepoints=nrow(z),A00_aligned_frequency=base,
     last_aligned_frequency=last$phase_aligned_frequency,
     change_from_A00=last$phase_aligned_frequency-base,
     endpoint_rate_per_generation=(last$phase_aligned_frequency-base)/last$generation,
     stringsAsFactors=FALSE)
 }
}
TT<-do.call(rbind,traj)
S<-if(length(perrep))do.call(rbind,perrep)else data.frame()
W(TT,'56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv')
if(nrow(S))W(S,'56_Per_population_endpoint_changes_DESCRIPTIVE.csv')

## Distinct source of biological variation: selected SNPs in the same R5 LD
## cluster, and across nearby clusters, may follow identical haplotypes.
## Do not report endpoint differences as independent site selection tests.
plot_regions<-function(primary=TRUE){
  pdf(file.path(out_dir,if(primary)
    '56_V3_Surviving_GA1_GA2_GA4_exact_trajectories.pdf' else
    '56_V3_GA3_extinct_exploratory_trajectories.pdf'),width=11,height=8.5)
  for(g in as.character(A$group)){
    ix<-key$group==g
    sites<-as.character(key$rs[ix]);if(!length(sites))next
    group_dat<-TT[TT$group==g,,drop=FALSE]
    allowed <- if(primary)c('A00','GA1','GA2','GA4') else c('A00','GA3')
    group_dat<-group_dat[group_dat$population%in%allowed,,drop=FALSE]
    ymax<-range(group_dat$phase_aligned_frequency,na.rm=TRUE)
    if(!all(is.finite(ymax)))next
    par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(0,0,3,0))
    for(sid in head(sites,4)){
      z<-group_dat[group_dat$SNP_ID==sid,,drop=FALSE]
      plot(NA,xlim=c(0,max(z$generation,na.rm=TRUE)),ylim=c(0,1),
        xlab='Generation',ylab='Phase-aligned frequency',
        main=paste('V:',g,' / ',sid))
      abline(h=c(0,.5,1),lty=3,col='grey85')
      b<-z[z$sample_id=='A00',,drop=FALSE]
      if(nrow(b)&&is.finite(b$phase_aligned_frequency[1]))
        points(0,b$phase_aligned_frequency,pch=4,cex=1.5,lwd=2)
      gs<-if(primary)c('GA1','GA2','GA4')else 'GA3'
      co<-c(GA1='#2C7BB6',GA2='#1A9850',GA3='#C51B7D',GA4='#D95F02')
      for(ga in gs){
        zz<-z[z$population==ga&z$generation>0&is.finite(z$phase_aligned_frequency),,drop=FALSE]
        if(!nrow(zz))next
        zz<-zz[order(zz$generation),,drop=FALSE]
        if(nrow(b)&&is.finite(b$phase_aligned_frequency[1]))
          zz<-rbind(transform(b,population=ga),zz)
        lines(zz$generation,zz$phase_aligned_frequency,col=co[ga],lwd=2,
          lty=if(ga=='GA3')2 else 1)
        points(zz$generation,zz$phase_aligned_frequency,col=co[ga],pch=16,cex=.9)
      }
      legend('topleft',legend=c('A00',gs),pch=c(4,rep(16,length(gs))),
             lty=c(NA,rep(1,length(gs))),col=c('black',co[gs]),cex=.75,bty='n')
    }
    if(length(sites)<4L){for(k in seq_len(4L-length(sites)))plot.new()}
    mtext(paste('Exact-SNP descriptive trajectory, R5 cluster',g,
         '; shared A00, not an independent replicate'),outer=TRUE,line=1,cex=.85)
  }
  dev.off()
}
plot_regions(TRUE)
if(any(TT$population=='GA3'))plot_regions(FALSE)
cat('SCRIPT 56 COMPLETED: ',length(unique(TT$SNP_ID)),
    ' exact SNPs; outputs: ',out_dir,'\n',sep='')
cat('Do not interpret within-replicate correlated trajectories, especially the',
    ' common A00 baseline, as independent selective evidence.\n')
