## 28_Regional_Poolseq_trajectories_and_selection.R
## Extract ALT read counts for historical SNPs, nominated variants and R5 tags.
## Reports FULL-vs-survivor drift-aware SNP tests where they exist; for low-p0
## loci reports descriptive per-population trajectories, NOT invalid q-values.
## No new fitness or selection coefficient is estimated by this script.

helper <- if(file.exists('09f_drift_functions.R')) '09f_drift_functions.R' else
          'drift_nullfit_update/09f_drift_functions.R'
source(helper)
script26_dir <- 'historical_QTL_R5/deep_dive'
out_dir <- file.path(script26_dir,'evolution')
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
script26_tags <- file.path(script26_dir,'26_R5_model_tag_SNPs.csv')
script26_focal <- '26_Focal_variants.csv'
legacy <- '23_Historical_region_SNPs.csv'
legacy_status <- 'historical_QTL_R5/23_Historical_SNP_status.csv'
selection_file <- 'selection_results/09_Selection_SNP_results.csv'
pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
full_snp_file <- 'selection_results/drift_oneG0_nullfit/09f_SNP_drift_Ne1000.csv'
survivor_snp_file <- 'selection_results/survivors_only_Ne1000/20_SNP_drift_Ne1000_survivors.csv'
full_cluster_file <- 'historical_QTL_R5/23_R5_all_cluster_boundaries_and_associations.csv'
survivor_cluster_file <- 'selection_results/survivors_only_Ne1000/20_Survivor_cluster_slope_summary.csv'
phase_file <- 'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv'
status_file <- '20_GA_population_status.csv'
include_all_historical_positions <- TRUE
exclude_G0_from_selection <- 'A6G140' # but EXPORT BOTH G0 samples
required<-c(script26_tags,script26_focal,selection_file,pool_file,full_snp_file,full_cluster_file)
if(any(!file.exists(required)))stop('Missing prerequisites: ',paste(required[!file.exists(required)],collapse=', '))

T<-read.csv(script26_tags,stringsAsFactors=FALSE)
F<-read.csv(script26_focal,stringsAsFactors=FALSE)
key<-rbind(data.frame(SNP_ID=T$tag_SNP_ID,region=T$region,source='R5_tag',
                       R5_group=as.character(T$group)),
 data.frame(SNP_ID=paste0(F$CHROM,'_',F$POS),region=F$region,
            source='focal',R5_group=NA_character_))
if(include_all_historical_positions&&file.exists(legacy)){
 zz<-read.csv(legacy,stringsAsFactors=FALSE)
 key<-rbind(key,data.frame(SNP_ID=paste0(zz$CHROM,'_',zz$POS),
    region=zz$historical_region,source='historical_exact',R5_group=NA_character_))
}
if(file.exists(legacy_status)){
 h<-read.csv(legacy_status,stringsAsFactors=FALSE)
 ix<-match(key$SNP_ID,h$SNP_ID)
 key$R5_group[is.na(key$R5_group)]<-as.character(h$R5_group_exact[ix][is.na(key$R5_group)])
}
key$R5_group[key$R5_group%in%c('NA','','NaN')]<-NA_character_
key<-key[!duplicated(paste(key$SNP_ID,key$region)),,drop=FALSE]

## Use all 13 aligned pools to reveal both G0 observations, including the
## bottlenecked pool excluded from the pre-specified selection baseline.
p<-read_drift_inputs(pool_file,selection_file,exclude_G0=character())
S<-p$sel;pop<-p$pop
if(!all(c('ID','rs')%in%names(S)))stop('Unexpected selection SNP manifest')
if(!identical(rownames(p$k_alt),as.character(S$rs)))
 warning('Check original SNP names; row-order was validated by read_drift_inputs')
full<-read.csv(full_snp_file,stringsAsFactors=FALSE)
if(anyDuplicated(full$ID))stop('Duplicated full SNP IDs')
surv<-if(file.exists(survivor_snp_file))read.csv(survivor_snp_file,stringsAsFactors=FALSE)else NULL
cl<-read.csv(full_cluster_file,stringsAsFactors=FALSE)
cl$group<-as.character(cl$group)
ph<-if(file.exists(phase_file))read.csv(phase_file,stringsAsFactors=FALSE)else NULL
svcl<-if(file.exists(survivor_cluster_file))read.csv(survivor_cluster_file,stringsAsFactors=FALSE)else NULL
if(!is.null(svcl))svcl$group<-as.character(svcl$group)
if(file.exists(status_file)){
 stat<-read.csv(status_file,stringsAsFactors=FALSE)
 if(!all(c('replicate','status')%in%names(stat)))stop('Malformed GA status CSV')
 pop$GA_status<-stat$status[match(pop$replicate,stat$replicate)]
 pop$GA_status[pop$generations==0]<-'baseline'
}else pop$GA_status<-ifelse(pop$generations==0,'baseline','not_provided')

u<-unique(key$SNP_ID);i<-match(u,S$ID)
allstatus<-data.frame(SNP_ID=u,Poolseq_in_corrected_08_and_09=!is.na(i))
write.csv(allstatus,file.path(out_dir,'28_SNP_coverage_status.csv'),row.names=FALSE)
selids<-which(!is.na(i))
if(length(selids)==0L)warning('No nominated SNP occurs in the corrected pooled GWAS-SNP universe')
traj<-list();repdelta<-list()
for(j in selids){
 sid<-u[j];ix<-i[j];depth<-as.numeric(p$depth[ix,]);alt<-as.numeric(p$k_alt[ix,])
 good<-is.finite(depth)&depth>0&is.finite(alt)&alt>=0&alt<=depth
 x<-pop
 x$SNP_ID<-sid;x$depth<-depth;x$ALT_reads<-alt
 x$ALT_frequency<-ifelse(good,alt/depth,NA_real_)
 x$selected_baseline_G0<-x$generations==0&x$sample_id!=exclude_G0_from_selection
 x$included_in_survivor_analysis<-x$generations==0|
      tolower(x$GA_status)=='survived'
 traj[[length(traj)+1L]]<-x
 bg<-x$ALT_frequency[x$generations==0&x$sample_id!=exclude_G0_from_selection]
 bg<-if(length(bg)==1L)bg else NA_real_
 for(r in unique(x$replicate[x$generations>0])){
    z<-x[x$replicate==r&x$generations>0&is.finite(x$ALT_frequency),,drop=FALSE]
    if(!nrow(z))next
    last<-z[which.max(z$generations),,drop=FALSE]
    b<-if(nrow(z)>=2&&length(unique(z$generations))>=2)
      unname(coef(lm(ALT_frequency~generations,data=z))[2])else NA_real_
    repdelta[[length(repdelta)+1L]]<-data.frame(SNP_ID=sid,replicate=r,
       GA_status=as.character(z$GA_status[1]),n_timepoints=nrow(z),
       final_generation=last$generations[1],G0_frequency_A00=bg,
       final_frequency=last$ALT_frequency[1],delta_from_A00=last$ALT_frequency[1]-bg,
       within_replicate_slope_per_generation=b,stringsAsFactors=FALSE)
 }
}
if(length(traj))write.csv(do.call(rbind,traj),file.path(out_dir,'28_ALT_read_counts_all_pools.csv'),row.names=FALSE)
if(length(repdelta))write.csv(do.call(rbind,repdelta),file.path(out_dir,'28_Per_GA_trajectory_summary.csv'),row.names=FALSE)

## Summarize EXACT sites, explicitly flagging whether the drift-q is defined.
site<-data.frame(SNP_ID=u,stringsAsFactors=FALSE)
m<-match(u,full$ID)
for(z in intersect(c('p0_hat','common_G0','slope_freq','SE_freq','p_drift','p_drift_BH_common'),names(full)))
 site[[paste0('all_GA_',z)]]<-full[[z]][m]
site$in_full_SNP_test<-!is.na(m)
if(!is.null(surv)){
 ix<-match(u,surv$ID)
 for(z in intersect(c('p0_hat','common_G0','slope_freq','SE_freq','p_drift','p_drift_BH_common'),names(surv)))
  site[[paste0('retained_GA_',z)]]<-surv[[z]][ix]
 site$in_retained_SNP_test<-!is.na(ix)
}
site<-merge(key,site,by='SNP_ID',all.x=TRUE,sort=FALSE)
site$selection_status<-ifelse(is.na(site$all_GA_p0_hat),
  'not_in_tested_SNP_universe',ifelse(!site$all_GA_common_G0,
  'low_or_high_ancestral_frequency_q_NOT_DEFINED',
  ifelse(is.finite(site$all_GA_p_drift_BH_common),
    'common_frequency_test_available','tested_but_no_finite_q')))
if(!is.null(ph)&&all(c('rs','group','phase')%in%names(ph))){
 k<-match(site$SNP_ID,ph$rs)
 site$R5_group_from_phase<-as.character(ph$group[k])
 site$R5_phase<-as.numeric(ph$phase[k])
 site$phase_aligned_full_slope<-site$R5_phase*site$all_GA_slope_freq
 if('retained_GA_slope_freq'%in%names(site))
  site$phase_aligned_retained_slope<-site$R5_phase*site$retained_GA_slope_freq
}
write.csv(site,file.path(out_dir,'28_Exact_variant_GWAS_and_selection_summary.csv'),row.names=FALSE)

## Compare R5 groups with FULL and retained-GA evidence, regardless of whether
## the focal historical SNP passes the common-frequency restriction.
groups<-unique(key[,c('region','R5_group')]);groups<-groups[!is.na(groups$R5_group),,drop=FALSE]
ind<-match(groups$R5_group,cl$group)
for(z in intersect(c('CHROM','start','end','n_slope_common','min_q_common',
                       'n_selected_common_q01','slope_median_common',
                       'effect_beta','q_beta','effect_SB','q_SB','effect_area.F','q_area.F'),names(cl)))
 groups[[paste0('full_',z)]]<-cl[[z]][ind]
if(!is.null(svcl)){
 ind<-match(groups$R5_group,svcl$group)
 for(z in intersect(c('n_survivor_common','n_selected_survivor','any_selected_survivor',
                     'slope_survivor_median_common'),names(svcl)))
  groups[[z]]<-svcl[[z]][ind]
}
for(trait in c('beta','SB','area.F')){
 eg<-paste0('full_effect_',trait)
 if(eg%in%names(groups)&&'full_slope_median_common'%in%names(groups)){
  groups[[paste0('direction_matches_',trait,'_all_GA')]]<-
    ifelse(is.finite(groups[[eg]])&is.finite(groups$full_slope_median_common),
      groups[[eg]]*groups$full_slope_median_common>0,NA)
 }
 if(eg%in%names(groups)&&'slope_survivor_median_common'%in%names(groups)){
  groups[[paste0('direction_matches_',trait,'_retained_GA')]]<-
    ifelse(is.finite(groups[[eg]])&is.finite(groups$slope_survivor_median_common),
      groups[[eg]]*groups$slope_survivor_median_common>0,NA)
 }
}
write.csv(groups,file.path(out_dir,'28_R5_region_full_vs_retained_GA_overlap.csv'),row.names=FALSE)
cat('Script 28 complete: Pool-seq allele counts and available drift-aware q values; ',
    'low-frequency SNPs are descriptive only.\n',sep='')

## Plot exact focal variants and their time courses; include the excluded G0
## pool as a distinct datum, not as a fitted baseline.
if(length(traj)){
 tr<-do.call(rbind,traj)
 focal_ids<-unique(paste0(F$CHROM,'_',F$POS))
 pdf(file.path(out_dir,'28_Focal_variant_ALT_frequency_trajectories.pdf'),width=8.5,height=5.5)
 for(sid in focal_ids){
   x<-tr[tr$SNP_ID==sid&is.finite(tr$ALT_frequency),,drop=FALSE]
   if(!nrow(x))next
   reps<-unique(x$replicate[x$generations>0]);cols<-setNames(rainbow(length(reps)),reps)
   par(mar=c(4,4,3,1))
   plot(NA,xlim=c(0,max(x$generations)),ylim=c(0,1),
        xlab='Generations in salt experiment',ylab='GWAS ALT frequency',
        main=paste('Exact variant',sid))
   for(r in reps){
     z<-x[x$replicate==r&x$generations>0,,drop=FALSE]
     z<-z[order(z$generations),,drop=FALSE]
     lines(z$generations,z$ALT_frequency,type='b',pch=16,col=cols[r],lwd=1.8)
   }
   G0<-x[x$generations==0,,drop=FALSE]
   if(nrow(G0)){
     points(G0$generations,G0$ALT_frequency,pch=ifelse(G0$sample_id=='A00',16,1),
            cex=1.4,col='black')
     aa<-G0$ALT_frequency[G0$sample_id=='A00']
     if(length(aa)==1&&is.finite(aa))abline(h=aa,lty=2,col='grey70')
   }
   legend('topleft',legend=c(reps,'G0 A00','other G0'),
      col=c(unname(cols), 'black','black'),pch=c(rep(16,length(reps)),16,1),
      bty='n',cex=.8)
 }
 dev.off()
}
