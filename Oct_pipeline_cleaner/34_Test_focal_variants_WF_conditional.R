## 34_Test_focal_variants_WF_conditional.R
## Neutral Wright-Fisher posterior-predictive tests of a PREDECLARED set
## of low/common-frequency nominated loci (not a new genome-wide FDR map).
## Fixed observed A00 counts; integrate true p0 uncertainty using Jeffreys
## Beta(.5,.5), WF drift and read sampling; optional all and retained GA sets.

helper <- if(file.exists('09f_drift_functions.R'))'09f_drift_functions.R' else
 'drift_nullfit_update/09f_drift_functions.R'
source(helper)
source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')

pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
orientation_file <- 'selection_results/09_Selection_SNP_results.csv'
status_file <- 'Oct_pipeline_cleaner/20_GA_population_status.csv'
raw_sync_file <- 'selection_results/frequency_strata/34a_Focal_raw_SYNC_ALT_counts.csv'
## Existing script-28b output is also accepted when 34a has not been run.
raw_sync_fallback <- 'historical_QTL_R5/deep_dive/evolution/raw_sync/28b_Exact_raw_SYNC_ALT_trajectories.csv'
focal_file <- 'Oct_pipeline_cleaner/26_Focal_variants.csv'
new_focal_file <- '30_Focal_variants_extended.csv'
old_tag_file <- 'historical_QTL_R5/deep_dive/26_R5_model_tag_SNPs.csv'
new_tag_file <- 'historical_QTL_R5/deep_dive_extended/26_R5_model_tag_SNPs.csv'
out_dir <- 'selection_results/frequency_strata'
excluded_G0 <- 'A6G140'
Ne_values <- c(500L,1000L,2000L)
n_sim_main <- 6000L
n_sim_sensitivity <- 2000L
## If the nominated set becomes >100 SNPs, reduce targets or start at 500 draws.
max_n_targets <- 150L
set.seed(20261029)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

req<-c(pool_file,orientation_file,status_file,focal_file)
if(any(!file.exists(req)))stop('Missing prerequisites: ',paste(req[!file.exists(req)],collapse=', '))
p<-read_drift_inputs(pool_file,orientation_file,exclude_G0=excluded_G0)
if(sum(p$pop$generations==0)!=1L ||
   p$pop$sample_id[p$pop$generations==0]!='A00')stop('Unexpected G0 baseline')
status<-read.csv(status_file,stringsAsFactors=FALSE)
status$status<-tolower(trimws(status$status))
if(!all(status$status%in%c('survived','extinct','unknown')))stop('Unknown GA statuses')
if(!setequal(as.character(status$replicate),
             as.character(unique(p$pop$replicate[p$pop$generations>0]))))
  stop('GA status manifest does not match Pool-seq populations')
keep_reps<-as.character(status$replicate[status$status=='survived'])
if(length(keep_reps)<2)stop('Need >=2 confirmed surviving GA replicates')

focal<-read.csv(focal_file,stringsAsFactors=FALSE)
target<-data.frame(SNP_ID=paste(focal$CHROM,focal$POS,sep='_'),
  region=as.character(focal$region),source='focal',stringsAsFactors=FALSE)
if(file.exists(new_focal_file)){
  a<-read.csv(new_focal_file,stringsAsFactors=FALSE)
  target<-rbind(target,data.frame(SNP_ID=paste(a$CHROM,a$POS,sep='_'),
    region=as.character(a$region),source='extended_focal',stringsAsFactors=FALSE))
}
for(f in c(old_tag_file,new_tag_file))if(file.exists(f)){
  a<-read.csv(f,stringsAsFactors=FALSE)
  target<-rbind(target,data.frame(SNP_ID=as.character(a$tag_SNP_ID),
    region=as.character(a$region),source='R5_tag',stringsAsFactors=FALSE))
}
## Optional user-edited predeclared targets, columns SNP_ID and region.
extra_file <- '34_Additional_target_SNPs.csv'
if(file.exists(extra_file)){
  a<-read.csv(extra_file,stringsAsFactors=FALSE)
  if(!all(c('SNP_ID','region')%in%names(a)))stop('Extra targets need SNP_ID and region')
  target<-rbind(target,data.frame(SNP_ID=a$SNP_ID,region=a$region,
                                      source='additional_predeclared',stringsAsFactors=FALSE))
}
target<-target[!is.na(target$SNP_ID)&grepl('^[IVX]+_[0-9]+$',target$SNP_ID),,drop=FALSE]
registry<-data.frame(SNP_ID=unique(target$SNP_ID),stringsAsFactors=FALSE)
registry$regions<-vapply(registry$SNP_ID,function(s)
  paste(unique(target$region[target$SNP_ID==s]),collapse=';'),character(1))
registry$sources<-vapply(registry$SNP_ID,function(s)
  paste(unique(target$source[target$SNP_ID==s]),collapse=';'),character(1))
if(nrow(registry)>max_n_targets)stop('More than ',max_n_targets,' candidate SNPs; curate targets first')
write.csv(registry,file.path(out_dir,'34_Predeclared_target_registry.csv'),row.names=FALSE)

raw <- NULL
for(f in c(raw_sync_file,raw_sync_fallback))if(file.exists(f)){
  raw<-read.csv(f,stringsAsFactors=FALSE)
  cat('Using supplemental raw SYNC counts from ',f,'\n');break
}
if(!is.null(raw)){
  if(!all(c('SNP_ID','sample_id','ALT_reads','depth','REF_agrees')%in%names(raw)))
    stop('Raw SYNC counts require SNP_ID/sample_id/ALT_reads/depth/REF_agrees')
  raw$SNP_ID<-as.character(raw$SNP_ID)
  if(anyDuplicated(paste(raw$SNP_ID,raw$sample_id)))
    stop('Duplicate raw SYNC SNP/sample rows')
}

## Cross-check focal loci occurring in BOTH sources for the old sample-column
## alignment bug. Small differences can arise through count reconstruction;
## substantial mismatch needs a direct sample/allele audit before inference.
if(!is.null(raw)){
  common <- intersect(unique(raw$SNP_ID),as.character(p$sel$rs))
  checks <- list()
  for(sid in common){
    a <- raw[raw$SNP_ID==sid,,drop=FALSE]
    j <- match(p$pop$sample_id,a$sample_id)
    if(anyNA(j))next
    ii <- match(sid,as.character(p$sel$rs))
    ck <- as.numeric(p$k_alt[ii,]);cd <- as.numeric(p$depth[ii,])
    rk <- as.numeric(a$ALT_reads[j]);rd <- as.numeric(a$depth[j])
    usable <- is.finite(ck)&is.finite(cd)&is.finite(rk)&is.finite(rd)
    checks[[length(checks)+1L]] <- data.frame(SNP_ID=sid,n_pools_compared=sum(usable),
       max_abs_ALT_read_difference=if(any(usable))max(abs(ck[usable]-rk[usable]))else NA_real_,
       max_abs_depth_difference=if(any(usable))max(abs(cd[usable]-rd[usable]))else NA_real_)
  }
  if(length(checks)){
    check_table<-do.call(rbind,checks)
    write.csv(check_table,file.path(out_dir,'34_Raw_vs_prepared_count_alignment_QC.csv'),row.names=FALSE)
    if(any(check_table$max_abs_ALT_read_difference>3 |
           check_table$max_abs_depth_difference>3,na.rm=TRUE))
      warning('Raw-versus-prepared counts differ at some sites. Inspect ',
              '34_Raw_vs_prepared_count_alignment_QC.csv before interpreting results')
  }
}

ix <- match(registry$SNP_ID,as.character(p$sel$rs))
counts <- list();diagnostics <- list()
for(i in seq_len(nrow(registry))){
  sid<-registry$SNP_ID[i]
  if(!is.na(ix[i])){
    k<-as.numeric(p$k_alt[ix[i],]);d<-as.numeric(p$depth[ix[i],]);
    src<-'script08_founder_filtered_GWAS_ALT'
  }else{
    if(is.null(raw)||!sid%in%raw$SNP_ID){
      diagnostics[[i]]<-data.frame(SNP_ID=sid,status='not_in_08_or_raw_sync',
        source='missing',n_pools=0,stringsAsFactors=FALSE)
      next
    }
    z<-raw[raw$SNP_ID==sid,,drop=FALSE]
    if(any(is.na(z$REF_agrees))||!all(z$REF_agrees)){
      diagnostics[[i]]<-data.frame(SNP_ID=sid,status='REF_mismatch',
        source='raw_SYNC',n_pools=nrow(z),stringsAsFactors=FALSE)
      next
    }
    zz<-match(p$pop$sample_id,z$sample_id)
    if(anyNA(zz)){
      diagnostics[[i]]<-data.frame(SNP_ID=sid,status='incomplete_raw_SYNC_samples',
        source='raw_SYNC',n_pools=nrow(z),stringsAsFactors=FALSE)
      next
    }
    k<-as.numeric(z$ALT_reads[zz]);d<-as.numeric(z$depth[zz]);
    src<-'raw_SYNC_VCF_ALT'
  }
  if(any(is.finite(k)&is.finite(d)&(k<0|k>d))){
    diagnostics[[i]]<-data.frame(SNP_ID=sid,status='invalid_counts',source=src,
      n_pools=length(d),stringsAsFactors=FALSE)
    next
  }
  diagnostics[[i]]<-data.frame(SNP_ID=sid,status='counts_available',source=src,
    n_pools=sum(is.finite(d)&d>0),stringsAsFactors=FALSE)
  counts[[sid]]<-data.frame(SNP_ID=sid,sample_id=p$pop$sample_id,
    generations=p$pop$generations,replicate=p$pop$replicate,
    ALT_reads=k,depth=d,source=src,stringsAsFactors=FALSE)
}
qc<-do.call(rbind,diagnostics)
write.csv(qc,file.path(out_dir,'34_Target_count_coverage_QC.csv'),row.names=FALSE)
if(!length(counts))stop('No nominated variants with usable counts; run 34a')
all_counts<-do.call(rbind,counts)
write.csv(all_counts,file.path(out_dir,'34_Targeted_ALT_counts_in_analysis_pools.csv'),row.names=FALSE)

fit_rows<-list();bincount<-0L
for(scenario in c('all_GA','retained_GA')){
  mask<-if(scenario=='all_GA')rep(TRUE,nrow(p$pop)) else
       p$pop$generations==0|as.character(p$pop$replicate)%in%keep_reps
  pop<-p$pop[mask,,drop=FALSE]
  for(Ne in Ne_values){
    B<-if(Ne==1000L)n_sim_main else n_sim_sensitivity
    for(sid in names(counts)){
      zz<-counts[[sid]][mask,,drop=FALSE]
      info<-conditional_neutral_info(zz$ALT_reads,zz$depth,pop,Ne=Ne)
      if(!info$ok){
        fit_rows[[length(fit_rows)+1L]]<-data.frame(SNP_ID=sid,
          scenario=scenario,Ne=Ne,status='not_testable',p0_ALT=info$p0,
          MAF_class=NA_character_,G0_ALT_reads=info$baseline_ALT_reads,
          G0_depth=info$baseline_depth,p0_lo=info$p0_lo,p0_hi=info$p0_hi,
          slope_freq=info$slope,mean_null_slope=info$mean_slope,
          centered_slope=NA_real_,conditional_SE=NA_real_,
          p_normal_approx=NA_real_,p_WF_conditional=NA_real_,
          n_sim=0L,extreme_draws=NA_integer_,
          posterior_predictive_normal_fp01=NA_real_,
          stringsAsFactors=FALSE)
        next
      }
      test<-mc_neutral_test(info,zz$ALT_reads,zz$depth,pop,Ne=Ne,n_sim=B)
      fit_rows[[length(fit_rows)+1L]]<-data.frame(SNP_ID=sid,
        scenario=scenario,Ne=Ne,status='tested',p0_ALT=info$p0,
        MAF_class=as.character(frequency_class(info$p0)),
        G0_ALT_reads=info$baseline_ALT_reads,G0_depth=info$baseline_depth,
        p0_lo=info$p0_lo,p0_hi=info$p0_hi,
        slope_freq=info$slope,mean_null_slope=info$mean_slope,
        centered_slope=test$obs,conditional_SE=info$conditional_SE,
        p_normal_approx=info$conditional_p_normal,
        p_WF_conditional=test$p,n_sim=test$B,
        extreme_draws=test$exceed,
        posterior_predictive_normal_fp01=test$normal_fp01,
        stringsAsFactors=FALSE)
      cat('Tested ',sid,' ',scenario,' Ne=',Ne,' p_MC=',format(test$p,digits=3),'\n',sep='')
    }
  }
}
res<-do.call(rbind,fit_rows)
res<-merge(res,registry,by='SNP_ID',all.x=TRUE,sort=FALSE)
res$source_counts<-qc$source[match(res$SNP_ID,qc$SNP_ID)]
res$p_BH_within_nominated_family_EXPLORATORY<-NA_real_
for(m in unique(paste(res$scenario,res$Ne,sep=':'))){
  ii<-which(paste(res$scenario,res$Ne,sep=':')==m&is.finite(res$p_WF_conditional))
  if(length(ii))res$p_BH_within_nominated_family_EXPLORATORY[ii]<-
    p.adjust(res$p_WF_conditional[ii],method='BH')
}
write.csv(res,file.path(out_dir,'34_Focal_WF_neutral_conditional_tests.csv'),row.names=FALSE)
write.csv(data.frame(scenario=c('all_GA','retained_GA'),
   included_GA=c(paste(unique(as.character(p$pop$replicate[p$pop$generations>0])),collapse=';'),
    paste(keep_reps,collapse=';')),baseline='A00',excluded_G0=excluded_G0),
   file.path(out_dir,'34_Analysis_sample_QC.csv'),row.names=FALSE)
cat('Script 34 completed. These are finite-resolution nominated-locus neutral tests,\n',
    'not genome-wide-calibrated selection discoveries.\n',sep='')
