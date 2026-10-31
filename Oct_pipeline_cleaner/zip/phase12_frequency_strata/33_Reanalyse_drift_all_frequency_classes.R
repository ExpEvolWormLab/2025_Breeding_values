## 33_Reanalyse_drift_all_frequency_classes.R
## Genome-wide reanalysis of existing (founder-filtered) GWAS SNP universe.
## Null: shared ancestral p0, independent neutral diploid WF populations,
## binomial Pool-seq reads. A00 baseline is conditioned on, not re-sampled.
## Normal p-values are APPROXIMATE, particularly at very low MAF.
## Script 34b tests their calibration; script 34 tests nominated sites by MC.

suppressPackageStartupMessages(library(parallel))
helper <- if(file.exists('09f_drift_functions.R')) '09f_drift_functions.R' else
  'drift_nullfit_update/09f_drift_functions.R'
source(helper)
source('33_Frequency_drift_helpers.R')

pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
orientation_file <- 'selection_results/09_Selection_SNP_results.csv'
status_file <- '20_GA_population_status.csv'
previous_full_file <- 'selection_results/drift_oneG0_nullfit/09f_SNP_drift_Ne1000.csv'
previous_survivor_file <- 'selection_results/survivors_only_Ne1000/20_SNP_drift_Ne1000_survivors.csv'
out_dir <- 'selection_results/frequency_strata'
Ne <- 1000L
excluded_G0 <- 'A6G140'
chunk_size <- 400L
mc_cores <- 6L
## Pre-specified MAF classes; never define them based on the resulting p-values.
if(.Platform$OS.type=='windows')mc_cores<-1L
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

if(!file.exists(pool_file)||!file.exists(orientation_file)||!file.exists(status_file))
  stop('Required files are absent. Check pool_file, orientation_file, status_file.')
p <- read_drift_inputs(pool_file,orientation_file,exclude_G0=excluded_G0)
if(sum(p$pop$generations==0)!=1L ||
   p$pop$sample_id[p$pop$generations==0]!='A00')
  stop('Expected verified A00 baseline and excluded A6G140')
if(!identical(colnames(p$k_alt),p$pop$sample_id))stop('Pool column order mismatch')
status<-read.csv(status_file,stringsAsFactors=FALSE)
if(!all(c('replicate','status')%in%names(status)))stop('Malformed status table')
status$status<-tolower(trimws(as.character(status$status)))
if(!all(status$status%in%c('survived','extinct','unknown')))stop('Unknown status category')
reps<-as.character(unique(p$pop$replicate[p$pop$generations>0]))
if(!setequal(reps,as.character(status$replicate)))stop('Status table and GA replicates differ')
kept_reps<-as.character(status$replicate[status$status=='survived'])
if(length(kept_reps)<2L)stop('Survivor analysis needs >=2 confirmed surviving replicates')
write.csv(status,file.path(out_dir,'33_GA_status_used.csv'),row.names=FALSE)

all_sets <- list(all_GA=rep(TRUE,nrow(p$pop)),
  retained_GA=p$pop$generations==0|as.character(p$pop$replicate)%in%kept_reps)
results <- list()
for(scenario in names(all_sets)){
  keep<-all_sets[[scenario]]
  pop<-p$pop[keep,,drop=FALSE]
  K<-p$k_alt[,keep,drop=FALSE]
  D<-p$depth[,keep,drop=FALSE]
  cat('Scenario ',scenario,': ',nrow(K),' SNPs, ',ncol(K),' pools\n',sep='')
  ix <- seq_len(nrow(K))
  chunks <- split(ix,ceiling(ix/chunk_size))
  a <- parallel::mclapply(chunks,function(ii){
    t(vapply(ii,function(i){
      z <- conditional_neutral_info(K[i,],D[i,],pop,Ne=Ne)
      c(p0=z$p0,p0_lo=z$p0_lo,p0_hi=z$p0_hi,
        ALT_reads_G0=z$baseline_ALT_reads,depth_G0=z$baseline_depth,
        slope=z$slope,expected_neutral_slope=z$mean_slope,
        conditional_SE=z$conditional_SE,conditional_Z=z$conditional_z,
        p_normal=z$conditional_p_normal,n_pools=z$n,
        n_reps=z$n_evolved_reps)
    },FUN.VALUE=setNames(rep(0,12),c('p0','p0_lo','p0_hi','ALT_reads_G0',
          'depth_G0','slope','expected_neutral_slope','conditional_SE',
          'conditional_Z','p_normal','n_pools','n_reps'))))
  },mc.cores=mc_cores,mc.preschedule=TRUE)
  a<-do.call(rbind,a)
  if(nrow(a)!=nrow(p$sel))stop('SNP order changed')
  res<-data.frame(SNP_ID=as.character(p$sel$rs),ID=as.character(p$sel$ID),
    CHROM=as.character(p$sel$CHROM),POS=as.integer(p$sel$POS),
    scenario=scenario,Ne=Ne,
    p0_ALT_posterior=a[,'p0'],
    p0_ALT_CI_lower=a[,'p0_lo'],p0_ALT_CI_upper=a[,'p0_hi'],
    G0_ALT_reads=a[,'ALT_reads_G0'],G0_depth=a[,'depth_G0'],
    MAF_G0=pmin(a[,'p0'],1-a[,'p0']),
    slope_freq_per_generation=a[,'slope'],
    neutral_expected_slope=a[,'expected_neutral_slope'],
    conditional_SE=a[,'conditional_SE'],
    z_normal_approx=a[,'conditional_Z'],
    p_normal_approx=a[,'p_normal'],
    n_pools=a[,'n_pools'],n_evolved_reps=a[,'n_reps'],
    stringsAsFactors=FALSE)
  res$MAF_class<-as.character(frequency_class(res$p0_ALT_posterior))
  res$G0_ALT_not_detected<-is.finite(res$G0_ALT_reads)&res$G0_ALT_reads==0
  ## Intentionally NO q-value from the uncalibrated normal approximation.
  ## This output is a hypothesis screen, not a new genome-wide selection map.
  fn<-file.path(out_dir,paste0('33_All_SNPs_',scenario,'_Ne',Ne,'.csv'))
  write.csv(res,fn,row.names=FALSE)
  write.csv(pop,file.path(out_dir,paste0('33_Samples_',scenario,'.csv')),row.names=FALSE)
  results[[scenario]]<-res
  cat('Wrote ',fn,'\n',sep='')
  oldfile<-if(scenario=='all_GA')previous_full_file else previous_survivor_file
  if(file.exists(oldfile)){
    old<-read.csv(oldfile,stringsAsFactors=FALSE)
    o<-match(res$SNP_ID,old$rs)
    if(anyNA(o))stop('Original 09f/20 SNPs differ from script 33')
    dd<-abs(res$slope_freq_per_generation-old$slope_freq[o])
    maxdiff<-if(any(is.finite(dd)))max(dd[is.finite(dd)])else NA_real_
    write.csv(data.frame(scenario=scenario,compared=sum(is.finite(dd)),
       max_abs_slope_difference=maxdiff),
       file.path(out_dir,paste0('33_QC_09f_slope_match_',scenario,'.csv')),row.names=FALSE)
    if(is.finite(maxdiff) && maxdiff>1e-7)
      stop('Script-09f slope mismatch: inspect corrected sample orientation, helper and Ne')
  }
}

summ <- list()
for(s in names(results)){
  z<-results[[s]]
  for(b in unique(c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50'))){
    x<-z[z$MAF_class==b & is.finite(z$p_normal_approx),,drop=FALSE]
    summ[[length(summ)+1L]]<-data.frame(scenario=s,MAF_class=b,
      n_tested=nrow(x),n_G0_zero_ALT=sum(x$G0_ALT_not_detected,na.rm=TRUE),
      median_G0_depth=if(nrow(x))median(x$G0_depth,na.rm=TRUE)else NA_real_,
      fraction_abs_Z_gt_2=if(nrow(x))mean(abs(x$z_normal_approx)>2)else NA_real_,
      n_p_approx_below_0.01=sum(x$p_normal_approx<.01),
      minimum_p_approx=if(nrow(x))min(x$p_normal_approx)else NA_real_)
  }
}
write.csv(do.call(rbind,summ),file.path(out_dir,'33_Frequency_class_summary.csv'),row.names=FALSE)

full<-results[['all_GA']];sr<-results[['retained_GA']]
if(!identical(full$SNP_ID,sr$SNP_ID))stop('Full and retained SNP orders differ')
paired<-data.frame(SNP_ID=full$SNP_ID,CHROM=full$CHROM,POS=full$POS,
  MAF_class=full$MAF_class,G0_ALT_reads=full$G0_ALT_reads,G0_depth=full$G0_depth,
  full_slope=full$slope_freq_per_generation,
  retained_slope=sr$slope_freq_per_generation,
  full_p_normal_approx=full$p_normal_approx,
  retained_p_normal_approx=sr$p_normal_approx)
write.csv(paired,file.path(out_dir,'33_Matched_full_vs_retained_all_frequencies.csv'),row.names=FALSE)
cat('Script 33 completed. Low-frequency approximate p-values are NOT calibrated discoveries.\n')
