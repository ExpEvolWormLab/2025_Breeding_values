## 43_Unified_config_and_helpers.R
## Frequency-INCLUSIVE, cluster-level, evolution-phase-aligned analysis.
## 43-46 scripts use the same immutable settings. Copy entire package to project root.
## Main: GA1/GA2/GA4. GA3 confirmed extinct: sensitivity only.

uf_config <- list(
  prefix = 'selection_results/unified_frequency_inclusive_R5',
  pool_file = 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData',
  orientation_file = 'selection_results/09_Selection_SNP_results.csv',
  phase_file = 'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv',
  snp_gwas_file = 'txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt',
  r5_gwas_file = 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
  weights_file = 'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Gradient_weight_QC.csv',
  status_file = '20_GA_population_status.csv',
  existing_snp_primary = 'selection_results/frequency_strata/33_All_SNPs_retained_GA_Ne1000.csv',
  existing_snp_supp = 'selection_results/frequency_strata/33_All_SNPs_all_GA_Ne1000.csv',
  raw_focal_file = 'selection_results/frequency_strata/34a_Focal_raw_SYNC_ALT_counts.csv',
  previous_focal_test_file = 'selection_results/frequency_strata/34_Focal_WF_neutral_conditional_tests.csv',
  Ne = 1000L,
  Ne_sensitivity = c(500L, 2000L),
  ## Coverage is a QUALITY threshold, never a minimum allele-frequency threshold.
  min_A00_depth = 30L,
  min_evolved_replicates = 2L,
  min_pools = 6L,
  block_bp = 2000000L,
  ## Main conditional null: A00 read-count posterior Beta(k+1/2,d-k+1/2)
  ## and diploid WF forward simulation. Adaptivity uses simultaneous
  ## Clopper-Pearson bounds; NOT naive optional-stop Monte Carlo p-values.
  mc_stages = c(3000L, 15000L, 75000L, 300000L, 1200000L, 3000000L),
  mc_chunk = 20000L,
  mc_family_error = 0.01,  # at least 99% simultaneous upper-bound coverage
  max_cores = 4L,
  batch_clusters = 20L,
  seed = 20261029L,
  bootstrap_B = 4000L,
  permutation_B = 3000L,
  run_all_GA_supp = TRUE,
  run_Ne_sensitivity = TRUE,
  Ne_sensitivity_B = 150000L
)

uf_check <- function(files){
  bad <- files[!file.exists(files)]
  if(length(bad))stop('MISSING INPUTS:\n',paste(' -',bad,collapse='\n'),call.=FALSE)
}
uf_read <- function(path, ...)read.csv(path,stringsAsFactors=FALSE,check.names=FALSE,...)
uf_write <- function(x,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);
  write.csv(x,path,row.names=FALSE,na='')}
uf_finite_median <- function(x){x<-x[is.finite(x)];if(length(x))median(x) else NA_real_}
uf_safe_cor <- function(a,b){z<-is.finite(a)&is.finite(b);
  if(sum(z)<8L || sd(a[z])==0 || sd(b[z])==0)return(NA_real_);
  cor(a[z],b[z])}
uf_status_and_samples <- function(p,cfg=uf_config){
  s<-uf_read(cfg$status_file)
  if(!all(c('replicate','population','status')%in%names(s)))stop('Status needs replicate,population,status')
  s$status<-tolower(trimws(s$status));s$population<-toupper(trimws(s$population))
  if(anyDuplicated(s$replicate)||any(!s$status%in%c('survived','extinct')))
    stop('Set explicit survived/extinct statuses; unknown is not accepted for the final run')
  expect<-c(GA1='survived',GA2='survived',GA3='extinct',GA4='survived')
  chk<-s$status[match(names(expect),s$population)]
  if(anyNA(chk)||any(chk!=expect))stop('Final manuscript requires GA1/GA2/GA4 survived, GA3 extinct')
  evolved<-unique(as.character(p$pop$replicate[p$pop$generations>0]))
  if(!setequal(evolved,as.character(s$replicate)))stop('Status IDs mismatch evolving Pool-seq replicates')
  if(!identical(colnames(p$k_alt),as.character(p$pop$sample_id)))stop('Corrected pool sample alignment failed')
  if(sum(p$pop$generations==0)!=1L || as.character(p$pop$sample_id[p$pop$generations==0])!='A00')
    stop('Require precisely the A00 baseline and exclude A6G140')
  keep_reps<-as.character(s$replicate[s$status=='survived'])
  keep<-p$pop$generations==0|as.character(p$pop$replicate)%in%keep_reps
  list(status=s,primary=keep,all=rep(TRUE,nrow(p$pop)))
}
uf_stop_duplicate <- function(x,what){if(anyDuplicated(x))stop('Duplicate ',what)}
uf_int <- function(x){as.integer(as.character(x))}
uf_pos_from_rs <- function(id)as.numeric(sub('^.*_','',as.character(id)))
uf_mc_bounds <- function(exceed,B,alpha){
  if(!is.finite(B)||B<1L)return(c(lower=NA_real_,upper=NA_real_))
  lower<-if(exceed==0)0 else qbeta(alpha,exceed,B-exceed+1L)
  upper<-if(exceed==B)1 else qbeta(1-alpha,exceed+1L,B-exceed)
  c(lower=lower,upper=upper)
}
uf_BH_BY <- function(x, method='BY'){
  ## p.adjust() automatically counts all nonmissing hypotheses.
  ## Caller must carry one entry per testable R5 cluster and zero other rows.
  p.adjust(x,method=method)
}
uf_block_boot <- function(dat, pred1, pred2, B=4000L, seed=1L){
  if(!all(c('slope','block',pred1,pred2)%in%names(dat)))stop('Bootstrap columns missing')
  ok<-is.finite(dat$slope)&is.finite(dat[[pred1]])&is.finite(dat[[pred2]])&!is.na(dat$block)
  z<-dat[ok,,drop=FALSE]
  blocks<-split(seq_len(nrow(z)),z$block)
  if(length(blocks)<8L)return(data.frame(delta_MSE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,
    bootstrap_sign_p=NA_real_,n_blocks=length(blocks),n_clusters=nrow(z)))
  set.seed(seed)
  obs<-mean((z$slope-z[[pred2]])^2-(z$slope-z[[pred1]])^2)
  draws<-replicate(B,{
    take<-unlist(blocks[sample.int(length(blocks),length(blocks),replace=TRUE)],use.names=FALSE)
    mean((z$slope[take]-z[[pred2]][take])^2-(z$slope[take]-z[[pred1]][take])^2)
  })
  ## Two-sided percentile-bootstrap sign check is exploratory, not a randomization p-value.
  data.frame(delta_MSE=obs,CI_low=unname(quantile(draws,.025)),
    CI_high=unname(quantile(draws,.975)),
    bootstrap_sign_p=min(1,2*min((1+sum(draws<=0))/(length(draws)+1),
                                     (1+sum(draws>=0))/(length(draws)+1))),
    n_blocks=length(blocks),n_clusters=nrow(z))
}
