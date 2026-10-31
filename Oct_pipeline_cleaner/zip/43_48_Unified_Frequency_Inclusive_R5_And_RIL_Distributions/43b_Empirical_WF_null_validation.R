## 43b_Empirical_WF_null_validation.R
## PRE-RUN VALIDATION: empirical type-I error of the whole A00-conditional
## neutral test under FIXED true ancestral frequencies, NOT only a prior
## predictive check. Needed because p0 is estimated from finite A00 reads.
## This tests approximate .05/.01 tails; it cannot certify the ~1e-5 tail.
## Computationally costly. Set simulation counts BEFORE any final run.
source('43_Unified_config_and_helpers.R')
cfg <- uf_config
source('09f_drift_functions.R')
source('33_Frequency_drift_helpers.R')
out <- cfg$prefix
uf_check(c(file.path(out,'43_FROZEN_selection_tags.csv'),cfg$pool_file,
           cfg$orientation_file,cfg$status_file))
tags <- uf_read(file.path(out,'43_FROZEN_selection_tags.csv'))
p <- read_drift_inputs(cfg$pool_file,cfg$orientation_file,exclude_G0='A6G140')
sets <- uf_status_and_samples(p,cfg)
keep <- sets$primary
pop <- p$pop[keep,,drop=FALSE]
idx <- match(as.character(tags$tag_SNP_ID),as.character(p$sel$ID))
if(anyNA(idx))stop('Tag IDs missing from oriented Pool-seq data')
D <- p$depth[idx,keep,drop=FALSE]
## Choose templates by ancestral coverage; preserve REAL evolved read depths.
valid <- which(rowSums(is.finite(D)&D>0)==ncol(D))
if(length(valid)<5L)stop('Not enough complete-depth tags for null validation')
d0 <- D[valid,which(pop$generations==0),drop=TRUE]
ord <- valid[order(d0)]
tmpl <- unique(ord[pmax(1L,round(seq(1,length(ord),length.out=7L)))])
true_frequencies <- c(.001,.005,.02,.07,.17,.35,.5)
n_outer_per_frequency <- 400L
n_inner_per_outer <- 2500L
set.seed(cfg$seed+703L)
res <- list();j <- 0L
for(p_true in true_frequencies){
  pp <- numeric(0)
  for(tidx in tmpl){
    depths <- as.numeric(D[tidx,])
    ## Simulate BOTH A00 read sampling and subsequent independent WF paths.
    nb <- ceiling(n_outer_per_frequency/length(tmpl))
    Y <- simulate_wf_site(p0=p_true,depths=depths,pop=pop,Ne=cfg$Ne,
                          n_sim=nb,pool_individuals=NA_real_)
    K <- round(sweep(Y,2,depths,'*'))
    for(rr in seq_len(nrow(K))){
      k <- as.numeric(K[rr,]);d <- depths
      info <- conditional_neutral_info(k,d,pop,Ne=cfg$Ne)
      if(!isTRUE(info$ok))next
      sims <- simulate_conditional_neutral(info,k,d,pop,Ne=cfg$Ne,
                                           n_sim=n_inner_per_outer)
      observed <- abs(info$slope-info$mean_slope)
      pp <- c(pp,(1+sum(abs(sims)>=observed))/(length(sims)+1))
    }
  }
  if(length(pp)>n_outer_per_frequency)pp <- pp[seq_len(n_outer_per_frequency)]
  if(length(pp)<100L)stop('Insufficient successful neutral datasets')
  for(alpha in c(.05,.01)){
    exceed <- sum(pp<=alpha);N <- length(pp)
    ## 95% one-sided binomial CI for the rejection rate.
    lower <- if(exceed==0L)0 else qbeta(.01/(length(true_frequencies)*2),exceed,N-exceed+1L)
    upper <- if(exceed==N)1 else qbeta(1-.01/(length(true_frequencies)*2),exceed+1L,N-exceed)
    j <- j+1L
    res[[j]] <- data.frame(true_ancestral_ALT_frequency=p_true,
      MAF=pmin(p_true,1-p_true),nominal_rejection=alpha,
      empirical_rejection=exceed/N,CI_99_familywise_one_sided_lower=lower,
      CI_99_familywise_one_sided_upper=upper,outer_datasets=N,
      inner_MC_per_dataset=n_inner_per_outer,
      clear_anti_conservatism=lower>alpha,Ne=cfg$Ne)
  }
  cat('p0=',p_true,' : ',length(pp),' valid null replicates\n',sep='')
}
result <- do.call(rbind,res)
uf_write(result,file.path(out,'43b_Empirical_null_calibration_FIXED_p0.csv'))
failed <- any(result$clear_anti_conservatism)
## This is a diagnostic and a gate for interpretation, not proof of correct
## genome-wide extreme-tail calibration. Flag even a single failed condition.
status <- if(failed)'FAIL: fixed-p0 empirical calibration showed anti-conservative p-values' else
  'NO_DETECTED_INFLATION: limited power; extreme genomic tail still unvalidated'
writeLines(c(status,paste('n_outer=',n_outer_per_frequency),
             paste('n_inner=',n_inner_per_outer)),
 file.path(out,'43b_Calibration_interpretation_gate.txt'))
if(failed)warning(status,'; do NOT claim calibrated genome-wide FDR until corrected.')
cat('43b DONE: ',status,'\n',sep='')
