## 44_Calibrated_WF_selection_all_R5_clusters.R
## Manuscript-grade, computationally intensive and RESUMABLE.
## One exact, independently GWAS-chosen tag per TESTABLE R5 cluster.
## Conditional on A00 read count; simulate discrete WF drift and sequencing
## at each observed time point. No minor-allele-frequency filter.
## Use simultaneously valid upper bounds on MC p-values with BY correction.
## Do NOT use script-33 normal p-values as q-values.
source('Oct_pipeline_cleaner/43_Unified_config_and_helpers.R')
cfg<-uf_config
source('Oct_pipeline_cleaner/09f_drift_functions.R');source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')
out<-cfg$prefix;uf_check(c(cfg$pool_file,cfg$orientation_file,cfg$status_file,
 file.path(out,'43_Unified_R5_cluster_registry.csv'),
 file.path(out,'43_FROZEN_selection_tags.csv')))
registry<-uf_read(file.path(out,'43_Unified_R5_cluster_registry.csv'))
validation_gate<-file.path(out,'43b_Calibration_interpretation_gate.txt')
if(!file.exists(validation_gate))stop('Run 43b empirical fixed-p0 null validation before 44')
signature_file<-file.path(out,'43b_calibration_input_signature.txt')
if(!file.exists(signature_file) ||
   !identical(readLines(signature_file,warn=FALSE)[1L],uf_calibration_signature(cfg)))
 stop('Stale calibration gate: rerun 43b after raw-sync integration or changed tags/config')
validation_status<-readLines(validation_gate,warn=FALSE)[1]
if(grepl('^FAIL',validation_status))stop('Fixed-p0 empirical null calibration FAILED. Do not run significance analysis until the null model is fixed.')
source('Oct_pipeline_cleaner/43x_Integrate_rawSYNC_R5_SNPs.R')
p<-uf_pool_with_R5_raw_recovery(cfg)
sets<-uf_status_and_samples(p,cfg)
if(anyDuplicated(p$sel$ID))stop('Duplicate oriented Pool-seq SNP IDs')
row_index<-setNames(seq_len(nrow(p$sel)),as.character(p$sel$ID))
usable<-registry[!is.na(registry$tag_SNP_ID)&nzchar(registry$tag_SNP_ID),,drop=FALSE]
if(!nrow(usable))stop('No eligible tags; run script 43')
if(anyNA(match(usable$tag_SNP_ID,names(row_index))))stop('Tag absent from corrected Pool-seq object')
if(anyDuplicated(usable$tag_SNP_ID))stop('Duplicate frozen tags')

## On a restart, checkpoint ensures the tested tags and parameters match.
## Changing any of these requires a new output directory; never mix nulls.
make_signature<-function(scen, m){
 paste(scen,uf_calibration_signature(cfg),cfg$Ne,
   paste(cfg$mc_stages,collapse=':'),cfg$mc_family_error,
   cfg$min_A00_depth,m,paste(usable$tag_SNP_ID,collapse=':'),sep='|')
}

mc_one_site<-function(j,scenario,idx_keep,n_tests,scenario_index){
  g<-usable[j,,drop=FALSE]
  ix<-row_index[as.character(g$tag_SNP_ID)]
  k<-as.numeric(p$k_alt[ix,idx_keep]);d<-as.numeric(p$depth[ix,idx_keep])
  pop<-p$pop[idx_keep,,drop=FALSE]
  info<-conditional_neutral_info(k,d,pop,Ne=cfg$Ne)
  if(!isTRUE(info$ok))return(data.frame(scenario=scenario,group=g$group,
    tag_SNP_ID=g$tag_SNP_ID,tag_phase=g$tag_phase,
    n_sim=0,exceed=NA_integer_,p_mc_raw=NA_real_,p_MC_lower=NA_real_,
    p_MC_upper=NA_real_,p_MC_upper_simultaneous=NA_real_,
    observed_ALT_slope=NA_real_,observed_phase_slope=NA_real_,
    posterior_A00_ALT=NA_real_,A00_ALT_reads=k[which(pop$generations==0)[1]],
    A00_depth=d[which(pop$generations==0)[1]],
    stop_reason='WF_fit_failed',Ne=cfg$Ne,stringsAsFactors=FALSE))
  ## Bonferroni across all sites, interim stages, scenarios, conditional
  ## MC intervals. Coverage is over stochastic MC uncertainty, NOT over the
  ## biological uncertainty in Ne or assumptions of the WF process.
  n_scen<-if(cfg$run_all_GA_supp)2L else 1L
  alpha<-cfg$mc_family_error/(n_tests*length(cfg$mc_stages)*n_scen)
  harmonic<-sum(1/seq_len(n_tests))
  max_BY_p<-0.05/harmonic
  min_BY_p<-0.05/(n_tests*harmonic)
  set.seed(as.integer(cfg$seed + scenario_index*100003L + j))
  B<-0L; exceed<-0L; reason<-'max_MC_budget'
  observed<-abs(info$slope-info$mean_slope)
  bounds<-c(lower=0,upper=1)
  if(!is.finite(observed)||!is.finite(info$mean_slope))stop('Nonfinite conditional statistic')
  if(observed==0){reason<-'zero_observed_deviation';bounds<-c(lower=1,upper=1)}
  else for(target in cfg$mc_stages){
    needed<-as.integer(target-B)
    while(needed>0L){
      take<-min(needed,cfg$mc_chunk)
      sim<-simulate_conditional_neutral(info,k,d,pop,Ne=cfg$Ne,n_sim=take)
      exceed<-exceed+sum(abs(sim)>=observed)
      B<-B+length(sim);needed<-needed-length(sim)
    }
    bounds<-uf_mc_bounds(exceed,B,alpha)
    ## Valid early stopping: ANY returned upper bound is simultaneous,
    ## including the adaptively chosen stopping stage.
    if(bounds['lower']>max_BY_p){reason<-'incompatible_with_BY_05';break}
    if(bounds['upper']<=min_BY_p){reason<-'resolved_Bonferroni_BY_rank1';break}
  }
  ## The raw MC estimator is recorded ONLY as a diagnostic: optional-stop
  ## raw p estimates cannot be used directly for genome-wide FDR.
  data.frame(scenario=scenario,group=as.character(g$group),
    tag_SNP_ID=as.character(g$tag_SNP_ID),tag_phase=as.integer(g$tag_phase),
    n_sim=B,exceed=exceed,
    p_mc_raw=(exceed+1)/(B+1),
    p_MC_lower=unname(bounds['lower']),p_MC_upper=unname(bounds['upper']),
    p_MC_upper_simultaneous=unname(bounds['upper']),
    observed_ALT_slope=info$slope,
    observed_phase_slope=as.numeric(g$tag_phase)*info$slope,
    posterior_A00_ALT=info$p0,
    A00_ALT_reads=info$baseline_ALT_reads,A00_depth=info$baseline_depth,
    stop_reason=reason,Ne=cfg$Ne,stringsAsFactors=FALSE)
}

run_one_scenario<-function(scenario,keep,scenario_index){
  main<-scenario=='primary_survivors'
  ## Precheck all tags under THIS scenario, then include all fit-valid tests
  ## in multiplicity; no selection on observed evolution magnitude.
  valid<-vapply(seq_len(nrow(usable)),function(j){
    ix<-row_index[as.character(usable$tag_SNP_ID[j])]
    z<-conditional_neutral_info(as.numeric(p$k_alt[ix,keep]),
       as.numeric(p$depth[ix,keep]),p$pop[keep,,drop=FALSE],cfg$Ne)
    isTRUE(z$ok)
  },logical(1))
  if(main && any(!valid))stop('Frozen primary tag not estimable. Check script33/43 QC')
  n_tests<-sum(valid)
  if(n_tests<100L)stop('Too few testable R5 cluster tags: ',n_tests)
  cat('\n',scenario,': ',n_tests,' WF tests, Ne=',cfg$Ne,'\n',sep='')
  signature<-make_signature(scenario,n_tests)
  checkpoint<-file.path(out,paste0('44_checkpoint_',scenario,'.rds'))
  results<-list()
  if(file.exists(checkpoint)){
    old<-readRDS(checkpoint)
    if(!identical(old$signature,signature))
      stop('Checkpoint settings mismatch. DO NOT mix analyses. Rename old checkpoint before changing config.')
    results<-old$results
    cat('Resuming ',length(results),' finished clusters\n')
  }
  pending<-which(valid & !usable$group%in%names(results))
  batches<-split(pending,ceiling(seq_along(pending)/cfg$batch_clusters))
  if(length(pending))for(b in seq_along(batches)){
    jj<-batches[[b]]
    nc<-if(.Platform$OS.type=='windows')1L else min(length(jj),cfg$max_cores)
    rr<-if(nc>1L)parallel::mclapply(jj,function(j)
      mc_one_site(j,scenario,keep,n_tests,scenario_index),
      mc.cores=nc,mc.preschedule=FALSE,mc.set.seed=FALSE) else
      lapply(jj,function(j)mc_one_site(j,scenario,keep,n_tests,scenario_index))
    for(z in rr){
      if(inherits(z,'try-error')||!is.data.frame(z))stop('A worker failed; inspect stdout and checkpoint')
      results[[as.character(z$group)]]<-z
    }
    tmp<-paste0(checkpoint,'.tmp');saveRDS(list(signature=signature,results=results),tmp)
    if(!file.rename(tmp,checkpoint))stop('Checkpoint rename failed')
    cat('  completed ',length(results),' / ',n_tests,'; batch ',b,'/',length(batches),'\n',sep='')
  }
  if(!length(results))stop('No WF tests completed')
  ans<-do.call(rbind,results);rownames(ans)<-NULL
  ans<-ans[match(usable$group,ans$group),,drop=FALSE]
  ## NA from supplementary untestable tags are displayed, never silently
  ## assigned p=1 or removed from primary membership.
  tested<-is.finite(ans$p_MC_upper_simultaneous)
  if(sum(tested)!=n_tests)stop('WF result/testable count mismatch')
  ans$q_BY_conservative<-NA_real_;ans$q_BH_conservative<-NA_real_
  ans$q_BY_conservative[tested]<-p.adjust(ans$p_MC_upper_simultaneous[tested],method='BY')
  ans$q_BH_conservative[tested]<-p.adjust(ans$p_MC_upper_simultaneous[tested],method='BH')
  ans$significant_BY_05<-tested&!is.na(ans$q_BY_conservative)&ans$q_BY_conservative<.05
  ans$significant_BH_05<-tested&!is.na(ans$q_BH_conservative)&ans$q_BH_conservative<.05
  ans$significant_BY_01<-tested&!is.na(ans$q_BY_conservative)&ans$q_BY_conservative<.01
  ans$MC_family_error_budget<-cfg$mc_family_error
  ans$Ne_primary<-cfg$Ne
  ans$valid_primary_clusters<-sum(registry$reason=='eligible')
  uf_write(ans,file.path(out,paste0('44_',scenario,'_WF_cluster_selection.csv')))
  ans
}

primary<-run_one_scenario('primary_survivors',sets$primary,1L)
if(cfg$run_all_GA_supp){
  supplement<-run_one_scenario('supp_all_GA_including_extinct_GA3',sets$all,2L)
}else supplement<-NULL

## Predeclared Ne sensitivity for only the primary BY hits plus ten
## most promising primary tags: never reinterpret as genome-wide q-values.
if(cfg$run_Ne_sensitivity){
  hits<-which(primary$significant_BY_05)
  top<-order(primary$p_MC_upper_simultaneous,na.last=NA)
  choose<-unique(c(hits,head(top,10L)))
  sensitivity<-list();zidx<-0L
  for(Ne in cfg$Ne_sensitivity){
    for(j in choose){
      ix<-row_index[as.character(usable$tag_SNP_ID[j])]
      k<-as.numeric(p$k_alt[ix,sets$primary]);d<-as.numeric(p$depth[ix,sets$primary])
      pop<-p$pop[sets$primary,,drop=FALSE]
      info<-conditional_neutral_info(k,d,pop,Ne=Ne)
      if(!isTRUE(info$ok))next
      set.seed(as.integer(cfg$seed+Ne*13L+j))
      B<-0L;ex<-0L;obs<-abs(info$slope-info$mean_slope)
      while(B<cfg$Ne_sensitivity_B){
        take<-min(cfg$mc_chunk,cfg$Ne_sensitivity_B-B)
        sim<-simulate_conditional_neutral(info,k,d,pop,Ne=Ne,n_sim=take)
        ex<-ex+sum(abs(sim)>=obs);B<-B+length(sim)
      }
      zidx<-zidx+1L
      sensitivity[[zidx]]<-data.frame(group=usable$group[j],tag_SNP_ID=usable$tag_SNP_ID[j],
        Ne=Ne,mc_exceed=ex,n_sim=B,p_MC_fixed_budget=(ex+1)/(B+1),
        p_one_sided_upper_99=if(ex==B)1 else qbeta(.99,ex+1,B-ex),
        stringsAsFactors=FALSE)
    }
  }
  if(length(sensitivity))uf_write(do.call(rbind,sensitivity),
    file.path(out,'44_Sensitivity_Ne_target_tags_NOT_genomewide_q.csv'))
}

summary_one<-function(z){data.frame(scenario=z$scenario[1],n_tested=sum(is.finite(z$p_MC_upper_simultaneous)),
  n_BY_q05=sum(z$significant_BY_05,na.rm=TRUE),n_BY_q01=sum(z$significant_BY_01,na.rm=TRUE),
  n_BH_q05=sum(z$significant_BH_05,na.rm=TRUE),
  total_MC_draws=sum(z$n_sim,na.rm=TRUE),
  n_baseline_ALT_zero=sum(z$A00_ALT_reads==0,na.rm=TRUE),
  n_max_budget=sum(z$stop_reason=='max_MC_budget',na.rm=TRUE))}
summary<-summary_one(primary)
if(!is.null(supplement))summary<-rbind(summary,summary_one(supplement))
uf_write(summary,file.path(out,'44_Selection_summary.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'44_sessionInfo.txt'))
cat('44 DONE. PRIMARY BY q<.05: ',sum(primary$significant_BY_05,na.rm=TRUE),
 ' out of ',sum(is.finite(primary$p_MC_upper_simultaneous)), ' testable clusters.\n',sep='')
