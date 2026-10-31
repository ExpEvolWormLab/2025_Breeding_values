## 38_Calibrate_rare_WF_by_site.R
## Calibrate the *conditional* neutral-score statistic by discrete-generation
## Wright-Fisher + sequencing-read simulation. Never interpret approximate
## rare-SNP p-values from script 33 as a calibrated genome-wide FDR screen.
## Fixed baseline A00; main scenario = GA1,GA2,GA4; GA3 supplementary.
## Requires the ORIGINAL 09f drift helper, 33 helper, corrected script-08 RData.

helper09 <- if(file.exists('09f_drift_functions.R')) '09f_drift_functions.R' else
  'drift_nullfit_update/09f_drift_functions.R'
source(helper09)
source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')

pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
orientation_file <- 'selection_results/09_Selection_SNP_results.csv'
screen_dir <- 'selection_results/frequency_strata'
status_file <- 'Oct_pipeline_cleaner/20_GA_population_status.csv'
predeclared_file <- file.path(screen_dir,'34_Predeclared_target_registry.csv')
output_dir <- file.path(screen_dir,'rare_calibration')
Ne_primary <- 1000L
Ne_sensitivity <- c(500L,2000L)
random_sites_per_class <- 16L    # independent of observed selection statistic
screened_candidates_per_rare_class <- 6L  # ranked by |Z|; EXPLORATORY
sim_random <- 3000L
sim_candidates <- 12000L
min_G0_depth <- 40L
seed <- 20261029L
run_all_GA_supplement <- TRUE
set.seed(seed)
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

if(!all(file.exists(c(pool_file,orientation_file,status_file))))stop('Missing Poolseq inputs or status table')
status <- read.csv(status_file,stringsAsFactors=FALSE)
if(!all(c('replicate','status') %in% names(status)))stop('Malformed GA status file')
status$status <- tolower(trimws(as.character(status$status)))
keep_reps <- as.character(status$replicate[status$status=='survived'])
if(length(keep_reps)!=3L)warning('Expected three surviving GA populations; found ',length(keep_reps))
if(!any(status$status=='extinct'))warning('No GA explicitly classified extinct; review the status file')

p <- read_drift_inputs(pool_file,orientation_file,exclude_G0='A6G140')
if(!identical(colnames(p$k_alt),as.character(p$pop$sample_id)))stop('Pool-column alignment failure')
if(sum(p$pop$generations==0)!=1L||p$pop$sample_id[p$pop$generations==0]!='A00')
  stop('Require one validated A00 baseline')

classes <- c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50')
scenarios <- list(retained_GA=p$pop$generations==0|as.character(p$pop$replicate)%in%keep_reps)
if(run_all_GA_supplement)scenarios$all_GA <- rep(TRUE,nrow(p$pop))
all_results <- list(); all_null <- list(); all_qc <- list()

for(scenario in names(scenarios)){
  screen_file <- file.path(screen_dir,paste0('33_All_SNPs_',scenario,'_Ne1000.csv'))
  if(!file.exists(screen_file))stop('Run script 33 first: ',screen_file)
  ss <- read.csv(screen_file,stringsAsFactors=FALSE)
  if(anyDuplicated(ss$ID) || anyDuplicated(p$sel$ID))stop('Duplicate SNP IDs')
  order_ix <- match(as.character(ss$ID),as.character(p$sel$ID))
  if(anyNA(order_ix))stop('Script-33 and Poolseq SNP universes differ')
  good <- is.finite(ss$z_normal_approx)&is.finite(ss$G0_depth)&
          ss$G0_depth>=min_G0_depth&!is.na(ss$MAF_class)
  ix_rand <- integer();ix_candidate <- integer()
  for(b in classes){
    ii <- which(good&ss$MAF_class==b)
    if(length(ii))ix_rand <- c(ix_rand,sample(ii,min(length(ii),random_sites_per_class)))
    ## Selection on large |Z| is a SCREEN; MC p-values are not GW-significant.
    if(b %in% classes[1:3]&&length(ii)){
      ii <- setdiff(ii,ix_rand)
      ord <- ii[order(abs(ss$z_normal_approx[ii]),decreasing=TRUE)]
      ix_candidate <- c(ix_candidate,head(ord,screened_candidates_per_rare_class))
    }
  }
  ix_declared <- integer()
  if(file.exists(predeclared_file)){
    decl <- read.csv(predeclared_file,stringsAsFactors=FALSE)
    if('SNP_ID'%in%names(decl))ix_declared <- which(ss$ID%in%decl$SNP_ID & good)
  }
  selected <- unique(c(ix_rand,ix_candidate,ix_declared))
  roles <- ifelse(selected%in%ix_rand,'random_calibration',
                  ifelse(selected%in%ix_declared,'predeclared_in_GWAS_universe',
                         'observed_extreme_screen'))
  write.csv(data.frame(scenario=scenario,SNP_ID=ss$ID[selected],
                       MAF_class=ss$MAF_class[selected],role=roles),
            file.path(output_dir,paste0('38_Sites_',scenario,'.csv')),row.names=FALSE)
  keep <- scenarios[[scenario]]
  pp <- p$pop[keep,,drop=FALSE]
  for(j in seq_along(selected)){
    ii <- selected[j]; row <- order_ix[ii]
    k <- as.numeric(p$k_alt[row,keep]);d <- as.numeric(p$depth[row,keep])
    role <- roles[j]
    Ne_grid <- if(role=='random_calibration')Ne_primary else c(Ne_primary,Ne_sensitivity)
    for(Ne in Ne_grid){
      info <- conditional_neutral_info(k,d,pp,Ne=Ne)
      if(!isTRUE(info$ok)){
        all_qc[[length(all_qc)+1L]] <- data.frame(scenario=scenario,ID=ss$ID[ii],
                                                   Ne=Ne,reason='fit_not_estimable')
        next
      }
      B <- if(role=='random_calibration')sim_random else sim_candidates
      mc <- mc_neutral_test(info,k,d,pp,Ne=Ne,n_sim=B)
      draw <- mc$sim
      approx_p_null <- 2*pnorm(abs(draw/info$conditional_SE),lower.tail=FALSE)
      all_results[[length(all_results)+1L]] <- data.frame(
        scenario=scenario,ID=as.character(ss$ID[ii]),
        MAF_class=as.character(ss$MAF_class[ii]),role=role,Ne=Ne,
        n_sim=mc$B,G0_reads=info$baseline_ALT_reads,G0_depth=info$baseline_depth,
        posterior_mean_p0=info$p0,
        observed_slope=info$slope,expected_neutral_slope=info$mean_slope,
        observed_Z=info$conditional_z,
        p_normal_approx=info$conditional_p_normal,
        p_conditional_WF_MC=mc$p,
        mc_floor=1/(mc$B+1),
        n_null_exceed=mc$exceed,
        null_stat_SD=mc$neutral_sd,
        null_fp_normal_p05=mean(approx_p_null<.05),
        null_fp_normal_p01=mean(approx_p_null<.01),
        null_fp_normal_p001=mean(approx_p_null<.001),
        null_fp_normal_p0001=mean(approx_p_null<1e-4),
        stringsAsFactors=FALSE)
      if(role=='random_calibration'&&Ne==Ne_primary)
        all_null[[length(all_null)+1L]] <- data.frame(
          scenario=scenario,ID=ss$ID[ii],MAF_class=ss$MAF_class[ii],
          normal_p_under_neutral=approx_p_null,stringsAsFactors=FALSE)
    }
    if(j%%10L==0L)cat('Scenario ',scenario,': ',j,'/',length(selected),' sites\n',sep='')
  }
}
if(!length(all_results))stop('No valid MC tests; inspect script 33 and inputs')
res <- do.call(rbind,all_results)
write.csv(res,file.path(output_dir,'38_Sitewise_conditional_WF_results.csv'),row.names=FALSE)
if(length(all_qc))write.csv(do.call(rbind,all_qc),file.path(output_dir,'38_Site_QC_failures.csv'),row.names=FALSE)
if(length(all_null)){
  nul <- do.call(rbind,all_null)
  ## Average per-SITE false-positive rates: site is the sampling unit,
  ## not each of its correlated WF simulations.
  nul$fp05 <- as.numeric(nul$normal_p_under_neutral<.05)
  nul$fp01 <- as.numeric(nul$normal_p_under_neutral<.01)
  nul$fp001 <- as.numeric(nul$normal_p_under_neutral<.001)
  nul$fp0001 <- as.numeric(nul$normal_p_under_neutral<1e-4)
  site <- aggregate(cbind(fp05,fp01,fp001,fp0001)~scenario+MAF_class+ID,
                    data=nul,FUN=mean)
  bybin <- aggregate(cbind(fp05,fp01,fp001,fp0001)~scenario+MAF_class,site,mean)
  n_sites <- aggregate(ID~scenario+MAF_class,site,length)
  names(n_sites)[3] <- 'n_independent_sites'
  bybin <- merge(bybin,n_sites,by=c('scenario','MAF_class'))
  write.csv(site,file.path(output_dir,'38_Null_false_positive_by_site.csv'),row.names=FALSE)
  write.csv(bybin,file.path(output_dir,'38_Null_calibration_by_frequency.csv'),row.names=FALSE)
  pdf(file.path(output_dir,'38_Calibration_vs_nominal.pdf'),width=9,height=5)
  for(s in unique(bybin$scenario)){
    y <- bybin[bybin$scenario==s,,drop=FALSE]
    y <- y[match(classes,y$MAF_class),,drop=FALSE]
    barplot(t(as.matrix(y[,c('fp05','fp01','fp001'),drop=FALSE])),beside=TRUE,
      names.arg=y$MAF_class,las=2,ylab='Neutral false-positive fraction',
      main=paste('Approximate normal-test calibration:',s),
      col=c('grey40','grey65','grey85'))
    legend('topright',legend=c('nominal 0.05','nominal 0.01','nominal 0.001'),
           fill=c('grey40','grey65','grey85'),bty='n',cex=.8)
  }
  dev.off()
}
cat('Script 38 complete. Targeted MC p-values are not genome-wide q-values.\n')
