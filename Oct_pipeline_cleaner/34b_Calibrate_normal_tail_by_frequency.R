## 34b_Calibrate_normal_tail_by_frequency.R
## Frequency-balanced, randomly sampled NULL simulations to check calibration
## of the baseline-conditional normal approximation from script 33.
## Deliberately sampled without reference to observed test p-values.
## This does NOT prove genome-wide FDR calibration at ~1e-7 or rarer.

helper <- if(file.exists('09f_drift_functions.R'))'09f_drift_functions.R' else
  'drift_nullfit_update/09f_drift_functions.R'
source(helper)
source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')
pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
orientation_file <- 'selection_results/09_Selection_SNP_results.csv'
all_loci_file <- 'selection_results/frequency_strata/33_All_SNPs_all_GA_Ne1000.csv'
out_dir <- 'selection_results/frequency_strata'
exclude_G0 <- 'A6G140'
Ne <- 1000L
sites_per_class <- 12L
n_sim_per_site <- 2500L
min_G0_depth <- 40L
set.seed(20261029)

if(!file.exists(all_loci_file))stop('Run 33 first')
screen<-read.csv(all_loci_file,stringsAsFactors=FALSE)
p<-read_drift_inputs(pool_file,orientation_file,exclude_G0=exclude_G0)
if(!identical(as.character(screen$SNP_ID),as.character(p$sel$rs)))
 stop('Cannot map script-33 results to the site-by-pool count matrices')
classes <- c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50')
selected <- integer()
for(b in classes){
  i<-which(screen$MAF_class==b & is.finite(screen$p_normal_approx) &
             screen$G0_depth>=min_G0_depth)
  if(length(i))selected<-c(selected,sample(i,min(sites_per_class,length(i))))
}
if(!length(selected))stop('No testable loci for neutral calibration')
write.csv(screen[selected,c('SNP_ID','MAF_class','G0_ALT_reads','G0_depth')],
  file.path(out_dir,'34b_Neutral_calibration_random_sites.csv'),row.names=FALSE)
rows<-list();sim_p<-list()
for(j in seq_along(selected)){
  i<-selected[j]
  info<-conditional_neutral_info(p$k_alt[i,],p$depth[i,],p$pop,Ne=Ne)
  if(!info$ok)next
  x<-simulate_conditional_neutral(info,p$k_alt[i,],p$depth[i,],p$pop,
                                  Ne=Ne,n_sim=n_sim_per_site)
  zp<-x/info$conditional_SE
  pp<-2*pnorm(abs(zp),lower.tail=FALSE)
  rows[[length(rows)+1L]]<-data.frame(
    SNP_ID=screen$SNP_ID[i],MAF_class=as.character(screen$MAF_class[i]),
    G0_ALT_reads=info$baseline_ALT_reads,G0_depth=info$baseline_depth,
    Ne=Ne,n_sim=length(pp),null_mean_z=mean(zp),null_sd_z=sd(zp),
    nominal_fp_05=mean(pp<.05),nominal_fp_01=mean(pp<.01),
    nominal_fp_001=mean(pp<.001),nominal_fp_0001=mean(pp<.0001),
    stringsAsFactors=FALSE)
  sim_p[[length(sim_p)+1L]]<-data.frame(MAF_class=as.character(screen$MAF_class[i]),
                     normal_p=pp,stringsAsFactors=FALSE)
  cat('Calibrated ',j,'/',length(selected),' ',screen$SNP_ID[i],'\n',sep='')
}
if(!length(rows))stop('No valid neutral simulations')
loc<-do.call(rbind,rows)
write.csv(loc,file.path(out_dir,'34b_Null_calibration_by_locus.csv'),row.names=FALSE)
all_p<-do.call(rbind,sim_p)
by_bin<-do.call(rbind,lapply(split(all_p$normal_p,all_p$MAF_class),function(x)
  data.frame(n_sim=length(x),nominal_fp_05=mean(x<.05),
    nominal_fp_01=mean(x<.01),nominal_fp_001=mean(x<.001),
    nominal_fp_0001=mean(x<.0001))))
by_bin$MAF_class<-rownames(by_bin)
rownames(by_bin)<-NULL
write.csv(by_bin,file.path(out_dir,'34b_Null_calibration_by_frequency.csv'),row.names=FALSE)
cat('Script 34b complete. Null simulations assess the analytic approximation;\n',
    'they do not validate the biological Ne model or extreme BH tail.\n',sep='')
