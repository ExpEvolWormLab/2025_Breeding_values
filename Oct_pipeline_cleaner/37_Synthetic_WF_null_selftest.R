## 37_Synthetic_WF_null_selftest.R (optional quick debugging check)
## Validates that the simulated null statistic has the analytic expected
## center and scale at boundary, rare, common and high-frequency alleles.
## This checks implementation only; it does not validate biological Ne.

helper <- if(file.exists('09f_drift_functions.R'))'09f_drift_functions.R' else
  'drift_nullfit_update/09f_drift_functions.R'
source(helper)
source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')
set.seed(20261029)
pop <- data.frame(sample_id=c('A00','GA110','GA136','GA150','GA210','GA236',
                              'GA250','GA310','GA336','GA410','GA436','GA450'),
 generations=c(0,10,36,50,10,36,50,10,36,10,36,50),
 replicate=c(1,3,3,3,4,4,4,5,5,6,6,6))
d <- rep(120,nrow(pop));Ne<-1000L
rows <- list()
for(k0 in c(0,1,5,30,60,115,120)){
  truep <- (k0+.5)/121
  k <- rep(round(d[1]*truep),length(d));k[1] <- k0
  info <- conditional_neutral_info(k,d,pop,Ne=Ne)
  if(!info$ok)stop('Invalid synthetic fit at k0=',k0)
  x <- simulate_conditional_neutral(info,k,d,pop,Ne=Ne,n_sim=3000)
  ratio<-sd(x)/info$conditional_SE
  center<-mean(x)/info$conditional_SE
  if(!is.finite(ratio)||ratio<.8||ratio>1.2||abs(center)>.15)
    stop('Synthetic conditional covariance mismatch at k0=',k0,
         ' ratio=',ratio,' mean_z=',center)
  rows[[length(rows)+1L]]<-data.frame(k0=k0,p0=info$p0,
      simulated_SD_over_analytic=ratio,simulated_center_z=center)
}
write.csv(do.call(rbind,rows),'37_Synthetic_WF_null_selftest.csv',row.names=FALSE)
cat('Synthetic WF null moment checks passed.\n')
