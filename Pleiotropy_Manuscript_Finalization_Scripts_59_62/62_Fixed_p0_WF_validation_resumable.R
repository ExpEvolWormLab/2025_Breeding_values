# 62_Fixed_p0_WF_validation_resumable.R
# PRIORITY: audit type-I error of the ORIGINAL fixed-A00-count posterior-
# conditional Wright-Fisher test using simulated data at known true p0.
# This is a more extensively powered version of old script 43b; it does NOT
# automatically approve 12 selection calls or validate ~1e-5 tails.
# Source identical functions as scripts 43b/44. Requires original Pool-seq
# RData and raw-SYNC recovery inputs. Run from project root.
# Warning: may be computationally expensive; pilot first, then final.

CONFIG <- list(
  mode = 'pilot', # change to 'final' for publication-stage diagnostics
  true_p0 = c(.001,.005,.020,.070,.170,.350,.500),
  alpha = c(.01,.05),
  Ne_grid = c(1000L), # optional sensitivity: c(500L,1000L,2000L)
  n_outer_pilot = 400L, n_inner_pilot = 2500L,
  n_outer_final = 2400L, n_inner_final = 5000L,
  seed = 20261103L,
  out = 'selection_results/manuscript_finalization/62_neutral_calibration'
)
if(!CONFIG$mode%in%c('pilot','final'))stop('mode must be pilot or final')
# Separate pilot/final and Ne settings so finished checkpoints cannot collide.
CONFIG$out <- file.path(CONFIG$out,CONFIG$mode,
                        paste0('Ne_',paste(CONFIG$Ne_grid,collapse='_')))
if(any(CONFIG$alpha<=0|CONFIG$alpha>=1))stop('Invalid nominal alpha')

source('Oct_pipeline_cleaner/43_Unified_config_and_helpers.R')
source('Oct_pipeline_cleaner/09f_drift_functions.R')
source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')
cfg<-uf_config
uf_check(c(file.path(cfg$prefix,'43_FROZEN_selection_tags.csv'),
           cfg$pool_file,cfg$orientation_file,cfg$status_file,
           cfg$phase_file,cfg$raw_focal_file))
source('Oct_pipeline_cleaner/43x_Integrate_rawSYNC_R5_SNPs.R')

p<-uf_pool_with_R5_raw_recovery(cfg)
sets<-uf_status_and_samples(p,cfg)
keep<-sets$primary
pop<-p$pop[keep,,drop=FALSE]
if(sum(pop$generations==0)!=1L)stop('Exactly one A00 baseline required')
tags<-uf_read(file.path(cfg$prefix,'43_FROZEN_selection_tags.csv'))
if(!'tag_SNP_ID'%in%names(tags))stop('Missing frozen tag_SNP_ID')
j<-match(tags$tag_SNP_ID,p$sel$ID)
if(anyNA(j))stop('Frozen tag absent in reconstructed Pool-seq object')
D<-p$depth[j,keep,drop=FALSE]
valid<-which(rowSums(is.finite(D)&D>0)==ncol(D))
if(length(valid)<10L)stop('Too few complete-depth SNP templates')
d0<-D[valid,which(pop$generations==0),drop=TRUE]
ord<-valid[order(d0)]
# Prespecified depth templates selected before seeing any simulated p-values.
tmpl<-unique(ord[pmax(1L,round(seq(1,length(ord),length.out=7L)))])
if(length(tmpl)<5L)stop('Insufficient distinct depth templates')
Nout<-if(CONFIG$mode=='pilot') CONFIG$n_outer_pilot else CONFIG$n_outer_final
Nin<-if(CONFIG$mode=='pilot') CONFIG$n_inner_pilot else CONFIG$n_inner_final
if(Nin<1000L || Nout<200L)warning('Small simulation size; expect limited power')
if(!all(is.finite(CONFIG$true_p0)&CONFIG$true_p0>0&CONFIG$true_p0<1))stop('Invalid p0 grid')
dir.create(CONFIG$out,recursive=TRUE,showWarnings=FALSE)
# Freeze and store file hashes. Checkpoint reuse never crosses simulation
# settings or different tag/pool files.
inputs<-c(file.path(cfg$prefix,'43_FROZEN_selection_tags.csv'),cfg$pool_file,
 cfg$orientation_file,cfg$status_file,cfg$phase_file,cfg$raw_focal_file)
hashes<-unname(tools::md5sum(inputs))
signature<-paste(c(hashes,CONFIG$mode,Nout,Nin,CONFIG$Ne_grid,
 CONFIG$true_p0,CONFIG$seed,p$pop$sample_id[keep]),collapse='|')
sig_file<-file.path(CONFIG$out,'62_input_signature.txt')
if(file.exists(sig_file)){
 if(!identical(readLines(sig_file,warn=FALSE),signature))
  stop('Input signature changed; use a new output directory (do not mix checkpoints).')
}else writeLines(signature,sig_file)

# Reserve a fixed number of complete-data simulations per depth template.
# Invalid fits receive p=1 (not dropped) for unconditional null rejection.
all_rows<-list()
for(Ne in CONFIG$Ne_grid){
 for(pp in CONFIG$true_p0){
   # A deterministic distribution of outer datasets over fixed depth templates.
   counts<-rep(Nout%/%length(tmpl),length(tmpl))
   if(Nout%%length(tmpl)>0)counts[seq_len(Nout%%length(tmpl))]<-counts[seq_len(Nout%%length(tmpl))]+1L
   for(tt in seq_along(tmpl)){
     key<-sprintf('Ne%d_p%.5f_template%02d_%s',Ne,pp,tt,CONFIG$mode)
     checkpoint<-file.path(CONFIG$out,paste0('62_',key,'.rds'))
     if(file.exists(checkpoint)){
       a<-readRDS(checkpoint)
       if(nrow(a)!=counts[tt])stop('Checkpoint dimension mismatch: ',checkpoint)
       all_rows[[length(all_rows)+1L]]<-a
       next
     }
     depths<-as.numeric(D[tmpl[tt],]); n<-counts[tt]
     # Fixed true p0 generates BOTH A00 reads and independent WF trajectories.
     # This is the crucial difference from a prior-predictive check.
     set.seed(as.integer(CONFIG$seed+Ne*17+round(pp*1e6)+tt*100003L))
     Y<-simulate_wf_site(p0=pp,depths=depths,pop=pop,Ne=Ne,n_sim=n,
                         pool_individuals=NA_real_)
     K <- round(sweep(Y, 2, depths, '*'))

depth_matrix <- matrix(
  rep(depths, each = n),
  nrow = n,
  ncol = length(depths)
)

K[K < 0] <- 0
K[K > depth_matrix] <- depth_matrix[K > depth_matrix]

stopifnot(
  is.matrix(K),
  nrow(K) == n,
  ncol(K) == length(depths)
)
     outp<-rep(1,n);valid_fit<-rep(FALSE,n)
     for(rr in seq_len(n)){
       k<-as.numeric(K[rr,]);d<-depths
       info<-conditional_neutral_info(k,d,pop,Ne=Ne)
       if(!isTRUE(info$ok))next
       sim<-simulate_conditional_neutral(info,k,d,pop,Ne=Ne,n_sim=Nin)
       if(length(sim)!=Nin || any(!is.finite(sim)))stop('Invalid conditional simulations')
       obs<-abs(info$slope-info$mean_slope)
       outp[rr]<-(1+sum(abs(sim)>=obs))/(Nin+1)
       valid_fit[rr]<-TRUE
     }
     a<-data.frame(Ne=Ne,true_p0=pp,MAF=pmin(pp,1-pp),
        template=tt,sim_index=seq_len(n),
        p_value=outp,fit_valid=valid_fit,
        n_inner=Nin,mode=CONFIG$mode)
     saveRDS(a,checkpoint)
     all_rows[[length(all_rows)+1L]]<-a
     cat('Completed ',key,' / ',n,' outer null datasets\n',sep='')
   }
 }
}
R<-do.call(rbind,all_rows)
if(anyDuplicated(paste(R$Ne,R$true_p0,R$template,R$sim_index)))
 stop('Duplicate simulation IDs')
write.csv(R,file.path(CONFIG$out,paste0('62_',CONFIG$mode,'_outer_pvalues.csv')),row.names=FALSE)

agg<-list()
for(Ne in CONFIG$Ne_grid)for(pp in CONFIG$true_p0){
 zz<-R[R$Ne==Ne & abs(R$true_p0-pp)<1e-12,,drop=FALSE]
 if(nrow(zz)!=Nout)stop('Incomplete simulations for Ne=',Ne,' p0=',pp)
 for(a in CONFIG$alpha){
   k<-sum(zz$p_value<=a)
   n<-nrow(zz)
   ci<-binom.test(k,n)$conf.int
   p_inflation<-pbinom(k-1L,size=n,prob=a,lower.tail=FALSE)
   agg[[length(agg)+1L]]<-data.frame(Ne=Ne,true_p0=pp,
    MAF=pmin(pp,1-pp),nominal=a,n_outer=n,n_inner=Nin,
    n_reject=k,rejection_rate=k/n,lower95=ci[1L],upper95=ci[2L],
    invalid_fits=sum(!zz$fit_valid),
    p_one_sided_inflation=p_inflation,
    Monte_Carlo_resolution=1/(Nin+1))
 }
}
A<-do.call(rbind,agg)
A$Holm_p_inflation<-p.adjust(A$p_one_sided_inflation,method='holm')
A$inflation_detected<-A$Holm_p_inflation<.05 & A$rejection_rate>A$nominal
write.csv(A,file.path(CONFIG$out,paste0('62_',CONFIG$mode,'_calibration_summary.csv')),row.names=FALSE)
gate<-if(any(A$inflation_detected)){
 'FAIL: detectable anti-conservative behavior at at least one tested p0/alpha. Do not present 12 BY hits as calibrated selection discoveries.'
}else{
 'INCONCLUSIVE: no detected inflation at tested levels; the extreme genome-wide significance tail remains UNVALIDATED.'
}
writeLines(c(gate,'This script does not replace the original 43b/44 testing gate.',
 'This test includes A00 read sampling at fixed true p0.',
 'Non-estimable outer simulations counted as non-rejections (p=1).',
 sprintf('N_outer=%d, N_inner=%d, MC minimum attainable p=%.6g',Nout,Nin,1/(Nin+1))),
 file.path(CONFIG$out,paste0('62_',CONFIG$mode,'_interpretation.txt')))
pdf(file.path(CONFIG$out,paste0('62_',CONFIG$mode,'_typeI_diagnostics.pdf')),width=9,height=5)
op<-par(mfrow=c(1,length(CONFIG$alpha)),mar=c(5,4.3,3,1))
for(aa in CONFIG$alpha){
 ss<-A[A$nominal==aa & A$Ne==cfg$Ne,,drop=FALSE]
 if(!nrow(ss))next
 x<-ss$true_p0
 plot(x,ss$rejection_rate,log='x',pch=19,xlab='True ancestral ALT frequency',
   ylab='Empirical rejection rate',main=paste('Nominal alpha =',aa),
   ylim=range(c(0,ss$upper95,aa)))
 arrows(x,ss$lower95,x,ss$upper95,angle=90,code=3,length=.03,col='gray40')
 abline(h=aa,lty=2,col='red3')
}
par(op);dev.off()
writeLines(capture.output(sessionInfo()),file.path(CONFIG$out,'62_sessionInfo.txt'))
print(A)
cat('62 COMPLETE: ',gate,'\n')
