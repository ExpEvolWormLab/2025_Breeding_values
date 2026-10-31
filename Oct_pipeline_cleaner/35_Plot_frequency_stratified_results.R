## 35_Plot_frequency_stratified_results.R
## Manuscript exploratory figures: observed GWAS-SNP responses by G0 MAF,
## conditional null-calibration diagnostics, nominated low-frequency examples.
## Does NOT turn exploratory normal p-values into genome-wide discoveries.

out_dir <- 'selection_results/frequency_strata'
fig_dir <- file.path(out_dir,'figures')
dir.create(fig_dir,recursive=TRUE,showWarnings=FALSE)
allfile<-file.path(out_dir,'33_Matched_full_vs_retained_all_frequencies.csv')
calfile<-file.path(out_dir,'34b_Null_calibration_by_frequency.csv')
testfile<-file.path(out_dir,'34_Focal_WF_neutral_conditional_tests.csv')
countfile<-file.path(out_dir,'34_Targeted_ALT_counts_in_analysis_pools.csv')
required<-c(allfile,calfile,testfile,countfile)
if(any(!file.exists(required)))
 stop('Run scripts 33, 34 and 34b before plotting. Missing: ',
      paste(required[!file.exists(required)],collapse=', '))

z<-read.csv(allfile,stringsAsFactors=FALSE)
cal<-read.csv(calfile,stringsAsFactors=FALSE)
test<-read.csv(testfile,stringsAsFactors=FALSE)
count<-read.csv(countfile,stringsAsFactors=FALSE)
classes<-c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50')
class_colors<-setNames(c('#bc4b51','#de8f05','#a9a431','#238b45','#31688e'),classes)

pdf(file.path(fig_dir,'35_Genomewide_frequency_strata.pdf'),width=10,height=7)
par(mfrow=c(2,2),mar=c(4.6,4.6,3.2,1))
for(s in c('full','retained')){
  sl<-z[[paste0(s,'_slope')]]
  boxplot(split(sl,z$MAF_class)[classes],outline=FALSE,
    ylab='ALT-frequency slope per generation',main=paste('GWAS SNPs:',s,'GA'),
    las=2,col=unname(class_colors),cex.axis=.8)
  abline(h=0,lty=2,col='grey50')
}
valid<-is.finite(z$full_slope)&is.finite(z$retained_slope)
if(sum(valid)>10000L){set.seed(20261029);valid[sample(which(valid),sum(valid)-10000L)]<-FALSE}
plot(z$full_slope[valid],z$retained_slope[valid],pch=16,cex=.3,
   col=adjustcolor('black',alpha.f=.14),
   xlab='Full GA frequency slope',ylab='Retained GA frequency slope',
   main='Matched SNP responses')
abline(a=0,b=1,col='firebrick',lty=2)
bars<-table(factor(z$MAF_class,levels=classes))
barplot(bars,col=unname(class_colors),las=2,ylab='Number of GWAS SNPs',
        main='G0 posterior-mean MAF classes',cex.names=.7)
dev.off()

pdf(file.path(fig_dir,'35_Neutral_calibration_by_frequency.pdf'),width=9,height=6)
par(mar=c(7,5,3,1))
cal<-cal[match(classes,cal$MAF_class),,drop=FALSE]
bar<-rbind(cal$nominal_fp_05,cal$nominal_fp_01,cal$nominal_fp_001)
barplot(bar,beside=TRUE,names.arg=classes,las=2,
        col=c('#31688e','#de8f05','#bc4b51'),
        ylab='Neutral simulated fraction below nominal cutoff',
        main='Frequency dependence of normal-tail calibration')
abline(h=.05,lty=2,col='#31688e');abline(h=.01,lty=2,col='#de8f05');
abline(h=.001,lty=2,col='#bc4b51')
legend('topright',legend=c('p<0.05','p<0.01','p<0.001'),
       fill=c('#31688e','#de8f05','#bc4b51'),bty='n')
dev.off()

## Focal locus diagnostics, 1000-Ne main analysis for each scenario.
a<-test[test$Ne==1000 & test$status=='tested',,drop=FALSE]
if(nrow(a)){
  pdf(file.path(fig_dir,'35_Nominated_locus_WF_evidence.pdf'),width=11,height=7)
  par(mar=c(5,8,3,2))
  ids<-unique(a$SNP_ID)
  aa<-a[a$scenario=='all_GA',,drop=FALSE];bb<-a[a$scenario=='retained_GA',,drop=FALSE]
  ids<-ids[ids%in%aa$SNP_ID & ids%in%bb$SNP_ID]
  if(length(ids)){
    oa<-match(ids,aa$SNP_ID);ob<-match(ids,bb$SNP_ID)
    aa<-aa[oa,,drop=FALSE];bb<-bb[ob,,drop=FALSE]
    pp<-rbind(pmax(aa$p_WF_conditional,1/(aa$n_sim+1)),
              pmax(bb$p_WF_conditional,1/(bb$n_sim+1)))
    mx<-max(-log10(pp),na.rm=TRUE)
    plot(NA,xlim=c(0,mx*1.1),ylim=c(.5,length(ids)+.5),yaxt='n',
         xlab='-log10(conditional WF neutral p; capped by MC resolution)',ylab='',
         main='Nominated sites: full GA vs retained GA')
    axis(2,at=seq_along(ids),labels=ids,las=2,cex.axis=.65)
    points(-log10(pp[1,]),seq_along(ids)+.13,pch=16,col='#31688e')
    points(-log10(pp[2,]),seq_along(ids)-.13,pch=17,col='#de8f05')
    abline(v=-log10(.05),lty=2,col='grey50')
    legend('topright',legend=c('All GA','Retained GA'),pch=c(16,17),
           col=c('#31688e','#de8f05'),bty='n')
  } else {plot.new();text(.5,.5,'No matched focal tests')}
  dev.off()
}

## Exact ALT frequency curves for selected historical loci.
want<-c('V_6322080','V_6637695','V_6826421','IV_8892663','IV_10499135',
        'I_7167738','I_7491656','I_7596359')
pdf(file.path(fig_dir,'35_Focal_low_frequency_trajectories.pdf'),width=9,height=6)
n_trajectory_pages <- 0L
for(sid in want){
  a<-count[count$SNP_ID==sid,,drop=FALSE]
  if(!nrow(a))next
  n_trajectory_pages <- n_trajectory_pages + 1L
  a$ALT_frequency<-ifelse(a$depth>0,a$ALT_reads/a$depth,NA_real_)
  repids<-unique(a$replicate[a$generations>0])
  colors<-setNames(c('#0072B2','#D55E00','#009E73','#CC79A7')[seq_along(repids)],repids)
  plot(NA,xlim=c(0,max(a$generations)),ylim=c(0,1),
       xlab='Salt-selection generation',ylab='Sequenced ALT frequency',
       main=paste('Low-frequency focal SNP',sid))
  g0<-a[a$generations==0,,drop=FALSE]
  if(nrow(g0))points(g0$generations,g0$ALT_frequency,pch=16,col='black',cex=1.2)
  for(r in repids){
    x<-a[a$replicate==r & a$generations>0,,drop=FALSE]
    x<-x[order(x$generations),,drop=FALSE]
    lines(x$generations,x$ALT_frequency,type='b',pch=16,lwd=2,col=colors[as.character(r)])
  }
  legend('topleft',legend=c(paste('replicate',repids),'A00'),
    col=c(unname(colors),'black'),pch=16,bty='n',cex=.85)
}
if(n_trajectory_pages==0L){plot.new();text(.5,.5,'No named low-frequency target counts found')}
dev.off()
cat('Script 35 complete: ',fig_dir,'\n',sep='')
