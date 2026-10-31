## 18_CHRIV_founder_and_Poolseq_variants.R
## Direct validation of candidate positions in raw Pool-seq SYNC and original VCF.
## Needs indexed .vcf.gz and raw .sync. Must supply 16 founder sample IDs.
## Does not assume that a zero read count establishes a true de novo mutation.

positions <- c(8892663L,10499135L)
chromosome <- 'IV'
vcf_file <- '/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz'
sync_file <- 'Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync'
founder_ids_file <- '18_founder_sample_ids.txt' # exactly 16 VCF sample IDs, one per line
manifest_file <- 'Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv'
out_dir <- 'selection_results/chrIV_candidate_mutations'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(vcf_file))stop('Set vcf_file in script 18')
if(!file.exists(sync_file))stop('Cannot locate raw Pool-seq SYNC')
if(!file.exists(founder_ids_file))stop('Supply the 16 VCF founder sample IDs in ',founder_ids_file)
founders<-trimws(readLines(founder_ids_file,warn=FALSE))
founders<-founders[nzchar(founders)]
if(length(founders)!=16L || anyDuplicated(founders))
  stop('Founder list must contain exactly 16 distinct VCF sample IDs')

vcf_samples<-system2('bcftools',c('query','-l',shQuote(vcf_file)),stdout=TRUE)
if(!all(founders %in% vcf_samples))stop('Founder IDs absent from VCF: ',
  paste(setdiff(founders,vcf_samples),collapse=', '))
## Query all genotype columns, explicitly match sample names from bcftools query -l.
variant_rows <- list()
for(pos in positions) {
  reg<-paste0(chromosome,':',pos,'-',pos)
  raw<-system2('bcftools',c('query','-r',shQuote(reg),'-f',
    shQuote('%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'),shQuote(vcf_file)),
    stdout=TRUE,stderr=FALSE)
  if(length(raw)!=1L) {
    warning('Expected one biallelic VCF entry at ',reg,'; got ',length(raw))
    next
  }
  x<-strsplit(raw,'\t',fixed=TRUE)[[1]]
  if(length(x)!=4+length(vcf_samples))stop('VCF sample count mismatch at ',pos)
  gt<-x[4+match(founders,vcf_samples)]
  if(any(grepl('1',gt,fixed=TRUE)))warning('ALT allele occurs in some founder genotypes at ',pos)
  if(any(grepl('\\.',gt)))warning('Missing founder genotypes at ',pos)
  variant_rows[[as.character(pos)]]<-data.frame(chr=x[1],position=as.integer(x[2]),
    REF=x[3],ALT=x[4],founder_id=founders,GT=gt,
    ALT_present=grepl('1',gt,fixed=TRUE),stringsAsFactors=FALSE)
}
if(!length(variant_rows))stop('Neither candidate position found in VCF; verify assembly/chromosome names')
founder_table<-do.call(rbind,variant_rows)
write.csv(founder_table,file.path(out_dir,'18_Founder_genotypes.csv'),row.names=FALSE)

## Raw SYNC column convention: CHR POS REF then one A:T:C:G:N:DEL block
## for each sample, 1-based sync_column as recorded in the manifest.
meta<-read.csv(manifest_file,stringsAsFactors=FALSE)
if(!all(c('sample_id','sync_column','generations','replicate')%in%names(meta)))
  stop('Check manifest field names')
con<-if(grepl('\\.gz$',sync_file))gzfile(sync_file,'rt') else file(sync_file,'rt')
open(con)
selected<-list();count_lines<-0L
repeat {
  ln<-readLines(con,n=50000L,warn=FALSE)
  if(!length(ln))break
  count_lines<-count_lines+length(ln)
  ## Efficient prefilter to exact chromosome and positions.
  take<-grepl(paste0('^',chromosome,'\\s+(',paste(positions,collapse='|'),')\\s'),ln)
  if(any(take)){
    for(line in ln[take]) {
      z<-strsplit(line,'\\s+')[[1]]
      selected[[paste(z[1],z[2],sep=':')]]<-z
    }
  }
}
close(con)
cat('Scanned',count_lines,'SYNC lines; found',length(selected),'positions\n')
base_names<-c('A','T','C','G','N','DEL')
results<-NULL
for(pos in positions){
  key<-paste(chromosome,pos,sep=':')
  z<-selected[[key]]
  if(is.null(z)){warning('Position not found in raw SYNC: ',key);next}
  v<-subset(founder_table,position==pos)
  if(!nrow(v))next
  ref<-unique(v$REF);alt<-unique(v$ALT)
  if(length(ref)!=1||length(alt)!=1|| nchar(ref)!=1||nchar(alt)!=1||
     !ref%in%base_names[1:4]||!alt%in%base_names[1:4]){
     warning('Non-simple SNP at ',pos);next
  }
  for(j in seq_len(nrow(meta))){
    col<-as.integer(meta$sync_column[j])
    if(col>length(z))stop('SYNC column missing for ',meta$sample_id[j])
    counts<-as.numeric(strsplit(z[col],':',fixed=TRUE)[[1]])
    if(length(counts)<4 || any(!is.finite(counts[1:4])))next
    names(counts)[1:4]<-base_names[1:4]
    cov<-sum(counts[base_names[1:4]])
    alt_n<-counts[alt]
    results<-rbind(results,data.frame(chr=chromosome,position=pos,
      sample_id=meta$sample_id[j],generation=meta$generations[j],
      replicate=meta$replicate[j],REF=ref,ALT=alt,
      alt_reads=alt_n,coverage=cov,alt_frequency=if(cov>0)alt_n/cov else NA_real_))
  }
}
if(is.null(results))stop('No usable candidate allele counts')
write.csv(results,file.path(out_dir,'18_Candidate_mutation_trajectories.csv'),row.names=FALSE)

pdf(file.path(out_dir,'18_CHRIV_mutation_trajectories.pdf'),width=9,height=4)
par(mfrow=c(1,2),mar=c(4.5,4.4,3,1))
for(pos in positions){
  z<-subset(results,position==pos)
  if(!nrow(z))next
  plot(NA,xlim=range(z$generation),ylim=c(0,max(.35,z$alt_frequency,na.rm=TRUE)),
       xlab='Generation',ylab='ALT allele frequency',
       main=paste0('IV:',pos),bty='n')
  g0<-subset(z,generation==0)
  if(nrow(g0))points(g0$generation,g0$alt_frequency,pch=1,cex=1.2)
  pops<-sort(unique(z$replicate[z$generation>0]))
  col<-grDevices::hcl.colors(length(pops),'Dark 3')
  for(k in seq_along(pops)){
    a<-subset(z,replicate==pops[k] & generation>0)
    a<-a[order(a$generation),]
    lines(a$generation,a$alt_frequency,type='b',pch=16,col=col[k],lwd=1.5)
  }
  legend('topleft',legend=paste('evolved',pops),col=col,pch=16,lty=1,bty='n',cex=.65)
}
par(mfrow=c(1,1));dev.off()
cat('Mutation output:',out_dir,'\n')
cat('Zero reads is a detection limit, not evidence for true absence.\n')

## Existing script-15 prediction errors: proximity sensitivity using the
## same two fixed mutation coordinates, without claiming causality.
err_file <- 'selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/final_chrIV/15_CHRIV_cluster_prediction_errors.csv'
if(file.exists(err_file)){
 d<-read.csv(err_file,stringsAsFactors=FALSE)
 if(all(c('position','delta_sqerr_full_vs_three')%in%names(d))){
  distance<-vapply(d$position,function(z)min(abs(z-positions)),numeric(1))
  summary<-do.call(rbind,lapply(c(250000,500000,1000000,2000000),function(w){
    near<-distance<=w
    data.frame(window_bp=w,n_near=sum(near),n_else=sum(!near),
      mean_reduction_near=if(any(near))mean(d$delta_sqerr_full_vs_three[near],na.rm=TRUE) else NA_real_,
      mean_reduction_else=if(any(!near))mean(d$delta_sqerr_full_vs_three[!near],na.rm=TRUE) else NA_real_)
  }))
  write.csv(summary,file.path(out_dir,'18_CHRIV_pred_error_proximity.csv'),row.names=FALSE)
  print(summary)
  cat('CAUTION: Proximity windows contain few clusters; descriptive comparison only.\n')
 }
}
