## 28b_Exact_SNP_raw_SYNC_recovery.R
## OPTIONAL but important for founder-private or post-founder mutations not in
## script-08's founder-SNP-filtered matrices. Reads the original SYNC directly,
## for ALL nominated exact focal variants and current R5 tag SNPs.
## Does not calculate GWAS or a drift-selection q-value.

sync_file <- 'Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync'
manifest_file <- 'Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv'
founder_signature_file <- 'historical_QTL_R5/deep_dive/26_Founder_allele_signature_all_sites.csv'
tag_file <- 'historical_QTL_R5/deep_dive/26_R5_model_tag_SNPs.csv'
focal_file <- 'Oct_pipeline_cleaner/26_Focal_variants.csv'
out_dir <- 'historical_QTL_R5/deep_dive/evolution/raw_sync'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
req<-c(sync_file,manifest_file,founder_signature_file,tag_file,focal_file)
if(any(!file.exists(req)))stop('Missing: ',paste(req[!file.exists(req)],collapse=', '))
meta<-read.csv(manifest_file,stringsAsFactors=FALSE)
if(!all(c('sample_id','sync_column','generations','replicate')%in%names(meta)))stop('Malformed verified sync manifest')
if(anyDuplicated(meta$sample_id)||anyDuplicated(meta$sync_column))stop('Duplicate sample/column in manifest')
sig<-read.csv(founder_signature_file,stringsAsFactors=FALSE)
tags<-read.csv(tag_file,stringsAsFactors=FALSE)
focal<-read.csv(focal_file,stringsAsFactors=FALSE)
targets<-unique(c(paste0(focal$CHROM,'_',focal$POS),tags$tag_SNP_ID))
sig<-sig[!duplicated(sig$SNP_ID),,drop=FALSE]
regions<-sig[match(targets,sig$SNP_ID),c('SNP_ID','CHROM','POS','REF','ALT'),drop=FALSE]
regions<-regions[is.finite(regions$POS)&nchar(regions$REF)==1L &nchar(regions$ALT)==1L &
          regions$ALT%in%c('A','C','G','T'),,drop=FALSE]
if(!nrow(regions))stop('No biallelic nominated SNP annotations from script 26')
key<-setNames(seq_len(nrow(regions)),paste(regions$CHROM,regions$POS,sep=':'))
bychr<-split(regions$POS,regions$CHROM)
regex<-lapply(names(bychr),function(chr){
 pos<-unique(bychr[[chr]])
 paste0('^',chr,'[[:space:]]+(',paste(pos,collapse='|'),')[[:space:]]')
})
names(regex)<-names(bychr)

conn<-if(grepl('\\.gz$',sync_file))gzfile(sync_file,'rt')else file(sync_file,'rt')
open(conn)
found<-list();n_lines<-0L
repeat{
 lines<-readLines(conn,n=50000L,warn=FALSE)
 if(!length(lines))break
 n_lines<-n_lines+length(lines)
 for(rx in regex){
   take<-grepl(rx,lines,perl=TRUE)
   if(any(take))for(ln in lines[take]){
      x<-strsplit(ln,'[[:space:]]+',perl=TRUE)[[1]]
      kk<-paste(x[1],x[2],sep=':')
      if(kk%in%names(key))found[[kk]]<-x
   }
 }
 if(n_lines%%1000000L<50000L)cat('SYNC lines scanned:',n_lines,'\n')
}
close(conn)
cat('SYNC rows scanned: ',n_lines,'; nominated loci found: ',length(found),'\n',sep='')

res<-list();qc<-list(); bases<-c('A','T','C','G','N','DEL')
for(j in seq_len(nrow(regions))){
 v<-regions[j,];id<-paste(v$CHROM,v$POS,sep=':')
 x<-found[[id]]
 qc[[j]]<-data.frame(SNP_ID=v$SNP_ID,found_in_raw_SYNC=!is.null(x),
   VCF_REF=v$REF,VCF_ALT=v$ALT,
   sync_REF=if(!is.null(x))x[3]else NA_character_,
   REF_agrees=if(!is.null(x))toupper(x[3])==toupper(v$REF)else NA)
 if(is.null(x))next
 if(length(x)<max(meta$sync_column))stop('Incomplete sync columns at ',id)
 for(i in seq_len(nrow(meta))){
   cc<-as.numeric(strsplit(x[meta$sync_column[i]],':',fixed=TRUE)[[1]])
   if(length(cc)<4L||any(!is.finite(cc[1:4])))next
   names(cc)[1:4]<-bases[1:4]
   dep<-sum(cc[bases[1:4]])
   a<-cc[v$ALT]
   res[[length(res)+1L]]<-data.frame(SNP_ID=v$SNP_ID,CHROM=v$CHROM,
     POS=v$POS,REF=v$REF,ALT=v$ALT,
     sample_id=meta$sample_id[i],generation=meta$generations[i],
     replicate=meta$replicate[i],ALT_reads=unname(a),depth=dep,
     ALT_frequency=if(dep>0)a/dep else NA_real_,
     sync_REF=x[3],REF_agrees= toupper(x[3])==toupper(v$REF),
     stringsAsFactors=FALSE)
 }
}
qc<-do.call(rbind,qc)
write.csv(qc,file.path(out_dir,'28b_Raw_SYNC_exact_locus_QC.csv'),row.names=FALSE)
if(!length(res))stop('None of the nominated loci is extractable from raw SYNC')
res<-do.call(rbind,res)
write.csv(res,file.path(out_dir,'28b_Exact_raw_SYNC_ALT_trajectories.csv'),row.names=FALSE)

## Per replicate descriptive summary. Biological ALT counts use base identity;
## REF mismatch is kept visible for assembly/orientation audit.
out<-list()
for(sid in unique(res$SNP_ID)){
 x<-res[res$SNP_ID==sid&is.finite(res$ALT_frequency),,drop=FALSE]
 bg<-x$ALT_frequency[x$generation==0 & x$sample_id=='A00']
 for(r in unique(x$replicate[x$generation>0])){
   z<-x[x$replicate==r & x$generation>0,,drop=FALSE]
   z<-z[order(z$generation),,drop=FALSE]
   if(!nrow(z))next
   out[[length(out)+1L]]<-data.frame(SNP_ID=sid,replicate=r,
       final_generation=tail(z$generation,1),final_ALT_frequency=tail(z$ALT_frequency,1),
       G0_A00=if(length(bg)==1L)bg else NA_real_,
       delta_from_A00=if(length(bg)==1L)tail(z$ALT_frequency,1)-bg else NA_real_,
       n_timepoints=nrow(z),REF_agrees=all(z$REF_agrees),stringsAsFactors=FALSE)
 }
}
if(length(out))write.csv(do.call(rbind,out),file.path(out_dir,'28b_Per_GA_raw_SYNC_changes.csv'),row.names=FALSE)

pdf(file.path(out_dir,'28b_Nominated_exact_variant_trajectories.pdf'),width=8,height=5)
for(sid in unique(res$SNP_ID)){
 x<-res[res$SNP_ID==sid&is.finite(res$ALT_frequency),,drop=FALSE]
 if(!nrow(x))next
 reps<-unique(x$replicate[x$generation>0]);cols<-setNames(rainbow(length(reps)),reps)
 plot(NA,xlim=c(0,max(x$generation)),ylim=c(0,1),
      xlab='Generations in salt experiment',ylab='GWAS ALT frequency',main=sid)
 for(r in reps){
   z<-x[x$replicate==r & x$generation>0,,drop=FALSE]
   z<-z[order(z$generation),,drop=FALSE]
   lines(z$generation,z$ALT_frequency,type='b',pch=16,col=cols[r],lwd=1.5)
 }
 z<-x[x$generation==0,,drop=FALSE]
 if(nrow(z))points(z$generation,z$ALT_frequency,pch=ifelse(z$sample_id=='A00',16,1),col='black')
 legend('topleft',legend=c(reps,'A00 baseline','other G0'),
        col=c(cols,'black','black'),pch=c(rep(16,length(reps)),16,1),bty='n',cex=.8)
}
dev.off()
cat('Raw SYNC check complete. No drift-adjusted p/q was estimated for low-frequency sites.\n')
