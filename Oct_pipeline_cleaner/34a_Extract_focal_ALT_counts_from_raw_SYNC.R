## 34a_Extract_focal_ALT_counts_from_raw_SYNC.R
## Read ONLY nominated SNPs from the original raw SYNC, preserving exact
## reference/ALT nucleotides supplied by the VCF-derived founder file.
## Does not create genome-wide selection p-values.

sync_file <- 'Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync'
manifest_file <- 'Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv'
founder_signature_file <- 'historical_QTL_R5/deep_dive_extended/26_Founder_allele_signature_all_sites.csv'
tag_file <- 'historical_QTL_R5/deep_dive_extended/26_R5_model_tag_SNPs.csv'
focal_file <- 'Oct_pipeline_cleaner/26_Focal_variants.csv'
out_dir <- 'selection_results/frequency_strata'
scan_chunk_lines <- 25000L
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(founder_signature_file))founder_signature_file <-
  'historical_QTL_R5/deep_dive/26_Founder_allele_signature_all_sites.csv'
if(!file.exists(tag_file))tag_file <-
  'historical_QTL_R5/deep_dive/26_R5_model_tag_SNPs.csv'
required <- c(sync_file,manifest_file,founder_signature_file,tag_file,focal_file)
if(any(!file.exists(required)))
  stop('Missing required files: ',paste(required[!file.exists(required)],collapse=', '))

meta <- read.csv(manifest_file,stringsAsFactors=FALSE)
if(nrow(meta)!=13 || !all(c('sample_id','sync_column','generations','replicate')%in%names(meta)) ||
   anyDuplicated(meta$sample_id) || anyDuplicated(meta$sync_column))
  stop('Unexpected sample manifest; use the verified original sync-column mapping')
if(!identical(sort(as.integer(meta$sync_column)),4:16))stop('13 selected columns must be 4:16')
sig <- read.csv(founder_signature_file,stringsAsFactors=FALSE)
tag <- read.csv(tag_file,stringsAsFactors=FALSE)
focal <- read.csv(focal_file,stringsAsFactors=FALSE)
req_sig <- c('SNP_ID','CHROM','POS','REF','ALT')
if(!all(req_sig%in%names(sig)))stop('Founder VCF signatures lack REF/ALT')
chosen <- unique(c(as.character(tag$tag_SNP_ID),paste(focal$CHROM,focal$POS,sep='_')))
## Remove duplicate positions from annotation and reject ambiguous multi-allelic sites.
sig <- sig[sig$SNP_ID%in%chosen,,drop=FALSE]
## A SNP may be listed under multiple R5 groups; retain one allele annotation
## if all duplicated records agree on CHROM/POS/REF/ALT.
for(sid in unique(sig$SNP_ID[duplicated(sig$SNP_ID)])){
  vv <- sig[sig$SNP_ID==sid,c('CHROM','POS','REF','ALT'),drop=FALSE]
  if(nrow(unique(vv)) != 1L)stop('Conflicting VCF alleles for ',sid)
}
sig <- sig[!duplicated(sig$SNP_ID),,drop=FALSE]
sig$REF <- toupper(as.character(sig$REF));sig$ALT<-toupper(as.character(sig$ALT))
sig <- sig[!is.na(sig$POS)&sig$REF%in%c('A','C','G','T')&
           sig$ALT%in%c('A','C','G','T')&sig$REF!=sig$ALT,,drop=FALSE]
if(!nrow(sig))stop('No target sites have unambiguous biallelic VCF annotation')
key <- paste(sig$CHROM,sig$POS,sep=':')
if(anyDuplicated(key))stop('Duplicate chromosome-position targets')
index <- setNames(seq_along(key),key)
## Scan streaming blocks, avoiding loading the giant sync into memory.
conn <- if(grepl('\\.gz$',sync_file))gzfile(sync_file,'rt')else file(sync_file,'rt')
open(conn)
found <- vector('list',length(key));names(found)<-key
by_chr <- split(sig$POS,as.character(sig$CHROM))
rx <- vapply(names(by_chr),function(ch){
   paste0('^',ch,'[[:space:]]+(',paste(unique(by_chr[[ch]]),collapse='|'),')[[:space:]]')
}, character(1))
count <- 0L
repeat {
  lines <- readLines(conn,n=scan_chunk_lines,warn=FALSE)
  if(!length(lines))break
  count <- count+length(lines)
  for(regex in rx){
    hit <- grep(regex,lines,perl=TRUE)
    if(!length(hit))next
    for(ln in lines[hit]){
      fields <- strsplit(ln,'[[:space:]]+',perl=TRUE)[[1]]
      if(length(fields)<3)next
      id <- paste(fields[1],fields[2],sep=':')
      ix <- unname(index[id])
      if(length(ix) && !is.na(ix))found[[ix]] <- fields
    }
  }
  if(count%%1000000L<scan_chunk_lines)cat('SYNC rows read: ',count,'\n')
}
close(conn)
qc<-list();res<-list();base_names<-c('A','T','C','G')
for(i in seq_len(nrow(sig))){
  row <- sig[i,];x <- found[[i]]
  matches <- !is.null(x) && length(x)>=max(meta$sync_column) &&
              toupper(x[3]) == row$REF
  qc[[i]]<-data.frame(SNP_ID=row$SNP_ID,REF=row$REF,ALT=row$ALT,
      in_raw_SYNC=!is.null(x),sync_REF=if(is.null(x))NA_character_ else x[3],
      REF_agrees=if(is.null(x))NA else toupper(x[3])==row$REF,
      usable=matches,stringsAsFactors=FALSE)
  if(!matches)next
  for(j in seq_len(nrow(meta))){
    cc<-suppressWarnings(as.numeric(strsplit(x[meta$sync_column[j]],':',fixed=TRUE)[[1]]))
    if(length(cc)<4L || any(!is.finite(cc[1:4])))next
    names(cc)[1:4] <- base_names
    dep <- sum(cc[base_names]);alt <- cc[row$ALT]
    res[[length(res)+1L]]<-data.frame(SNP_ID=row$SNP_ID,CHROM=row$CHROM,
       POS=row$POS,REF=row$REF,ALT=row$ALT,sample_id=meta$sample_id[j],
       generation=meta$generations[j],replicate=meta$replicate[j],
       ALT_reads=as.integer(alt),depth=as.integer(dep),
       ALT_frequency=if(dep>0)unname(alt/dep) else NA_real_,
       REF_agrees=TRUE,stringsAsFactors=FALSE)
  }
}
qc<-do.call(rbind,qc)
write.csv(qc,file.path(out_dir,'34a_Raw_SYNC_target_QC.csv'),row.names=FALSE)
if(!length(res))stop('No usable nominated SNPs in raw SYNC; inspect 34a QC')
res<-do.call(rbind,res)
write.csv(res,file.path(out_dir,'34a_Focal_raw_SYNC_ALT_counts.csv'),row.names=FALSE)
cat('Extracted ',length(unique(res$SNP_ID)),' of ',nrow(qc),
    ' unambiguous targeted loci; scanned ',count,' SYNC lines.\n',sep='')
