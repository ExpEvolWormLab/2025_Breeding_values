# 67b_Correct_WS283_curated_features.R
# Correct script 67's functional classification by using only curated
# WormBase-source gene/transcript/CDS/exon/UTR features; no VEP needed.
# Run from the CeMEE project root (base R only).
CFG <- list(
  input_dir = 'selection_results/manuscript_finalization/67_WS283_local_annotations',
  out_dir = 'selection_results/manuscript_finalization/67_WS283_local_annotations'
)
O <- read.csv(file.path(CFG$input_dir,'67_exact_overlap_GFF3_features.csv'),
              stringsAsFactors=FALSE,check.names=FALSE)
B <- read.csv(file.path(CFG$input_dir,'67_WS283_broad_variant_annotation.csv'),
              stringsAsFactors=FALSE,check.names=FALSE)
required <- c('SNP_ID','source','type','ID','Parent','Name','locus','Alias','strand')
stopifnot(all(required %in% names(O)),all(c('SNP_ID','REF','ALT') %in% names(B)))
O <- O[!is.na(O$source) & O$source=='WormBase',,drop=FALSE]
priority <- c('CDS','five_prime_UTR','three_prime_UTR','intron',
              'exon','mRNA','nc_primary_transcript','gene')
labels <- c('CDS'='coding sequence (CDS)',
            'five_prime_UTR'='5-prime UTR',
            'three_prime_UTR'='3-prime UTR',
            'intron'='intron',
            'exon'='exon (no specified UTR/CDS)',
            'mRNA'='transcribed gene region',
            'nc_primary_transcript'='nc primary transcript',
            'gene'='gene body')
choose_class <- function(types){
    z <- priority[priority %in% types]
    if(!length(z))return('no curated WormBase feature')
    unname(labels[z[1]])
}
rows <- list(); transcript_rows <- list()
for(i in seq_len(nrow(B))){
  v<-B[i,]; w<-O[O$SNP_ID == v$SNP_ID,,drop=FALSE]
  g<-w[w$type=='gene',,drop=FALSE]
  if(nrow(g)!=1L)stop('Expected exactly one WormBase gene for ',v$SNP_ID,
                      '; found ',nrow(g),' (inspect manually)')
  gene<-g[1,]
  gene_name<-if(!is.na(gene$locus) && nzchar(gene$locus))gene$locus else
    if(!is.na(gene$Alias) && nzchar(gene$Alias))strsplit(gene$Alias,',',fixed=TRUE)[[1]][1] else gene$Name
  tr<-w[w$type %in% c('mRNA','nc_primary_transcript'),,drop=FALSE]
  rows[[i]]<-data.frame(SNP_ID=v$SNP_ID, REF=v$REF, ALT=v$ALT,
    gene=gene_name, gene_id=gene$ID, gene_strand=gene$strand,
    WS283_WormBase_feature=choose_class(w$type),
    n_gene_transcript_models_overlapping=nrow(tr),
    transcript_names=paste(tr$Name,collapse=';'),
    source_only='WormBase GFF3 source (excluding auxiliary gene predictors)',
    coding_change_computed=FALSE,
    reporting_caveat='Transcript overlap, not functional validation',
    stringsAsFactors=FALSE)
  if(nrow(tr))for(k in seq_len(nrow(tr))){
    tx<-tr[k,]
    sub<-w[!is.na(w$Parent) & vapply(w$Parent,function(p){
      if(is.na(p))return(FALSE)
      tx$ID %in% strsplit(p,',',fixed=TRUE)[[1]]
    },logical(1)),,drop=FALSE]
    ftypes<-unique(sub$type)
    cls<-choose_class(ftypes)
    if(cls %in% c('transcribed gene region','nc primary transcript'))cls<-'transcribed gene region'
    transcript_rows[[length(transcript_rows)+1L]]<-data.frame(
      SNP_ID=v$SNP_ID,REF=v$REF,ALT=v$ALT,gene=gene_name,gene_id=gene$ID,
      transcript_id=tx$ID,transcript_name=tx$Name,strand=gene$strand,
      consequence_from_gff=cls,feature_types=paste(sort(ftypes),collapse=';'),
      stringsAsFactors=FALSE)
  }
}
dir.create(CFG$out_dir,recursive=TRUE,showWarnings=FALSE)
A<-do.call(rbind,rows)
T<-do.call(rbind,transcript_rows)
write.csv(A,file.path(CFG$out_dir,'67b_WS283_curated_variant_summary.csv'),row.names=FALSE)
write.csv(T,file.path(CFG$out_dir,'67b_WS283_transcript_feature_summary.csv'),row.names=FALSE)
print(A[,c('SNP_ID','gene','WS283_WormBase_feature')])
cat('Output:',normalizePath(CFG$out_dir),'\n')
