# 68_Offline_WS283_CDS_codon_consequences.R
# Calculate codon and putative amino-acid changes at CDS-overlapping variants,
# using the user's exact WS283 genomic FASTA and GFF3 (no VEP; base R only).
# IMPORTANT: Sequence- and frame-QC are mandatory before interpreting effects.
# Works with plain or gzip-compressed FASTA and GFF3.

CFG <- list(
  fasta = 'WS283/20220216_c_elegans_WS283.genome.fa', # EDIT: the exact FASTA used in script 64
  gff = 'WS283/c_elegans.WS283.VEP.gff3.gz', # EDIT
  annotation_dir = 'selection_results/manuscript_finalization/67_WS283_local_annotations',
  out = 'selection_results/manuscript_finalization/68_WS283_CDS_codon_consequences',
  chunk_lines = 50000L
)

need <- c('fasta','gff')
for(nm in need) if(!file.exists(CFG[[nm]]))stop('Set CFG$',nm,' to a file that exists: ',CFG[[nm]])
variants_file <- file.path(CFG$annotation_dir,'67b_WS283_curated_variant_summary.csv')
feature_file <- file.path(CFG$annotation_dir,'67_exact_overlap_GFF3_features.csv')
if(!file.exists(variants_file) || !file.exists(feature_file))
  stop('Run script 67b and keep the original script 67 exact-overlap CSV in CFG$annotation_dir')
V <- read.csv(variants_file,stringsAsFactors=FALSE)
F <- read.csv(feature_file,stringsAsFactors=FALSE)
F <- F[!is.na(F$source) & F$source=='WormBase',,drop=FALSE]
V <- V[grepl('^coding sequence',V$WS283_WormBase_feature),,drop=FALSE]
if(!nrow(V))stop('No CDS-overlapping variants')
V$chr <- sub('_[0-9]+$','',V$SNP_ID)
V$pos <- as.integer(sub('^.+_','',V$SNP_ID))
stopifnot(all(grepl('^[ACGT]$',V$REF)),all(grepl('^[ACGT]$',V$ALT)))

# GFF3 attribute extractor, handling comma-separated Parent fields.
attr_get <- function(s,key) {
  a <- strsplit(s,';',fixed=TRUE)[[1]]
  t <- a[startsWith(a,paste0(key,'='))]
  if(length(t))substring(t[1],nchar(key)+2L) else NA_character_
}

# For each coding variant, identify CDS-overlapping transcript IDs from the
# curated WormBase-source rows (not from auxiliary gene prediction tracks).
selected <- list()
for(i in seq_len(nrow(V))){
 w <- F[F$SNP_ID==V$SNP_ID[i] & F$type=='CDS',,drop=FALSE]
 for(p in w$Parent)if(!is.na(p))for(tid in strsplit(p,',',fixed=TRUE)[[1]]){
   selected[[length(selected)+1L]] <- data.frame(
      SNP_ID=V$SNP_ID[i],gene=V$gene[i], transcript_id=tid,stringsAsFactors=FALSE)
 }
}
S <- unique(do.call(rbind,selected))
if(!nrow(S))stop('Missing curated WormBase CDS Parent annotations')
target_tx <- unique(S$transcript_id)
target_chr <- unique(V$chr)

# Stream the GFF and collect EVERY CDS segment for each transcript, not just
# segments overlapping the candidate SNP.
con <- if(grepl('\\.gz$',CFG$gff,ignore.case=TRUE))gzfile(CFG$gff,'rt') else file(CFG$gff,'rt')
segments <- list()
repeat {
  lines <- readLines(con,n=CFG$chunk_lines,warn=FALSE)
  if(!length(lines))break
  lines <- lines[grepl('\tWormBase\tCDS\t',lines,fixed=TRUE)]
  if(!length(lines))next
  for(s in lines){
    z<-strsplit(s,'\t',fixed=TRUE)[[1]]
    if(length(z)<9L)next
    ch<-sub('^chr','',z[1],ignore.case=TRUE)
    ch<-sub('^CHROMOSOME_','',ch,ignore.case=TRUE)
    if(!(ch %in% target_chr))next
    parents<-attr_get(z[9],'Parent')
    if(is.na(parents))next
    tids<-intersect(strsplit(parents,',',fixed=TRUE)[[1]],target_tx)
    if(!length(tids))next
    for(tid in tids){
      segments[[length(segments)+1L]]<-data.frame(
        chr=ch,start=as.integer(z[4]),end=as.integer(z[5]),strand=z[7],
        phase=suppressWarnings(as.integer(z[8])),transcript_id=tid,
        stringsAsFactors=FALSE)
    }
  }
}
close(con)
if(!length(segments))stop('No matching WormBase CDS segments; check GFF assembly and source')
CDS<-unique(do.call(rbind,segments))

# FASTA reader: holds only I/IV/V sequences, not all the genome.
# Accumulate ONE string per 50,000-line chunk (avoid repeated per-line
# vector growth on large FASTA files).
normalize_chr <- function(h){
  y<-strsplit(sub('^>','',h),'[[:space:]]+')[[1]][1]
  y<-sub('^chr','',y,ignore.case=TRUE)
  sub('^CHROMOSOME_','',y,ignore.case=TRUE)
}
con<-if(grepl('\\.gz$',CFG$fasta,ignore.case=TRUE))gzfile(CFG$fasta,'rt') else file(CFG$fasta,'rt')
seq_parts<-setNames(lapply(target_chr,function(z)list()),target_chr)
current<-NA_character_
append_part<-function(ch,lines) {
  if(is.na(ch)||!(ch%in%target_chr)||!length(lines))return(invisible(NULL))
  x<-toupper(paste0(trimws(lines),collapse=''))
  if(nzchar(x))seq_parts[[ch]][[length(seq_parts[[ch]])+1L]] <<- x
  invisible(NULL)
}
repeat {
  zz<-readLines(con,n=CFG$chunk_lines,warn=FALSE)
  if(!length(zz))break
  marks<-which(startsWith(zz,'>'))
  start<-1L
  if(length(marks))for(k in marks){
    if(k>start)append_part(current,zz[start:(k-1L)])
    current<-normalize_chr(zz[k])
    start<-k+1L
  }
  if(start<=length(zz))append_part(current,zz[start:length(zz)])
}
close(con)
SREF<-lapply(seq_parts,function(z)paste0(unlist(z,use.names=FALSE),collapse=''))
if(any(!nzchar(unlist(SREF,use.names=FALSE))))
 stop('Missing chromosome(s) in FASTA: ',paste(target_chr[!nzchar(unlist(SREF,use.names=FALSE))],collapse=','))
get_seq<-function(ch,s,e) {
  if(s<1L||e>nchar(SREF[[ch]])||e<s)stop('FASTA coordinate out of range for ',ch,':',s,'-',e)
  substring(SREF[[ch]],s,e)
}
complement <- function(x){
  chartr('ACGTN','TGCAN',x)
}
revcomp <- function(x){
  paste0(rev(strsplit(complement(x),'',fixed=TRUE)[[1]]),collapse='')
}
code<-c(
 TTT='F',TTC='F',TTA='L',TTG='L',TCT='S',TCC='S',TCA='S',TCG='S',
 TAT='Y',TAC='Y',TAA='*',TAG='*',TGT='C',TGC='C',TGA='*',TGG='W',
 CTT='L',CTC='L',CTA='L',CTG='L',CCT='P',CCC='P',CCA='P',CCG='P',
 CAT='H',CAC='H',CAA='Q',CAG='Q',CGT='R',CGC='R',CGA='R',CGG='R',
 ATT='I',ATC='I',ATA='I',ATG='M',ACT='T',ACC='T',ACA='T',ACG='T',
 AAT='N',AAC='N',AAA='K',AAG='K',AGT='S',AGC='S',AGA='R',AGG='R',
 GTT='V',GTC='V',GTA='V',GTG='V',GCT='A',GCC='A',GCA='A',GCG='A',
 GAT='D',GAC='D',GAA='E',GAG='E',GGT='G',GGC='G',GGA='G',GGG='G')
translate<-function(codon){
  if(nchar(codon)!=3L||grepl('[^ACGT]',codon))return(NA_character_)
  unname(code[codon])
}
results <- list()
for(i in seq_len(nrow(S))){
 v <- V[V$SNP_ID==S$SNP_ID[i],,drop=FALSE][1,]
 d <- CDS[CDS$transcript_id==S$transcript_id[i],,drop=FALSE]
 if(!nrow(d))next
 if(length(unique(d$strand))!=1L||length(unique(d$chr))!=1L)stop('Inconsistent CDS strand/chr for ',S$transcript_id[i])
 strand<-d$strand[1]
 if(!strand %in% c('+','-'))stop('Unexpected strand: ',strand)
 d<-d[if(strand=='+')order(d$start)else order(-d$end),,drop=FALSE]
 d$len<-d$end-d$start+1L
 expected_phase <- (3L - c(0L, head(cumsum(d$len),-1L)) %% 3L) %% 3L
 phase_consistent <- all(!is.na(d$phase)) && all(d$phase==expected_phase)
 if(!phase_consistent) warning('Nonstandard CDS phase for ',S$transcript_id[i],': check before use')
 coding_chunks<-mapply(function(ch,st,en){
    x<-get_seq(ch,st,en)
    if(strand=='-')revcomp(x)else x
 },d$chr,d$start,d$end,USE.NAMES=FALSE)
 coding<-paste0(coding_chunks,collapse='')
 hits<-which(d$start<=v$pos & d$end>=v$pos)
 if(length(hits)!=1L)stop('Variant does not uniquely overlap a CDS segment: ',v$SNP_ID,' / ',S$transcript_id[i])
 j<-hits[1]
 local_idx<-if(strand=='+')v$pos-d$start[j]+1L else d$end[j]-v$pos+1L
 cds_idx<-sum(d$len[seq_len(j-1L)])+local_idx
 ref_base<-get_seq(v$chr,v$pos,v$pos)
 if(ref_base!=v$REF)stop('REF does not match FASTA for ',v$SNP_ID, ': ',ref_base, ' vs ',v$REF)
 coding_ref<-if(strand=='+')v$REF else complement(v$REF)
 coding_alt<-if(strand=='+')v$ALT else complement(v$ALT)
 stopifnot(substr(coding,cds_idx,cds_idx)==coding_ref)
 codon_first<-as.integer(((cds_idx-1L)%/%3L)*3L+1L)
 codon_last<-codon_first+2L
 ref_codon<-if(codon_last<=nchar(coding))substring(coding,codon_first,codon_last)else NA_character_
 alt_codon<-ref_codon
 if(!is.na(ref_codon)) substr(alt_codon,cds_idx-codon_first+1L,cds_idx-codon_first+1L)<-coding_alt
 ref_aa<-if(!is.na(ref_codon))translate(ref_codon)else NA_character_
 alt_aa<-if(!is.na(alt_codon))translate(alt_codon)else NA_character_
 effect<-if(!phase_consistent||is.na(ref_aa)||is.na(alt_aa))'needs_manual_review' else
   if(ref_aa==alt_aa)'synonymous_variant' else
   if(ref_aa!='*'&&alt_aa=='*')'stop_gained' else
   if(ref_aa=='*'&&alt_aa!='*')'stop_lost' else
   if(cds_idx<=3L&&ref_aa=='M'&&alt_aa!='M')'start_lost' else 'missense_variant'
 results[[length(results)+1L]]<-data.frame(
  SNP_ID=v$SNP_ID,gene=v$gene,transcript_id=S$transcript_id[i],
  strand=strand,genomic_REF=v$REF,genomic_ALT=v$ALT,
  fasta_ref_validated=TRUE,CDS_len=nchar(coding), CDS_position=cds_idx,
  reference_codon=ref_codon,alternate_codon=alt_codon,
  reference_AA=ref_aa,alternate_AA=alt_aa,
  reference_codon_index=as.integer((cds_idx-1L)%/%3L+1L),
  CDS_phase_consistent=phase_consistent,
  predicted_effect=effect,
  stringsAsFactors=FALSE)
}
A<-if(length(results))do.call(rbind,results)else data.frame()
dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)
write.csv(A,file.path(CFG$out,'68_WormBase_WS283_coding_consequences.csv'),row.names=FALSE)
print(A)
cat('NOTE: predicted coding consequences only; not biochemical or causal validation.\n')
