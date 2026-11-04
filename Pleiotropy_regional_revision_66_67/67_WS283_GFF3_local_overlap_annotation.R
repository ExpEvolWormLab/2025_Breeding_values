# 67_WS283_GFF3_local_overlap_annotation.R
# Approximate local sequence-ontology evidence from your OWN WS283 GFF3.
# No VEP installation, no bgzip, no index: works with ordinary .gff3.gz.
# This DOES NOT compute missense, synonymous or splice consequences.
# Strictly local feature-overlap and proximity annotations, even when web VEP is down.
# Run from CeMEE project root. Base R only.

CFG <- list(
  gff = 'WS283/c_elegans.WS283.VEP.gff3.gz', # edit with full path
  variants = 'selection_results/manuscript_finalization/65_focal_variant_evidence/65_Regional_variant_evidence_for_manuscript.csv',
  out = 'selection_results/manuscript_finalization/67_WS283_local_annotations',
  flank_bp = 20000L,
  chunk_lines = 30000L
)

if(!file.exists(CFG$gff))stop('Set CFG$gff to actual WS283 .gff3.gz path')
if(!file.exists(CFG$variants))stop('Missing script 65 eight-variant table: ',CFG$variants)
V<-read.csv(CFG$variants,stringsAsFactors=FALSE)
stopifnot(all(c('chr','pos','SNP_ID','REF','ALT')%in%names(V)))
V<-V[!duplicated(V$SNP_ID),,drop=FALSE]
if(anyNA(V$REF)|anyNA(V$ALT))stop('REF and ALT alleles must be exact and reference-validated')
if('QC'%in%names(V) && any(V$QC!='EXACT_REF_VALIDATED'))stop('Variant reference QC failed')
V$pos<-as.integer(V$pos)
# Stream GFF3, filter by location first to avoid holding whole annotation in memory.
con<-if(grepl('\\.gz$',CFG$gff,ignore.case=TRUE))gzfile(CFG$gff,'rt') else file(CFG$gff,'rt')
collect<-list();n<-0L
repeat {
  ss<-readLines(con,n=CFG$chunk_lines,warn=FALSE)
  if(!length(ss))break
  ss<-ss[!grepl('^#|^$',ss)]
  if(!length(ss))next
  chr_names<-unique(V$chr)
  ss<-ss[vapply(strsplit(ss,'\t',fixed=TRUE),function(z){
    if(length(z)<5)return(FALSE)
    ch<-sub('^chr','',z[1]);if(!ch%in%chr_names)return(FALSE)
    st<-suppressWarnings(as.integer(z[4]));en<-suppressWarnings(as.integer(z[5]))
    if(!is.finite(st)||!is.finite(en))return(FALSE)
    any(V$chr==ch & V$pos>=st-CFG$flank_bp & V$pos<=en+CFG$flank_bp)
  },logical(1))]
  if(!length(ss))next
  spl<-strsplit(ss,'\t',fixed=TRUE)
  spl<-spl[lengths(spl)>=9L]
  if(!length(spl))next
  d<-as.data.frame(do.call(rbind,lapply(spl,function(z)z[1:9])),stringsAsFactors=FALSE)
  names(d)<-c('chr','source','type','start','end','score','strand','phase','attributes')
  d$chr<-sub('^chr','',d$chr)
  d$start<-as.integer(d$start);d$end<-as.integer(d$end)
  collect[[length(collect)+1L]]<-d
  n<-n+nrow(d)
}
close(con)
if(!length(collect))stop('No overlapping GFF3 features. Check WS283 assembly and chromosome labels')
F<-do.call(rbind,collect);rownames(F)<-NULL
# Extract GFF3 tags without guessing transcript-to-gene identifiers.
get_attr<-function(z,key){
 out<-rep(NA_character_,length(z))
 for(i in seq_along(z)){
  x<-strsplit(z[i],';',fixed=TRUE)[[1]]
  y<-x[startsWith(x,paste0(key,'='))]
  if(length(y))out[i]<-sub(paste0('^',key,'='),'',y[1])
 }
 out
}
for(k in c('ID','Parent','Name','gene','gene_id','locus','Alias','description'))F[[k]]<-get_attr(F$attributes,k)
# Include each exact overlapping feature and nearest annotated gene intervals.
OV<-list();NEAR<-list()
for(i in seq_len(nrow(V))){
 v<-V[i,]
 f<-F[F$chr==v$chr & F$start<=v$pos & F$end>=v$pos,,drop=FALSE]
 if(nrow(f)){
   f$SNP_ID<-v$SNP_ID
   f$ref<-v$REF;f$alt<-v$ALT
   OV[[length(OV)+1L]]<-f
 }
 # GFF3 gene type may vary: include gene, pseudogene and annotated RNA genes.
 g<-F[F$chr==v$chr & grepl('gene|pseudogene',F$type,ignore.case=TRUE),,drop=FALSE]
 if(nrow(g)){
   g$distance_bp<-pmax(g$start-v$pos,v$pos-g$end,0L)
   g<-g[order(g$distance_bp,g$start),,drop=FALSE]
   g<-g[seq_len(min(5L,nrow(g))),,drop=FALSE]
   g$SNP_ID<-v$SNP_ID
   NEAR[[length(NEAR)+1L]]<-g
 }
}
O<-if(length(OV))do.call(rbind,OV)else F[FALSE,]
G<-if(length(NEAR))do.call(rbind,NEAR)else F[FALSE,]
if(!'SNP_ID'%in%names(O))O$SNP_ID<-character(nrow(O))
if(!'SNP_ID'%in%names(G))G$SNP_ID<-character(nrow(G))
if(!'distance_bp'%in%names(G))G$distance_bp<-numeric(nrow(G))
dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)
write.csv(O,file.path(CFG$out,'67_exact_overlap_GFF3_features.csv'),row.names=FALSE)
write.csv(G,file.path(CFG$out,'67_nearest_annotated_genes.csv'),row.names=FALSE)
rows<-lapply(seq_len(nrow(V)),function(i){
 v<-V[i,];f<-O[O$SNP_ID==v$SNP_ID,,drop=FALSE]
 types<-unique(f$type)
 classification <- if('CDS'%in%types)'CDS-overlap (protein change not computed)' else
   if('exon'%in%types)'exonic/non-CDS-overlap' else
   if(any(grepl('transcript|mRNA|RNA',types,ignore.case=TRUE)))'transcribed region, outside annotated exon' else
   if(any(grepl('gene',types,ignore.case=TRUE)))'genic; transcript structure not resolved' else
   'no GFF3 feature directly overlaps'
 g<-G[G$SNP_ID==v$SNP_ID,,drop=FALSE]
 data.frame(SNP_ID=v$SNP_ID,REF=v$REF,ALT=v$ALT,
  overlap_types=paste(types,collapse=';'),broad_class=classification,
  nearest_gene_id=if(nrow(g))g$ID[1]else NA_character_,
  nearest_gene_name=if(nrow(g))g$Name[1]else NA_character_,
  nearest_gene_distance_bp=if(nrow(g))g$distance_bp[1]else NA_integer_,
  nearest_gene_strand=if(nrow(g))g$strand[1]else NA_character_,
  stringsAsFactors=FALSE)
})
Z<-do.call(rbind,rows)
write.csv(Z,file.path(CFG$out,'67_WS283_broad_variant_annotation.csv'),row.names=FALSE)
writeLines(c(
 'This tool uses the actual WS283 GFF3 annotation supplied locally.',
 'It does not resolve alternative transcript-specific consequences, CDS codon changes, synonymous/missense calls or strand-specific splice effects.',
 'A variant inside a gene is not proof of causality. Do not infer functional impact from nearest gene.',
 'When VEP returns, compare against the identical WS283 assembly/transcript versions.'
),file.path(CFG$out,'67_interpretation_limits.txt'))
print(Z)
