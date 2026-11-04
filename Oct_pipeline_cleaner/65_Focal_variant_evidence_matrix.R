# 65_Focal_variant_evidence_matrix.R
# Integrates frozen SNP manifest, allele QC, VEP (if available), exact R5
# membership (if available), and current genome-wide GWAS-evolution outcomes.
# It deliberately distinguishes a tagged R5 region from an identified causal variant.
# Run after scripts 63/64. Base R only.

CFG <- list(
  out='selection_results/manuscript_finalization/65_focal_variant_evidence',
  variant_dir='selection_results/manuscript_finalization/64_variant_consequences',
  phase_csv='data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv',
  prediction_candidates=c('selection_results/drift_oneG0_nullfit/R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
                          '45_All_R5_clusters_GWAS_evolution_predictions.csv'),
  gwas_candidates=c('selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
                    '11_Clusters_selection_and_GWAS_09f.csv')
)
first_file <- function(x){f<-x[file.exists(x)];if(length(f))f[1]else NA_character_}
D<-CFG$variant_dir
manifest<-file.path(D,'64_Frozen_focal_variant_manifest.csv')
alleles<-file.path(D,'64_Focal_alleles_and_QC.csv')
if(!file.exists(manifest))stop('Run 64 first; missing ',manifest)
a<-read.csv(manifest,stringsAsFactors=FALSE)
a <- a[!duplicated(a$SNP_ID),,drop=FALSE]
if(file.exists(alleles)){
  b<-read.csv(alleles,stringsAsFactors=FALSE)
  if(anyDuplicated(b$SNP_ID))stop('Duplicate allele QC SNP_ID')
  ind<-match(a$SNP_ID,b$SNP_ID)
  for(k in c('REF','ALT','QC'))if(k%in%names(b))a[[k]]<-b[[k]][ind]
}

# Exact group membership MUST derive from the frozen R5 phase registry, never
# physical-nearest-cluster matching or a guessed historical cluster number.
a$R5_group<-NA_character_;a$membership_status<-'UNVERIFIED_NO_PHASE_REGISTRY'
if(file.exists(CFG$phase_csv)){
  p<-read.csv(CFG$phase_csv,stringsAsFactors=FALSE)
  if(!all(c('rs','group')%in%names(p)))stop('Unexpected phase registry format')
  if(anyDuplicated(p$rs))stop('Duplicate SNP rs identifier in phase registry')
  j<-match(a$SNP_ID,p$rs)
  a$R5_group<-as.character(p$group[j])
  a$membership_status<-ifelse(!is.na(j),'EXACT_MEMBER','NOT_IN_R5')
}

pfile<-first_file(CFG$prediction_candidates)
if(!is.na(pfile)){
  pred<-read.csv(pfile,stringsAsFactors=FALSE)
  pred<-pred[pred$scenario=='primary_survivors',,drop=FALSE]
  if(anyDuplicated(pred$group))stop('Duplicated main R5 prediction groups')
  ind<-match(a$R5_group,as.character(pred$group))
  kcols<-intersect(c('primary_slope_median_phase','pred_full','LOCO_full',
                     'n_eligible','selection_p_MC_upper_simultaneous',
                     'selection_significant_BY_05'),names(pred))
  for(k in kcols)a[[paste0('R5_',k)]]<-pred[[k]][ind]
  a$R5_prediction_group_resolved <- !is.na(ind)
}

# VEP transcript-level table is often many-to-one by variant. Preserve all
# transcript consequences separately rather than collapsing to a best impact.
vepfile<-file.path(D,'64_VEP_transcript_consequences.tsv')
vepdetails<-NULL
if(file.exists(vepfile)){
 lines<-readLines(vepfile,warn=FALSE)
 j<-grep('^#Uploaded_variation\\t|^Uploaded_variation\\t',lines)[1L]
 if(is.na(j)) warning('Could not find VEP header; inspect ',vepfile)
 else{
  vepdetails<-read.delim(text=paste(sub('^#','',lines[j]),paste(lines[-seq_len(j)],collapse='\n'),sep='\n'),
                         stringsAsFactors=FALSE,check.names=FALSE)
  if(!'Location'%in%names(vepdetails))stop('VEP lacks Location field')
  loc<-sub('^(chr)?([^:]+):([0-9]+)(-.*)?$','\\2_\\3',vepdetails$Location)
  vepdetails$SNP_ID<-loc
  vepdetails<-vepdetails[vepdetails$SNP_ID%in%a$SNP_ID,,drop=FALSE]
  if(nrow(vepdetails)){
    # VEP uses multiple effects separated by commas or ampersands; this is a
    # descriptive collapsed overview and NEVER supersedes transcript table.
    a$n_VEP_transcripts <- as.integer(table(factor(vepdetails$SNP_ID,levels=a$SNP_ID)))
    a$VEP_effects_all_transcripts <- vapply(a$SNP_ID,function(id){
        s<-unique(as.character(vepdetails$Consequence[vepdetails$SNP_ID==id]))
        paste(s[nzchar(s)&s!='-'],collapse=';')},character(1))
    a$VEP_genes_all_transcripts <- vapply(a$SNP_ID,function(id){
        s<-unique(as.character(vepdetails$SYMBOL[vepdetails$SNP_ID==id]))
        paste(s[nzchar(s)&s!='-'],collapse=';')},character(1))
  }
 }
}
# Normalize NA to explicit unavailable; important for downstream manuscript copy.
if(!'n_VEP_transcripts'%in%names(a))a$n_VEP_transcripts<-NA_integer_
if(!'VEP_effects_all_transcripts'%in%names(a))a$VEP_effects_all_transcripts<-NA_character_
if(!'VEP_genes_all_transcripts'%in%names(a))a$VEP_genes_all_transcripts<-NA_character_

dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)
write.csv(a,file.path(CFG$out,'65_Regional_variant_evidence_for_manuscript.csv'),row.names=FALSE)
if(!is.null(vepdetails))write.csv(vepdetails,file.path(CFG$out,'65_All_transcript_consequences.csv'),row.names=FALSE)
writeLines(c('Read this before using the evidence table:',
  'REF/ALT alleles MUST pass script 64 reference verification.',
  'Protein coding consequences depend on the transcript and isoform.',
  'Trait associations belong to linked R5 clusters, not automatically the exact mutation.',
  'A mutation inside a known gene is not necessarily the cause of the observed phenotype.',
  'Post-domestication origin of IV:10499135 and V:6322080 is supplied by prior independent history.',
  'Do not infer a mutation origin from non-detection in only 16 sampled founders.',
  'If the phase registry is absent, exact membership is not established.',
  capture.output(sessionInfo())),file.path(CFG$out,'65_Interpretation_and_QC.txt'))
print(a[,intersect(c('SNP_ID','candidate_gene','QC','R5_group','membership_status',
                    'VEP_effects_all_transcripts'),names(a)),drop=FALSE])
