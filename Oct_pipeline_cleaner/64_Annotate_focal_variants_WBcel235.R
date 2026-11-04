# 64_Annotate_focal_variants_WBcel235.R
# Transcript-level effects require REF/ALT and matched WBcel235 annotation.
# Uses bcftools + (if available) VEP; never infers missense/splicing from gene proximity.
# Supported input is the SAME indexed CeMEE genotypes VCF used in scripts 24/26/30a.
# Edit the VCF and FASTA paths below. Requires bcftools, optional VEP.

CFG <- list(
 vcf = '/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz',
 fasta = 'WS283/20220216_c_elegans_WS283.genome.fa',  # REQUIRED for reference-allele verification
 vep_cmd = 'vep', # or full executable path. VEP with matching cache is required for consequence calls.
 use_vep = FALSE,
 out = 'selection_results/manuscript_finalization/64_variant_consequences',
 custom_gff = '/Volumes/nematode_genetics/reference_genomes/c_elegans.PRJNA13758.WS283.annotations.gff3.gz',  # optional bgzipped/tabix-indexed WBcel235 GFF3 for VEP --gff
 min_nearby_bp = 10000L
)
# All positions are WBcel235, 1-based. Origins and candidate genes are from historical analysis notes,
# not from a functional validation experiment. The impact column will be filled by VEP only.
sites <- data.frame(
 chr=c('I','I','I','IV','IV','V','V','V'),
 pos=c(7167738L,7491656L,7596359L,8892663L,10499135L,6322080L,6637695L,6826421L),
 region=c(rep('I_RC301_PX174',3),rep('IV_CB4856',2),rep('V_multilocus',3)),
 hypothesis=c('RC301_PX174_shared','RC301_private','RC301_private',
              'post_founder_candidate','confirmed_domestication_mutation',
              'confirmed_domestication_mutation','CB4507_linked','CB4507_linked'),
 candidate_gene=c('','','','','tax-6','cest-24','',''),
 stringsAsFactors=FALSE
)
site_path <- function(z) paste0(z$chr,'_',z$pos)
sites$SNP_ID <- site_path(sites)
dir.create(CFG$out,recursive=TRUE,showWarnings=FALSE)
write.csv(sites,file.path(CFG$out,'64_Frozen_focal_variant_manifest.csv'),row.names=FALSE)
if(!nzchar(Sys.which('bcftools'))) stop('bcftools unavailable. Install it before running script 64.')
if(!file.exists(CFG$vcf))stop('Set CFG$vcf to indexed CeMEE VCF. Manifest was written to ',CFG$out)
if(!nzchar(CFG$fasta)||!file.exists(CFG$fasta))
 stop('Set CFG$fasta to the EXACT WBcel235 reference used for this VCF (required REF QC).')

# Resolve chromosome label from the VCF header. No silent coordinate liftover.
headers <- suppressWarnings(system2('bcftools',c('view','-h',shQuote(CFG$vcf)),stdout=TRUE,stderr=FALSE))
contigs <- sub('^##contig=<ID=([^,>]+).*$', '\\1',grep('^##contig=<ID=',headers,value=TRUE))
resolve_contig <- function(chr){
  tests<-c(chr,paste0('chr',chr))
  h<-tests[tests%in%contigs]
  if(length(h)!=1L)stop('Cannot resolve VCF contig for ',chr,' in ',CFG$vcf)
  h
}
sites$VCF_CHROM <- vapply(sites$chr,resolve_contig,character(1))
positions <- unique(sites[c('VCF_CHROM','pos')])
intervals <- file.path(CFG$out,'64_exact_positions.tsv')
write.table(data.frame(positions$VCF_CHROM,positions$pos,positions$pos),
 intervals,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
selected <- file.path(CFG$out,'64_selected_exact_variants.vcf')
err <- file.path(CFG$out,'64_bcftools_stderr.txt')
status<-system2('bcftools',c('view','-R',shQuote(intervals),'-Ov','-o',shQuote(selected),shQuote(CFG$vcf)),
 stdout=FALSE,stderr=err)
if(status!=0L)stop('bcftools extraction failed; inspect ',err)
# Fail on incorrect reference, without rewriting the alleles.
validated <- file.path(CFG$out,'64_reference_checked.vcf')
status <- system2('bcftools',c('norm','-f',shQuote(CFG$fasta),'-c','e','-Ov','-o',
                                    shQuote(validated),shQuote(selected)),
                  stdout=FALSE,stderr=err)
if(status!=0L)stop('FASTA/VCF REF mismatch, wrong assembly or contig aliases. See ',err)
# Extract exact SNP alleles and report absent / multiallelic coordinates. Never guess REF/ALT.
fmt <- '%CHROM\\t%POS\\t%REF\\t%ALT\\n'
ln <- suppressWarnings(system2('bcftools',c('query','-f',shQuote(fmt),shQuote(validated)),
                                  stdout=TRUE,stderr=err))
if(length(ln)){
  v<-read.table(text=ln,sep='\t',stringsAsFactors=FALSE,header=FALSE,fill=TRUE)
  names(v)<-c('VCF_CHROM','pos','REF','ALT')
  v$pos<-as.integer(v$pos)
  v<-v[nchar(v$REF)==1L & nchar(v$ALT)==1L & !grepl(',',v$ALT),,drop=FALSE]
  v<-v[!duplicated(v[,c('VCF_CHROM','pos')]) &
       !duplicated(v[,c('VCF_CHROM','pos')],fromLast=TRUE),,drop=FALSE]
}else v<-data.frame(VCF_CHROM=character(),pos=integer(),REF=character(),ALT=character())
A<-merge(sites,v,by=c('VCF_CHROM','pos'),all.x=TRUE,sort=FALSE)
A<-A[match(sites$SNP_ID,A$SNP_ID),,drop=FALSE]
A$QC <- ifelse(is.na(A$REF),'MISSING_OR_MULTIALLELIC_IN_VCF','EXACT_REF_VALIDATED')
write.csv(A,file.path(CFG$out,'64_Focal_alleles_and_QC.csv'),row.names=FALSE)
if(any(A$QC!='EXACT_REF_VALIDATED'))
 warning('Some focal sites are missing/ambiguous. Do not infer their consequence from the manuscript notes.')

# Transcript consequence annotation via Ensembl VEP. Do not claim that a VEP
# consequence proves a phenotypic or fitness effect. Transcript isoforms matter.
vep<-Sys.which(CFG$vep_cmd)
if(CFG$use_vep && nzchar(vep)){
  if(!nzchar(CFG$custom_gff)){
    # Cache mode; install a caenorhabditis_elegans cache matched to local VEP.
    annotation_arg <- c('--cache','--offline','--species','caenorhabditis_elegans',
                        '--assembly','WBcel235')
  }else{
    if(!file.exists(CFG$custom_gff))stop('custom GFF missing')
    annotation_arg<-c('--offline','--gff',shQuote(CFG$custom_gff))
  }
  fields<-paste(c('Uploaded_variation','Location','Allele','Gene','Feature',
                  'Feature_type','Consequence','IMPACT','SYMBOL','BIOTYPE',
                  'EXON','INTRON','HGVSc','HGVSp','Amino_acids','Codons'),collapse=',')
  out_tsv<-file.path(CFG$out,'64_VEP_transcript_consequences.tsv')
  vep_log<-file.path(CFG$out,'64_VEP_stderr.txt')
  arg<-c('--input_file',shQuote(validated),'--output_file',shQuote(out_tsv),
         '--tab','--force_overwrite','--no_stats','--symbol','--biotype',
         '--canonical','--hgvs','--fasta',shQuote(CFG$fasta),
         '--distance',as.character(CFG$min_nearby_bp),'--fields',shQuote(fields),
         annotation_arg)
  status<-system2(vep,arg,stdout=FALSE,stderr=vep_log)
  if(status!=0L)warning('VEP failed: check ',vep_log,'; VEP cache version, FASTA and GFF/assembly must match.')
}else warning('VEP not installed or disabled. Allele QC only; no predicted protein/splice effects.')
writeLines(c('WBcel235 variant consequences — source documentation:',
 'https://static.ensembl.org/info/docs/tools/vep/index.html',
 'https://static.ensembl.org/info/docs/tools/vep/script/vep_cache.html',
 'https://parasite.wormbase.org/Caenorhabditis_elegans_prjna13758/Info/Index',
 'Variant consequence != causal phenotype/fitness effect.',
 'Review ALL isoforms, especially alternative exons in tax-6.',
 'If intronic, inspect splice proximity and transcript-specific consequences.',
 'Context-dependent upstream/downstream annotation needs regulatory data for claims of regulation.',
 'Gene in notes is a research hypothesis, NOT evidence from this script.',
 'Check whether VCF has reference-matching alleles and whether samples are the correct RIL/founder universe.',
 capture.output(sessionInfo())),file.path(CFG$out,'64_functional_annotation_interpretation.txt'))
cat('64 complete. See ',CFG$out,'\n',sep='')
