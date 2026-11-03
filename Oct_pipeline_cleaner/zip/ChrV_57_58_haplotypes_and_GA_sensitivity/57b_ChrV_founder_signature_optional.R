## 57b — Optional four-anchor founder genotype signatures for chr V.
## Requires the SAME imputed, biallelic reference VCF/assembly as used to build
## the RIL genotypes, plus the sixteen founder sample IDs.
## This yields four-marker allele signatures, NOT rigorous founder ancestry calls.
## Run after 57; the 57 and 58 analyses do not require this script.
options(stringsAsFactors=FALSE)
out <- 'selection_results/chrV_deep_dive/57_haplotypes'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
## Modify these two paths if necessary (from the original script 26).
vcf <- '/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz'
founder_file <- '18_founder_sample_ids.txt'
## Set TRUE only after confirming VCF ALT at all four positions corresponds
## to the genotype_matrix ALT coding and Pool-seq GWAS-ALT orientation.
confirm_identical_ALT_allele_coding <- FALSE
chrom_prefix <- ''  # set to 'chr' only if VCF uses chrV rather than V
write_out <- function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na='')
read_tab <- function(f)read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
tags_path <- file.path(out,'57_Four_anchor_marker_and_ancestral_orientation_QC.csv')
if(!file.exists(tags_path))stop('Run 57 first; missing ',tags_path)
if(!file.exists(vcf))stop('VCF not found: edit `vcf` at top of 57b')
if(!file.exists(founder_file))stop('Missing founder sample ID file: ',founder_file)
if(!nzchar(Sys.which('bcftools')))stop('bcftools required on PATH')
tags <- read_tab(tags_path)
if(nrow(tags)!=4L||anyDuplicated(tags$SNP_ID))stop('Expected four unambiguous frozen anchor tags')
founders <- trimws(readLines(founder_file,warn=FALSE));founders<-founders[nzchar(founders)]
if(length(founders)!=16L||anyDuplicated(founders))stop('Expect exactly 16 unique founders')
vcf_ids <- suppressWarnings(system2('bcftools',c('query','-l',shQuote(vcf)),stdout=TRUE,stderr=FALSE))
if(!all(founders%in%vcf_ids))stop('Founders missing from VCF: ',paste(setdiff(founders,vcf_ids),collapse=', '))
founders <- vcf_ids[vcf_ids%in%founders] # actual VCF order
positions <- as.integer(sub('^.*_','',tags$SNP_ID))
posfile<-file.path(out,'57b_Four_exact_VCF_positions.tsv')
write.table(data.frame(CHROM=paste0(chrom_prefix,'V'),POS=sort(positions),END=sort(positions)),
  posfile,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
fmt <- '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'
err <- file.path(out,'57b_bcftools_stderr.txt')
cmd <- c('query','-s',shQuote(paste(founders,collapse=',')),
         '-f',shQuote(fmt),'-R',shQuote(posfile),shQuote(vcf))
lines<-suppressWarnings(system2('bcftools',cmd,stdout=TRUE,stderr=err))
status <- attr(lines,'status')
if((!is.null(status)&&status!=0L)||!length(lines))
  stop('bcftools failed or no records; inspect ',err,' (contigs/VCF/index)')
rows <- strsplit(lines,'\t',fixed=TRUE)
rows <- rows[lengths(rows)==(4L+length(founders))]
if(!length(rows))stop('No complete biallelic founder VCF records')
Z <- do.call(rbind,rows)
ids <- paste(toupper(sub('^chr','',Z[,1],ignore.case=TRUE)),as.integer(Z[,2]),sep='_')
if(anyDuplicated(ids)||!setequal(ids,tags$SNP_ID))
  stop('VCF has duplicate, absent or unexpected records for anchor tag positions')
Z <- Z[match(tags$SNP_ID,ids),,drop=FALSE]
if(any(nchar(Z[,3])!=1L|nchar(Z[,4])!=1L|grepl(',',Z[,4],fixed=TRUE)))
  stop('Multiallelic or non-SNP VCF variant among focal four positions')
gt <- sub(':.*$','',Z[,-(1:4),drop=FALSE]);dim(gt) <- c(4L,length(founders))
parse_gt<-function(x){d<-rep(NA_integer_,length(x));d[x%in%c('0|0','0/0','0')]<-0L
 d[x%in%c('1|1','1/1','1')]<-2L;d[x%in%c('0|1','1|0','0/1','1/0')]<-1L;d}
G <- matrix(parse_gt(as.vector(gt)),4L,length(founders),dimnames=list(tags$SNP_ID,founders))
write_out(data.frame(tags[,c('group','SNP_ID','rare_allele','A00_ALT_freq')],
                     VCF_REF=Z[,3],VCF_ALT=Z[,4],n_founder_called=rowSums(!is.na(G))),
          '57b_Founder_exact_marker_alleles_and_coverage.csv')
calls <- data.frame(founder=founders,stringsAsFactors=FALSE)
for(i in seq_len(nrow(tags)))calls[[paste0('ALT_dosage_',tags$group[i])]] <- as.integer(G[i,])
write_out(calls,'57b_Founder_raw_VCF_ALT_dosages.csv')
if(!confirm_identical_ALT_allele_coding){
  writeLines(c('FOUNDER SIGNATURES NOT JOINED TO RIL HAPLOTYPES:',
    'confirm_identical_ALT_allele_coding is FALSE.',
    'Check REF/ALT identities for each SNP in genotype_matrix and Pool-seq input',
    'against VCF_REF/VCF_ALT, then explicitly set TRUE if they are consistent.'),
    file.path(out,'57b_Orientation_verification_required.txt'))
  message('Raw founder ALT dosages saved; no RIL pattern comparison without allele verification.')
} else {
  b <- matrix(NA_integer_,nrow=length(founders),ncol=nrow(tags),
              dimnames=list(founders,as.character(tags$group)))
  for(i in seq_len(nrow(tags))){
    v <- as.numeric(G[i,]);ok <- is.finite(v)&v%in%c(0,2)
    b[ok,i] <- if(tags$rare_allele[i]=='ALT')as.integer(v[ok]==2)
                else as.integer(v[ok]==0)
  }
  pattern <- apply(b,1,function(z)if(anyNA(z))NA_character_ else paste0(z,collapse=''))
  found <- data.frame(founder=founders,rare_allele_4tag_signature=pattern,
                     fully_called=!is.na(pattern),stringsAsFactors=FALSE)
  write_out(found,'57b_Founder_four_marker_signature_NOT_ancestry.csv')
  ril <- read_tab(file.path(out,'57_RIL_four_marker_haplotype_frequencies.csv'))
  matching<-merge(ril,aggregate(founder~rare_allele_4tag_signature,found,FUN=function(z)paste(z,collapse=';')),
                  by.x='rare_haplotype_4tag',by.y='rare_allele_4tag_signature',all.x=TRUE,sort=FALSE)
  names(matching)[names(matching)=='founder']<-'matching_founder_signatures_not_ancestry'
  write_out(matching,'57b_RIL_patterns_matching_16_founder_signatures_PROXY.csv')
  warning('Founder four-tag matches are NOT phased local ancestry estimates and ',
          'cannot distinguish shared ancestral backgrounds with identical four marker calls.')
}
cat('57b complete. Founder VCF signatures are descriptive, not ancestry calls.\n')
