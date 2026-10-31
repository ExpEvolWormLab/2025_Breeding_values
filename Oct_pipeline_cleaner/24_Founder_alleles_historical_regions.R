## 24_Founder_alleles_historical_regions.R
## Extract exact historical-region SNP genotypes from the indexed VCF,
## classify which founders carry each allele (especially PX174 and RC301),
## and link founder allele patterns to CURRENT R5 clusters.
## Run 23 first. This is an allele-sharing screen, NOT phased ancestry inference.

## ------------------------------- Settings -----------------------------
vcf_file <- '/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz'
founder_ids_file <- '18_founder_sample_ids.txt' # 16 VCF sample IDs, one per line
legacy_file <- '23_Historical_region_SNPs.csv'
match_file <- 'historical_QTL_R5/23_Historical_SNP_status.csv' # optional
out_dir <- 'historical_QTL_R5/founder_alleles'
px_name <- 'PX174'  # matches one founder ID substring
rc_name <- 'RC301'  # matches one founder ID substring
chrom_order <- c('I','II','III','IV','V','X')
## If VCF uses chrI, chrII, etc., set this to 'chr' instead of ''.
vcf_chr_prefix <- ''
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(vcf_file))stop('Configure vcf_file to your indexed VCF')
if(!file.exists(founder_ids_file))stop('Supply the same 16 founders used in script 18')
if(!nzchar(Sys.which('bcftools')))stop('bcftools must be on PATH')
founders <- trimws(readLines(founder_ids_file,warn=FALSE))
founders <- founders[nzchar(founders)]
if(length(founders)!=16L||anyDuplicated(founders))stop('Expected 16 distinct founder IDs')
legacy <- read.csv(legacy_file,stringsAsFactors=FALSE)
legacy$CHROM <- toupper(sub('^chr','',legacy$CHROM,ignore.case=TRUE))
legacy$POS <- as.integer(legacy$POS)
legacy$SNP_ID <- paste(legacy$CHROM,legacy$POS,sep='_')
orig_legacy <- legacy
legacy <- unique(legacy[,c('CHROM','POS','SNP_ID')])
legacy <- legacy[order(match(legacy$CHROM,chrom_order),legacy$POS),,drop=FALSE]

## Read the VCF sample-order header; all-sample query retains this order.
vs <- suppressWarnings(system2('bcftools',c('query','-l',shQuote(vcf_file)),
                                  stdout=TRUE,stderr=FALSE))
vs <- vs[nzchar(vs)]
if(!all(founders%in%vs))stop('Founders not found in VCF: ',
                            paste(setdiff(founders,vs),collapse=', '))
founder_col <- match(founders,vs)
sel <- which(vapply(c(px_name,rc_name),function(nm)sum(grepl(nm,founders,fixed=TRUE))==1L,logical(1)))
if(length(sel)!=2L)warning('Check PX174 and RC301 spelling in founder list')
px <- grep(px_name,founders,fixed=TRUE,value=TRUE)
rc <- grep(rc_name,founders,fixed=TRUE,value=TRUE)

## bcftools -R needs CHROM,POS,POS and requires the VCF to be indexed.
regions <- file.path(out_dir,'24_bcf_query_positions.tsv')
write.table(data.frame(CHROM=paste0(vcf_chr_prefix,legacy$CHROM),
                       START=legacy$POS,END=legacy$POS),regions,
            quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
fmt <- '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'
err_file <- file.path(out_dir,'24_bcftools_stderr.txt')
q <- suppressWarnings(system2('bcftools',c('query','-R',shQuote(regions),
           '-f',shQuote(fmt),shQuote(vcf_file)),stdout=TRUE,stderr=err_file))
if(!length(q))stop('No VCF records returned. Inspect ',err_file,
  '; check VCF path/index, assembly and chromosome prefix')
if(!is.null(attr(q,'status'))&&attr(q,'status')!=0L)
  stop('bcftools query failed. Inspect ',err_file)

parse_gt <- function(z) {
  z <- sub(':.*$','',z)
  out <- rep(NA_real_,length(z))
  out[z%in%c('0/0','0|0','0')]<-0
  out[z%in%c('0/1','1/0','0|1','1|0')]<-1
  out[z%in%c('1/1','1|1','1')]<-2
  out
}
records <- list(); i<-0L
for(line in q){
  x <- strsplit(line,'\t',fixed=TRUE)[[1]]
  if(length(x)!=4L+length(vs)){
    warning('Skipped VCF row with unexpected number of sample genotypes')
    next
  }
  chr <- toupper(sub('^chr','',x[1],ignore.case=TRUE))
  pos <- suppressWarnings(as.integer(x[2]))
  id <- paste(chr,pos,sep='_')
  if(!id%in%legacy$SNP_ID)next
  valid_snp <- nchar(x[3])==1L && nchar(x[4])==1L &&
    x[3]%in%c('A','C','G','T') && x[4]%in%c('A','C','G','T')
  dosage <- if(valid_snp)parse_gt(x[4L+founder_col]) else rep(NA_real_,length(founders))
  if(valid_snp && all(is.na(dosage)))warning('No parsable founder genotypes at ',id)
  i<-i+1L
  records[[i]]<-data.frame(SNP_ID=id,CHROM=chr,POS=pos,REF=x[3],ALT=x[4],
    founder_id=founders,alt_dosage=dosage,
    genotype_status=if(valid_snp)'biallelic_SNP' else 'not_simple_SNP',
    stringsAsFactors=FALSE)
}
if(length(records))geno <- do.call(rbind,records) else geno <- data.frame(
  SNP_ID=character(),CHROM=character(),POS=integer(),REF=character(),ALT=character(),
  founder_id=character(),alt_dosage=numeric(),genotype_status=character())
write.csv(geno,file.path(out_dir,'24_Historical_SNP_founder_genotypes_long.csv'),row.names=FALSE)

classes <- list();n <-0L
for(id in unique(legacy$SNP_ID)){
  z <- geno[geno$SNP_ID==id,,drop=FALSE]
  if(!nrow(z)){
    n<-n+1L; classes[[n]]<-data.frame(SNP_ID=id,REF=NA_character_,ALT=NA_character_,
      n_founder_ALT=NA_integer_,founder_ALT_names=NA_character_,PX174_dosage=NA_real_,
      RC301_dosage=NA_real_,PX_RC_class='not_in_VCF',n_VCF_records=0L)
    next
  }
  if(length(unique(paste(z$REF,z$ALT)))>1L){
    warning('Multiple REF/ALT alleles for ',id,'; marking ambiguous')
  }
  zz <- z[seq_len(min(nrow(z),length(founders))),,drop=FALSE]
  g <- setNames(zz$alt_dosage,zz$founder_id)
  a<-if(length(px)==1L)g[px] else NA_real_
  b<-if(length(rc)==1L)g[rc] else NA_real_
  cls<-if(any(!is.finite(c(a,b))))'missing_or_non_biallelic' else
    if(a==2&&b==2)'PX174_AND_RC301_ALT' else
    if(a==2&&b==0)'PX174_only_ALT' else
    if(a==0&&b==2)'RC301_only_ALT' else
    if(a==0&&b==0)'neither_ALT' else 'heterozygous_or_other'
  n<-n+1L
  classes[[n]]<-data.frame(SNP_ID=id,REF=zz$REF[1],ALT=zz$ALT[1],
    n_founder_ALT=sum(g==2,na.rm=TRUE),
    founder_ALT_names=paste(names(g)[which(g==2)],collapse=';'),
    PX174_dosage=unname(a),RC301_dosage=unname(b),PX_RC_class=cls,
    n_VCF_records=nrow(z)/length(founders),stringsAsFactors=FALSE)
}
summary <- do.call(rbind,classes)
meta <- orig_legacy[,c('SNP_ID','historical_region','historical_trait')]
summary <- merge(meta,summary,by='SNP_ID',all.x=TRUE,sort=FALSE)
if(file.exists(match_file)){
  mm <- read.csv(match_file,stringsAsFactors=FALSE)
  keep <- intersect(c('SNP_ID','historical_region','R5_group_exact','mapping_status',
                      'GWAS_q_beta','GWAS_q_area','selection_q_common',
                      'selection_freq_slope_ALT'),names(mm))
  if(all(c('SNP_ID','historical_region')%in%keep)){
    summary <- merge(summary,mm[,keep,drop=FALSE],
      by=c('SNP_ID','historical_region'),all.x=TRUE,sort=FALSE)
  }
}
write.csv(summary,file.path(out_dir,'24_Historical_SNP_founder_signature_classes.csv'),row.names=FALSE)
cc <- aggregate(rep.int(1L,nrow(summary)),by=summary[c('historical_region','PX_RC_class')],FUN=sum)
names(cc)[3] <- 'n_historical_SNPs'
write.csv(cc,file.path(out_dir,'24_Founder_signature_counts_by_region.csv'),row.names=FALSE)
cat('\nFounders PX174:',paste(px,collapse=';'),' RC301:',paste(rc,collapse=';'),'\n')
cat('Historical positions:',nrow(legacy),' ; found simple VCF entries:',
    length(unique(geno$SNP_ID)), '\n')
print(cc,row.names=FALSE)
cat('All founder allele signatures are descriptive; shared alleles do not prove founder ancestry.\n')
