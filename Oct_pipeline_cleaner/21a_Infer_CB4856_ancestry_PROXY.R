## 21a_Infer_CB4856_ancestry_PROXY.R
## OPTIONAL approximation of local CB4856 founder ancestry from founder VCF
## genotypes and RIL genotypes already saved by script 00.
## Use documented founder-origin probabilities if available: these are superior.
## This script DOES NOT establish phased haplotype ancestry or mutation origin.

vcf_file <- '/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz'
founder_ids_file <- '18_founder_sample_ids.txt'
cb_founder_name <- 'CB4856'  # Set to the exact VCF sample name if different.
preparation_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
positions <- c(8892663L,10499135L)
flanking_window_bp <- 150000L
max_markers_per_site <- 800L
min_markers_per_RIL <- 25L
min_score_difference <- 0.07  # Difference in mismatch fraction to CB vs closest other
max_CB_mismatch <- 0.25
out_file <- '21_CB4856_local_ancestry.csv'
out_qc <- '21_CB4856_local_ancestry_PROXY_QC.csv'

if (!file.exists(vcf_file)) stop('Set vcf_file to your indexed founder/RIL VCF')
if (!file.exists(founder_ids_file)) stop('Supply 18_founder_sample_ids.txt')
if(!nzchar(Sys.which('bcftools')))stop('bcftools is required in PATH')
founders <- trimws(readLines(founder_ids_file,warn=FALSE))
founders <- founders[nzchar(founders)]
if(length(founders)!=16L||anyDuplicated(founders))stop('Expected 16 unique founder sample IDs')
if(!cb_founder_name %in% founders)
  stop('Set cb_founder_name to the exact CB4856 sample name in the founder list')
if(!file.exists(preparation_file)&&file.exists(paste0(preparation_file,'.RData')))
  preparation_file <- paste0(preparation_file,'.RData')
e <- new.env(parent=emptyenv())
load(preparation_file,envir=e)
if(!'genotype_matrix'%in%ls(e))stop('Missing RIL genotype_matrix')
G <- e$genotype_matrix
ids <- rownames(G)
if(is.null(ids)||anyDuplicated(ids))stop('Missing or duplicate RIL genotype IDs')

## Genotype parsing: alt dosage 0 / 1 / 2, missing otherwise.
gt2dosage <- function(x){
  x<-sub(':.*$','',x)
  out<-rep(NA_real_,length(x))
  out[x%in%c('0/0','0|0','0')]<-0
  out[x%in%c('0/1','1/0','0|1','1|0')]<-1
  out[x%in%c('1/1','1|1','1')]<-2
  out
}
call_bcftools <- function(pos) {
  reg<-paste0('IV:',max(1,pos-flanking_window_bp),'-',pos+flanking_window_bp)
  fmt <- '%POS\t%REF\t%ALT[\t%GT]\n'
  z <- suppressWarnings(system2('bcftools',c('query','-s',
    shQuote(paste(founders,collapse=',')),'-r',shQuote(reg),'-f',shQuote(fmt),
    shQuote(vcf_file)),stdout=TRUE,stderr=FALSE))
  if(!length(z))stop('bcftools returned no records for ',reg)
  split <- strsplit(z,'\t',fixed=TRUE)
  split <- split[lengths(split)==(3L+length(founders))]
  if(!length(split))stop('No parsable founder records at ',reg)
  rr<-do.call(rbind,split)
  list(pos=as.integer(rr[,1]),ref=rr[,2],alt=rr[,3],
       GT=matrix(gt2dosage(rr[,-(1:3),drop=FALSE]),ncol=length(founders),
                 dimnames=list(NULL,founders)))
}

anc <- matrix(NA_real_,length(ids),length(positions),
              dimnames=list(ids,paste0('CB4856_',positions)))
qcs <- list()
for(k in seq_along(positions)) {
  pos<-positions[k]
  v<-call_bcftools(pos)
  take <- v$pos!=pos & nchar(v$ref)==1 & nchar(v$alt)==1 &
    !grepl(',',v$alt,fixed=TRUE) &
    rowSums(is.na(v$GT))==0 &
    rowSums(v$GT==0 | v$GT==2)==length(founders)
  marker <- paste0('IV_',v$pos)
  take <- take & marker%in%colnames(G) &
    rowSums(v$GT != v$GT[,cb_founder_name])>0
  ii<-which(take)
  if(length(ii)>max_markers_per_site)
    ii<-ii[unique(round(seq(1,length(ii),length.out=max_markers_per_site)))]
  if(length(ii)<min_markers_per_RIL){
    warning('Only ',length(ii),' informative markers around IV:',pos)
    next
  }
  Xm <- G[,marker[ii],drop=FALSE]
  fGT <- v$GT[ii,,drop=FALSE]
  ## Each RIL's local similarity to all founders, including CB4856. Different
  ## founder haplotypes can have almost identical patterns: mark ambiguous NA.
  for(j in seq_along(ids)) {
    x <- Xm[j,]
    ok <- is.finite(x) & x>=0 & x<=2
    n <- sum(ok)
    if(n<min_markers_per_RIL){
      qcs[[length(qcs)+1L]] <- data.frame(pop_label=ids[j],position=pos,
        n_informative=n,CB_mismatch=NA_real_,best_other_mismatch=NA_real_,
        CB_minus_other=NA_real_,CB4856_status=NA_real_)
      next
    }
    diff<-colMeans(abs(sweep(fGT[ok,,drop=FALSE],1,x[ok],'-'))/2)
    d_cb<-diff[cb_founder_name]
    d_other<-min(diff[names(diff)!=cb_founder_name])
    margin<-d_cb-d_other
    val<-if(margin<=-min_score_difference && d_cb<=max_CB_mismatch) 1 else
         if(margin>=min_score_difference) 0 else NA_real_
    anc[j,k]<-val
    qcs[[length(qcs)+1L]]<-data.frame(pop_label=ids[j],position=pos,
      n_informative=n,CB_mismatch=d_cb,best_other_mismatch=d_other,
      CB_minus_other=margin,CB4856_status=val)
  }
  cat('IV:',pos,' informative markers=',length(ii),
      '; confident CB4856=',sum(anc[,k]==1,na.rm=TRUE),
      '; confident other=',sum(anc[,k]==0,na.rm=TRUE),
      '; ambiguous=',sum(is.na(anc[,k])),'\n',sep='')
}
write.csv(data.frame(pop_label=ids,anc,check.names=FALSE),out_file,
          row.names=FALSE)
if(length(qcs))write.csv(do.call(rbind,qcs),out_qc,row.names=FALSE)
cat('Created local genotype-similarity PROXY: ',out_file,'\n',
    'Inspect QC before using in script 21. Do NOT interpret the proxy as\n',
    'verified founder ancestry when nearby founder haplotypes are identical.\n',sep='')
