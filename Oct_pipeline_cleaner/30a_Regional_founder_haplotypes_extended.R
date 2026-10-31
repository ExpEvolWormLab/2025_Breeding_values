## Auto-generated extended version for phase 11 local-QTL competition analyses.
## 26_Regional_founder_haplotypes.R
## Founder allele signatures and RIL local-genotype similarity for nominated R5 QTLs.
## Founder similarity is a PROXY, not a phased ancestry call.
## Run from project root; copy 26_*.csv beside this script.

vcf_file <- '/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz'
founder_ids_file <- '18_founder_sample_ids.txt'
preparation_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
group_file <- '30_Regions_and_R5_groups_extended.csv'
focal_file <- '30_Focal_variants_extended.csv'
r5_file <- 'txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt'
qtl_file <- 'historical_QTL_R5/23_R5_all_cluster_boundaries_and_associations.csv'
legacy_file <- 'historical_QTL_R5/23_Historical_SNP_status.csv'
output_dir <- 'historical_QTL_R5/deep_dive_extended'
chrom_prefix <- '' # use 'chr' only if VCF contigs are chrI, chrII, ...
window_bp <- 100000L
max_panel_markers <- 200L
min_panel_markers <- 12L
min_margin <- 0.06
max_mismatch <- 0.25
max_exact_per_group <- 80L
set.seed(20261029)
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

if(!file.exists(vcf_file))stop('Edit vcf_file to point to the indexed VCF')
if(!nzchar(Sys.which('bcftools')))stop('bcftools is required on PATH')
reqfiles<-c(founder_ids_file,preparation_file,group_file,focal_file,r5_file,qtl_file)
if(any(!file.exists(reqfiles)))stop('Missing project files: ',paste(reqfiles[!file.exists(reqfiles)],collapse=', '))
founders <- trimws(readLines(founder_ids_file,warn=FALSE));founders<-founders[nzchar(founders)]
if(length(founders)!=16L||anyDuplicated(founders))stop('Expected exactly 16 unique VCF founder sample IDs')
vcfs <- suppressWarnings(system2('bcftools',c('query','-l',shQuote(vcf_file)),stdout=TRUE,stderr=FALSE))
if(!all(founders%in%vcfs))stop('Some founder sample IDs are missing from the VCF: ',paste(setdiff(founders,vcfs),collapse=','))
## bcftools emits samples in VCF order when -s is supplied with this ordering;
## use the actual VCF header order to avoid any hidden column swaps.
founders <- vcfs[vcfs%in%founders]
if(length(founders)!=16L)stop('VCF founder-order resolution failed')
groups<-read.csv(group_file,stringsAsFactors=FALSE)
groups<-groups[toupper(as.character(groups$include))=='TRUE',,drop=FALSE]
groups$group<-as.character(groups$group)
focal<-read.csv(focal_file,stringsAsFactors=FALSE)
r5<-read.table(r5_file,header=TRUE,sep='\t',check.names=FALSE,stringsAsFactors=FALSE,quote='',comment.char='')
r5$group<-as.character(r5$group);r5$POS<-as.integer(r5$POS)
qtl<-read.csv(qtl_file,stringsAsFactors=FALSE,check.names=FALSE)
qtl$group<-as.character(qtl$group)
if(anyDuplicated(qtl$group))stop('Duplicate R5 group metadata')
if(!all(groups$group%in%qtl$group))stop('Nominated R5 groups missing from script 23 output')
if(!all(groups$group%in%r5$group))stop('Nominated R5 groups missing from current R5 membership')
legacy<-if(file.exists(legacy_file))read.csv(legacy_file,stringsAsFactors=FALSE) else NULL

z<-new.env(parent=baseenv());load(preparation_file,envir=z)
if(!'genotype_matrix'%in%ls(z))stop('Missing genotype_matrix in preparation RData')
G<-z$genotype_matrix
if(!length(rownames(G))||!length(colnames(G)))stop('RIL genotype matrix needs RIL IDs and SNP IDs')

parse_gt<-function(x){x<-sub(':.*$','',x);v<-rep(NA_real_,length(x));
 v[x%in%c('0/0','0|0','0')]<-0;v[x%in%c('0/1','1/0','0|1','1|0')]<-1;
 v[x%in%c('1/1','1|1','1')]<-2;v}
query_vcf <- function(region=NULL,positions_file=NULL,section='unknown'){
 fmt<-'%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'
 args<-c('query','-s',shQuote(paste(founders,collapse=',')))
 if(!is.null(region))args<-c(args,'-r',shQuote(region))
 if(!is.null(positions_file))args<-c(args,'-R',shQuote(positions_file))
 err<-file.path(output_dir,paste0('26_bcf_stderr_',section,'.txt'))
 lines<-suppressWarnings(system2('bcftools',c(args,'-f',shQuote(fmt),shQuote(vcf_file)),
                                      stdout=TRUE,stderr=err))
 st<-attr(lines,'status');if(!is.null(st)&&st!=0L)
   stop('bcftools failed in ',section,'. Read ',err,'; verify contigs, index and VCF path.')
 if(!length(lines))return(data.frame())
 rows<-strsplit(lines,'\t',fixed=TRUE)
 rows<-rows[lengths(rows)==4L+length(founders)]
 if(!length(rows))return(data.frame())
 dat<-do.call(rbind,rows)
 chr<-toupper(sub('^chr','',dat[,1],ignore.case=TRUE))
 pos<-as.integer(dat[,2]);ref<-dat[,3];alt<-dat[,4]
 valid<-!is.na(pos)&nchar(ref)==1L&nchar(alt)==1L&
     ref%in%c('A','C','G','T')&alt%in%c('A','C','G','T')
 if(!any(valid))return(data.frame())
 dat<-dat[valid,,drop=FALSE]
 GT<-matrix(parse_gt(dat[,-(1:4),drop=FALSE]),ncol=length(founders),
             dimnames=list(NULL,founders))
 out<-data.frame(CHROM=chr[valid],POS=pos[valid],REF=ref[valid],ALT=alt[valid],
                 SNP_ID=paste(chr[valid],pos[valid],sep='_'),stringsAsFactors=FALSE)
 for(k in seq_along(founders))out[[founders[k]]]<-GT[,k]
 out
}

## Select the explicitly supplied historical SNPs, requested focal mutations,
## and a small representative panel of current R5 SNPs for each nominated group.
selected<-list();tagrows<-list();ii<-0L
for(j in seq_len(nrow(groups))){
 x<-groups[j,];ss<-r5[r5$group==x$group,,drop=FALSE]
 ## Complete mapping from the eight phenotype labels to the current R5 SNP GWAS columns.
 pcols <- c(SF='pvalue_T1',SB='pvalue_T2',FS='pvalue_T3',
            FB='pvalue_T4',BS='pvalue_T5',BF='pvalue_T6',
            'area.F'='pvalue_T7',beta='pvalue_T8')
 target <- as.character(x$target_trait)[1]
 if(!target %in% names(pcols))
   stop('Unknown target_trait: ',target,' for R5 group ',x$group,
        '. Expected: ',paste(names(pcols),collapse=', '))
 pvcol <- unname(pcols[target])
 if(!pvcol %in% names(ss))
   stop('GWAS p-value column ',pvcol,' is absent for group ',x$group)
 if(!nrow(ss)) stop('R5 group ',x$group,' has no member SNPs')
 pval <- suppressWarnings(as.numeric(ss[[pvcol]]))
 pval[!is.finite(pval)] <- Inf
 ## Prefer exact old SNPs when available, then minimum SNP p-value.
 oldrs<-if(!is.null(legacy))legacy$SNP_ID[legacy$R5_group_exact==x$group & !is.na(legacy$R5_group_exact)]else character()
 eligible<-which(ss$rs%in%oldrs)
 cand<-if(length(eligible))eligible[order(pval[eligible],na.last=TRUE)] else order(pval,na.last=TRUE)
 tag<-as.character(ss$rs[cand[1]])
 tagrows[[j]]<-data.frame(region=x$region,CHROM=x$CHROM,group=x$group,
   target_trait=x$target_trait,tag_SNP_ID=tag,
   tag_SNP_p=as.numeric(pval[match(tag,ss$rs)]),
   R5_cluster_ID=qtl$R5_cluster_ID[match(x$group,qtl$group)],
   R5_start=qtl$start[match(x$group,qtl$group)],
   R5_end=qtl$end[match(x$group,qtl$group)],stringsAsFactors=FALSE)
 old<-which(ss$rs%in%oldrs)
 more<-unique(round(seq(1,nrow(ss),length.out=min(nrow(ss),25L))))
 take<-unique(c(old,more,cand[seq_len(min(10L,length(cand)))]))
 if(length(take)>max_exact_per_group){
   take<-unique(c(cand[1],head(old,40L),head(more,20L)))
   take<-head(take,max_exact_per_group)
 }
 ii<-ii+1L
 selected[[ii]]<-data.frame(region=x$region,group=x$group,CHROM=x$CHROM,
   POS=ss$POS[take],SNP_ID=ss$rs[take],from_R5=TRUE,stringsAsFactors=FALSE)
}
selected<-do.call(rbind,selected)
for(j in seq_len(nrow(focal)))selected<-rbind(selected,data.frame(
  region=focal$region[j],group=NA_character_,CHROM=focal$CHROM[j],
  POS=as.integer(focal$POS[j]),SNP_ID=paste(focal$CHROM[j],focal$POS[j],sep='_'),
  from_R5=FALSE,stringsAsFactors=FALSE))
selected<-unique(selected)
tags<-do.call(rbind,tagrows)
write.csv(tags,file.path(output_dir,'26_R5_model_tag_SNPs.csv'),row.names=FALSE)
write.csv(selected,file.path(output_dir,'26_Selected_founder_panel_positions.csv'),row.names=FALSE)
qpositions<-unique(selected[,c('CHROM','POS')]);qpositions<-qpositions[order(qpositions$CHROM,qpositions$POS),]
positions_path<-file.path(output_dir,'26_exact_VCF_positions.tsv')
write.table(data.frame(CHROM=paste0(chrom_prefix,qpositions$CHROM),
 POS=qpositions$POS,END=qpositions$POS),positions_path,
 quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
exact<-query_vcf(positions_file=positions_path,section='exact')
if(!nrow(exact))stop('No exact records. Check chromosome prefix or VCF assembly')
if(anyDuplicated(exact$SNP_ID))warning('Multiple VCF alleles at same position; omitting ambiguous sites')
exact<-exact[!duplicated(exact$SNP_ID)&!duplicated(exact$SNP_ID,fromLast=TRUE),]
all_exact<-merge(selected,exact,by=c('CHROM','POS','SNP_ID'),all.x=TRUE,sort=FALSE)
fn<-founders
all_exact$n_founders_ALT<-rowSums(all_exact[,fn,drop=FALSE]==2,na.rm=TRUE)
all_exact$n_founders_called<-rowSums(!is.na(all_exact[,fn,drop=FALSE]))
all_exact$ALT_founder_names<-apply(all_exact[,fn,drop=FALSE],1,function(v)
 paste(names(v)[which(!is.na(v)&v==2)],collapse=';'))
all_exact$founder_signature<-ifelse(all_exact$n_founders_called==0,'not_in_VCF',
 ifelse(all_exact$n_founders_ALT==0,'ALT_absent_from_founders',
 ifelse(all_exact$n_founders_ALT==1,'ALT_founder_private','ALT_shared_founders')))
write.csv(all_exact,file.path(output_dir,'26_Founder_allele_signature_all_sites.csv'),row.names=FALSE)

## Marker-rich local windows for RIL-to-founder similarity. Exclude all focal
## mutation SNPs AND candidate tag SNPs from ancestry classification to avoid
## conditioning on the exact variant being tested.
block<-list();qc<-list();scores<-list()
excluded<-unique(c(focal$POS,paste0('_',tags$tag_SNP_ID)))
for(j in seq_len(nrow(tags))){
 g<-tags[j,];center<-round((g$R5_start+g$R5_end)/2)
 ## Use current cluster centroid where possible; for broad clusters use the
 ## nominated SNP tag so a distant bounding box cannot define local ancestry.
 center<-as.integer(sub('^.*_','',g$tag_SNP_ID))
 reg<-paste0(chrom_prefix,g$CHROM,':',max(1L,center-window_bp),'-',center+window_bp)
 v<-query_vcf(region=reg,section=paste0('panel_',g$group))
 if(nrow(v)){
   v<-v[!duplicated(v$SNP_ID)&!duplicated(v$SNP_ID,fromLast=TRUE),,drop=FALSE]
   v<-v[v$SNP_ID%in%colnames(G)&!v$SNP_ID%in%c(focal$CHROM,paste0(focal$CHROM,'_',focal$POS),tags$tag_SNP_ID),,drop=FALSE]
   gt<-as.matrix(v[,fn,drop=FALSE]);storage.mode(gt)<-'double'
   valid<-if(nrow(v))rowSums(is.na(gt))==0L &
     rowSums(gt==0 | gt==2)==length(fn) &
     apply(gt,1,function(z)length(unique(z))>1) else logical(0)
   v<-v[valid,,drop=FALSE]
 }
 if(nrow(v)>max_panel_markers){
   ## evenly spaced; preserve positional coverage of the genomic background
   v<-v[unique(round(seq(1,nrow(v),length.out=max_panel_markers))),,drop=FALSE]
 }
 n<-nrow(v)
 status<-if(n<min_panel_markers)'too_few_informative_markers'else'proxy_only'
 qc[[j]]<-data.frame(region=g$region,group=g$group,CHROM=g$CHROM,
    panel_center=center,window_bp=window_bp,n_informative_markers=n,status=status)
 if(n<min_panel_markers){warning('Insufficient founder markers for group ',g$group,' (',n,')');next}
 GG<-as.matrix(G[,v$SNP_ID,drop=FALSE]);storage.mode(GG)<-'double'
 ff<-as.matrix(v[,fn,drop=FALSE]);storage.mode(ff)<-'double'
 best<-rep(NA_character_,nrow(GG));second<-best;margin<-rep(NA_real_,nrow(GG))
 ncall<-integer(nrow(GG));bestdiff<-rep(NA_real_,nrow(GG))
 for(i in seq_len(nrow(GG))){
    val<-GG[i,];ok<-is.finite(val)&val>=0&val<=2
    ncall[i]<-sum(ok)
    if(ncall[i]<min_panel_markers)next
    dd<-colMeans(abs(sweep(ff[ok,,drop=FALSE],1,val[ok],'-'))/2)
    ord<-order(dd);bestdiff[i]<-dd[ord[1]]
    margin[i]<-dd[ord[2]]-dd[ord[1]]
    if(bestdiff[i]<=max_mismatch && margin[i]>=min_margin){
      best[i]<-names(dd)[ord[1]];second[i]<-names(dd)[ord[2]]
    }
 }
 scores[[length(scores)+1L]]<-data.frame(pop_label=rownames(GG),region=g$region,
    group=g$group,panel_center=center,n_markers=ncall,
    best_founder_proxy=best,second_founder_proxy=second,
    mismatch_best=bestdiff,margin_best_over_second=margin,
    confidently_classified=!is.na(best),stringsAsFactors=FALSE)
 cat('group ',g$group,': informative markers=',n,
     ', confident RILs=',sum(!is.na(best)),'/',length(best),'\n',sep='')
}
write.csv(do.call(rbind,qc),file.path(output_dir,'26_Local_founder_panel_QC.csv'),row.names=FALSE)
if(length(scores)){
 score<-do.call(rbind,scores)
 write.csv(score,file.path(output_dir,'26_RIL_founder_similarity_PROXIES.csv'),row.names=FALSE)
 counts<-as.data.frame(table(score$region,score$group,score$best_founder_proxy,useNA='ifany'))
 names(counts)<-c('region','group','best_founder_proxy','n_RIL')
 write.csv(counts,file.path(output_dir,'26_Local_founder_similarity_counts.csv'),row.names=FALSE)
}
## Per nominated group, summarize founding ALT across selected variants.
fm<-aggregate(cbind(n_founders_ALT,n_founders_called)~region+group,
  data=all_exact[!is.na(all_exact$group),,drop=FALSE],FUN=median,na.rm=TRUE)
write.csv(fm,file.path(output_dir,'26_Selected_marker_founder_summary.csv'),row.names=FALSE)
cat('Script 26 complete. PROXY founder identities require inspection before any causal claims.\n')
## Release the full RIL genotype matrix before the next script is sourced.
rm(G,z);invisible(gc())
