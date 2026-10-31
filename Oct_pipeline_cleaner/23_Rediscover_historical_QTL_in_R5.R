## 23_Rediscover_historical_QTL_in_R5.R
## Locate explicitly supplied historical SNP lists in the CURRENT R5 pipeline.
## Never force a historical region to be one current QTL: R5 clusters are
## defined independently by LD and signed GWAS effects.
## Run from project root after scripts 02, 05, 09f and 11.
## No GWAS or SNP-selection model is rerun; all outputs are descriptive.

## ------------------------------- Settings -----------------------------
legacy_file <- '23_Historical_region_SNPs.csv'
r5_snp_file <- 'txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt'
r5_cluster_file <- 'RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R5_with_residuals.RData'
all_snp_gwas_file <- 'GridLMM_results_NaCl/All_Gwas.RData' # optional
phase_file <- 'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv' # optional
selection_snp_file <- 'selection_results/drift_oneG0_nullfit/09f_SNP_drift_Ne1000.csv' # optional
selection_cluster_file <- 'selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv' # optional
out_dir <- 'historical_QTL_R5'
nearest_window_bp <- 100000L # nearest SNP is a PROXY, never an exact match
pdf_width <- 10
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
chrom_order <- c('I','II','III','IV','V','X')
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

norm_chr <- function(x) toupper(sub('^CHR', '', as.character(x), ignore.case=TRUE))
nval <- function(x) as.numeric(as.character(x))
safemin <- function(x) {x <- x[is.finite(x)]; if(length(x))min(x) else NA_real_}
safemax <- function(x) {x <- x[is.finite(x)]; if(length(x))max(x) else NA_real_}

legacy <- read.csv(legacy_file, stringsAsFactors=FALSE, check.names=FALSE)
if(!all(c('historical_region','CHROM','historical_trait','POS') %in% names(legacy)))
  stop('Legacy file needs historical_region,CHROM,historical_trait,POS')
legacy$CHROM <- norm_chr(legacy$CHROM)
legacy$POS <- as.integer(legacy$POS)
legacy$SNP_ID <- paste(legacy$CHROM,legacy$POS,sep='_')
legacy <- unique(legacy)
if(anyNA(legacy$POS) || anyDuplicated(paste(legacy$historical_region,legacy$SNP_ID)))
  stop('Invalid or duplicate historical SNP assignments')
if(any(!legacy$CHROM %in% chrom_order))stop('Unexpected historical chromosome')
if(any(vapply(split(legacy$CHROM,legacy$historical_region),function(z)length(unique(z)),integer(1))!=1))
  stop('Historical region labels must describe one chromosome each')

snp <- read.table(r5_snp_file,header=TRUE,sep='\t',check.names=FALSE,
                  stringsAsFactors=FALSE,quote='',comment.char='')
req <- c('rs','CHROM','POS','group')
if(!all(req %in% names(snp)))stop('R5 SNP clustering file lacks: ',paste(setdiff(req,names(snp)),collapse=', '))
snp$CHROM <- norm_chr(snp$CHROM)
snp$POS <- as.integer(snp$POS)
snp$rs <- as.character(snp$rs)
snp$group <- as.character(snp$group)
if(anyDuplicated(snp$rs) || anyNA(snp$group))stop('Duplicated/unassigned R5 SNP IDs')

ce <- new.env(parent=globalenv())
load(r5_cluster_file,envir=ce)
if(!'all_GWAS_clusters' %in% ls(ce))stop('R5 cluster object not found')
qtl <- ce$all_GWAS_clusters
if(!all(c('group','chr','ps') %in% names(qtl)))stop('Unexpected script-05 QTL fields')
qtl$group <- as.character(qtl$group)
qtl$chr <- norm_chr(qtl$chr)
if(anyDuplicated(qtl$group))stop('Duplicate R5 cluster group IDs')
if(!setequal(unique(snp$group),qtl$group))
  warning('SNP and cluster GWAS group universes differ; unmatched rows will be flagged')

## 1. Exact legacy SNP -> R5 mapping; NO implicit nearest assignment.
legacy$in_R5 <- legacy$SNP_ID %in% snp$rs
idx <- match(legacy$SNP_ID,snp$rs)
legacy$R5_group_exact <- snp$group[idx]
legacy$R5_raw_component <- if('cluster'%in%names(snp))as.character(snp$cluster[idx]) else NA_character_
legacy$R5_SNP_p_beta <- if('pvalue_T8'%in%names(snp))nval(snp$pvalue_T8[idx]) else NA_real_
legacy$R5_SNP_p_area <- if('pvalue_T7'%in%names(snp))nval(snp$pvalue_T7[idx]) else NA_real_

## 2. Optional current SNP-wide mapping: a legacy SNP can be tested but not
## included in the R5 significant-SNP clustering universe.
legacy$in_all_GWAS <- NA
legacy$GWAS_p_beta <- NA_real_
legacy$GWAS_q_beta <- NA_real_
legacy$GWAS_p_area <- NA_real_
legacy$GWAS_q_area <- NA_real_
if(file.exists(all_snp_gwas_file)) {
  ee <- new.env(parent=globalenv())
  load(all_snp_gwas_file,envir=ee)
  if('all_GWAS'%in%ls(ee) && all(c('rs','pvalue_T7','pvalue_T8')%in%names(ee$all_GWAS))) {
    aa <- ee$all_GWAS
    m <- match(legacy$SNP_ID, as.character(aa$rs))
    legacy$in_all_GWAS <- !is.na(m)
    for (pair in list(c('area','7'),c('beta','8'))) {
      pv <- nval(aa[[paste0('pvalue_T',pair[2])]])
      qv <- p.adjust(pv,method='BH')
      legacy[[paste0('GWAS_p_',pair[1])]] <- pv[m]
      legacy[[paste0('GWAS_q_',pair[1])]] <- qv[m]
    }
  } else warning('All_Gwas.RData is present but does not contain expected all_GWAS fields')
}

## Nearest R5 SNP within the SAME chromosome; LABEL as proximity only.
legacy$nearest_R5_rs <- NA_character_
legacy$nearest_R5_distance_bp <- NA_real_
legacy$nearest_R5_group_PROXY <- NA_character_
for(chr in unique(legacy$CHROM)) {
  j <- which(legacy$CHROM==chr)
  z <- snp[snp$CHROM==chr,,drop=FALSE]
  if(!nrow(z))next
  z <- z[order(z$POS),,drop=FALSE]
  ii <- findInterval(legacy$POS[j],z$POS)
  left <- pmax(1L,pmin(nrow(z),ii))
  right <- pmax(1L,pmin(nrow(z),ii+1L))
  choose <- ifelse(abs(z$POS[left]-legacy$POS[j])<=abs(z$POS[right]-legacy$POS[j]),left,right)
  dd <- abs(z$POS[choose]-legacy$POS[j])
  legacy$nearest_R5_rs[j] <- z$rs[choose]
  legacy$nearest_R5_distance_bp[j] <- dd
  legacy$nearest_R5_group_PROXY[j] <- ifelse(dd<=nearest_window_bp,z$group[choose],NA_character_)
}
legacy$mapping_status <- ifelse(legacy$in_R5,'exact_R5_SNP',
                         ifelse(!is.na(legacy$nearest_R5_group_PROXY),
                                'nearby_R5_SNP_only_NOT_member','no_R5_SNP_within_window'))
write.csv(legacy,file.path(out_dir,'23_Historical_SNP_status.csv'),row.names=FALSE)

## 3. Cluster boundaries calculated from TRUE member SNPs, not from the
## centroid reported by script 03. A cluster may have disjoint membership.
split_group <- split(snp,snp$group)
group_bounds <- do.call(rbind,lapply(split_group,function(z) data.frame(
  group=as.character(z$group[1]),CHROM=z$CHROM[1],
  start=min(z$POS),end=max(z$POS),n_R5_member_SNPs=nrow(z),
  stringsAsFactors=FALSE)))
rownames(group_bounds) <- NULL
j <- match(group_bounds$group,qtl$group)
group_bounds$R5_cluster_ID <- if('rs'%in%names(qtl))as.character(qtl$rs[j]) else NA_character_
group_bounds$R5_centroid <- nval(qtl$ps[j])
for(tt in traits){
  key <- paste0(tt,'_NaCl')
  eff <- paste0('beta_',key)
  p <- paste0('p_value_REML_',key)
  q <- paste0('p_value_REML.cor_',key)
  group_bounds[[paste0('effect_',tt)]] <- if(eff%in%names(qtl))nval(qtl[[eff]][j]) else NA_real_
  group_bounds[[paste0('p_',tt)]] <- if(p%in%names(qtl))nval(qtl[[p]][j]) else NA_real_
  group_bounds[[paste0('q_',tt)]] <- if(q%in%names(qtl))nval(qtl[[q]][j]) else NA_real_
}

## 4. Optional SNP- and cluster-level Pool-seq selection evidence.
## Historical positions absent from the common-frequency analysis are retained
## and explicitly flagged as not tested, not as selection-negative.
legacy$selection_p0 <- legacy$selection_freq_slope_ALT <- legacy$selection_q_common <- NA_real_
legacy$selection_tested_common <- FALSE
if(file.exists(selection_snp_file)){
  sel <- read.csv(selection_snp_file,stringsAsFactors=FALSE)
  if(all(c('ID','p0_hat','slope_freq','p_drift_BH_common','common_G0')%in%names(sel))){
    m <- match(legacy$SNP_ID,sel$ID)
    legacy$selection_p0 <- nval(sel$p0_hat[m])
    legacy$selection_freq_slope_ALT <- nval(sel$slope_freq[m])
    legacy$selection_q_common <- nval(sel$p_drift_BH_common[m])
    legacy$selection_tested_common <- !is.na(m)&!is.na(sel$common_G0[m])&sel$common_G0[m]&
      is.finite(legacy$selection_q_common)
    write.csv(legacy,file.path(out_dir,'23_Historical_SNP_status.csv'),row.names=FALSE)
  }
}
if(file.exists(selection_cluster_file)){
  sc <- read.csv(selection_cluster_file,stringsAsFactors=FALSE)
  if('group'%in%names(sc) && !anyDuplicated(as.character(sc$group))){
    k <- match(group_bounds$group,as.character(sc$group))
    for(nm in intersect(c('n_slope_common','n_selected_common_q01',
                          'any_selected_common_q01','min_q_common',
                          'slope_median_common','slope_mean_common'),names(sc)))
      group_bounds[[nm]] <- sc[[nm]][k]
  }
}
write.csv(group_bounds,file.path(out_dir,'23_R5_all_cluster_boundaries_and_associations.csv'),row.names=FALSE)

## 5. Region-to-R5 cluster candidates: preserve THREE evidence tiers.
## exact: one or more old SNPs are R5 members.
## interval_only: no exact old SNP, but a CURRENT R5 member falls between the
## extreme positions of a historical region. This does not mean same causal QTL.
## nearest_only: no R5 SNP in historical interval, best R5 SNP within 100 kb.
regions <- split(legacy,legacy$historical_region)
entries <- list();reports <- list(); cidx <- 0L
for(rr in names(regions)){
  z <- regions[[rr]]
  chr <- unique(z$CHROM)
  a <- min(z$POS);b <- max(z$POS)
  exact_groups <- unique(na.omit(z$R5_group_exact))
  within <- snp[snp$CHROM==chr&snp$POS>=a&snp$POS<=b,,drop=FALSE]
  within_groups <- unique(as.character(within$group))
  proposed <- union(exact_groups,within_groups)
  nearest_groups <- character()
  if(!length(proposed)) nearest_groups <- unique(na.omit(z$nearest_R5_group_PROXY))
  proposed <- union(proposed,nearest_groups)
  reports[[rr]] <- data.frame(historical_region=rr,CHROM=chr,
    historical_trait=unique(z$historical_trait)[1],
    legacy_start=a,legacy_end=b,n_legacy_SNPs=nrow(z),
    n_exact_R5_SNPs=sum(z$in_R5),n_in_current_all_GWAS=sum(z$in_all_GWAS,na.rm=TRUE),
    n_R5_SNPs_inside_interval=nrow(within),
    n_R5_clusters_exact=length(exact_groups),
    n_R5_clusters_interval=length(within_groups),
    n_candidate_R5_clusters=length(proposed),
    region_mapping_status=if(length(exact_groups))'exact_old_SNP_overlap' else
      if(length(within_groups))'R5_SNPs_inside_legacy_span_only' else
      if(length(nearest_groups))'nearest_R5_proxy_only' else 'unmapped',
    stringsAsFactors=FALSE)
  for(g in proposed){
    q <- group_bounds[group_bounds$group==g,,drop=FALSE]
    if(nrow(q)!=1L){warning('R5 group missing from metadata: ',g);next}
    old_exact <- z$POS[!is.na(z$R5_group_exact)&z$R5_group_exact==g]
    member_pos <- snp$POS[snp$group==g]
    min_dist <- min(vapply(z$POS,function(p) min(abs(p-member_pos)),numeric(1)))
    cidx <- cidx+1L
    q$historical_region <- rr
    q$historical_trait <- unique(z$historical_trait)[1]
    q$legacy_start <- a;q$legacy_end <- b
    q$n_exact_old_SNPs <- length(old_exact)
    q$n_current_SNPs_within_legacy_span <- sum(snp$group==g & snp$POS>=a & snp$POS<=b)
    q$nearest_old_to_member_bp <- min_dist
    q$match_tier <- if(g%in%exact_groups)'exact_legacy_SNP_member' else
      if(g%in%within_groups)'current_R5_SNP_in_legacy_span_only' else
      'nearest_only_NOT_overlap'
    q$target_trait_q <- q[[paste0('q_',unique(z$historical_trait)[1])]]
    q$target_trait_effect <- q[[paste0('effect_',unique(z$historical_trait)[1])]]
    q$target_trait_q05 <- is.finite(q$target_trait_q)&q$target_trait_q<0.05
    entries[[cidx]]<-q
  }
}
summary <- do.call(rbind,reports)
rownames(summary)<-NULL
write.csv(summary,file.path(out_dir,'23_Historical_region_summary.csv'),row.names=FALSE)
if(length(entries)) {
  candidates <- do.call(rbind,entries)
  candidates <- candidates[order(match(candidates$historical_region,names(regions)),
    match(candidates$match_tier,c('exact_legacy_SNP_member',
      'current_R5_SNP_in_legacy_span_only','nearest_only_NOT_overlap')),
    candidates$target_trait_q,candidates$nearest_old_to_member_bp),,drop=FALSE]
  rownames(candidates)<-NULL
  write.csv(candidates,file.path(out_dir,'23_Historical_region_to_R5_QTL_candidates.csv'),row.names=FALSE)
  members <- snp[snp$group %in% candidates$group,,drop=FALSE]
  members$matching_historical_regions <- vapply(members$group,function(g)
    paste(unique(candidates$historical_region[candidates$group==g]),collapse=';'),character(1))
  if(file.exists(phase_file)){
    pha <- read.csv(phase_file,stringsAsFactors=FALSE)
    if(all(c('rs','phase')%in%names(pha)) && !anyDuplicated(pha$rs))
      members$R5_phase <- pha$phase[match(members$rs,pha$rs)]
  }
  write.csv(members,file.path(out_dir,'23_Candidate_R5_cluster_member_SNPs.csv'),row.names=FALSE)
} else {
  candidates <- data.frame(historical_region=character(),group=character())
  write.csv(candidates,file.path(out_dir,'23_Historical_region_to_R5_QTL_candidates.csv'),row.names=FALSE)
}

## 6. One page per historical region: visual audit of the redefinition.
pdf(file.path(out_dir,'23_Historical_to_R5_region_maps.pdf'),width=pdf_width,height=6)
for(rr in names(regions)){
  z <- regions[[rr]];chr <- unique(z$CHROM)
  a<-min(z$POS);b<-max(z$POS)
  xlim <- c(max(0,a-100000),b+100000)
  cands <- candidates[candidates$historical_region==rr,,drop=FALSE]
  n<-nrow(cands)
  op<-par(mar=c(4.3,12,3,2))
  plot(NA,xlim=xlim,ylim=c(0,max(3,n+1)),xlab=paste0('Chromosome ',chr,' position (bp)'),
    ylab='',yaxt='n',bty='n',main=paste0(rr,' | ',n,' candidate R5 clusters'))
  abline(v=c(a,b),col='gray85',lty=3)
  rug(z$POS,col='gray35',lwd=1)
  if(n){
    for(i in seq_len(n)){
      d<-cands[i,,drop=FALSE]
      col<-if(d$match_tier=='exact_legacy_SNP_member')'#277e77' else
        if(d$match_tier=='current_R5_SNP_in_legacy_span_only')'#b79054' else '#999999'
      segments(d$start,i,d$end,i,lwd=4,col=col)
      points(d$R5_centroid,i,pch=16,col=col,cex=.75)
    }
    labs<-paste0('R5 group ',cands$group,' | q(',unique(z$historical_trait)[1],')=',
                  ifelse(is.finite(cands$target_trait_q),signif(cands$target_trait_q,2),'NA'))
    axis(2,at=seq_len(n),labels=labs,las=1,cex.axis=.57)
  }
  mtext('Ticks = historical SNP positions; green = exact overlap; gold = interval only',
        side=1,line=3,cex=.67)
  par(op)
}
dev.off()

cat('\nHistorical regions:',nrow(summary),' (including 2 chromosome IV candidate mutations)\n')
print(summary[,c('historical_region','n_legacy_SNPs','n_exact_R5_SNPs',
  'n_candidate_R5_clusters','region_mapping_status')],row.names=FALSE)
cat('Outputs in: ',out_dir,'\n',sep='')
cat('CAUTION: exact overlap, span overlap and nearest proximity are distinct evidence tiers.\n')
