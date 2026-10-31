## 29_Rank_additional_QTL_case_studies.R
## Rank NEW R5 QTL case studies for follow-up WITHOUT refitting GWAS or
## treating multiple selected SNPs within one LD cluster as independent.
## Explicitly separate (a) GWAS strength, (b) drift-aware selection,
## (c) direction agreement, (d) full versus retained-GA robustness,
## (e) multiple traits, (f) novelty relative to historical loci.

qtl_file <- 'historical_QTL_R5/23_R5_all_cluster_boundaries_and_associations.csv'
survivor_file <- 'selection_results/survivors_only_Ne1000/20_Survivor_cluster_slope_summary.csv'
legacy_file <- '23_Historical_region_SNPs.csv'
regions_file <- 'Oct_pipeline_cleaner/26_Regions_and_R5_groups.csv'
historical_match_file <- 'historical_QTL_R5/23_Historical_region_to_R5_QTL_candidates.csv'
out_dir <- 'historical_QTL_R5/deep_dive/discovery'
q_primary <- .05
q_strong <- .01
q_selection <- .01
min_common_members <- 3L
min_distance_to_historical_bp <- 200000L
min_distance_between_novel_centers_bp <- 400000L
n_shortlist <- 25L
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
req<-c(qtl_file,legacy_file,regions_file)
if(any(!file.exists(req)))stop('Missing: ',paste(req[!file.exists(req)],collapse=','))
z<-read.csv(qtl_file,stringsAsFactors=FALSE,check.names=FALSE)
z$group<-as.character(z$group)
if(anyDuplicated(z$group))stop('Duplicate R5 groups')
traits<-c('SF','SB','FS','FB','BS','BF','area.F','beta')
if(!all(paste0('q_',traits)%in%names(z)))stop('R5 GWAS columns missing')
if(!all(c('CHROM','start','end','R5_centroid','n_slope_common','min_q_common',
          'slope_median_common')%in%names(z)))stop('Missing cluster selection summaries; run scripts 23, 09f, 11')
qmat<-as.matrix(z[,paste0('q_',traits),drop=FALSE]);storage.mode(qmat)<-'double'
z$n_GWAS_traits_q05<-rowSums(qmat<q_primary,na.rm=TRUE)
z$n_GWAS_traits_q01<-rowSums(qmat<q_strong,na.rm=TRUE)
z$GWAS_traits_q05<-apply(qmat,1,function(v)
 paste(traits[is.finite(v)&v<q_primary],collapse=';'))
z$minimum_GWAS_q<-apply(qmat,1,function(v)
  if(any(is.finite(v)))min(v,na.rm=TRUE)else NA_real_)
z$passes_GWAS<-(is.finite(z$q_beta)&z$q_beta<q_primary)|
 (is.finite(z[['q_area.F']])&z[['q_area.F']]<q_primary)|
 z$n_GWAS_traits_q05>=2L
z$selected_full<-is.finite(z$min_q_common)&z$min_q_common<q_selection
z$adequate_common_cov<-is.finite(z$n_slope_common)&z$n_slope_common>=min_common_members
z$beta_direction_full<-ifelse(is.finite(z$effect_beta)&is.finite(z$slope_median_common),
 sign(z$effect_beta*z$slope_median_common),NA_real_)
z$area_direction_full<-ifelse(is.finite(z[['effect_area.F']])&is.finite(z$slope_median_common),
 sign(z[['effect_area.F']]*z$slope_median_common),NA_real_)
if(file.exists(survivor_file)){
 s<-read.csv(survivor_file,stringsAsFactors=FALSE)
 if(anyDuplicated(as.character(s$group)))stop('Duplicate survivor groups')
 i<-match(z$group,as.character(s$group))
 for(nm in intersect(c('n_survivor_common','n_selected_survivor',
                      'slope_survivor_median_common','any_selected_survivor'),names(s)))
    z[[nm]]<-s[[nm]][i]
 z$selected_retained<-is.finite(z$n_selected_survivor)&z$n_selected_survivor>0
 z$direction_stable_full_vs_retained<-ifelse(
    is.finite(z$slope_median_common)&is.finite(z$slope_survivor_median_common),
    sign(z$slope_median_common)==sign(z$slope_survivor_median_common),NA)
}else{
 z$selected_retained<-NA
 z$direction_stable_full_vs_retained<-NA
 warning('No script-20 survivor results: robustness score unavailable')
}
## Nominated and historical locus exclusion based on the supplied original SNP
## spans, not on broad R5 cluster bounding boxes that can be discontinuous.
old<-read.csv(legacy_file,stringsAsFactors=FALSE)
old$CHROM<-toupper(sub('^chr','',old$CHROM,ignore.case=TRUE))
intervals<-do.call(rbind,lapply(split(old,old$historical_region),function(v)
 data.frame(historical_region=v$historical_region[1],CHROM=v$CHROM[1],
            start=min(v$POS),end=max(v$POS))))
intervals$start<-pmax(1,intervals$start-min_distance_to_historical_bp)
intervals$end<-intervals$end+min_distance_to_historical_bp
z$near_historical_region<-FALSE
z$historical_region_nearby<-''
for(i in seq_len(nrow(z))){
    m<-intervals$CHROM==z$CHROM[i] &
      z$R5_centroid[i]>=intervals$start & z$R5_centroid[i]<=intervals$end
    z$near_historical_region[i]<-any(m)
    if(any(m))z$historical_region_nearby[i]<-paste(intervals$historical_region[m],collapse=';')
}
focus<-read.csv(regions_file,stringsAsFactors=FALSE)
focus<-focus[toupper(as.character(focus$include))=='TRUE',,drop=FALSE]
z$already_nominated_group<-z$group%in%as.character(focus$group)
z$already_historical_R5_group<-FALSE
if(file.exists(historical_match_file)){
 hh<-read.csv(historical_match_file,stringsAsFactors=FALSE)
 z$already_historical_R5_group<-z$group%in%as.character(hh$group)
}
## Finite evidence. Do not award significance if no SNP passed common filter.
z$score<-0L
z$score<-z$score+ifelse(is.finite(z$q_beta)&z$q_beta<.01,3L,
                      ifelse(is.finite(z$q_beta)&z$q_beta<.05,2L,0L))
z$score<-z$score+ifelse(is.finite(z[['q_area.F']])&z[['q_area.F']]<.01,2L,
                      ifelse(is.finite(z[['q_area.F']])&z[['q_area.F']]<.05,1L,0L))
z$score<-z$score+ifelse(z$n_GWAS_traits_q05>=3,2L,
                      ifelse(z$n_GWAS_traits_q05>=2,1L,0L))
z$score<-z$score+ifelse(z$selected_full,3L,
                       ifelse(is.finite(z$min_q_common)&z$min_q_common<.05,1L,0L))
z$score<-z$score+ifelse(z$selected_retained%in%TRUE,2L,0L)
z$score<-z$score+ifelse(z$direction_stable_full_vs_retained%in%TRUE,1L,0L)
z$score<-z$score+ifelse(z$beta_direction_full==1 &
                        is.finite(z$q_beta)&z$q_beta<.05,1L,0L)
z$score<-z$score+ifelse(z$area_direction_full==1 &
                        is.finite(z[['q_area.F']])&z[['q_area.F']]<.05,1L,0L)
z$score<-z$score+ifelse(z$adequate_common_cov,1L,0L)
z$strong_joint_case <- z$passes_GWAS & z$adequate_common_cov & z$selected_full &
 (!is.na(z$direction_stable_full_vs_retained)&z$direction_stable_full_vs_retained)
z$reason<-paste0('GWAS traits q<0.05: ',z$GWAS_traits_q05,
 '; selected full SNP: ',z$selected_full,
 '; selected retained SNP: ',z$selected_retained,
 '; slope sign stable: ',z$direction_stable_full_vs_retained,
 '; n common SNPs: ',z$n_slope_common)
z$selection_calibration_note<-'Drift-aware BH at SNP level is provisional in the extreme WF tail; overlapping SNPs are not independent'
z$novel_unstudied<-!z$near_historical_region & !z$already_nominated_group &
   !z$already_historical_R5_group
allrank<-z[order(!z$novel_unstudied,-z$score,z$minimum_GWAS_q,na.last=TRUE),]
write.csv(allrank,file.path(out_dir,'29_All_R5_QTL_evidence_score.csv'),row.names=FALSE)
## Focus on novel loci with BOTH a phenotype association and measurable
## common-frequency response; do not force a candidate on each chromosome.
new<-z[z$novel_unstudied & z$passes_GWAS & z$adequate_common_cov,,drop=FALSE]
new<-new[order(-new$score,new$min_q_common,new$minimum_GWAS_q,na.last=TRUE),]
keep<-logical(nrow(new));seen<-list()
for(i in seq_len(nrow(new))){
 chr<-new$CHROM[i];pos<-new$R5_centroid[i]
 pp<-seen[[chr]]
 if(is.null(pp)||all(abs(pos-pp)>=min_distance_between_novel_centers_bp)){
   keep[i]<-TRUE;seen[[chr]]<-c(pp,pos)
 }
}
short<-head(new[keep,,drop=FALSE],n_shortlist)
write.csv(short,file.path(out_dir,'29_Novel_QTL_shortlist_spatially_thinned.csv'),row.names=FALSE)
selected<-z[z$strong_joint_case & z$passes_GWAS,,drop=FALSE]
selected<-selected[order(-selected$score,selected$minimum_GWAS_q),]
write.csv(selected,file.path(out_dir,'29_GWAS_selection_direction_consistent_cases.csv'),row.names=FALSE)
cat('Total R5 clusters: ',nrow(z),'\nNovel candidates after spatial thinning: ',nrow(short),
    '\nGenome-wide strong joint cases: ',nrow(selected),'\n',sep='')
cat('Selection is inferred from common SNPs, not necessarily the QTL causal SNP.\n')
