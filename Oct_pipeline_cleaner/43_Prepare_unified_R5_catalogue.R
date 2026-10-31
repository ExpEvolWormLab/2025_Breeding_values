## 43_Prepare_unified_R5_catalogue.R
## Run FROM THE PROJECT ROOT after corrected scripts 08, 09f, 11, 13, 33.
## Every R5 GWAS cluster is retained in the catalogue. No MAF restriction.
## ONE pre-evolution, phenotype-based tag per testable cluster will be used
## for calibrated selection testing. ALL eligible SNPs form prediction slope.
source('Oct_pipeline_cleaner/43_Unified_config_and_helpers.R')
cfg<-uf_config
uf_check(c(cfg$pool_file,cfg$orientation_file,cfg$phase_file,cfg$snp_gwas_file,
           cfg$r5_gwas_file,cfg$weights_file,cfg$status_file,cfg$existing_snp_primary))
out<-cfg$prefix;dir.create(out,recursive=TRUE,showWarnings=FALSE)
source('Oct_pipeline_cleaner/09f_drift_functions.R')
source('Oct_pipeline_cleaner/33_Frequency_drift_helpers.R')
source('Oct_pipeline_cleaner/43x_Integrate_rawSYNC_R5_SNPs.R')
p<-uf_pool_with_R5_raw_recovery(cfg)
sets<-uf_status_and_samples(p,cfg)
uf_write(sets$status,file.path(out,'43_GA_status_used.csv'))
uf_write(p$pop[sets$primary,,drop=FALSE],file.path(out,'43_Primary_sample_manifest.csv'))

phase<-uf_read(cfg$phase_file)
rs_gwas<-read.table(cfg$snp_gwas_file,header=TRUE,sep='\t',stringsAsFactors=FALSE,
                    check.names=FALSE,quote='',comment.char='')
qtl<-uf_read(cfg$r5_gwas_file)
ss<-uf_read(cfg$existing_snp_primary)
ss<-uf_add_missing_script33_R5_rows(ss,p,sets$primary,cfg,'retained_GA')
if(!all(c('rs','group','phase')%in%names(phase)))stop('Malformed R5 phase table')
if(!all(c('rs','group')%in%names(rs_gwas)))stop('Malformed R5 SNP GWAS table')
if(!all(c('group','chr','ps')%in%names(qtl)))stop('Malformed 11 cluster GWAS table')
if(!all(c('ID','SNP_ID','CHROM','POS','G0_depth','G0_ALT_reads','MAF_G0',
          'slope_freq_per_generation','n_pools','n_evolved_reps')%in%names(ss)))
  stop('Script 33 missing columns. Rerun corrected phase12 script33')
uf_stop_duplicate(phase$rs,'R5 phase IDs');uf_stop_duplicate(rs_gwas$rs,'GWAS SNP IDs')
uf_stop_duplicate(ss$ID,'script33 IDs');uf_stop_duplicate(qtl$group,'cluster IDs')
if(anyNA(phase$phase)||!all(phase$phase%in%c(-1,1)))stop('R5 phase must be +/-1')
if(!setequal(as.character(phase$group),as.character(qtl$group)))
  stop('R5 phase and GWAS cluster group universes do not agree')
mp<-match(phase$rs,as.character(rs_gwas$rs))
if(anyNA(mp)||any(as.character(phase$group)!=as.character(rs_gwas$group[mp])))
  stop('R5 phase and original GWAS SNP group IDs disagree')

## Rank each SNP using the strongest GWAS signal ACROSS EIGHT TRAITS.
## This ranking involves RIL phenotypes only, not Pool-seq frequency changes.
pcols<-grep('^pvalue_T[1-8]$',names(rs_gwas),value=TRUE)
if(length(pcols)!=8L)stop('Expected GWAS pvalue_T1..pvalue_T8 in script-02 table')
P<-as.matrix(rs_gwas[mp,pcols,drop=FALSE]);storage.mode(P)<-'double'
P[!is.finite(P)|P<0|P>1]<-NA_real_
rankp<-apply(P,1,function(z){z<-z[is.finite(z)];if(length(z))min(z)else Inf})

ix<-match(as.character(phase$rs),as.character(ss$ID))
phase$poolseq_present<-!is.na(ix)
phase$gwas_min_p<-rankp
phase$G0_depth<-ss$G0_depth[ix]
phase$G0_ALT_reads<-ss$G0_ALT_reads[ix]
phase$G0_MAF<-ss$MAF_G0[ix]
phase$slope_ALT_main<-ss$slope_freq_per_generation[ix]
phase$slope_phase_main<-as.numeric(phase$phase)*phase$slope_ALT_main
phase$passes_coverage<-phase$poolseq_present & is.finite(phase$G0_depth)&
  phase$G0_depth>=cfg$min_A00_depth &
  is.finite(ss$n_pools[ix])&ss$n_pools[ix]>=cfg$min_pools &
  is.finite(ss$n_evolved_reps[ix]) & ss$n_evolved_reps[ix]>=cfg$min_evolved_replicates &
  is.finite(phase$slope_ALT_main)
phase$CHROM<-sub('_.*$','',as.character(phase$rs))
phase$POS<-uf_pos_from_rs(phase$rs)
if(any(phase$passes_coverage & (!is.finite(phase$POS)|is.na(phase$CHROM))))
  stop('SNP positions could not be parsed')
phase$tag_for_neutral_test<-FALSE

## Cluster responses are based on all quality-eligible phase-aligned SNPs;
## never choose SNPs using drift significance or observed slope direction.
bygrp<-split(seq_len(nrow(phase)),as.character(phase$group))
records<-vector('list',length(bygrp));ii<-0L
for(g in names(bygrp)){
  ii<-ii+1L
  ixg<-bygrp[[g]];eligible<-ixg[phase$passes_coverage[ixg]]
  qq<-qtl[match(g,as.character(qtl$group)),,drop=FALSE]
  if(length(eligible)){
    ## SNP GWAS ranking; tie-break by ancestral depth; then fixed SNP ID.
    ord<-eligible[order(phase$gwas_min_p[eligible],-phase$G0_depth[eligible],phase$rs[eligible])]
    tagix<-ord[1];phase$tag_for_neutral_test[tagix]<-TRUE
    tagged<-as.character(phase$rs[tagix])
    med<-uf_finite_median(phase$slope_phase_main[eligible])
    tag_phase<-as.integer(phase$phase[tagix])
    median_MAF<-median(phase$G0_MAF[eligible])
    nzero<-sum(phase$G0_ALT_reads[eligible]==0,na.rm=TRUE)
  }else{
    tagged<-NA_character_;med<-NA_real_;tag_phase<-NA_integer_
    median_MAF<-NA_real_;nzero<-NA_integer_
  }
  records[[ii]]<-data.frame(group=g,chr=as.character(qq$chr[1]),
    ps=as.numeric(qq$ps[1]),n_R5_SNPs=length(ixg),
    n_present=sum(phase$poolseq_present[ixg]),n_eligible=length(eligible),
    n_A00_zero_ALT_eligible=nzero,median_A00_MAF=median_MAF,
    primary_slope_median_phase=med,tag_SNP_ID=tagged,
    tag_phase=tag_phase,tag_GWAS_min_p=if(length(eligible))phase$gwas_min_p[tagix]else NA_real_,
    tag_A00_depth=if(length(eligible))phase$G0_depth[tagix]else NA_real_,
    tag_A00_ALT_reads=if(length(eligible))phase$G0_ALT_reads[tagix]else NA_real_,
    reason=if(length(eligible))'eligible' else 'no_phase_SNP_passing_Poolseq_coverage',
    stringsAsFactors=FALSE)
}
registry<-do.call(rbind,records)
registry<-registry[match(as.character(qtl$group),registry$group),,drop=FALSE]
if(anyNA(registry$group)||anyDuplicated(registry$group)||nrow(registry)!=nrow(qtl))
  stop('Catalogue does not contain precisely the original R5 groups')
if(anyDuplicated(registry$tag_SNP_ID[!is.na(registry$tag_SNP_ID)]))stop('Tag assigned twice')
if(sum(phase$tag_for_neutral_test)!=sum(registry$reason=='eligible'))stop('Tag count mismatch')

## Freeze the primary-eligible SNP set for the GA3 supplementary comparison.
if(file.exists(cfg$existing_snp_supp)){
  supp<-uf_read(cfg$existing_snp_supp)
  supp<-uf_add_missing_script33_R5_rows(supp,p,sets$all,cfg,'all_GA')
  uf_stop_duplicate(supp$ID,'script33 supplementary SNP IDs')
  kk<-match(phase$rs,as.character(supp$ID))
  phase$slope_phase_supp<-as.numeric(phase$phase)*supp$slope_freq_per_generation[kk]
  sup<-vapply(bygrp,function(ind){
    v<-ind[phase$passes_coverage[ind]]
    uf_finite_median(phase$slope_phase_supp[v])},numeric(1))
  registry$supp_slope_median_phase<-unname(sup[registry$group])
}else registry$supp_slope_median_phase<-NA_real_

## Hard-stop if historical R5 focal loci have failed to join the unified catalogue.
for(fid in cfg$must_include_R5_focal_ids){
 fgrp<-as.character(phase$group[match(fid,phase$rs)])
 rr<-registry[as.character(registry$group)==fgrp,,drop=FALSE]
 sr<-phase[phase$rs==fid,,drop=FALSE]
 if(nrow(rr)!=1L||nrow(sr)!=1L||!isTRUE(sr$passes_coverage) ||
    !is.finite(sr$slope_phase_main)||rr$n_eligible<1L)
   stop('Previously mapped R5 focal SNP missing from primary genome-wide analysis: ',fid,
       ' (group ',fgrp,'). Inspect 43x raw integration QC.')
}
uf_write(registry,file.path(out,'43_Unified_R5_cluster_registry.csv'))
uf_write(phase,file.path(out,'43_All_SNP_phase_slope_and_inclusion_QC.csv'))
uf_write(registry[!is.na(registry$tag_SNP_ID),c('group','chr','ps','tag_SNP_ID',
   'tag_phase','tag_GWAS_min_p','tag_A00_depth','tag_A00_ALT_reads')],
   file.path(out,'43_FROZEN_selection_tags.csv'))
summary<-data.frame(n_catalogue_clusters=nrow(registry),
 n_with_poolseq_data=sum(registry$n_present>0),
 n_with_eligible_response=sum(registry$reason=='eligible'),
 n_previous_common_frequency=675L,
 n_zero_ALT_tags=sum(registry$tag_A00_ALT_reads==0,na.rm=TRUE),
 min_A00_depth=cfg$min_A00_depth,
 n_assigned_SNPs=nrow(phase),
 n_eligible_SNPs=sum(phase$passes_coverage),
 n_primary_replicates=3L,
 stringsAsFactors=FALSE)
uf_write(summary,file.path(out,'43_Catalogue_size_and_QC.csv'))
writeLines(capture.output(sessionInfo()),file.path(out,'43_sessionInfo.txt'))
cat('43 DONE: ',summary$n_with_eligible_response,' / ',nrow(registry),
    ' clusters with usable Pool-seq evolution. No ancestral-MAF exclusions.\n',sep='')
