## 43x_Integrate_rawSYNC_R5_SNPs.R
## Recovery of exact ALT counts from raw sync for SNPs already in R5 GWAS.
## NOT new LD groups, NOT a separate hypothesis family.
## Required before the final unified frequency-inclusive main-text run.
## This file is sourced by 43, 43b and 44; it does not run by itself.

uf_pool_with_R5_raw_recovery <- function(cfg) {
  if (!exists('read_drift_inputs',mode='function'))
    stop('Source 09f_drift_functions.R before invoking R5 recovery')
  check <- c(cfg$pool_file,cfg$orientation_file,cfg$phase_file,cfg$raw_focal_file)
  uf_check(check)
  p <- read_drift_inputs(cfg$pool_file,cfg$orientation_file,
                         exclude_G0='A6G140')
  raw <- uf_read(cfg$raw_focal_file)
  phase <- uf_read(cfg$phase_file)
  req <- c('SNP_ID','sample_id','CHROM','POS','REF','ALT','generation',
           'replicate','ALT_reads','depth','REF_agrees')
  if (!all(req %in% names(raw))) stop('Raw sync counts missing: ',
       paste(setdiff(req,names(raw)),collapse=', '))
  if (!all(c('rs','group','phase') %in% names(phase)))
    stop('R5 phase file is missing rs/group/phase columns')
  raw$SNP_ID <- as.character(raw$SNP_ID)
  raw$sample_id <- as.character(raw$sample_id)
  phase$rs <- as.character(phase$rs)
  if (anyDuplicated(phase$rs)) stop('Duplicate SNP rs in R5 phase file')
  if (any(raw$REF_agrees %in% FALSE,na.rm=TRUE))
    stop('Raw sync REF mismatch among recovered focal SNPs')
  if (anyNA(raw$REF_agrees) || anyNA(raw$ALT_reads) || anyNA(raw$depth))
    stop('Missing REF validation, read counts, or depth in raw sync table')
  if (any(raw$ALT_reads<0 | raw$depth<0 | raw$ALT_reads>raw$depth))
    stop('Invalid ALT reads or coverage in raw sync table')
  focal <- cfg$must_include_R5_focal_ids
  if (!length(focal)) stop('The prespecified R5 focal-site audit list is empty')
  missing_phase <- setdiff(focal,phase$rs)
  if (length(missing_phase)) stop('Expected focal sites absent from R5 GWAS phase file: ',
       paste(missing_phase,collapse=', '))
  missing_raw <- setdiff(focal,raw$SNP_ID)
  if (length(missing_raw)) stop('Expected focal sites absent from recovered raw sync: ',
       paste(missing_raw,collapse=', '))
  raw <- raw[raw$SNP_ID %in% phase$rs, , drop=FALSE]
  ## Verify the mapping of *every* supplied raw SNP before using any counts.
  key <- paste(raw$SNP_ID,raw$sample_id,sep='||')
  if (anyDuplicated(key)) stop('Duplicate raw-sync SNP/sample combinations')
  raw$CHROM <- as.character(raw$CHROM)
  raw$POS <- as.numeric(raw$POS)
  if (any(paste(raw$CHROM,raw$POS,sep='_') != raw$SNP_ID))
    stop('Raw-sync SNP ID disagrees with its genomic coordinates')
  if (!all(p$pop$sample_id %in% unique(raw$sample_id)))
    stop('Missing expected Pool-seq samples in raw-SYNC recovery')
  ref <- p$pop
  ix_pop <- match(raw$sample_id,as.character(ref$sample_id))
  ## A6G140 is the only G0 deliberately excluded from the current p object.
  kept <- !is.na(ix_pop)
  if (any(kept & (as.numeric(raw$generation) != as.numeric(ref$generations[ix_pop]) |
                     as.character(raw$replicate) != as.character(ref$replicate[ix_pop]))))
    stop('Raw-sync sample manifest disagrees with corrected population mapping')
  raw <- raw[kept,,drop=FALSE]
  if (anyDuplicated(p$sel$ID))stop('Duplicate original GWAS-oriented Pool-seq IDs')
  ids_all <- unique(raw$SNP_ID)
  ids_missing <- setdiff(ids_all,as.character(p$sel$ID))
  if (anyDuplicated(as.character(p$sel$rs)))stop('Duplicate original Pool-seq rs IDs')
  new_counts <- matrix(NA_real_,nrow=length(ids_missing),ncol=nrow(ref),
    dimnames=list(ids_missing,as.character(ref$sample_id)))
  new_depth <- new_counts
  new_sel <- p$sel[rep(1L,length(ids_missing)),,drop=FALSE]
  if(length(ids_missing)) {
    new_sel[] <- NA
    for (i in seq_along(ids_missing)) {
      id <- ids_missing[i]
      z <- raw[raw$SNP_ID==id,,drop=FALSE]
      ii <- match(ref$sample_id,z$sample_id)
      if (anyNA(ii) || nrow(z)!=nrow(ref))
        stop('Incomplete 12-sample raw-sync trajectory for ',id)
      if (length(unique(z$REF))!=1 || length(unique(z$ALT))!=1 ||
          any(z$REF==z$ALT))stop('Ambiguous REF/ALT designation at ',id)
      new_counts[i,] <- as.numeric(z$ALT_reads[ii])
      new_depth[i,] <- as.numeric(z$depth[ii])
      new_sel$rs[i] <- id
      new_sel$ID[i] <- id
      new_sel$CHROM[i] <- as.character(z$CHROM[1])
      new_sel$POS[i] <- as.numeric(z$POS[1])
      new_sel$orientation[i] <- 'minor_is_alt'
    }
    ## For the appended rows, the stored 'minor' matrix is now explicitly ALT.
    ## No additional flip is applied at read_drift_inputs time.
    p$sel <- rbind(p$sel,new_sel)
    p$k_alt <- rbind(p$k_alt,new_counts)
    p$depth <- rbind(p$depth,new_depth)
  }
  if(nrow(p$sel)!=nrow(p$k_alt)||!identical(dim(p$k_alt),dim(p$depth)))
    stop('Raw-sync integration produced misaligned Pool-seq matrices')
  if(!identical(colnames(p$k_alt),colnames(p$depth)) ||
     !identical(colnames(p$k_alt),as.character(ref$sample_id)))
    stop('Raw-sync integration changed Pool-seq sample column order')
  if(anyDuplicated(p$sel$ID))stop('Duplicate after raw-sync R5 integration')
  if(!all(focal %in% p$sel$ID))stop('Five prespecified R5 loci are not in the integrated pool')
  qc <- data.frame(SNP_ID=ids_all,
    group=as.character(phase$group[match(ids_all,phase$rs)]),
    from_raw_SYNC=TRUE,
    already_present_in_original_oriented_Poolseq=ids_all %in% setdiff(as.character(p$sel$ID),ids_missing),
    integrated_into_unified_R5_pool=ids_all %in% p$sel$ID,
    newly_recovered=ids_all %in% ids_missing,
    stringsAsFactors=FALSE)
  uf_write(qc,file.path(cfg$prefix,'43x_R5_rawSYNC_integration_QC.csv'))
  p$raw_recovered_ids <- ids_missing
  p$raw_R5_ids <- ids_all
  p
}

uf_add_missing_script33_R5_rows <- function(ss,p,keep,cfg,scenario) {
  if(!exists('conditional_neutral_info',mode='function'))
    stop('Source 33_Frequency_drift_helpers.R before R5 slope recovery')
  if(!all(c('ID','SNP_ID','G0_depth','G0_ALT_reads','MAF_G0',
       'slope_freq_per_generation','n_pools','n_evolved_reps') %in% names(ss)))
    stop('Script 33 table lacks required columns')
  if(anyDuplicated(ss$ID))stop('Script 33 contains duplicate site IDs')
  ## Only ever inject rows explicitly checked against the raw sync and R5.
  missing <- setdiff(p$raw_R5_ids,as.character(ss$ID))
  if(!length(missing))return(ss)
  pop <- p$pop[keep,,drop=FALSE]
  rows <- vector('list',length(missing))
  for(i in seq_along(missing)){
    id <- missing[i]
    j <- match(id,as.character(p$sel$ID))
    if(is.na(j))stop('Raw-R5 SNP has no Pool-seq counts: ',id)
    info <- conditional_neutral_info(as.numeric(p$k_alt[j,keep]),
       as.numeric(p$depth[j,keep]),pop,Ne=cfg$Ne)
    if(!isTRUE(info$ok))
      warning('Recovered site not estimable under drift null: ',id)
    z <- ss[1L,,drop=FALSE];z[] <- NA
    z$ID <- id;z$SNP_ID <- id
    chr <- sub('_.*$','',id)
    z$CHROM <- chr;z$POS <- as.numeric(sub('^.*_','',id))
    if('scenario'%in%names(z))z$scenario<-scenario
    if('Ne'%in%names(z))z$Ne<-cfg$Ne
    if('p0_ALT_posterior'%in%names(z))z$p0_ALT_posterior<-info$p0
    if('p0_ALT_CI_lower'%in%names(z))z$p0_ALT_CI_lower<-info$p0_lo
    if('p0_ALT_CI_upper'%in%names(z))z$p0_ALT_CI_upper<-info$p0_hi
    z$G0_ALT_reads<-info$baseline_ALT_reads
    z$G0_depth<-info$baseline_depth
    z$MAF_G0<-if(is.finite(info$p0))min(info$p0,1-info$p0)else NA_real_
    z$slope_freq_per_generation<-info$slope
    z$n_pools<-info$n
    z$n_evolved_reps<-info$n_evolved_reps
    if('neutral_expected_slope'%in%names(z))z$neutral_expected_slope<-info$mean_slope
    if('conditional_SE'%in%names(z))z$conditional_SE<-info$conditional_SE
    if('z_normal_approx'%in%names(z))z$z_normal_approx<-info$conditional_z
    if('p_normal_approx'%in%names(z))z$p_normal_approx<-info$conditional_p_normal
    if('G0_ALT_not_detected'%in%names(z))
      z$G0_ALT_not_detected<-is.finite(info$baseline_ALT_reads)&&info$baseline_ALT_reads==0
    if('MAF_class'%in%names(z)) z$MAF_class<-if(is.finite(z$MAF_G0))
      as.character(frequency_class(info$p0))else NA_character_
    rows[[i]]<-z
  }
  rbind(ss,do.call(rbind,rows))
}
