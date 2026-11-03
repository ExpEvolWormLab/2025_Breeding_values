## 55_ChrV_RIL_LD_structure.R
## Reconstruct RIL-genotype linkage across ALL R5 groups with representative
## positions at V:6--8 Mb (including groups without Pool-seq response).
## Requires original project's RIL genotype_matrix and R5 phase/group membership.
## This is an LD / distinguishability analysis, NOT a causal fine-map.
## Run from the project root, after script 54.

options(stringsAsFactors=FALSE)
set.seed(20261102L)
root_out <- 'selection_results/chrV_deep_dive'
dir.create(root_out,recursive=TRUE,showWarnings=FALSE)
max_snps_per_group <- 8L
min_paired_RIL <- 40L
high_LD_r2 <- 0.80
first_path <- function(x) { x<-x[file.exists(x)]; if(length(x))x[1L] else NA_character_ }
read_req <- function(paths,desc) {
 f<-first_path(paths);if(is.na(f))stop('Missing ',desc,':\n',paste(paths,collapse='\n'),call.=FALSE)
 message(desc,': ',f);read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
}
w <- function(x,n)write.csv(x,file.path(root_out,n),row.names=FALSE,na='')
A <- read_req(c(file.path(root_out,'54_V3_12_cluster_prediction_audit.csv')),'script 54 V_3 audit')
Q <- read_req(c('selection_results/drift_oneG0_nullfit/R5/11_Clusters_selection_and_GWAS_09f.csv',
                '11_Clusters_selection_and_GWAS_09f.csv'),'complete R5 GWAS cluster table')
R <- Q[Q$chr=='V' & Q$ps>=6000000 & Q$ps<8000000,c('group','ps'),drop=FALSE]
A <- unique(rbind(A[,c('group','ps'),drop=FALSE],R))
A <- A[order(A$ps),,drop=FALSE]
if(anyDuplicated(A$group))stop('Duplicated R5 groups in local 6--8 Mb membership')
message('Including ',nrow(A),' R5 groups, with or without measurable Pool-seq response')
phase <- read_req(c('selection_results/unified_frequency_inclusive_R5/43_All_SNP_phase_slope_and_inclusion_QC.csv',
                    'data_BLUE_GridLMM_NaCl/SNP_phase_clusters_R5.csv'),'R5 SNP-to-cluster phase table')
if(!all(c('rs','group','phase')%in%names(phase)))stop('Phase mapping must have rs, group, phase')
if(anyDuplicated(phase$rs))stop('Duplicated SNP IDs in phase table')
phase$group <- as.character(phase$group)
phase$phase <- as.numeric(phase$phase)
if(any(!phase$phase%in%c(-1,1)))stop('Invalid phase signs')
phase$position <- suppressWarnings(as.numeric(sub('^.*_','',as.character(phase$rs))))
phase$chromosome <- sub('_.*$','',as.character(phase$rs))
phase <- phase[phase$group%in%as.character(A$group)&phase$chromosome=='V',,drop=FALSE]
if(!nrow(phase))stop('No phase entries for V:6--8 Mb R5 groups')
full_phase_counts <- table(phase$group)

f <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
if(!file.exists(f))stop('Missing original RIL input: ',f)
e <- new.env(parent=globalenv()); load(f,envir=e)
if(!exists('genotype_matrix',envir=e,inherits=FALSE))stop('No genotype_matrix in ',f)
G <- e$genotype_matrix
if(is.data.frame(G))G <- as.matrix(G)
if(!is.matrix(G)||is.null(rownames(G))||is.null(colnames(G)))
 stop('RIL genotype_matrix must be RIL-by-SNP with row/column names')
if(anyDuplicated(colnames(G)))stop('Duplicate RIL genotype SNP names')
phase <- phase[phase$rs%in%colnames(G),,drop=FALSE]
if(!nrow(phase))stop('No R5 phase markers found in genotype_matrix')
## Convert only selected columns to double below; converting all GWAS SNPs
## would duplicate the large genotype matrix in memory.

## Candidate tags are deliberately NOT claimed to be causal variants.
Lpath <- first_path(c('pleiotropy_results/R5/linkage_diagnostics/17_Within_cluster_lead_comparisons.csv',
                      '17_Within_cluster_lead_comparisons.csv'))
L <- if(is.na(Lpath))data.frame() else read.csv(Lpath,stringsAsFactors=FALSE)
Epath <- first_path(c('selection_results/unified_frequency_inclusive_R5/45_All_R5_clusters_GWAS_evolution_predictions.csv',
                      '45_All_R5_clusters_GWAS_evolution_predictions.csv'))
E <- if(is.na(Epath))data.frame() else read.csv(Epath,stringsAsFactors=FALSE)
if(nrow(E))E<-E[E$scenario=='primary_survivors' & E$reason=='eligible',,drop=FALSE]

selected <- list();group_qc <- list()
for(i in seq_len(nrow(A))){
  group <- as.character(A$group[i]); members <- phase[phase$group==group,,drop=FALSE]
  if(!nrow(members)) {
    group_qc[[length(group_qc)+1L]]<-data.frame(group=group,
      n_phase_member=as.integer(if(group%in%names(full_phase_counts))full_phase_counts[[group]]else 0L),
      n_genotype_member=0L,n_tag_selected=0L,start=NA,end=NA)
    next
  }
  selected_here <- character()
  # Primary neutral tag is different from a GWAS cluster effect representative.
  if(nrow(E)){
    j<-match(group,as.character(E$group))
    if(!is.na(j))selected_here<-c(selected_here,as.character(E$selection_tag_SNP_ID[j]))
  }
  if(nrow(L)) {
    z<-L[as.character(L$group)==group,,drop=FALSE]
    selected_here<-c(selected_here,as.character(z$lead_A),as.character(z$lead_B))
  }
  selected_here<-c(selected_here,intersect(c('V_6322080','V_6300718',
    'V_6637695','V_6826421','V_6959758'),members$rs))
  selected_here<-unique(selected_here[selected_here%in%members$rs])
  members<-members[order(members$position,members$rs),,drop=FALSE]
  centroid<-as.numeric(A$ps[i])
  near<-members$rs[order(abs(members$position-centroid))][1]
  selected_here<-unique(c(selected_here,near))
  # Add broad coverage of large clusters instead of treating all SNPs as independent.
  choose <- unique(round(seq(1,nrow(members),length.out=max_snps_per_group)))
  for(snp in members$rs[choose]){
    if(length(selected_here)>=max_snps_per_group)break
    selected_here<-unique(c(selected_here,snp))
  }
  if(length(selected_here)>max_snps_per_group)
     warning('Group ',group,' retains ',length(selected_here),
      ' prioritized markers, above per-group target; all are retained')
  selected_here<-selected_here[selected_here%in%members$rs]
  pick<-members[match(selected_here,members$rs),,drop=FALSE]
  pick$source <- ifelse(pick$rs==near,'nearest_centroid','regional_or_spaced_marker')
  selected[[length(selected)+1L]] <- pick
  group_qc[[length(group_qc)+1L]]<-data.frame(group=group,
    n_phase_member=as.integer(full_phase_counts[[group]]),
    n_genotype_member=nrow(members),n_tag_selected=nrow(pick),
    start=min(members$position),end=max(members$position))
}
T <- if(length(selected))do.call(rbind,selected) else stop('No selected markers')
T <- T[!duplicated(T$rs),,drop=FALSE]
if(nrow(T)<2L)stop('Need at least two informative selected RIL markers to compute LD')
if(anyNA(T$phase))stop('Missing phase in selected tags')
Qc <- do.call(rbind,group_qc)
w(Qc,'55_Focal_cluster_physical_span_and_marker_QC.csv')

Z <- G[,as.character(T$rs),drop=FALSE]
if(!is.numeric(Z))storage.mode(Z)<-'double'
if(any(is.finite(Z)&(Z<0|Z>2)))stop('Unexpected RIL genotype dosages outside [0,2]')
for(k in seq_len(ncol(Z)))if(T$phase[k]==-1)Z[,k] <- 2-Z[,k]
T$n_genotyped_RIL <- colSums(is.finite(Z))
T$mean_aligned_dosage <- colMeans(Z,na.rm=TRUE)
T$sd_aligned_dosage <- apply(Z,2,sd,na.rm=TRUE)
T$is_polymorphic <- is.finite(T$sd_aligned_dosage)&T$sd_aligned_dosage>0
T$position <- as.numeric(T$position)
w(T[,intersect(c('rs','group','position','phase','source','n_genotyped_RIL',
                       'mean_aligned_dosage','sd_aligned_dosage','is_polymorphic'),names(T))],
  '55_Selected_RIL_markers.csv')

## Pairwise signed correlations reflect phase orientation. R2 ignores phase.
C <- suppressWarnings(cor(Z,use='pairwise.complete.obs'))
n <- ncol(Z)
rs <- as.character(T$rs)
rows<-vector('list',max(0L,n*(n-1L)/2L));a<-0L
for(i in seq_len(n-1L))for(j in (i+1L):n){
  a<-a+1L
  good<-is.finite(Z[,i])&is.finite(Z[,j]);nn<-sum(good)
  r<-if(nn>=min_paired_RIL&&is.finite(C[i,j]))C[i,j] else NA_real_
  rows[[a]]<-data.frame(rs_A=rs[i],rs_B=rs[j],group_A=T$group[i],group_B=T$group[j],
    pos_A=T$position[i],pos_B=T$position[j],distance_bp=abs(T$position[i]-T$position[j]),
    n_RIL=nn,aligned_r=r,r2=r*r,across_R5_clusters=T$group[i]!=T$group[j],
    high_LD=is.finite(r)&r*r>=high_LD_r2)
}
pairs<-do.call(rbind,rows)
w(pairs,'55_Pairwise_RIL_LD_across_chrV_markers.csv')

## Inter-cluster LD evidence. Max tag R2 is suggestive; no claim of independence
## can be based on absence of R2 at only eight sampled markers.
pairs2<-pairs[pairs$across_R5_clusters,,drop=FALSE]
if(nrow(pairs2)){
  pairs2$pair<-vapply(seq_len(nrow(pairs2)),function(i)
    paste(sort(c(as.character(pairs2$group_A[i]),as.character(pairs2$group_B[i]))),collapse='__'),character(1))
  keep<-split(seq_len(nrow(pairs2)),pairs2$pair)
  cross<-do.call(rbind,lapply(names(keep),function(key){
    ii<-keep[[key]];rr<-pairs2$r2[ii]
    data.frame(group_pair=key,n_tag_pairs=length(ii),n_estimable=sum(is.finite(rr)),
      max_RIL_tag_r2=if(any(is.finite(rr)))max(rr,na.rm=TRUE)else NA_real_,
      median_RIL_tag_r2=if(any(is.finite(rr)))median(rr,na.rm=TRUE)else NA_real_,
      n_tag_pairs_r2_ge_0.8=sum(rr>=high_LD_r2,na.rm=TRUE))
  }))
  cross<-cross[order(-cross$max_RIL_tag_r2),,drop=FALSE]
  w(cross,'55_Inter_cluster_LD_summary.csv')
}

## Descriptive effective dimension among selected, polymorphic, mean-imputed
## marker dosages; this is NOT the number of independent selected loci.
retain<-which(T$is_polymorphic & T$n_genotyped_RIL>=min_paired_RIL)
if(length(retain)>=2L){
  M<-Z[,retain,drop=FALSE]
  for(j in seq_len(ncol(M)))M[!is.finite(M[,j]),j]<-mean(M[,j],na.rm=TRUE)
  R<-cor(M)
  eg<-pmax(eigen(R,symmetric=TRUE,only.values=TRUE)$values,0)
  neff<-(sum(eg)^2)/sum(eg^2)
  writeLines(c(paste('n_RIL',nrow(Z)),paste('n_markers_total',ncol(Z)),
    paste('n_polymorphic_markers',length(retain)),
    paste('effective_dimension_of_selected_RIL_dosages',signif(neff,5)),
    'This dimensionality depends on marker sampling, NOT causal independent sites.'),
    file.path(root_out,'55_Selected_RIL_marker_dimension.txt'))
}

## Plot inter-marker LD in coordinate order; long-range LD can cross R5 clusters.
pdf(file.path(root_out,'55_Regional_RIL_LD_heatmap.pdf'),width=11,height=10)
ord<-order(T$position,T$group)
Aheat<-C[ord,ord,drop=FALSE]^2
Aheat[!is.finite(Aheat)]<-NA_real_
par(mar=c(10,10,4,2))
image(seq_len(nrow(Aheat)),seq_len(ncol(Aheat)),Aheat,axes=FALSE,
      col=colorRampPalette(c('white','#88BEE8','#164A85'))(100),zlim=c(0,1),
      xlab='Selected R5 marker',ylab='Selected R5 marker',
      main='chr V: RIL r² (phase invariant; selected informative SNPs)')
axis(1,at=seq_along(ord),labels=paste0(T$group[ord],':',round(T$position[ord]/1e6,2)),las=2,cex.axis=.45)
axis(2,at=seq_along(ord),labels=paste0(T$group[ord],':',round(T$position[ord]/1e6,2)),las=2,cex.axis=.45)
box();dev.off()
cat('55 FINISHED: RIL LD results in ',root_out,'\n',sep='')
cat('CAUTION: LD r2 and founder similarity are not proof of molecular pleiotropy,',
    ' epistasis, or independent selected mutations.\n')
