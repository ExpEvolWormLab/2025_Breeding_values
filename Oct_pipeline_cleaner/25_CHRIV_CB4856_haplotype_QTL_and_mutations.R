## 25_CHRIV_CB4856_haplotype_QTL_and_mutations.R
## Three-stage AUDIT of chrIV: (i) current R5 GWAS QTLs near the two
## mutations; (ii) overlap with drift-aware selection evidence; (iii) effects
## of CB4856-like background (H1,H2) *before* adding mutation terms (M1,M2).
## Reuses the already kinship-adjusted fitted models of script 21; does not
## claim to have estimated interactions that lack genotype support.

## ---------------------------- SETTINGS --------------------------------
region_file <- 'historical_QTL_R5/23_R5_all_cluster_boundaries_and_associations.csv'
r5_snp_file <- 'txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R5.txt'
cluster_genotype_file <- 'data_BLUE_GridLMM_NaCl/genotypes_matrix_clusters_R5.txt'
mutation_QC_file <- 'mutation_effects_chrIV/CB4856_all8/21_Genotype_and_ancestry_QC.csv'
mutation_effect_file <- 'mutation_effects_chrIV/CB4856_all8/21_All_eight_trait_effects.csv'
mutation_model_QC_file <- 'mutation_effects_chrIV/CB4856_all8/21_Model_identifiability_and_power_QC.csv'
survivor_clusters <- 'selection_results/survivors_only_Ne1000/20_Survivor_cluster_slope_summary.csv' # optional
out_dir <- 'historical_QTL_R5/chrIV_CB4856'
mutation_positions <- c(8892663L,10499135L)
window_bp <- 750000L
traits <- c('SF','SB','FS','FB','BS','BF','area.F','beta')
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

req_files <- c(region_file,mutation_QC_file,mutation_effect_file)
if(any(!file.exists(req_files)))stop('Run script 23 and 21 first; missing: ',
                                   paste(req_files[!file.exists(req_files)],collapse=', '))
r5 <- read.csv(region_file,stringsAsFactors=FALSE,check.names=FALSE)
m <- read.csv(mutation_QC_file,stringsAsFactors=FALSE,check.names=FALSE)
eff <- read.csv(mutation_effect_file,stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c('group','CHROM','start','end','q_SB','q_beta')%in%names(r5)))
  stop('Missing R5 region fields (run script 23)')
if(!all(c('pop_label','mut1_dosage','mut2_dosage','CB4856_8892663',
         'CB4856_10499135')%in%names(m)))stop('Script-21 genotype/ancestry fields missing')
if(!all(c('trait','analysis','term','beta_as_homozygote_contrast','SE',
         'p_wald','q_BH_eight_traits')%in%names(eff)))stop('Script-21 effect fields missing')

## 1. R5 regional clusters with ACTUAL MEMBER SNPs in mutation windows.
## A cluster can span a long interval without any members near a mutation;
## bounding-box intersection alone is not sufficient evidence.
ss <- read.table(r5_snp_file,header=TRUE,sep='\t',stringsAsFactors=FALSE,
                 check.names=FALSE,quote='',comment.char='')
ss$group <- as.character(ss$group)
ss$POS <- as.numeric(ss$POS)
ss$CHROM <- as.character(ss$CHROM)
ss <- ss[ss$CHROM=='IV',,drop=FALSE]
if(!nrow(ss))stop('No R5 chromosome IV SNPs found')
near_member <- do.call(rbind,lapply(split(ss,ss$group),function(z){
  data.frame(group=as.character(z$group[1]),
    distance_to_mut1_bp=min(abs(z$POS-mutation_positions[1])),
    distance_to_mut2_bp=min(abs(z$POS-mutation_positions[2])),
    n_member_SNPs_near_mut1=sum(abs(z$POS-mutation_positions[1])<=window_bp),
    n_member_SNPs_near_mut2=sum(abs(z$POS-mutation_positions[2])<=window_bp))
}))
near_member <- near_member[near_member$n_member_SNPs_near_mut1>0 |
                             near_member$n_member_SNPs_near_mut2>0,,drop=FALSE]
near_iv <- merge(r5,near_member,by='group',sort=FALSE)
near_iv <- near_iv[near_iv$CHROM=='IV',,drop=FALSE]

near_iv$mutation_1_is_R5_member <- FALSE
near_iv$mutation_2_is_R5_member <- FALSE
if(file.exists('historical_QTL_R5/23_Historical_SNP_status.csv')){
  h <- read.csv('historical_QTL_R5/23_Historical_SNP_status.csv',stringsAsFactors=FALSE)
  for(j in 1:2){
    k <- which(h$SNP_ID==paste0('IV_',mutation_positions[j]))
    if(length(k)==1L&&!is.na(h$R5_group_exact[k]))
      near_iv[[paste0('mutation_',j,'_is_R5_member')]] <-
        near_iv$group==as.character(h$R5_group_exact[k])
  }
}
if(file.exists(survivor_clusters)){
  s <- read.csv(survivor_clusters,stringsAsFactors=FALSE)
  k <- match(as.character(near_iv$group),as.character(s$group))
  for(nm in intersect(c('n_survivor_common','n_selected_survivor',
                        'any_selected_survivor','slope_survivor_median_common'),names(s)))
    near_iv[[nm]]<-s[[nm]][k]
}
near_iv$R5_SB_q05 <- is.finite(near_iv$q_SB)&near_iv$q_SB<.05
near_iv$R5_beta_q05 <- is.finite(near_iv$q_beta)&near_iv$q_beta<.05
near_iv$R5_selected_full_q01 <- if('any_selected_common_q01'%in%names(near_iv))
  !is.na(near_iv$any_selected_common_q01)&near_iv$any_selected_common_q01 else NA
near_iv$R5_selected_survivor_q01 <- if('any_selected_survivor'%in%names(near_iv))
  !is.na(near_iv$any_selected_survivor)&near_iv$any_selected_survivor else NA

## 2. Does the cluster-tagging haplotype segregate without either mutation?
## Cluster dosage is already PHASE-ALIGNED within the R5 LD cluster.
## Correlations to H are proxies for a CB4856-associated region, NOT proof
## that the cluster SNPs belong to a single phased CB4856 haplotype.
helper_cor <- function(x,y){v <- is.finite(x)&is.finite(y)
  if(sum(v)<10L||sd(x[v])==0||sd(y[v])==0)return(NA_real_)
  cor(x[v],y[v])}

for(nm in c('r_cluster_M1','r_cluster_M2','r_cluster_H1','r_cluster_H2',
            'r_cluster_H1_no_mutations','r_cluster_H2_no_mutations',
            'n_RILs_overlap','n_nonmutant_RILs'))near_iv[[nm]] <- NA_real_
if(file.exists(cluster_genotype_file)&&nrow(near_iv)) {
  GG <- as.matrix(read.table(cluster_genotype_file,header=TRUE,row.names=1,
                             check.names=FALSE,stringsAsFactors=FALSE))
  storage.mode(GG) <- 'double'
  id <- match(colnames(GG),m$pop_label)
  v <- !is.na(id)
  if(sum(v)<30L)warning('Few RIL IDs matched between cluster genotype and script 21')
  m <- m[id[v],,drop=FALSE]
  GG <- GG[,v,drop=FALSE]
  for(i in seq_len(nrow(near_iv))){
    id <- as.character(near_iv$R5_cluster_ID[i])
    if(!id%in%rownames(GG))next
    g <- as.numeric(GG[id,])
    M1 <- m$mut1_dosage/2;M2 <- m$mut2_dosage/2
    H1 <- m$CB4856_8892663;H2<-m$CB4856_10499135
    nonmut <- is.finite(M1)&is.finite(M2)&M1==0&M2==0
    near_iv$r_cluster_M1[i] <- helper_cor(g,M1)
    near_iv$r_cluster_M2[i] <- helper_cor(g,M2)
    near_iv$r_cluster_H1[i] <- helper_cor(g,H1)
    near_iv$r_cluster_H2[i] <- helper_cor(g,H2)
    near_iv$r_cluster_H1_no_mutations[i] <- helper_cor(g[nonmut],H1[nonmut])
    near_iv$r_cluster_H2_no_mutations[i] <- helper_cor(g[nonmut],H2[nonmut])
    near_iv$n_RILs_overlap[i] <- sum(is.finite(g))
    near_iv$n_nonmutant_RILs[i] <- sum(nonmut&is.finite(g))
  }
} else warning('Cluster genotype matrix unavailable: regional GWAS/selection annotations still exported')
near_iv$candidate_CB4856_background_QTL_without_mutations <-
  !near_iv$mutation_1_is_R5_member & !near_iv$mutation_2_is_R5_member &
  (is.finite(near_iv$r_cluster_H1_no_mutations)&
     abs(near_iv$r_cluster_H1_no_mutations)>.5 |
   is.finite(near_iv$r_cluster_H2_no_mutations)&
     abs(near_iv$r_cluster_H2_no_mutations)>.5)
## The 0.5 correlation cutoff is only a descriptive screen, NOT a p-value.
write.csv(near_iv,file.path(out_dir,'25_R5_chrIV_QTL_selection_CB4856_vs_mutations.csv'),
          row.names=FALSE)

## 3. Organize the EXISTING staged, kinship-adjusted script-21 fits:
## Stage A: H1+H2 only (plus source-population and LOCO kinship).
## Stage B: H1+H2+M1+M2. Mutation terms are conditional on H.
## This makes the two estimands transparent; it is NOT a new fit.
old <- eff[eff$analysis=='with_CB4856'&
             eff$term%in%c('H1','H2','M1','M2')&eff$trait%in%traits,,drop=FALSE]
old$stage <- ifelse(old$term%in%c('H1','H2'),
             'A_CB4856_background_only_without_mutation_covariates',
             'B_mutations_conditional_on_CB4856_and_other_mutation')
old$interpretation <- ifelse(old$term%in%c('H1','H2'),
 'Haplotype association; may include the effect of linked mutation carriers',
 'Mutation association within observed haplotype backgrounds; not causal proof')
old <- old[order(match(old$trait,traits),match(old$term,c('H1','H2','M1','M2'))),,drop=FALSE]
write.csv(old,file.path(out_dir,'25_StageA_haplotype_StageB_conditional_mutations_all8.csv'),row.names=FALSE)

if(file.exists(mutation_model_QC_file)){
  qc <- read.csv(mutation_model_QC_file,stringsAsFactors=FALSE)
  qc <- qc[qc$analysis=='with_CB4856',,drop=FALSE]
  kk <- intersect(c('trait','n_RIL','p_H_given_origin',
      'p_mutations_given_H','p_M1M2_given_additive',
      'H_model_estimable','adjusted_two_mutations_estimable'),names(qc))
  write.csv(qc[,kk,drop=FALSE],file.path(out_dir,'25_Nested_haplotype_then_mutation_tests.csv'),row.names=FALSE)
}

## 4. Cell table: observed mutation contrast within CB4856-like backgrounds.
M1 <- m$mut1_dosage/2; M2 <- m$mut2_dosage/2
H1 <- m$CB4856_8892663;H2<-m$CB4856_10499135
cell <- as.data.frame(table(M1=factor(round(M1),levels=0:1),
 M2=factor(round(M2),levels=0:1),
 H1=factor(round(H1),levels=0:1),H2=factor(round(H2),levels=0:1)))
names(cell)[5] <- 'n_RILs'
write.csv(cell,file.path(out_dir,'25_CB4856_and_two_mutation_genotype_cells.csv'),row.names=FALSE)

if(nrow(near_iv)) {
  print(near_iv[,intersect(c('group','start','end','q_SB','q_beta',
       'R5_selected_full_q01','R5_selected_survivor_q01',
       'r_cluster_H2_no_mutations',
       'candidate_CB4856_background_QTL_without_mutations'),names(near_iv))],row.names=FALSE)
}
cat('\nOutputs: ',out_dir,'\n',sep='')
cat('Interpretation: cluster-haplotype correlation in nonmutant RILs is descriptive.\n')
cat('A single H+M fit does not identify mutation x H when carriers occur only on H=1.\n')
