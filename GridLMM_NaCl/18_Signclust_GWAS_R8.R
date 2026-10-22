rm(list=ls());gc()
library(qqman)
library(data.table)
library(dplyr)
library(igraph)
library(parallel)
threshold_clust=.8

selection_clusters = read.table(paste0("txt_GridLMM_NaCl/all_sign_SELECTION_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
dim(selection_clusters) # 3k selected SNPs

diff_signs = tapply(selection_clusters$s,selection_clusters$group,function(x){length(unique(sign(x)))})
groups_diff_signs = names(diff_signs)[as.numeric(diff_signs)==2]

load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
GWAS_results = read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
table(GWAS_results$sign_overall) # 33k sign SNPs

## Overlap between the two ?

dim(subset(selection_clusters,ID%in%GWAS_results$rs)) # 895 out of 3k / More than just by chance ?

### We can export a list of GWAS clusters that are significant - here meaning at least one SNP selected
dim(subset(GWAS_results,rs%in%selection_clusters$ID)) # 895 SNPs
length(unique(subset(GWAS_results,rs%in%selection_clusters$ID)$group)) # 94 GWAS clusters
write.table(unique(subset(GWAS_results,rs%in%selection_clusters$ID)$group),file='Selected_GWAS_clusters_R8.txt')

