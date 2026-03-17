rm(list=ls());gc()
library(data.table)

threshold_clust = 0.8
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
#load(file="RData/GWAS_results_concatenated.RData")
gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)
#Kmat <- as.matrix(read.table("data_BLUE/kinship.txt"))

genotypes <- read.table(file=paste0("data_BLUE_clusters/genotypes_clusters_R",round(10*threshold_clust),".txt"))
#genotypes2=read.table(file="data_BLUE/genotypes_for_gemma.txt")

table(genotypes$V1==all_GWAS_clusters$rs) # Not in the good order
genotypes=genotypes[order(genotypes$V1),]
all_GWAS_clusters=all_GWAS_clusters[order(all_GWAS_clusters$rs),]
table(genotypes$V1==all_GWAS_clusters$rs) # Fixed

pred_T1 = unlist(apply(genotypes[,4:969],2,function(x){sum(x[all_GWAS_clusters$is_sign_T1]*all_GWAS_clusters$beta_T1[all_GWAS_clusters$is_sign_T1])}))
plot(pred_T1,gemma_pheno$V1)

pred_T7 = unlist(apply(genotypes[,4:969],2,function(x){sum(x[all_GWAS_clusters$is_sign_T7]*all_GWAS_clusters$beta_T7[all_GWAS_clusters$is_sign_T7])}))
pred_T14 = unlist(apply(genotypes[,4:969],2,function(x){sum(x[all_GWAS_clusters$is_sign_T14]*all_GWAS_clusters$beta_T14[all_GWAS_clusters$is_sign_T14])}))


plot(pred_T14,gemma_pheno$V14)
abline(a=0,b=1/5)
plot(pred_T7,gemma_pheno$V7)

cor(pred_T7,pred_T14)
cor(na.omit(gemma_pheno))[7,14]

cor(pred_T1,pred_T7)
cor(na.omit(gemma_pheno))[1,7]







