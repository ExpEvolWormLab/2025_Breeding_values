rm(list=ls());gc()
library(data.table)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")
vect_P_traits_bothEnv <- paste0(vect_P_traits,"_NaCl")

threshold_clust = 0.2
load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
genotypes <- read.table(file=paste0("data_BLUE_GridLMM_NaCl//genotypes_matrix_clusters_R",round(10*threshold_clust),".txt"))

table(row.names(genotypes)==all_GWAS_clusters$rs) # OK
genotypes=genotypes[order(row.names(genotypes)),]
all_GWAS_clusters=all_GWAS_clusters[order(all_GWAS_clusters$rs),]
table(row.names(genotypes)==all_GWAS_clusters$rs)  # Fixed

pred_traits=NULL
for(i in 1:8){
pred_traits = cbind(pred_traits,unlist(apply(genotypes,2,function(x){
  sum(x[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])]]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])][all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])})))
}

load('RData/Gmatrices_P.RData')
all_traits_corrected=all_traits_corrected[order(all_traits_corrected$pop_label),]
pred_traits=pred_traits[order(rownames(pred_traits)),]
table(rownames(pred_traits)==all_traits_corrected$pop_label)
gemma_pheno=all_traits_corrected[,vect_P_traits_bothEnv]

i=7
pdf("Predicted_phen__NaCl_R2.pdf",h=4,w=4)
plot(gemma_pheno[,i],pred_traits[,i],pch=16,xaxt="n",yaxt="n",xlab="Mean RILs phenotype",ylab="Predicted phenotype from sign. clusters",
     bty="n")
abline(lm(pred_traits[,i]~gemma_pheno[,i]))
dev.off()

df.lm=NULL
for(i in 1:8){
  temp.lm=lm(pred_traits[,i]~gemma_pheno[,i])
  df.lm=rbind(df.lm,c(coef(temp.lm),summary(temp.lm)$r.squared))
}
summary(df.lm[,2])
summary(df.lm[,3])

all_cor_explained = NULL
for(i in 1:8){
  isNA= is.na(pred_traits[,i]) | is.na(gemma_pheno[,i])
  all_cor_explained=rbind(all_cor_explained,(cbind(cor(pred_traits[!isNA,i], gemma_pheno[!isNA,i])^2,
        summary(lm(pred_traits[,i]~gemma_pheno[,i]))$adj.r.squared)))
}

pred_traits_noP=NULL
pred_traits_onlyP=NULL

all_GWAS_clusters$sign_overall = rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv)])
sum(all_GWAS_clusters$sign_overall>=1)/nrow(all_GWAS_clusters) # 99% significant in r=0.2

### Sign ~ Length Cluster ?
clust_length = data.frame(clusterID=names(table(all_GWAS_clusters$clusterID)),clust_length=as.numeric(table(all_GWAS_clusters$clusterID)))

all_GWAS_clusters_onlyP = subset(all_GWAS_clusters,sign_overall>=2)
genotypes_onlyP=genotypes[all_GWAS_clusters$sign_overall>=2,]

for(i in 1:8){
  pred_traits_onlyP = cbind(pred_traits_onlyP,unlist(apply(genotypes_onlyP,2,function(x){
    sum(x[all_GWAS_clusters_onlyP[,paste0("is_sign_",vect_P_traits_bothEnv[i])]]*all_GWAS_clusters_onlyP[,paste0("beta_",vect_P_traits_bothEnv[i])][all_GWAS_clusters_onlyP[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])})))
}

all_GWAS_clusters_noP = subset(all_GWAS_clusters,sign_overall==1)
genotypes_noP=genotypes[all_GWAS_clusters$sign_overall==1,]
colSums(all_GWAS_clusters_noP[,paste0("is_sign_",vect_P_traits_bothEnv)])
for(i in 1:8){
  pred_traits_noP = cbind(pred_traits_noP,unlist(apply(genotypes_noP,2,function(x){
    sum(x[all_GWAS_clusters_noP[,paste0("is_sign_",vect_P_traits_bothEnv[i])]]*all_GWAS_clusters_noP[,paste0("beta_",vect_P_traits_bothEnv[i])][all_GWAS_clusters_noP[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])})))
}

cor_df=NULL
for(i in 1:8){
  for(j in 1:8){
    if(i!=j){
      cor_df=rbind(cor_df,data.frame(
                   trait1=i,trait2=j,
                   cor_pheno=cor(na.omit(gemma_pheno))[i,j],
                   cor_all=cor(pred_traits)[i,j],
                   cor_noP=cor(pred_traits_noP)[i,j],
                   cor_onlyP=cor(pred_traits_onlyP)[i,j]
                     ))
    }
  }
}
head(cor_df)

pdf("Phenotypic_correlations_pleiotropy_NaCl_R2.pdf",h=5,w=5)
plot(cor_df$cor_pheno,cor_df$cor_all,pch=16,xlab="Phenotypic correlations",ylab="Correlation between predicted phenotypes")
abline(lm(cor_df$cor_all~cor_df$cor_pheno),col="black")
points(cor_df$cor_pheno,cor_df$cor_noP,pch=16,col="green")
abline(lm(cor_df$cor_noP~cor_df$cor_pheno),col="green")
points(cor_df$cor_pheno,cor_df$cor_onlyP,pch=16,col="red")
abline(lm(cor_df$cor_onlyP~cor_df$cor_pheno),col="red")
abline(a=0,b=1,lty=2)
legend(0,-.4,c("All clusters","Pleiotropic clusters","Non pleiotropic clusters"),lwd=4,col=c("black","red","green"),bty="n",cex=.8)
dev.off()


cor_df_betNaCl = subset(cor_df,trait2==16)
dim(cor_df_betNaCl)

cor_df_betNaCl$trait1=vect_P_traits_bothEnv[cor_df_betNaCl$trait1]
cor_df_betNaCl$trait2=vect_P_traits_bothEnv[cor_df_betNaCl$trait2]

cor_df_betNaCl$ratio = cor_df_betNaCl$cor_onlyP/cor_df_betNaCl$cor_all

all_GWAS_clusters_betaNaCl=subset(all_GWAS_clusters,is_sign_beta_NaCl)
dim(all_GWAS_clusters_betaNaCl) # 37

# library(gplots)
# library(limma) 
# vennDiagram(vennCounts(
#   cbind(as.numeric(all_GWAS_clusters$is_sign_beta_NaCl),
#         as.numeric(all_GWAS_clusters$is_sign_area.F_NGM),
#         as.numeric(all_GWAS_clusters$is_sign_area.F_NaCl))),
#   names=c("beta","Trait NGM","Trait NaCl"),main='All')
# 
# vennDiagram(vennCounts(
#   cbind(as.numeric(all_GWAS_clusters$is_sign_beta_NaCl),
#         as.numeric(all_GWAS_clusters$is_sign_SF_NGM),
#         as.numeric(all_GWAS_clusters$is_sign_SF_NaCl))),
#   names=c("beta","Trait NGM","Trait NaCl"),main='All')
# 
# sum(all_GWAS_clusters$is_sign_area.F_NGM)
# sum(all_GWAS_clusters$is_sign_area.F_NaCl)
# 
# sum(all_GWAS_clusters$is_sign_area.F_NGM & all_GWAS_clusters$is_sign_beta_NaCl)/sum(all_GWAS_clusters$is_sign_area.F_NGM)
# sum(all_GWAS_clusters$is_sign_area.F_NaCl & all_GWAS_clusters$is_sign_beta_NaCl)/sum(all_GWAS_clusters$is_sign_area.F_NaCl)
# 
# sum(all_GWAS_clusters$is_sign_area.F_NGM & all_GWAS_clusters$is_sign_beta_NGM)/sum(all_GWAS_clusters$is_sign_area.F_NGM)
# sum(all_GWAS_clusters$is_sign_area.F_NaCl & all_GWAS_clusters$is_sign_beta_NGM)/sum(all_GWAS_clusters$is_sign_area.F_NaCl)
# 
# ### SF
# sum(all_GWAS_clusters$is_sign_SF_NGM & all_GWAS_clusters$is_sign_beta_NaCl)/sum(all_GWAS_clusters$is_sign_SF_NGM)
# sum(all_GWAS_clusters$is_sign_SF_NaCl & all_GWAS_clusters$is_sign_beta_NaCl)/sum(all_GWAS_clusters$is_sign_SF_NaCl)
# 
# sum(all_GWAS_clusters$is_sign_SF_NGM & all_GWAS_clusters$is_sign_beta_NGM)/sum(all_GWAS_clusters$is_sign_SF_NGM)
# sum(all_GWAS_clusters$is_sign_SF_NaCl & all_GWAS_clusters$is_sign_beta_NGM)/sum(all_GWAS_clusters$is_sign_SF_NaCl)
# 
# ### SB
# sum(all_GWAS_clusters$is_sign_SB_NGM & all_GWAS_clusters$is_sign_beta_NaCl)/sum(all_GWAS_clusters$is_sign_SB_NGM)
# sum(all_GWAS_clusters$is_sign_SB_NaCl & all_GWAS_clusters$is_sign_beta_NaCl)/sum(all_GWAS_clusters$is_sign_SB_NaCl)
# 
# sum(all_GWAS_clusters$is_sign_SB_NGM & all_GWAS_clusters$is_sign_beta_NGM)/sum(all_GWAS_clusters$is_sign_SB_NGM)
# sum(all_GWAS_clusters$is_sign_SB_NaCl & all_GWAS_clusters$is_sign_beta_NGM)/sum(all_GWAS_clusters$is_sign_SB_NaCl)
# 
# 
# plot(density(all_GWAS_clusters$beta_beta_NGM[all_GWAS_clusters$is_sign_beta_NGM]))
# lines(density(all_GWAS_clusters$beta_beta_NaCl[all_GWAS_clusters$is_sign_beta_NaCl]),col="red")
# 
# 
# plot(density(all_GWAS_clusters$beta_area.F_NGM[all_GWAS_clusters$is_sign_area.F_NGM]))
# lines(density(all_GWAS_clusters$beta_area.F_NaCl[all_GWAS_clusters$is_sign_area.F_NaCl]),col="red")
# 
# plot(density(all_GWAS_clusters$beta_SF_NGM[all_GWAS_clusters$is_sign_SF_NGM]))
# lines(density(all_GWAS_clusters$beta_SF_NaCl[all_GWAS_clusters$is_sign_SF_NaCl]),col="red")



