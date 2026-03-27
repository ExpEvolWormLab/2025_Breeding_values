rm(list=ls());gc()
library(data.table)
gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)

threshold_clust = 0.6
load(paste0("RData/GWAS_cluster_results_BLUE_Qvalues_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
genotypes <- read.table(file=paste0("data_BLUE_clusters/genotypes_clusters_R",round(10*threshold_clust),".txt"))

table(genotypes$V1==all_GWAS_clusters$rs) # Not in the good order
genotypes=genotypes[order(genotypes$V1),]
all_GWAS_clusters=all_GWAS_clusters[order(all_GWAS_clusters$rs),]
table(genotypes$V1==all_GWAS_clusters$rs) # Fixed

pred_traits=NULL
for(i in 1:14){
pred_traits = cbind(pred_traits,unlist(apply(genotypes[,4:969],2,function(x){
  sum(x[all_GWAS_clusters[,paste0("is_sign_T",i)]]*all_GWAS_clusters[,paste0("beta_T",i)][all_GWAS_clusters[,paste0("is_sign_T",i)]])})))
}

plot(gemma_pheno[,1],pred_traits[,1],pch=16)
abline(lm(pred_traits[,1]~gemma_pheno[,1]))
df.lm=NULL
for(i in 1:14){
  df.lm=rbind(df.lm,coef(lm(pred_traits[,i]~gemma_pheno[,i])))
}
summary(df.lm[,2])

pred_traits_noP=NULL
pred_traits_onlyP=NULL

all_GWAS_clusters$sign_overall = rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:14)])
sum(all_GWAS_clusters$sign_overall>=1)/nrow(all_GWAS_clusters) # 38% significant with R=0.8 / 44% with R=0.6/ 89% in r=0.2

### Sign ~ Length Cluster ?
clust_length = data.frame(clusterID=names(table(all_GWAS_clusters$clusterID)),clust_length=as.numeric(table(all_GWAS_clusters$clusterID)))

all_GWAS_clusters_onlyP = subset(all_GWAS_clusters,sign_overall>=2)
genotypes_onlyP=genotypes[all_GWAS_clusters$sign_overall>=2,]

for(i in 1:14){
  pred_traits_onlyP = cbind(pred_traits_onlyP,unlist(apply(genotypes_onlyP[,4:969],2,function(x){
    sum(x[all_GWAS_clusters_onlyP[,paste0("is_sign_T",i)]]*all_GWAS_clusters_onlyP[,paste0("beta_T",i)][all_GWAS_clusters_onlyP[,paste0("is_sign_T",i)]])})))
}

all_GWAS_clusters_noP = subset(all_GWAS_clusters,sign_overall==1)
genotypes_noP=genotypes[all_GWAS_clusters$sign_overall==1,]

for(i in 1:14){
  pred_traits_noP = cbind(pred_traits_noP,unlist(apply(genotypes_noP[,4:969],2,function(x){
    sum(x[all_GWAS_clusters_noP[,paste0("is_sign_T",i)]]*all_GWAS_clusters_noP[,paste0("beta_T",i)][all_GWAS_clusters_noP[,paste0("is_sign_T",i)]])})))
}

#cor_df=data.frame(trait1=paste0("T",1:14),trait2=paste0("T",1:14),cor_pheno=NA,cor_all=NA,cor_noP=NA,cor_onlyP=NA)
cor_df=NULL
for(i in 1:14){
  for(j in 1:14){
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

plot(cor_df$cor_pheno,cor_df$cor_all,pch=16,xlab="Phenotypic correlations",ylab="Correlation between predicted phenotypes")
abline(lm(cor_df$cor_all~cor_df$cor_pheno),col="black")
points(cor_df$cor_pheno,cor_df$cor_noP,pch=16,col="green")
abline(lm(cor_df$cor_noP~cor_df$cor_pheno),col="green")
points(cor_df$cor_pheno,cor_df$cor_onlyP,pch=16,col="red")
abline(lm(cor_df$cor_onlyP~cor_df$cor_pheno),col="red")
abline(a=0,b=1,lty=2)
