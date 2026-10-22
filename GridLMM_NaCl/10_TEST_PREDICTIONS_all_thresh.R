rm(list=ls());gc()
library(data.table)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")
vect_P_traits_bothEnv <- paste0(vect_P_traits,"_NaCl")

load('RData/Gmatrices_P.RData')
all_traits_corrected=all_traits_corrected[order(all_traits_corrected$pop_label),]
#pred_traits=pred_traits[order(rownames(pred_traits)),]
#table(rownames(pred_traits)==all_traits_corrected$pop_label)
gemma_pheno=all_traits_corrected[,vect_P_traits_bothEnv]
all.trait.pred=NULL
for(threshold_clust in c(0.2,0.8,0.4,0.6,0.9,1)){
  #threshold_clust = 0.2
  load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
  genotypes <- read.table(file=paste0("data_BLUE_GridLMM_NaCl//genotypes_matrix_clusters_R",round(10*threshold_clust),".txt"))
  
  table(row.names(genotypes)==all_GWAS_clusters$rs) # OK ?
  genotypes=genotypes[order(row.names(genotypes)),]
  all_GWAS_clusters=all_GWAS_clusters[order(all_GWAS_clusters$rs),]
  table(row.names(genotypes)==all_GWAS_clusters$rs)  # Fixed
  
  pred_traits=NULL
  for(i in 1:8){
  pred_traits = cbind(pred_traits,unlist(apply(genotypes,2,function(x){
    sum(x[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])]]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])][all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])})))
  }
  df.lm=NULL
  for(i in 1:8){
    temp.lm=lm(pred_traits[,i]~gemma_pheno[,i])
    df.lm=rbind(df.lm,c(coef(temp.lm),summary(temp.lm)$r.squared))
  }
  df.lm=data.frame(df.lm)
  names(df.lm)=c("Intercept","Slope")
  df.lm$thresh=threshold_clust
  df.lm$trait=paste0("T",1:8)
  
  all.trait.pred=rbind(all.trait.pred,df.lm)
}

boxplot(all.trait.pred$Slope~all.trait.pred$thresh)
abline(h=1,lwd=2,col="red")

area.pred = subset(all.trait.pred,trait=="T7")
area.pred=area.pred[order(area.pred$thresh),]
plot(area.pred$thresh,area.pred$Slope,type='b')
abline(h=1,lwd=2,col="red")
tapply(all.trait.pred$Slope,all.trait.pred$thresh,median)
tapply(all.trait.pred$Slope,all.trait.pred$thresh,mean)

