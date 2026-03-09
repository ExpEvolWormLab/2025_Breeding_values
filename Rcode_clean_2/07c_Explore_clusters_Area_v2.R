rm(list=ls());gc()
library(gplots)
library(limma) 

threshold_clust=0.6
all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters)

all_cluster_merged$sign_overall = rowSums(all_cluster_merged[,paste0("is_sign_T",1:14)])>0
all_GWAS_clusters$sign_overall = rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:14)])>0
table(all_cluster_merged$sign_overall)
table(all_GWAS_clusters$sign_overall) # only 1911 truly significant with R=0.6 / 3549 with R=0.8

all_GWAS_clusters$ps_min=NA
all_GWAS_clusters$ps_max=NA
## clust_minMax ?
for(i in unique(all_GWAS_clusters$group)){
  temp=subset(all_cluster_merged,group==i)
  
  all_GWAS_clusters[all_GWAS_clusters$group==i,]$ps_min=min(temp$POS)
  all_GWAS_clusters[all_GWAS_clusters$group==i,]$ps_max=max(temp$POS)
  
}

i_Trait = 7
i_chr = "II"

#pdf("Test_CHRII_Area_NGM_R06_v2.pdf",w=12)

all_cluster_merged_CHRI = subset(all_cluster_merged,CHROM==i_chr)
all_cluster_merged_CHRI_clust = subset(all_GWAS_clusters,CHROM==i_chr)

plot(all_cluster_merged_CHRI_clust$ps,all_cluster_merged_CHRI_clust[,paste0("beta_T",i_Trait)],pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3)
for(i in unique(all_cluster_merged_CHRI_clust$group)){
  temp=subset(all_cluster_merged_CHRI_clust,group==i)
  if(temp$is_sign_T7){
  lines(c(temp$ps_min,temp$ps_max),rep(temp[,paste0("beta_T",i_Trait)],2),col="red",lwd=7)
  }
  if(temp$is_sign_T14){
    lines(c(temp$ps_min,temp$ps_max),rep(temp[,paste0("beta_T",i_Trait+7)],2),col="blue",lwd=7)
  }
}

