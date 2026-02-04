rm(list=ls());gc()
library(gplots)
library(limma) 

threshold_clust=0.6
all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)

dim(all_cluster_merged)
dim(all_GWAS_clusters)

plot(density(subset(all_cluster_merged,is_sign_T7)$beta_T7))
lines(density(subset(all_cluster_merged,is_sign_T14)$beta_T14),col="orange")

boxplot(c(abs(subset(all_cluster_merged,is_sign_T7)$beta_T7),
        abs(subset(all_cluster_merged,is_sign_T14)$beta_T14))~rep(c("NGM","NaCl"),c(sum(all_cluster_merged$is_sign_T7),sum(all_cluster_merged$is_sign_T14))))

summary(abs(subset(all_cluster_merged,is_sign_T7)$beta_T7))
summary(abs(subset(all_cluster_merged,is_sign_T14)$beta_T14))

### Should be done at the cluster level ###

plot(density(subset(all_GWAS_clusters,is_sign_T7)$beta_T7))
lines(density(subset(all_GWAS_clusters,is_sign_T14)$beta_T14),col="orange")

boxplot(c(abs(subset(all_GWAS_clusters,is_sign_T7)$beta_T7),
          abs(subset(all_GWAS_clusters,is_sign_T14)$beta_T14))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_T7),sum(all_GWAS_clusters$is_sign_T14))))


boxplot(c((subset(all_GWAS_clusters,is_sign_T7 & beta_T7<0)$beta_T7),
          (subset(all_GWAS_clusters,is_sign_T14 & beta_T14<0)$beta_T14))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_T7 & all_GWAS_clusters$beta_T7<0),sum(all_GWAS_clusters$is_sign_T14 & all_GWAS_clusters$beta_T14<0))))

boxplot(c((subset(all_GWAS_clusters,is_sign_T7 & beta_T7>0)$beta_T7),
          (subset(all_GWAS_clusters,is_sign_T14 & beta_T14>0)$beta_T14))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_T7 & all_GWAS_clusters$beta_T7>0),sum(all_GWAS_clusters$is_sign_T14 & all_GWAS_clusters$beta_T14>0))))

summary(abs(subset(all_GWAS_clusters,is_sign_T7)$beta_T7))
summary(abs(subset(all_GWAS_clusters,is_sign_T14)$beta_T14))


### Allele frequencies ?
all_GWAS_clusters_T7=subset(all_GWAS_clusters,is_sign_T7)
all_GWAS_clusters_T7$pq.NGM = all_GWAS_clusters_T7$af.NGM * (1-all_GWAS_clusters_T7$af.NGM)
plot(all_GWAS_clusters_T7$pq.NGM,abs(all_GWAS_clusters_T7$beta_T7),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T7$beta_T7)+2])
plot(all_GWAS_clusters_T7$pq.NGM,1/(all_GWAS_clusters_T7$beta_T7^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T7$beta_T7)+2])


all_GWAS_clusters_T14=subset(all_GWAS_clusters,is_sign_T14)
plot(all_GWAS_clusters_T14$af.NGM,all_GWAS_clusters_T14$beta_T14,pch=16)
all_GWAS_clusters_T14$pq.NGM = all_GWAS_clusters_T14$af.NGM * (1-all_GWAS_clusters_T14$af.NGM)
plot(all_GWAS_clusters_T14$pq.NGM,abs(all_GWAS_clusters_T14$beta_T14),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T14$beta_T14)+2])
plot(all_GWAS_clusters_T14$pq.NGM,1/(all_GWAS_clusters_T14$beta_T14^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T14$beta_T14)+2])


test.mod <- lm(I(1/all_GWAS_clusters_T14$beta_T14^2)~pq.NGM + as.factor(sign(all_GWAS_clusters_T14$beta_T14)),data=all_GWAS_clusters_T14)
summary(test.mod)

### All SNPs ?
all_GWAS_clusters$pq.NGM=all_GWAS_clusters$af.NGM * (1-all_GWAS_clusters$af.NGM)
plot(all_GWAS_clusters$pq.NGM,1/(all_GWAS_clusters$beta_T7^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters$beta_T7)+2],ylim=c(0.1,2000000),log="y")


### Now plot both T7/T14

all_GWAS_clusters_Area=subset(all_GWAS_clusters,is_sign_T7 | is_sign_T14)
all_GWAS_clusters_Area$pq.NGM = all_GWAS_clusters_Area$af.NGM * (1-all_GWAS_clusters_Area$af.NGM)

all_GWAS_clusters_Area$sign_effect_retained = all_GWAS_clusters_Area$beta_T7
all_GWAS_clusters_Area$sign_effect_retained[all_GWAS_clusters_Area$is_sign_T14 & !all_GWAS_clusters_Area$is_sign_T7] = all_GWAS_clusters_Area$beta_T14[all_GWAS_clusters_Area$is_sign_T14 & !all_GWAS_clusters_Area$is_sign_T7]

plot(all_GWAS_clusters_Area$pq.NGM,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T7)+2])
plot(all_GWAS_clusters_Area$pq.NGM,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","red","red")[1+all_GWAS_clusters_Area$is_sign_T7])
plot(all_GWAS_clusters_Area$pq.NGM,1/(all_GWAS_clusters_Area$sign_effect_retained^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T7)+2])




