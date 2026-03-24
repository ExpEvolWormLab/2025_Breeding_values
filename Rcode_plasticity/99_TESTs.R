rm(list=ls());gc()
library(gplots)
library(limma) 

threshold_clust=0.8
#all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
#load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))

## Plasticity ?
all_cluster_merged=read.table(paste0("txt_BLUEs_plasticity//all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
load(paste0("RData/GWAS_cluster_results_BLUE_plasticity_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))


all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)

dim(all_cluster_merged)
dim(all_GWAS_clusters)

plot(density(subset(all_cluster_merged,is_sign_T1)$beta_T1))
lines(density(subset(all_cluster_merged,is_sign_T2)$beta_T2),col="orange")

boxplot(c(abs(subset(all_cluster_merged,is_sign_T1)$beta_T1),
        abs(subset(all_cluster_merged,is_sign_T2)$beta_T2))~rep(c("NGM","NaCl"),c(sum(all_cluster_merged$is_sign_T1),sum(all_cluster_merged$is_sign_T2))))

summary(abs(subset(all_cluster_merged,is_sign_T1)$beta_T1))
summary(abs(subset(all_cluster_merged,is_sign_T2)$beta_T2))

### Should be done at the cluster level ###

plot(density(subset(all_GWAS_clusters,is_sign_T1)$beta_T1))
lines(density(subset(all_GWAS_clusters,is_sign_T2)$beta_T2),col="orange")

boxplot(c(abs(subset(all_GWAS_clusters,is_sign_T1)$beta_T1),
          abs(subset(all_GWAS_clusters,is_sign_T2)$beta_T2))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_T1),sum(all_GWAS_clusters$is_sign_T2))))


boxplot(c((subset(all_GWAS_clusters,is_sign_T1 & beta_T1<0)$beta_T1),
          (subset(all_GWAS_clusters,is_sign_T2 & beta_T2<0)$beta_T2))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_T1 & all_GWAS_clusters$beta_T1<0),sum(all_GWAS_clusters$is_sign_T2 & all_GWAS_clusters$beta_T2<0))))

boxplot(c((subset(all_GWAS_clusters,is_sign_T1 & beta_T1>0)$beta_T1),
          (subset(all_GWAS_clusters,is_sign_T2 & beta_T2>0)$beta_T2))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_T1 & all_GWAS_clusters$beta_T1>0),sum(all_GWAS_clusters$is_sign_T2 & all_GWAS_clusters$beta_T2>0))))

summary(abs(subset(all_GWAS_clusters,is_sign_T1)$beta_T1))
summary(abs(subset(all_GWAS_clusters,is_sign_T2)$beta_T2))

### Allele frequencies ?
all_GWAS_clusters_T1=subset(all_GWAS_clusters,is_sign_T1)
all_GWAS_clusters_T1$pq.NGM = all_GWAS_clusters_T1$af_T1 * (1-all_GWAS_clusters_T1$af_T1)
plot(all_GWAS_clusters_T1$pq.NGM,abs(all_GWAS_clusters_T1$beta_T1),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T1$beta_T1)+2])
plot(all_GWAS_clusters_T1$pq.NGM,1/(all_GWAS_clusters_T1$beta_T1^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T1$beta_T1)+2])

all_GWAS_clusters_T2=subset(all_GWAS_clusters,is_sign_T2)
plot(all_GWAS_clusters_T2$af_T1,all_GWAS_clusters_T2$beta_T2,pch=16)
all_GWAS_clusters_T2$pq.NGM = all_GWAS_clusters_T2$af_T1 * (1-all_GWAS_clusters_T2$af_T1)
plot(all_GWAS_clusters_T2$pq.NGM,abs(all_GWAS_clusters_T2$beta_T2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T2$beta_T2)+2])
plot(all_GWAS_clusters_T2$pq.NGM,1/(all_GWAS_clusters_T2$beta_T2^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_T2$beta_T2)+2])


test.mod <- lm(I(1/all_GWAS_clusters_T2$beta_T2^2)~pq.NGM + as.factor(sign(all_GWAS_clusters_T2$beta_T2)),data=all_GWAS_clusters_T2)
summary(test.mod)

### All SNPs ?
all_GWAS_clusters$pq.NGM=all_GWAS_clusters$af_T1 * (1-all_GWAS_clusters$af_T1)
plot(all_GWAS_clusters$pq.NGM,1/(all_GWAS_clusters$beta_T1^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters$beta_T1)+2],ylim=c(0.1,2000000),log="y")


### Now plot both T1/T2

all_GWAS_clusters_Area=subset(all_GWAS_clusters,is_sign_T1 | is_sign_T2)
all_GWAS_clusters_Area$pq.NGM = all_GWAS_clusters_Area$af_T1 * (1-all_GWAS_clusters_Area$af_T1)

all_GWAS_clusters_Area$sign_effect_retained = all_GWAS_clusters_Area$beta_T1
all_GWAS_clusters_Area$sign_effect_retained[all_GWAS_clusters_Area$is_sign_T2 & !all_GWAS_clusters_Area$is_sign_T1] = all_GWAS_clusters_Area$beta_T2[all_GWAS_clusters_Area$is_sign_T2 & !all_GWAS_clusters_Area$is_sign_T1]

plot(all_GWAS_clusters_Area$pq.NGM,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T1)+2])
plot(all_GWAS_clusters_Area$pq.NGM,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","red","red")[1+all_GWAS_clusters_Area$is_sign_T1])
plot(all_GWAS_clusters_Area$pq.NGM,1/(all_GWAS_clusters_Area$sign_effect_retained^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T1)+2])

abline(lm(1/(all_GWAS_clusters_Area$sign_effect_retained^2)~all_GWAS_clusters_Area$pq.NGM),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T1)+2])




