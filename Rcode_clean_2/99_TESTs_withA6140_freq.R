rm(list=ls());gc()
library(gplots)
library(limma) 

threshold_clust=0.6
all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))

#######################################
### Load the A6140 allele frequencies ?
#######################################

a6140_freq = read.table("data/A6140_frequencies.txt",sep='\t',h=TRUE)[,c(1:4,6,7)]
names(a6140_freq)[5:6]=c("a6140_freq","a6140_cov")
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_cluster_merged)
dim(all_GWAS_clusters)

## Polarization ##

all_cluster_merged_with_freq <- merge(all_cluster_merged,a6140_freq)
head(all_cluster_merged_with_freq)

### Some are not the same alleles ?
table(all_cluster_merged_with_freq$allele0==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$allele0==all_cluster_merged_with_freq$minor)
table(all_cluster_merged_with_freq$allele1==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$allele1==all_cluster_merged_with_freq$minor)

to_remove = !(all_cluster_merged_with_freq$allele0==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$allele0==all_cluster_merged_with_freq$minor) |
  ! (all_cluster_merged_with_freq$allele1==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$allele1==all_cluster_merged_with_freq$minor)

all_cluster_merged_with_freq=subset(all_cluster_merged_with_freq,!to_remove)
dim(all_cluster_merged_with_freq) # 92.867

## Now re-polarize ##

all_cluster_merged_with_freq[all_cluster_merged_with_freq$allele1==all_cluster_merged_with_freq$minor,]$a6140_freq = 
  1- (all_cluster_merged_with_freq[all_cluster_merged_with_freq$allele1==all_cluster_merged_with_freq$minor,]$a6140_freq)

plot(all_cluster_merged_with_freq$a6140_freq,all_cluster_merged_with_freq$af.NGM,xlab="Freq in A6140",ylab="Freq in CeMEE",
     pch=16,col=adjustcolor("grey",alpha=.3))

## Group frequencies
group_Freq <- tapply(all_cluster_merged_with_freq$a6140_freq,all_cluster_merged_with_freq$group,mean,na.rm=TRUE)
group_Freq=data.frame(group=names(group_Freq),a6140_freq=1-as.numeric(group_Freq))
dim(merge(all_GWAS_clusters,group_Freq)) # 4201 with frequencies
all_GWAS_clusters_with_freq = merge(all_GWAS_clusters,group_Freq)
plot(all_GWAS_clusters_with_freq$a6140_freq,all_GWAS_clusters_with_freq$af.NGM)


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
all_GWAS_clusters_with_freq_T7=subset(all_GWAS_clusters_with_freq,is_sign_T7)
all_GWAS_clusters_with_freq_T7$pq.A6140 = all_GWAS_clusters_with_freq_T7$a6140_freq * (1-all_GWAS_clusters_with_freq_T7$a6140_freq)
plot(all_GWAS_clusters_with_freq_T7$pq.A6140,abs(all_GWAS_clusters_with_freq_T7$beta_T7),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_T7$beta_T7)+2])
plot(all_GWAS_clusters_with_freq_T7$pq.A6140,1/(all_GWAS_clusters_with_freq_T7$beta_T7^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_T7$beta_T7)+2])


all_GWAS_clusters_with_freq_T14=subset(all_GWAS_clusters_with_freq,is_sign_T14)
plot(all_GWAS_clusters_with_freq_T14$a6140_freq,all_GWAS_clusters_with_freq_T14$beta_T14,pch=16)
all_GWAS_clusters_with_freq_T14$pq.A6140 = all_GWAS_clusters_with_freq_T14$a6140_freq * (1-all_GWAS_clusters_with_freq_T14$a6140_freq)
plot(all_GWAS_clusters_with_freq_T14$pq.A6140,abs(all_GWAS_clusters_with_freq_T14$beta_T14),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_T14$beta_T14)+2])
plot(all_GWAS_clusters_with_freq_T14$pq.A6140,1/(all_GWAS_clusters_with_freq_T14$beta_T14^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_T14$beta_T14)+2])

test.mod <- lm(I(1/all_GWAS_clusters_with_freq_T14$beta_T14^2) ~ pq.A6140 + as.factor(sign(all_GWAS_clusters_with_freq_T14$beta_T14)),data=all_GWAS_clusters_with_freq_T14)
summary(test.mod)

### All SNPs ?
all_GWAS_clusters_with_freq$pq.A6140=all_GWAS_clusters_with_freq$a6140_freq * (1-all_GWAS_clusters_with_freq$a6140_freq)
plot(all_GWAS_clusters_with_freq$pq.A6140,1/(all_GWAS_clusters_with_freq$beta_T7^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq$beta_T7)+2],ylim=c(0.1,200000000),log="y")

### Now plot both T7/T14

all_GWAS_clusters_Area=subset(all_GWAS_clusters_with_freq,is_sign_T7 | is_sign_T14)
all_GWAS_clusters_Area$pq.A6140 = all_GWAS_clusters_Area$a6140_freq * (1-all_GWAS_clusters_Area$a6140_freq)

all_GWAS_clusters_Area$sign_effect_retained = all_GWAS_clusters_Area$beta_T7
all_GWAS_clusters_Area$sign_effect_retained[all_GWAS_clusters_Area$is_sign_T14 & !all_GWAS_clusters_Area$is_sign_T7] = all_GWAS_clusters_Area$beta_T14[all_GWAS_clusters_Area$is_sign_T14 & !all_GWAS_clusters_Area$is_sign_T7]

plot(all_GWAS_clusters_Area$pq.A6140,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T7)+2])
plot(all_GWAS_clusters_Area$pq.A6140,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","red","red")[1+all_GWAS_clusters_Area$is_sign_T7])
plot(all_GWAS_clusters_Area$pq.A6140,1/(all_GWAS_clusters_Area$sign_effect_retained^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_T7)+2])




