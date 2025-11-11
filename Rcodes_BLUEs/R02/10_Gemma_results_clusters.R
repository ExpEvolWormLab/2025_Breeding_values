rm(list=ls())
gc()
library(qqman)


t1_clust=read.table("gemma_results_clusters_BLUE_R02///gwas_results_clusters_BLUE_R02_T7.assoc.txt",h=TRUE)
all_cluster_merged = read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters_R02.txt'),sep='\t',h=TRUE)
t1_mean_no_clust_p <- tapply(all_cluster_merged$p_wald_T7,all_cluster_merged$group,mean)
t1_mean_no_clust_beta <- tapply(all_cluster_merged$beta_T7,all_cluster_merged$group,mean)


df_clusters = data.frame(rs=paste0("rs",as.numeric(names(t1_mean_no_clust_p))),
                         beta_t1_mean=as.numeric(t1_mean_no_clust_beta),
                         p_t1_mean=as.numeric(t1_mean_no_clust_p))
t1_merged <- (merge(t1_clust,df_clusters))
plot(t1_merged$beta,t1_merged$beta_t1_mean,pch=16)
plot(-log10(t1_merged$p_wald),-log10(t1_merged$p_t1_mean),pch=16)

### OK, it looks great, wait for the permutations ###

