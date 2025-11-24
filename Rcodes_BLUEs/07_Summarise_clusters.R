rm(list=ls());gc()

all_cluster_merged=read.table('txt_BLUEs/all_sign_GWAS_with_clusters_R02.txt',sep='\t',h=TRUE)

summary_clust <- data.frame(
nb_clust=length(unique(all_cluster_merged$group)),
nb_singletons=sum(table(all_cluster_merged$group)==1),
r2=0.2)

all_cluster_merged=read.table('txt_BLUEs/all_sign_GWAS_with_clusters_R04.txt',sep='\t',h=TRUE)
summary_clust <- rbind(summary_clust,data.frame(
  nb_clust=length(unique(all_cluster_merged$group)),
  nb_singletons=sum(table(all_cluster_merged$group)==1),
  r2=0.4))

all_cluster_merged=read.table('txt_BLUEs/all_sign_GWAS_with_clusters_R06.txt',sep='\t',h=TRUE)
summary_clust <- rbind(summary_clust,data.frame(
  nb_clust=length(unique(all_cluster_merged$group)),
  nb_singletons=sum(table(all_cluster_merged$group)==1),
  r2=0.6))
  
all_cluster_merged=read.table('txt_BLUEs/all_sign_GWAS_with_clusters.txt',sep='\t',h=TRUE)
summary_clust <- rbind(summary_clust,data.frame(
  nb_clust=length(unique(all_cluster_merged$group)),
  nb_singletons=sum(table(all_cluster_merged$group)==1),
  r2=0.8))
  
all_cluster_merged=read.table('txt_BLUEs/all_sign_GWAS_with_clusters_R09.txt',sep='\t',h=TRUE)
summary_clust <- rbind(summary_clust,data.frame(
  nb_clust=length(unique(all_cluster_merged$group)),
  nb_singletons=sum(table(all_cluster_merged$group)==1),
  r2=0.9))

all_cluster_merged=read.table('txt_BLUEs/all_sign_GWAS_with_clusters_R10.txt',sep='\t',h=TRUE)
summary_clust <- rbind(summary_clust,data.frame(
  nb_clust=length(unique(all_cluster_merged$group)),
  nb_singletons=sum(table(all_cluster_merged$group)==1),
  r2=0.99))

summary_clust <- rbind(summary_clust,data.frame(
  nb_clust=nrow(all_GWAS_sign),
  nb_singletons=nrow(all_GWAS_sign),
  r2=1))


plot(summary_clust$r2,summary_clust$nb_clust,type="b",xlab="r2 threshold",ylab="Nb clusters/singletons",ylim=c(0,130000))
points(summary_clust$r2,summary_clust$nb_singletons,type="b",col="red")

plot(summary_clust$r2,summary_clust$nb_clust,type="b",xlab="r2 threshold",ylab="Nb clusters/singletons",ylim=c(0,45500))
points(summary_clust$r2,summary_clust$nb_singletons,type="b",col="red")








