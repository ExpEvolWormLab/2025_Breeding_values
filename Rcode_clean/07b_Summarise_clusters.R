rm(list=ls());gc()

all_cluster_merged=NULL
summary_clust=NULL
for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,1)){
print(threshold_clust)
all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)


if(!is.null(summary_clust)){
  summary_clust_temp <- data.frame(
    nb_clust=length(unique(all_cluster_merged$group)),
    nb_singletons=sum(table(all_cluster_merged$group)==1),r2=threshold_clust)
  summary_clust=rbind(summary_clust,summary_clust_temp)
}else{
  summary_clust <- data.frame(
    nb_clust=length(unique(all_cluster_merged$group)),
    nb_singletons=sum(table(all_cluster_merged$group)==1),r2=threshold_clust)
    
}
}

pdf("plots_BLUEs/Summary_clusters_R2.pdf")
plot(summary_clust$r2,summary_clust$nb_clust,type="b",xlab="r2 threshold",ylab="Nb clusters/singletons",ylim=c(0,45000))
points(summary_clust$r2,summary_clust$nb_singletons,type="b",col="red")
legend(0.2,40000,c("Total number of clusters (including singletons)","Singletons"),col=c("black","red"),lwd=2)
dev.off()

