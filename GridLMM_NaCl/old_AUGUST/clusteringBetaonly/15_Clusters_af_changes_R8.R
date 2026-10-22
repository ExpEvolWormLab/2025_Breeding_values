rm(list=ls());gc()
threshold_clust=0.8

load(paste0("RData_GridLMM_NaCl/Beta_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)
vect_P_traits=c("area.F")

## Load the selection coefficients

clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
clusters_selection <- read.table('RData_GridLMM_NaCl/Beta_compute_S_parallel_group_R08.txt',h=TRUE)

clusters_selection_all=subset(clusters_selection,type=='all_gen')
clusters_selection_no0=subset(clusters_selection,type=='no0')
rm(clusters_selection)

all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters[,c('group',"clusterID")])

## Add group to the merged_pos
clusters_selection_all=merge(all_GWAS_clusters,clusters_selection_all)
clusters_selection_no0=merge(all_GWAS_clusters,clusters_selection_no0)

plot(clusters_selection_all$beta_beta_NaCl,clusters_selection_all$s,pch=16,cex=.8)
abline(h=0,lty=2,col="red",lwd=2)
abline(lm(s~beta_beta_NaCl-1,data=clusters_selection_all),lwd=3,col="red")
summary(lm(s~beta_beta_NaCl-1,data=clusters_selection_all))

plot(clusters_selection_no0$beta_beta_NaCl,clusters_selection_no0$s,pch=16,cex=.8)
abline(h=0,lty=2,col="red",lwd=2)
abline(lm(s~beta_beta_NaCl-1,data=clusters_selection_no0),lwd=3,col="red")
summary(lm(s~beta_beta_NaCl-1,data=clusters_selection_no0))

clusters_selection_no0$nb_sign = rowSums(clusters_selection_no0[,paste0("P_GA",c(1,2,4))]<0.05)
clusters_selection_all$nb_sign = rowSums(clusters_selection_all[,paste0("P_GA",c(1,2,4))]<0.05)

table(clusters_selection_no0$nb_sign)
table(clusters_selection_all$nb_sign)

subset(clusters_selection_all,group%in%c(411,414,415,418,419,420,421,422,423))
subset(clusters_selection_no0,group%in%c(411,414,415,418,419,420,421,422,423))

sort(c(unique(subset(clusters_selection_all,s<0 & beta_beta_NaCl<(-0.15))$group),
unique(subset(clusters_selection_all,s>0 & beta_beta_NaCl>(0.15))$group)))

subset(clusters_selection_no0,group%in%c(23:37))
subset(clusters_selection_all,group%in%c(23:37))
