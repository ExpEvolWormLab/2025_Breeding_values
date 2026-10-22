rm(list=ls());gc()
threshold_clust=0.8

load(paste0("RData_GridLMM_NaCl/Area_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)
vect_P_traits=c("area.F")

## Load the selection coefficients

clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Area_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
clusters_selection <- read.table('RData_GridLMM_NaCl/Area_compute_S_parallel_group_R08.txt',h=TRUE)

clusters_selection_all=subset(clusters_selection,type=='all_gen')
clusters_selection_no0=subset(clusters_selection,type=='no0')
rm(clusters_selection)

all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters[,c('group',"clusterID")])

## Add group to the merged_pos
clusters_selection_all=merge(all_GWAS_clusters,clusters_selection_all)
clusters_selection_no0=merge(all_GWAS_clusters,clusters_selection_no0)

plot(clusters_selection_all$s~clusters_selection_all$beta_area.F_NaCl,pch=16)

plot(clusters_selection_all$beta_area.F_NaCl,clusters_selection_all$s,pch=16,cex=.8)
abline(h=0,lty=2,col="red",lwd=2)
abline(lm(s~beta_area.F_NaCl-1,data=clusters_selection_all),lwd=3,col="red")
summary(lm(s~beta_area.F_NaCl-1,data=clusters_selection_all))


plot(clusters_selection_no0$beta_area.F_NaCl,clusters_selection_no0$s,pch=16,cex=.8)
abline(h=0,lty=2,col="red",lwd=2)
abline(lm(s~beta_area.F_NaCl-1,data=clusters_selection_no0),lwd=3,col="red")
summary(lm(s~beta_area.F_NaCl-1,data=clusters_selection_no0))
