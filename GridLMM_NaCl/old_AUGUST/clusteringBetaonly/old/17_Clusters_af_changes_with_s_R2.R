rm(list=ls());gc()
library(data.table)
threshold_clust=0.2
load(paste0("RData_GridLMM_NaCl/Beta_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)
vect_P_traits=c("beta")

s_est.env=new.env()
load('RData_GridLMM_NaCl/Compute_S_parallel.RData',envir=s_est.env)
ls(s_est.env)
merged_pos=s_est.env$merged_pos
s_estimates=s_est.env$results
rm(s_est.env);gc()
clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)

## Fix the rs label
s_estimates$rs_dot=s_estimates$rs
s_estimates$rs=paste(tstrsplit(s_estimates$rs_dot,"[.]")[[1]],tstrsplit(s_estimates$rs_dot,"[.]")[[2]],sep='_')
## Add group to the merged_pos
merged_pos_clust=merge(merged_pos,unique(clusters_positions[,c("rs","group")]))
s_estimates_clust=merge(s_estimates,unique(clusters_positions[,c("rs","group")]))
#s_estimates_clust = s_estimates[merged_pos$rs%in%merged_pos_clust$rs,]
dim(merged_pos_clust) # 31k
dim(s_estimates_clust)
vect_s_est=c("P","s","P_GA1","s_GA1","P_GA2","s_GA2","P_GA3","s_GA3","P_GA4","s_GA4")

### Can we compute the s of the clusters ? -> Yes but we need to re-polarize
clust_pol=read.table("data_clustPolarisation_NaCl/Beta_genotypes_clusters_pol_R2.txt",h=TRUE)
names(clust_pol)[1]="rs"
names(clust_pol)[2:3]=c("ref.GWAS","alt.GWAS")

merged_pos=merged_pos[,!names(merged_pos)%in%c("ref.GWAS","alt.GWAS")]
merged_pos=merge(merged_pos,clust_pol)

pol_consistent = ((merged_pos$major.SYNC==merged_pos$ref.GWAS) & (merged_pos$minor.SYNC==merged_pos$alt.GWAS))
pol_reversed = ((merged_pos$minor.SYNC==merged_pos$ref.GWAS) & (merged_pos$major.SYNC==merged_pos$alt.GWAS))
incompatible = !(pol_consistent | pol_reversed)
table(incompatible) # all FALSE

merged_pos$reverse=pol_reversed

s_estimates_clust=merge(s_estimates_clust,merged_pos[,c("rs","reverse")])
for(i in c("s",paste0("s_GA",1:4))){
  s_estimates_clust[,i]=as.numeric(s_estimates_clust[,i])
}
## Change reverse
s_estimates_clust$s[s_estimates_clust$reverse]= -s_estimates_clust$s[s_estimates_clust$reverse]
s_estimates_clust$s_GA1[s_estimates_clust$reverse]= -s_estimates_clust$s_GA1[s_estimates_clust$reverse]
s_estimates_clust$s_GA2[s_estimates_clust$reverse]= -s_estimates_clust$s_GA2[s_estimates_clust$reverse]
s_estimates_clust$s_GA3[s_estimates_clust$reverse]= -s_estimates_clust$s_GA3[s_estimates_clust$reverse]
s_estimates_clust$s_GA4[s_estimates_clust$reverse]= -s_estimates_clust$s_GA4[s_estimates_clust$reverse]

### Can we compute the s of the clusters ?
clust.s.est=matrix(NA,nrow=length(unique(merged_pos_clust$group)),ncol=length(vect_s_est))
k=0
for(i in vect_s_est){
  k=k+1
  s_estimates_clust[,i]=as.numeric(s_estimates_clust[,i])
  clust.s.est[,k]=as.numeric(tapply(s_estimates_clust[,i],s_estimates_clust$group,mean,na.rm=TRUE))
}
clust.s.est=as.data.frame(clust.s.est)
names(clust.s.est)=vect_s_est
clust.s.est$group=names(tapply(s_estimates_clust[,i],s_estimates_clust$group,mean,na.rm=TRUE))
head(clust.s.est)

### Compute the s estimates from the allele frequency changes ?
plot(clust.s.est$s,-log10(clust.s.est$P),pch=16)
write.table(clust.s.est,file='Beta_clust_s_est_02.txt',sep="\t",row.names=FALSE,quote=FALSE)

