rm(list=ls());gc()

#threshold_clust=0.2
#load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
#rm(chrID,traitID)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

# load the poolseq data now
poolseq_env=new.env()
load(file='SaveGLM.RData',envir=poolseq_env)
populations=poolseq_env$populations
all.af.caendr=poolseq_env$all.af
all.cov.caendr=poolseq_env$all.cov
all.counts.caendr=poolseq_env$all.counts
all.counts.alt.caendr=poolseq_env$all.counts.alt
snps_positions.sync.caendr=poolseq_env$snps_positions.sync


### 
head(all.af.caendr[,1:5])
snp_block=c(6637695, 6322080,6826421,6046042, 6073853,6138448,6992729,7014705,7278103)
snp_all=all.af.caendr[row.names(all.af.caendr)%in%paste0("V.",snp_block),]
dim(snp_all) # 9 SNPs

plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  tmp=populations$replicate_fact%in%c(1,i)
  lines(snp_all[1,tmp]~populations$generations[tmp])
}

for(p in 1:9){
  plot(snp_all[p,]~populations$generations,type='n',ylim=c(0,1))
  for(i in 3:6){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all[p,tmp]~populations$generations[tmp])
  }
  
}
plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  for(p in 1:9){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}

clust_s_est=read.table("Beta_clust_s_est_08.txt",h=TRUE,sep='\t')
subset(clust_s_est,group%in%c(418,415,419,411,421,423))



