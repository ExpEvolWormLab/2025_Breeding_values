rm(list=ls())
gc()
library(qqman)
library(data.table)

for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,0.99)){
  print(threshold_clust)
  #for(env_label in c("NGM","NaCl")){
  env_label="NaCl"
    for(traitID in c("beta","SF",'SB','FS','FB','BS','BF',"area.F")){
      traitJ=NULL
      for(chrID in  c("I","II","III","IV","V","X")){
        traitJ = rbind(traitJ,read.table(paste0("GridLMM_clusters_NaCl/R",round(10*threshold_clust),"/gwas_results_R",round(10*threshold_clust),"_",traitID,"_CHR_",chrID,"_",env_label,".txt"),h=TRUE)[,c("X_ID","beta.5","p_value_REML")])
      }
      names(traitJ)[1]="rs"
      traitJ$group = tstrsplit(traitJ$rs,"_")[[1]]
      traitJ$chr = tstrsplit(traitJ$rs,"_")[[2]]
      traitJ$ps = tstrsplit(traitJ$rs,"_")[[3]]
      
      ### Should be removed
      df_clusters=read.table(paste0("data_BLUE_GridLMM_NaCl/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
      df_clusters$rs=df_clusters$clusterID
      dim(traitJ)
      traitJ=merge(traitJ,df_clusters)
      dim(traitJ)
      table(traitJ$chr==traitJ$CHROM)
      table(traitJ$ps==traitJ$POS)
      
      traitJ$p_value_REML.cor <- p.adjust(traitJ$p_value_REML,method="BH")
      
      ## We will use the top 5% as sign threshold - keep this even if less col
      traitJ$is_sign <- traitJ$p_value_REML.cor<0.05
      table(traitJ$is_sign)
      
      concat_GWAS.temp=traitJ[,c("rs","chr","ps","beta.5","p_value_REML","p_value_REML.cor","is_sign","clusterID")] 
      names(concat_GWAS.temp)[4]="beta"
      
      for(p in c(4:7)){
        names(concat_GWAS.temp)[p]=paste0(names(concat_GWAS.temp)[p],"_",traitID,"_",env_label)
      }
      
      if(!exists("all_GWAS")){
        all_GWAS=concat_GWAS.temp
      }else{
        all_GWAS=merge(all_GWAS,concat_GWAS.temp,all=TRUE)
      }
      #print(table(traitJ$rs%in%all_GWAS$rs))
      print(dim(all_GWAS))
    }
    rm(concat_GWAS.temp,p,traitJ)

  all_GWAS_clusters = all_GWAS
  rm(all_GWAS)
  
  save(list=ls(),file=paste0("RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
}




