rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
library(parallel)
library(qqman)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

trait_ID=0

env_label='NaCl'
### Load all results and concatenate them
  #for(traitJ in vect_P_traits){
traitJ="area.F"
  print(traitJ)
    trait_ID=trait_ID+1
    # concatenate all the CHR
    all_GWAS.temp=NULL
    for(i in c("I","II","III","IV","V","X")){
      all_GWAS.temp=rbind(all_GWAS.temp,read.table(paste0("GridLMM_results_NaCl/gwas_results_",traitJ,"_CHR_",i,"_",env_label,".txt"),h=TRUE)[,c("X_ID","beta.5","p_value_REML")])
      
    }
    
    names(all_GWAS.temp)=c("rs",paste0("beta_T",trait_ID),paste0("pvalue_T",trait_ID))

    all_GWAS.temp$chr = factor(tstrsplit(all_GWAS.temp$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
    all_GWAS.temp$ps = as.numeric(tstrsplit(all_GWAS.temp$rs,"_")[[2]])
    all_GWAS.temp=all_GWAS.temp[order(all_GWAS.temp$chr,all_GWAS.temp$ps),]
    
    ### Load the pruning
    