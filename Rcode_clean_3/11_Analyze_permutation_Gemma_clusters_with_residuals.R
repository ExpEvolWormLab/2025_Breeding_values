rm(list=ls())
gc()
library(qqman)
library(data.table)

for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,1)){
  print(threshold_clust)
  for(j in 1:14){
    
    traitJ = read.table(paste0("gemma_results_BLUE_clusters/gwas_results_clusters_BLUE_R",round(10*threshold_clust),"_T",j,".assoc.txt"),h=TRUE)
    traitJ$group = tstrsplit(traitJ$rs,"_")[[1]]
    traitJ$chr = tstrsplit(traitJ$rs,"_")[[2]]
    traitJ$ps = tstrsplit(traitJ$rs,"_")[[3]]
    
    ### Should be removed
    df_clusters=read.table(paste0("data_BLUE_clusters/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
    df_clusters$rs=df_clusters$clusterID
    dim(traitJ)
    traitJ=merge(traitJ,df_clusters)
    dim(traitJ)
    table(traitJ$chr==traitJ$CHROM)
    table(traitJ$ps==traitJ$POS)
    
    # Genome wide significant threshold after fdr correction
    traitJ$p_wald_cor=p.adjust(traitJ$p_wald,method='BH')
    ## Take 10%
    #plot(-log10(traitJ$p_wald_cor))
    #abline(h=-log10(0.3))
    traitJ$is_sign <- traitJ$p_wald_cor<0.2
    table(traitJ$is_sign)
    
    concat_GWAS.temp=traitJ[,c("rs","chr","ps","n_miss","allele1","allele0","af","beta","se","p_wald","is_sign","clusterID")] 
    
    for(p in c(4,7:11)){
      names(concat_GWAS.temp)[p]=paste0(names(concat_GWAS.temp)[p],"_T",j)
    }
    
    #if(j<=7){
    #  names(concat_GWAS.temp)[4]="n_miss.NGM"
    #  names(concat_GWAS.temp)[7]="af.NGM"
    #}else{
    #  names(concat_GWAS.temp)[4]="n_miss.NaCl"
    #  names(concat_GWAS.temp)[7]="af.NaCl"
    #}
    
    if(!exists("all_GWAS")){
      all_GWAS=concat_GWAS.temp
    }else{
      all_GWAS=merge(all_GWAS,concat_GWAS.temp,all=TRUE)
    }
    #print(table(traitJ$rs%in%all_GWAS$rs))
    print(dim(all_GWAS))
  }
  rm(concat_GWAS.temp,concat.perm.sorted,j,p,traitJ,nb_NA,list_length)
  
  all_GWAS_clusters = all_GWAS
  rm(all_GWAS)
  #head(all_GWAS_clusters)
  
  
  save(list=ls(),file=paste0("RData/GWAS_cluster_results_BLUE_Qvalues_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
}





