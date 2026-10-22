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
  for(traitJ in vect_P_traits){
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
    
    if(trait_ID==1){
      all_GWAS=all_GWAS.temp
    }else{
      all_GWAS=merge(all_GWAS,all_GWAS.temp)
    }
    
  }

### Load the pruning
library("vcfR")
genotype_pruned <- read.vcfR('plink_pruning/pruning08/CeMEE.pruned.0.8.vcf.gz')
## 22 k SNPs kept
head(genotype_pruned@fix)

all_GWAS_08=subset(all_GWAS,rs%in%genotype_pruned@fix[,3])

for(i in 1:8){
  all_GWAS_08[,paste0("is_sign_T",i)] = p.adjust(all_GWAS_08[,paste0("pvalue_T",i)],method="BH")<0.05
  print(table(all_GWAS_08[,paste0("is_sign_T",i)]))
}
sign_all_df <- all_GWAS_08[,names(all_GWAS_08)%in%paste0("is_sign_T",1:7)]
table(rowSums(sign_all_df)>=1) #  200 SNPS at least sign for one trait
table(rowSums(sign_all_df)) # up to 5

genotype_pruned <- read.vcfR('plink_pruning/pruning02/CeMEE.pruned.0.2.vcf.gz')
all_GWAS_02=subset(all_GWAS,rs%in%genotype_pruned@fix[,3])

for(i in 1:8){
  all_GWAS_02[,paste0("is_sign_T",i)] = p.adjust(all_GWAS_02[,paste0("pvalue_T",i)],method="BH")<0.05
  print(table(all_GWAS_02[,paste0("is_sign_T",i)]))
}
sign_all_df <- all_GWAS_02[,names(all_GWAS_02)%in%paste0("is_sign_T",1:7)]
table(rowSums(sign_all_df)>=1) #  43 SNPS only at least sign for one trait
table(rowSums(sign_all_df)) # up to 3

    
    
    
    
    
    
    