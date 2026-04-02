rm(list=ls());gc()
library(vcfR)
library(ggplot2)
#install.packages("~/Downloads/rutilstimflutre-0.170.0.tar.gz")
library(rutilstimflutre)
library(ASRgenomics)
#devtools::install_github('deruncie/GridLMM')
library(GridLMM)
library(data.table)
library(matrixStats)
library(performance)
library(lme4)
library(AGHmatrix)
library(qqman)

temp.env=new.env()
load('RData/Preparatory_GLMM.RData',envir=temp.env)
genotype_pruned_matrix=temp.env$genotype_pruned_matrix
all_BLUEs_genotyped_NGM=temp.env$all_BLUEs_genotyped_NGM
all_BLUEs_genotyped_NaCl=temp.env$all_BLUEs_genotyped_NaCl
rm(temp.env);gc()

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")
CHR_genotype_matrix_pruned = tstrsplit(colnames(genotype_pruned_matrix),"_")[[1]]

## Reduced the pruned genotype matrix
for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,1)){
  
  # Load the genotype matrix
  genotyped_matrix <- read.table(paste0("data_BLUE_clusters_CHR/genotypes_matrix_clusters_R",round(10*threshold_clust),".txt"))
  CHR_genotype_matrix = tstrsplit(rownames(genotyped_matrix),"_")[[2]]
  
  for(traitID in c("beta","SF",'SB','FS','FB','BS','BF',"area.F")){
    print(traitID)
    all_BLUEs_genotyped_NGM$y=all_BLUEs_genotyped_NGM[,traitID]
    all_BLUEs_genotyped_NaCl$y=all_BLUEs_genotyped_NaCl[,traitID]
    
    all_BLUEs_genotyped_NGM$y= scale(all_BLUEs_genotyped_NGM$y)
    all_BLUEs_genotyped_NaCl$y= scale(all_BLUEs_genotyped_NaCl$y)
    
    for(chrID in c("I","II","III","IV","V","X")){
      
      genotype_matrix_CHR = t(genotyped_matrix[CHR_genotype_matrix==chrID,])
      CHR_genotype_matrix_pruned_CHR=genotype_pruned_matrix[,!CHR_genotype_matrix_pruned==chrID]
      GRM_CHR <- Gmatrix(CHR_genotype_matrix_pruned_CHR, method="VanRaden") 
      
      is.NA=colSums(is.na(genotype_matrix_CHR))
      
      gwas.NGM.CHR = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                                  test_formula = ~1, reduced_formula = ~0,
                                  data = all_BLUEs_genotyped_NGM, X = genotype_matrix_CHR[,is.NA==0],
                                  X_ID = 'pop_label', relmat = list(pop_label = GRM_CHR), mc.cores = 10 ,method = 'REML')
      
      NGM_results <- (gwas.NGM.CHR$results)
      write.table(NGM_results,paste0("GridLMM_clusters_perCHR/gwas_results_",traitID,"_CHR_",chrID,"_NGM.txt"),quote=FALSE,row.names=FALSE)
      
      gwas.NaCl.CHR = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                                   test_formula = ~1, reduced_formula = ~0,
                                   data = all_BLUEs_genotyped_NaCl, X = genotype_matrix_CHR[,is.NA==0],
                                   X_ID = 'pop_label', relmat = list(pop_label = GRM_CHR), mc.cores = 10 ,method = 'REML')
      
      NaCl_results <- (gwas.NaCl.CHR$results)
      write.table(NaCl_results,paste0("GridLMM_clusters_perCHR/gwas_results_",traitID,"_CHR_",chrID,"_NaCl.txt"),quote=FALSE,row.names=FALSE)
      
    }
    
  }
}