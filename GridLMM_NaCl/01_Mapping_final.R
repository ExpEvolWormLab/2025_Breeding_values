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

load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData')

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

CHR_genotype_matrix = tstrsplit(colnames(genotype_matrix),"_")[[1]]
CHR_genotype_matrix_pruned = tstrsplit(colnames(genotype_pruned_matrix),"_")[[1]]

## Reduced the pruned genotype matrix

for(traitID in c("beta","SF",'SB','FS','FB','BS','BF',"area.F")){
  print(traitID)
  
  all_BLUEs_genotyped_NaCl$y=all_BLUEs_genotyped_NaCl[,traitID]
  all_BLUEs_genotyped_NaCl$y= scale(all_BLUEs_genotyped_NaCl$y)
  
  for(chrID in c("I","II","III","IV","V","X")){
    
    genotype_matrix_CHR = genotype_matrix[,CHR_genotype_matrix==chrID]
    CHR_genotype_matrix_pruned_CHR=genotype_pruned_matrix[,!CHR_genotype_matrix_pruned==chrID]
    GRM_CHR <- Gmatrix(CHR_genotype_matrix_pruned_CHR, method="VanRaden") 

    is.NA=colSums(is.na(genotype_matrix_CHR))
  
    gwas.NaCl.CHR = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                             test_formula = ~1, reduced_formula = ~0,
                             data = all_BLUEs_genotyped_NaCl, X = genotype_matrix_CHR[,is.NA==0],
                             X_ID = 'pop_label', relmat = list(pop_label = GRM_CHR), mc.cores = 10 ,method = 'REML')
    
    NaCl_results <- (gwas.NaCl.CHR$results)
      write.table(NaCl_results,paste0("GridLMM_results_NaCl/gwas_results_",traitID,"_CHR_",chrID,"_NaCl.txt"),quote=FALSE,row.names=FALSE)
    
  }
  
}