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

load('RData/Preparatory_GLMM.RData')

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

for(traitID in vect_P_traits){
print(traitID)
all_BLUEs_genotyped_NGM$y=all_BLUEs_genotyped_NGM[,traitID]
all_BLUEs_genotyped_NaCl$y=all_BLUEs_genotyped_NaCl[,traitID]

all_BLUEs_genotyped_NGM$y= scale(all_BLUEs_genotyped_NGM$y)
all_BLUEs_genotyped_NaCl$y= scale(all_BLUEs_genotyped_NaCl$y)


gwas.NGM = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                         data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                        X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 10 ,method = 'REML')

NGM_results <- (gwas.NGM$results)
write.table(NGM_results,paste0("GridLMM_results_mixed_scaled//gwas_results_",traitID,"_NGM.txt"),quote=FALSE,row.names=FALSE)

gwas.NaCl = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str)  + (1|pop_label)+ (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                        data = all_BLUEs_genotyped_NaCl, X = genotype_matrix[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 10 ,method = 'REML')
NaCl_results <- (gwas.NaCl$results)
write.table(NGM_results,paste0("GridLMM_results_mixed_scaled//gwas_results_",traitID,"_NaCl.txt"),quote=FALSE,row.names=FALSE)

}
