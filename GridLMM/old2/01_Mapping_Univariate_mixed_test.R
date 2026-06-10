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


traitID="area.F"
print(traitID)
all_BLUEs_genotyped_NGM$y=all_BLUEs_genotyped_NGM[,traitID]
all_BLUEs_genotyped_NaCl$y=all_BLUEs_genotyped_NaCl[,traitID]

all_BLUEs_genotyped_NGM$y= scale(all_BLUEs_genotyped_NGM$y)
all_BLUEs_genotyped_NaCl$y= scale(all_BLUEs_genotyped_NaCl$y)

# Prune for faster results
genotype_pruned_02 <- read.vcfR('plink_pruning/pruning02/CeMEE.pruned.0.2.vcf.gz')

table(colnames(genotype_matrix)%in%genotype_pruned_02@fix[,3])
genotype_matrix=genotype_matrix[,colnames(genotype_matrix)%in%genotype_pruned_02@fix[,3]]
dim(genotype_matrix)
is.NA=colSums(is.na(genotype_matrix))

gwas.NGM = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                         data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                        X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

test.glm <- lmer(y~1+logD + rel_humidity + temperature + (1|date_str)  + (1|session_id),data=all_BLUEs_genotyped_NGM)
all_BLUEs_genotyped_NGM$y_corrected_residuals <- resid(test.glm)

gwas.NGM.2 = GridLMM_GWAS(formula = y_corrected_residuals ~ 1 + (1|pop_label),
                        test_formula = ~1, reduced_formula = ~0,
                        data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                        X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM.3 = GridLMM_GWAS(formula = y_corrected_residuals ~ 1 + (1|pop_label),
                          test_formula = ~1, reduced_formula = ~1,
                          data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                          X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM.4 = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                          test_formula = ~ 1 + logD + rel_humidity + temperature,
                          reduced_formula = ~0,
                          data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                          X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

#gwas.NGM.5 = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
#                          test_formula = ~ 1 + logD + rel_humidity + temperature,
#                          reduced_formula = ~1 +logD + rel_humidity + temperature + (1|date_str) + (1|session_id),
#                          data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
#                          X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

plot(gwas.NGM$results$beta.1,gwas.NGM.4$results$beta.1,pch=16)
plot(gwas.NGM$results$beta.1,gwas.NGM.2$results$beta.1,pch=16)
plot(gwas.NGM$results$beta.2,gwas.NGM.2$results$beta.2,pch=16)

plot(-log10(gwas.NGM$results$p_value_REML),-log10(gwas.NGM.4$results$p_value_REML.1),pch=16)
abline(a=0,b=1)
# Load the Gemma results
gemma_T7=(read.table("gemma_R02/gwas_results_BLUE_R02_T7.assoc.txt",h=TRUE))
gemma_T7=subset(gemma_T7,rs%in%gwas.NGM$results$X_ID)
table(gemma_T7$rs==gwas.NGM$results$X_ID)

plot(gwas.NGM$results$beta.2,gemma_T7$beta,pch=16)

plot(gwas.NGM$results$beta.5,gemma_T7$beta,pch=16)
plot(gwas.NGM.2$results$beta.2,gemma_T7$beta,pch=16)

plot(-log10(gwas.NGM$results$p_value_REML),-log10(gemma_T7$p_wald),pch=16)
abline(a=0,b=1)
plot(-log10(gwas.NGM.2$results$p_value_REML),-log10(gemma_T7$p_wald),pch=16)





