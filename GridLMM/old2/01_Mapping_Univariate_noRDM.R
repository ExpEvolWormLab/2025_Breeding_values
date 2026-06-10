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

load('RData/Preparatory_GLMM_scaled.RData')

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F")

traitID="SF"

all_BLUEs_genotyped_NGM$y=all_BLUEs_genotyped_NGM[,traitID]
all_BLUEs_genotyped_NaCl$y=all_BLUEs_genotyped_NaCl[,traitID]

all_BLUEs_genotyped_NGM$y= scale(all_BLUEs_genotyped_NGM$y)
all_BLUEs_genotyped_NaCl$y= scale(all_BLUEs_genotyped_NaCl$y)


gwas.NGM_no_rdm = GridLMM_GWAS(formula = y ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                         data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')
NGM_results_no_rdm <- (gwas.NGM$gwas.NGM_no_rdm)
write.table(NGM_results,paste0("GridLMM_results_nordm/gwas_results_",traitID,"_NGM.txt"),quote=FALSE,row.names=FALSE)

gwas.NGM_rdm = GridLMM_GWAS(formula = y ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                        test_formula = ~1, reduced_formula = ~0,
                        data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')
NGM_results_rdm <- (gwas.NGM_rdm$results)
write.table(NGM_results_rdm,paste0("GridLMM_results_nordm/gwas_results_withRDM_",traitID,"_NGM.txt"),quote=FALSE,row.names=FALSE)

rm(list=ls())
gc()


ngm1=read.table("GridLMM_results_nordm/gwas_results_SF_NGM.txt",h=TRUE)
ngm2=read.table("GridLMM_results_nordm/gwas_results_withRDM_SF_NGM.txt",h=TRUE)

all_BLUEs_genotyped_NGM$y=sample(all_BLUEs_genotyped_NGM$y)
gwas.NGM_rdm.perm1 = GridLMM_GWAS(formula = y ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                            test_formula = ~1, reduced_formula = ~0,
                            data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')
min(gwas.NGM_rdm.perm1$results$p_value_REML)




dim(ngm1)
dim(ngm2)

hist(ngm1$p_value_REML)
hist(p.adjust(ngm1$p_value_REML,method="fdr"))
table(p.adjust(ngm1$p_value_REML,method="fdr")<0.05)
table(p.adjust(ngm2$p_value_REML,method="fdr")<0.05)


genotype_pruned_08 <- read.vcfR('plink_pruning/pruning08/CeMEE.pruned.0.8.vcf.gz')
rs_kept_08=(genotype_pruned_08@fix)[,3]
ngm2$pcor=p.adjust(ngm2$p_value_REML)

ngm2_08=subset(ngm2,X_ID%in%rs_kept_08)
plot(ngm2_08$pcor,p.adjust(ngm2_08$p_value_REML))

table(p.adjust(ngm2$p_value_REML,method="fdr")<0.05)
table(p.adjust(ngm2_08$p_value_REML,method="fdr")<0.05)


plot(density(ngm1$p_value_REML))
lines(density(ngm2$p_value_REML),col="red")

genotype_matrix_pruned=genotype_matrix[,colnames(genotype_matrix)%in%rs_kept_08]
is.NA=colSums(is.na(genotype_matrix_pruned))
head(genotype_matrix[,1:4])

gwas.NGM_rdm = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                            test_formula = ~1, reduced_formula = ~0,
                            data = all_BLUEs_genotyped_NGM, X = genotype_matrix_pruned[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')


### New GRM - not pruned
library("vcfR")
genotype_pruned <- read.vcfR('plink_pruning/pruning099//CeMEE.pruned.0.99.vcf.gz')
genotype_pruned_matrix <- apply(genotype_pruned@gt[,2:ncol(genotype_pruned@gt)],1,function(x){
  as.numeric(substring(x,1,1))
})
genotype_pruned_matrix=2*genotype_pruned_matrix
strainID <- tstrsplit(colnames(genotype_pruned@gt)[2:ncol(genotype_pruned@gt)],'_')[[1]]
GRM99 <- Gmatrix((genotype_pruned_matrix), method="VanRaden") 
dim(GRM99) # 1001 / 1001

GRM99_sub = GRM99[strainID %in% rownames(genotype_matrix_pruned),strainID %in% rownames(genotype_matrix_pruned)]
strainID_sub = strainID[strainID %in% rownames(genotype_matrix_pruned)]
table(strainID_sub==rownames(genotype_matrix_pruned)) # Same order

plot(GRM99_sub,GRM_sub,pch=16)

colnames(GRM99_sub)=strainID_sub
rownames(GRM99_sub)=strainID_sub

gwas.NGM_rdm_p99 = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                            test_formula = ~1, reduced_formula = ~0,
                            data = all_BLUEs_genotyped_NGM, X = genotype_matrix_pruned[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM99_sub), mc.cores = 2 ,method = 'REML')

table(gwas.NGM_rdm_p99$results$p_value_REML<0.05)
table(p.adjust(gwas.NGM_rdm_p99$results$p_value_REML,method="fdr")<0.05)


#### Now do only CHR 4 - exclude it from the GRM or not

genotype_pruned <- read.vcfR('plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz')
genotype_pruned_matrix <- apply(genotype_pruned@gt[,2:1002],1,function(x){
  as.numeric(substring(x,1,1))
})
genotype_pruned_matrix=2*genotype_pruned_matrix
strainID <- tstrsplit(colnames(genotype_pruned@gt)[2:1002],'_')[[1]]

GRM <- Gmatrix((genotype_pruned_matrix), method="VanRaden") # 2112 / 2112
## Re-order
GRM_sub = GRM[strainID %in% rownames(genotype_matrix_pruned),strainID %in% rownames(genotype_matrix_pruned)]
strainID_sub = strainID[strainID %in% rownames(genotype_matrix_pruned)]
table(strainID_sub==rownames(genotype_matrix)) # Same order

CHR_GRM=tstrsplit(genotype_pruned@fix[,3],"_")[[1]]
GRM_no4 <- Gmatrix((genotype_pruned_matrix[,CHR_GRM!="IV"]), method="VanRaden") # 2112 / 2112
## Re-order
GRM_no4_sub = GRM_no4[strainID %in% rownames(genotype_matrix_pruned),strainID %in% rownames(genotype_matrix_pruned)]
strainID_sub = strainID[strainID %in% rownames(genotype_matrix_pruned)]
table(strainID_sub==rownames(genotype_matrix)) # Same order

SNP_chr= tstrsplit(colnames(genotype_matrix_pruned),"_")[[1]]

colnames(GRM_no4_sub)=strainID_sub
rownames(GRM_no4_sub)=strainID_sub

colnames(GRM_sub)=strainID_sub
rownames(GRM_sub)=strainID_sub


gwas.NGM_rdm_with4 = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                            test_formula = ~1, reduced_formula = ~0,
                            data = all_BLUEs_genotyped_NGM, X = genotype_matrix_pruned[,is.NA==0 & SNP_chr=="IV"], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM_rdm_without4 = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                                  test_formula = ~1, reduced_formula = ~0,
                                  data = all_BLUEs_genotyped_NGM, X = genotype_matrix_pruned[,is.NA==0 & SNP_chr=="IV"], X_ID = 'pop_label', relmat = list(pop_label = GRM_no4_sub), mc.cores = 2 ,method = 'REML')


plot(-log10(gwas.NGM_rdm_with4$results$p_value_REML),-log10(gwas.NGM_rdm_without4$results$p_value_REML))
abline(a=0,b=1)

plot(gwas.NGM_rdm_without4$results$beta.1,-log10(gwas.NGM_rdm_without4$results$p_value_REML),pch=16)
plot(gwas.NGM_rdm_without4$results$beta.1,-log10(p.adjust(gwas.NGM_rdm_without4$results$p_value_REML,method="fdr")),pch=16)
abline(h=(-log10(0.05)))

pos=as.numeric(tstrsplit(gwas.NGM_rdm_without4$results$X_ID,"_")[[2]])
plot(pos,-log10(p.adjust(gwas.NGM_rdm_without4$results$p_value_REML,method="fdr")))


## Only A6/CAs ?

gwas.NGM_rdm_A6140 = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                                  test_formula = ~1, reduced_formula = ~0,
                                  data = subset(all_BLUEs_genotyped_NGM,substring(all_BLUEs_genotyped_NGM$pop_label,1,2)%in%c("A6","CA")), X = genotype_matrix_pruned[,is.NA==0 & SNP_chr=="IV"], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM_rdm_A6140_sessionlate = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                                  test_formula = ~1, reduced_formula = ~0,
                                  data = subset(all_BLUEs_genotyped_NGM,session_id%in%c("2018","F","SB") & substring(all_BLUEs_genotyped_NGM$pop_label,1,2)%in%c("A6","CA")), X = genotype_matrix_pruned[,is.NA==0 & SNP_chr=="IV"], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM_rdm_A6140_sessionlateR = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                                              test_formula = ~1, reduced_formula = ~0,
                                              data = subset(all_BLUEs_genotyped_NGM,session_id%in%c("F","SB") & substring(all_BLUEs_genotyped_NGM$pop_label,1,2)%in%c("A6","CA")), X = genotype_matrix_pruned[,is.NA==0 & SNP_chr=="IV"], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM_rdm_A6140_sessionEarly = GridLMM_GWAS(formula = SF ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|session_id)+(1|pop_label),
                                               test_formula = ~1, reduced_formula = ~0,
                                               data = subset(all_BLUEs_genotyped_NGM,!session_id%in%c("F") & substring(all_BLUEs_genotyped_NGM$pop_label,1,2)%in%c("A6","CA")), X = genotype_matrix_pruned[,is.NA==0 & SNP_chr=="IV"], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

plot(-log10(gwas.NGM_rdm_A6140_sessionEarly$results$p_value_REML))


save(list=ls(),file="GridLMM_results_mixed/Test.RData")