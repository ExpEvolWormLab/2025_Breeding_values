rm(list=ls())
gc()
library(AGHmatrix)
library(data.table)
library("vcfR")

bimbam=read.table("data_BLUE/genotypes_for_gemma.txt")
head(bimbam[,1:5])

genotype_pruned_02 <- read.vcfR('plink_pruning/pruning02/CeMEE.pruned.0.2.vcf.gz')
genotype_pruned_08 <- read.vcfR('plink_pruning/pruning08/CeMEE.pruned.0.8.vcf.gz')
kept_pos = genotype_pruned_02@fix[,3]
kept_pos8 = genotype_pruned_08@fix[,3]


table(bimbam$V1%in%kept_pos) # 3502 only out of 314k
bimbam_02 = bimbam[bimbam$V1%in%kept_pos,]
bimbam_08 = bimbam[bimbam$V1%in%kept_pos8,]
dim(bimbam_02) # 3.5k
dim(bimbam_08) # 22k


write.table(GRM_sub,file="data_BLUE/kinship.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

### Run Gemma
###
# Use the latest installed gemma version - use BIMBAM format
#conda activate gemma-env 

### Only takes one column at a time
#for i in 2 3 4 5 6 7 8 9 10 11 12 13 14
#do
#cut -f $i data_BLUE/gemma_pheno.txt > data_BLUE/split/gemma_pheno_T${i}.txt
#nohup gemma -g data_BLUE/genotypes_for_gemma.txt -k data_BLUE/kinship.txt -lmm 4 -p data_BLUE/split/gemma_pheno_T${i}.txt -o gwas_results_BLUE_T${i} > gemma_output_BLUE_txt/output_gemma_t${i}.txt &
#done

