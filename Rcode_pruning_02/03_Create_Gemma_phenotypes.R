rm(list=ls())
gc()
library(AGHmatrix)
library(data.table)
library("vcfR")
library(qqman)

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

write.table(bimbam_02,"data_BLUE/genotypes_for_gemma_R02.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(bimbam_08,"data_BLUE/genotypes_for_gemma_R08.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


### Run Gemma
###
# Use the latest installed gemma version - use BIMBAM format
#conda activate gemma-env 

### Only takes one column at a time
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
#do
#cut -f $i data_BLUE/gemma_pheno.txt > data_BLUE/split/gemma_pheno_T${i}.txt
#nohup gemma -g data_BLUE/genotypes_for_gemma_R02.txt -k data_BLUE/kinship.txt -lmm 4 -p data_BLUE/split/gemma_pheno_T${i}.txt -o gwas_results_BLUE_R02_T${i} > gemma_output_BLUE_txt_R02/output_gemma_R02_t${i}.txt &
#nohup gemma -g data_BLUE/genotypes_for_gemma_R08.txt -k data_BLUE/kinship.txt -lmm 4 -p data_BLUE/split/gemma_pheno_T${i}.txt -o gwas_results_BLUE_R08_T${i} > gemma_output_BLUE_txt_R08/output_gemma_R08_t${i}.txt &
#done


### Read result

gemma_T1_R02=read.table("gemma_R02/gwas_results_BLUE_R02_T1.assoc.txt",h=TRUE)
gemma_T1_R02$chr = factor(tstrsplit(gemma_T1_R02$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
gemma_T1_R02$ps = as.numeric(tstrsplit(gemma_T1_R02$rs,"_")[[2]])
gemma_T1_R02=gemma_T1_R02[order(gemma_T1_R02$chr,gemma_T1_R02$ps),]

gemma_T1_R02$p_wald.cor=p.adjust(gemma_T1_R02$p_wald)
data_manhattan = data.frame(CHR=as.numeric(gemma_T1_R02$chr),BP=gemma_T1_R02$ps,SNP=c(1:nrow(gemma_T1_R02)),
                            P=-log10(gemma_T1_R02$p_wald.cor))

qqnorm(data_manhattan$P)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          ,cex=.6,suggestiveline=FALSE,genomewideline=FALSE)

###########

gemma_T1_R08=read.table("gemma_R08/gwas_results_BLUE_R08_T1.assoc.txt",h=TRUE)
gemma_T1_R08$chr = factor(tstrsplit(gemma_T1_R08$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
gemma_T1_R08$ps = as.numeric(tstrsplit(gemma_T1_R08$rs,"_")[[2]])
gemma_T1_R08=gemma_T1_R08[order(gemma_T1_R08$chr,gemma_T1_R08$ps),]

gemma_T1_R08$p_wald.cor=p.adjust(gemma_T1_R08$p_wald)
data_manhattan = data.frame(CHR=as.numeric(gemma_T1_R08$chr),BP=gemma_T1_R08$ps,SNP=c(1:nrow(gemma_T1_R08)),
                            P=-log10(gemma_T1_R08$p_wald))

qqnorm(data_manhattan$P)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          ,cex=.6,suggestiveline=FALSE,genomewideline=FALSE)



