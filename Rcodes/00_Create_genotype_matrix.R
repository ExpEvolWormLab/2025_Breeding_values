rm(list=ls())
gc()

library(vcfR)
Cemee_vcf <- read.vcfR('/Volumes/nematode_genetics/Cemee_data/vcf_2025/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz')

snp_positions <- Cemee_vcf@fix[,1:5]
genotype_info <- Cemee_vcf@gt
rm(Cemee_vcf);gc()

## remove format column. # colnames(genotype_info)[1]
genotype_info <- genotype_info[,!colnames(genotype_info)=="FORMAT"]
genotype_info=substring(genotype_info,1,3)
genotype_info=apply(genotype_info,2,function(x){
  x[!x%in%c("1|1","0|0")]=NA
  x[x=="0|0"]=0
  x[x=="1|1"]=2
  as.numeric(x)})
gc()

snp_positions=data.frame(snp_positions)
snp_positions$POS=as.numeric(snp_positions$POS)
save(list=ls(),file='RData/Genotype_info_Cemee.RData')
