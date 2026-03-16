rm(list=ls())
gc()
library(qqman)
library(data.table)
library(matrixStats)
for(j in 1:14){
traitJ = read.table(paste0("gemma_results_BLUE/gwas_results_BLUE_T",j,".assoc.txt"),h=TRUE)
# Load permutation results
perm.temp <- new.env()
load(paste0("Permutation_Gemma_BLUE_with_residuals/Permutation_gemma_trait_BLUE_",j,".RData"),envir=perm.temp)
List_output=perm.temp$List_output
rm(perm.temp);gc()
## Some errors - we will keep only complete lists
list_length <- unlist(lapply(List_output,nrow))
## Concatenate all the complete lists
List_output=List_output[list_length==max(list_length)]
print("Kept permutations:")
print(length(List_output))

# Extract the p values threshold - genome-wide
pwald_list_5p = min(unlist(lapply(List_output, function(x){sort(x[,3])[15880]})))


traitJ$is_sign <- traitJ$p_wald<pwald_list_5p
#table(traitJ$is_sign)
print(sum(traitJ$is_sign))

traitJ$chr = factor(tstrsplit(traitJ$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
traitJ$ps = as.numeric(tstrsplit(traitJ$rs,"_")[[2]])
traitJ=traitJ[order(traitJ$chr,traitJ$ps),]
dim(traitJ)

data_manhattan = data.frame(CHR=as.numeric(traitJ$chr),BP=traitJ$ps,SNP=c(1:nrow(traitJ)),
                            P=-log10(traitJ$p_wald),is_sign=traitJ$is_sign)

data_manhattan=subset(data_manhattan,CHR!=7)

pdf(paste0("plots_BLUEs/Manhattan_GWAS_Trait_Global_",j,".pdf"),w=10)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          highlight=data_manhattan$SNP[data_manhattan$is_sign],cex=.6,suggestiveline=FALSE,genomewideline=FALSE)
dev.off()

concat_GWAS.temp=traitJ[,c(1:9,13,16)]
for(p in c(4,7,8:11)){
names(concat_GWAS.temp)[p]=paste0(names(concat_GWAS.temp)[p],"_T",j)
}

if(!exists("all_GWAS")){
  all_GWAS=concat_GWAS.temp
}else{
all_GWAS=merge(all_GWAS,concat_GWAS.temp,all=TRUE)
}
#print(table(traitJ$rs%in%all_GWAS$rs))
print(dim(all_GWAS))

}
rm(concat_GWAS.temp,data_manhattan,j,p,traitJ,nb_NA,i,list_length)

# Simply replace NA values by FALSE
for(i in 1:14){
  if(sum(is.na(all_GWAS[,names(all_GWAS)==paste0("is_sign_T",i)]))>0){
all_GWAS[is.na(all_GWAS[,names(all_GWAS)==paste0("is_sign_T",i)]),paste0("is_sign_T",i)]=FALSE
  }
}
rm(i)
save(list=ls(),file="RData/GWAS_results_concatenated_Global.RData")

