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
#for(i in 1:length(List_output)){
#  names(List_output[[i]])[3]=paste("X",i)
#}

## Concatenate all the complete lists
List_output=List_output[list_length==max(list_length)]
print("Kept permutations:")
print(length(List_output))

# Extract the p values
pwald_list = lapply(List_output, function(x){x[,3]})
pwald_list=do.call(cbind,pwald_list)
concat.perm.sorted=List_output[[1]][,1:2]
rm(List_output);gc()

nb_NA <- rowSums(is.na(pwald_list))
## Remove those that have Nas>=5
pwald_list=pwald_list[nb_NA<5,]

pwald_list_sorted <- t(apply(pwald_list, 1, sort,na.last = TRUE))
pwald_list_sorted_columns <- t(apply(pwald_list, 2, sort,na.last = TRUE))


min_per_permutations <- pwald_list_sorted_columns[,1]
min_per_permutations=sort(min_per_permutations)
plot(min_per_permutations,type='l')

table(traitJ$p_wald<min_per_permutations[50]) # Q=0.05
table(traitJ$p_wald<min_per_permutations[100]) # Q=0.1
table(traitJ$p_wald<min_per_permutations[200]) # Q=0.2
table(traitJ$p_wald<min_per_permutations[200]) # Q=0.3

table(p.adjust(traitJ$p_wald,method="BH")<0.1)
table(p.adjust(traitJ$p_wald,method="BH")<0.2)



sign_thres <- pwald_list_sorted[,floor(ncol(pwald_list_sorted)*0.05)]

concat.perm.sorted$sign_thres=sign_thres

## We will use the top 5% as sign threshold - keep this even if less col
traitJ = merge(traitJ,concat.perm.sorted[,c("rs","sign_thres")])
traitJ$is_sign <- traitJ$p_wald<traitJ$sign_thres

#rm(sign_thres,pwald_list,pwald_list_sorted)
#rm(concat.perm.sorted);gc()
rm(List_output);gc()

traitJ$chr = factor(tstrsplit(traitJ$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
traitJ$ps = as.numeric(tstrsplit(traitJ$rs,"_")[[2]])
traitJ=traitJ[order(traitJ$chr,traitJ$ps),]
dim(traitJ)

data_manhattan = data.frame(CHR=as.numeric(traitJ$chr),BP=traitJ$ps,SNP=c(1:nrow(traitJ)),
                            P=-log10(traitJ$p_wald),is_sign=traitJ$is_sign)
data_manhattan=subset(data_manhattan,CHR!=7)

pdf(paste0("plots_BLUEs/Manhattan_GWAS_Trait_Q_",j,".pdf"),w=10)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          highlight=data_manhattan$SNP[data_manhattan$is_sign],cex=.6,suggestiveline=FALSE,genomewideline=FALSE)
dev.off()

concat_GWAS.temp=traitJ[,c(1:9,13,17)]
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
save(list=ls(),file="RData/GWAS_results_concatenated.RData")

