rm(list=ls())
gc()
library(qqman)
library(data.table)

for(j in 1:14){
traitJ = read.table(paste0("gemma_results_BLUE/gwas_results_BLUE_T",j,".assoc.txt"),h=TRUE)
# Load permutation results
perm.temp <- new.env()
load(paste0("Permutation_Gemma_BLUE/Permutation_gemma_trait_BLUE_",j,".RData"),envir=perm.temp)
List_output=perm.temp$List_output
## Some errors - we will keep only complete lists
list_length <- unlist(lapply(List_output,nrow))
for(i in 1:length(List_output)){
  names(List_output[[i]])[3]=paste("X",i)
}
List_output <- lapply(List_output, as.data.table)
concat.perm <- Reduce(function(x, y) merge(x, y), List_output)
rm(List_output);rm(perm.temp);gc()

nb_NA <- rowSums(is.na(concat.perm))
## Remove those that have Nas>=5
concat.perm=concat.perm[nb_NA<5,]
mat <- as.matrix(concat.perm[,3:ncol(concat.perm)])
mat_sorted <- t(apply(mat, 1, sort,na.last = TRUE))
dim(mat_sorted)
#concat.perm.sorted=data.frame((concat.perm.sorted))
sign_thres <- mat_sorted[,floor(ncol(mat_sorted)*0.05)]
concat.perm.sorted=cbind(data.frame(concat.perm[,1:2]),mat_sorted)
concat.perm.sorted$sign_thres=sign_thres
rm(concat.perm,sign_thres,mat,mat_sorted);gc()

## We will use the top 5% as sign threshold - keep this even if less col
traitJ = merge(traitJ,concat.perm.sorted[,c("chr","rs","sign_thres")])
traitJ$is_sign <- traitJ$p_wald<traitJ$sign_thres

#plot((1:length(traitJ$beta))[!traitJ$is_sign],traitJ$beta[!traitJ$is_sign],col=adjustcolor("grey",alpha=.6),pch=16,cex=.6,
#     ylim=c(min(traitJ$beta),max(traitJ$beta)))
#points((1:length(traitJ$beta))[traitJ$is_sign],traitJ$beta[traitJ$is_sign],col=adjustcolor("red",alpha=.6),pch=16,cex=.6)

data_manhattan = data.frame(CHR=as.numeric(as.factor(traitJ$chr)),BP=as.numeric(traitJ$ps),SNP=c(1:nrow(traitJ)),
                            P=-log10(traitJ$p_wald),is_sign=traitJ$is_sign)
data_manhattan=subset(data_manhattan,CHR!=6)
data_manhattan$CHR=as.numeric(factor(data_manhattan$CHR,levels = c(2:5,7,1)))
pdf(paste0("plots_BLUEs/Manhattan_GWAS_Trait_",j,".pdf"),w=10)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          highlight=data_manhattan$SNP[data_manhattan$is_sign],cex=.6,suggestiveline=FALSE,genomewideline=FALSE)
dev.off()

concat_GWAS.temp=traitJ[,c(1:9,13,17)]
for(p in 8:11){
names(concat_GWAS.temp)[p]=paste0(names(concat_GWAS.temp)[p],"_T",j)
}

if(j<=7){
  names(concat_GWAS.temp)[4]="n_miss.NGM"
  names(concat_GWAS.temp)[7]="af.NGM"
}else{
  names(concat_GWAS.temp)[4]="n_miss.NaCl"
  names(concat_GWAS.temp)[7]="af.NaCl"
}

if(!exists("all_GWAS")){
  all_GWAS=concat_GWAS.temp
}else{
all_GWAS=merge(all_GWAS,concat_GWAS.temp,all=TRUE)
}
#print(table(traitJ$rs%in%all_GWAS$rs))
print(dim(all_GWAS))

}
rm(concat_GWAS.temp,concat.perm.sorted,data_manhattan,j,p,traitJ,nb_NA,i,list_length)

# Simply replace NA values by FALSE
for(i in 1:14){
  if(sum(is.na(all_GWAS[,names(all_GWAS)==paste0("is_sign_T",i)]))>0){
all_GWAS[is.na(all_GWAS[,names(all_GWAS)==paste0("is_sign_T",i)]),paste0("is_sign_T",i)]=FALSE
  }
}
rm(i)
save(list=ls(),file="RData_BLUEs//GWAS_results_concatenated.RData")

