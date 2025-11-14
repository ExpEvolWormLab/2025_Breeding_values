rm(list=ls())
gc()
library(qqman)
library(data.table)
for(j in 1:14){

traitJ = read.table(paste0("gemma_results_clusters_BLUE_R02/gwas_results_clusters_BLUE_R02_T",j,".assoc.txt"),h=TRUE)
traitJ$group = tstrsplit(traitJ$rs,"_")[[1]]
traitJ$chr = tstrsplit(traitJ$rs,"_")[[2]]
traitJ$ps = tstrsplit(traitJ$rs,"_")[[3]]

### Should be removed
df_clusters=read.table("data_BLUE/data_clusters_BLUE_R02.txt",h=TRUE)
df_clusters$rs=df_clusters$clusterID
dim(traitJ)
traitJ=merge(traitJ,df_clusters)
dim(traitJ)
table(traitJ$chr==traitJ$CHROM)
table(traitJ$ps==traitJ$POS)

# Load permutation results
perm.temp <- new.env()
load(paste0("Permutation_Gemma_clusters_BLUE_R02/Permutation_gemma_trait_clusters_BLUE_R02_",j,".RData"),envir=perm.temp)
List_output=perm.temp$List_output
## Some errors - we will keep only complete lists
list_length <- unlist(lapply(List_output,nrow))
for(i in 1:length(List_output)){
  names(List_output[[i]])[3]=paste("X",i)
}
List_output <- lapply(List_output, as.data.table)
concat.perm <- Reduce(function(x, y) merge(x, y,all=TRUE), List_output)
dim(concat.perm)
rm(List_output);rm(perm.temp);gc()
nb_NA <- rowSums(is.na(concat.perm))
## Remove those that have Nas>=1%
concat.perm=concat.perm[nb_NA<floor(ncol(concat.perm)*0.01),]
concat.perm.sorted=apply(concat.perm[,3:ncol(concat.perm)],1,sort)
concat.perm.sorted=data.frame(t(concat.perm.sorted))
sign_thres <- concat.perm.sorted[,floor(ncol(concat.perm.sorted)*0.05)]
concat.perm.sorted=cbind(data.frame(concat.perm[,1:2]),concat.perm.sorted)
concat.perm.sorted$sign_thres=sign_thres
rm(concat.perm,sign_thres);gc()

## We will use the top 5% as sign threshold - keep this even if less col
traitJ = merge(traitJ,concat.perm.sorted[,c("rs","sign_thres")])
traitJ$is_sign <- traitJ$p_wald<traitJ$sign_thres
table(traitJ$is_sign)

concat_GWAS.temp=traitJ[,c("rs","chr","ps","n_miss","allele1","allele0","af","beta","se","p_wald","is_sign","clusterID")]
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
rm(concat_GWAS.temp,concat.perm.sorted,j,p,traitJ,nb_NA,i,list_length)

all_GWAS_clusters = all_GWAS
rm(all_GWAS)
head(all_GWAS_clusters)


save(list=ls(),file="RData_BLUEs/GWAS_cluster_results_BLUE_concatenated_R02.RData")

