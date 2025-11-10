rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
library(parallel)
library(bigsnpr)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

## Load cluster data
all_cluster_merged = read.table(file=paste0('txt_BLUEs/all_sign_GWAS_with_clusters_R02.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged=all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
## Load Cemee data
load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))

## Re-order
genotype_info=genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
snp_positions=snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

## Subset
genotype_info=genotype_info[snp_positions$ID%in%all_cluster_merged$rs,]
snp_positions=snp_positions[snp_positions$ID%in%all_cluster_merged$rs,]
dim(snp_positions)
### We are going to create a new matrix for gemma mapping - with a line per cluster
## Order groups/clusters
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
all_cluster_merged$group = as.numeric(factor(as.character(all_cluster_merged$group),levels=unique(all_cluster_merged$group)))

table(snp_positions$ID==all_cluster_merged$rs)
table(snp_positions$REF==all_cluster_merged$allele0)
geno_matrix_groups = matrix(data=NA,nrow=max(all_cluster_merged$group),ncol=1001)

for(i in 1:max(all_cluster_merged$group)){
  if(i==1000*round(i/1000)){
    print(i)
  }
  if(sum(all_cluster_merged$group==i)==1){
    geno_matrix_groups[i,]=as.numeric(genotype_info[all_cluster_merged$group==i,])
  }else{
  geno_group <- genotype_info[all_cluster_merged$group==i,]
  #geno_clusters <- rep(3,1001)
  #n_snp = sum(all_cluster_merged$group==i)
  #n_alt = colSums(geno_group)/2
  #geno_clusters[n_alt<(n_snp/5)]=0
  #geno_clusters[n_alt>(n_snp-(n_snp/5))]=2
  #geno_matrix_groups[i,]=geno_clusters
  
  ## Instead, just take the mean
  geno_matrix_groups[i,] = colMeans(geno_group,na.rm=TRUE)
}
}

geno_matrix_groups[is.na(geno_matrix_groups)] <- 3
table(geno_matrix_groups==3)
table(rowSums(geno_matrix_groups==3)>50)

### Create CHR/POS info for the clusters
df_clusters = data.frame(clusterID=1:max(all_cluster_merged$group),
  CHROM= tapply(as.character(all_cluster_merged$CHROM),all_cluster_merged$group,unique),
  POS= tapply((all_cluster_merged$POS),all_cluster_merged$group,mean))

### Change this - BIMBAM FORMAT

head(df_clusters)

bimbam <- data.frame(id=paste0("rs",1:nrow(geno_matrix_groups)),major="A",minor="a")
bimbam <- cbind(bimbam,geno_matrix_groups)
head(bimbam[,1:10])

df_clusters$id=paste0("rs",1:nrow(geno_matrix_groups))

write.table(bimbam,"data_BLUE/genotypes_clusters_R02.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(df_clusters,"data_BLUE/data_clusters_BLUE_R02.txt",row.names=FALSE,quote=FALSE)

### Run Gemma
###
# Use the latest installed gemma version
#conda activate gemma-env 

### Only takes one column at a time
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
#do
#cut -f $i data_BLUE/gemma_pheno.txt > phenotypic_data/gemma_pheno_T${i}.txt
#nohup gemma -g data_BLUE/genotypes_clusters_R02.txt -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p data_BLUE/gemma_pheno_T${i}.txt -o gwas_results_clusters_BLUE_R02_T${i} > gemma_output_BLUE_clusters_txt/output_gemma_clusters_R02_t${i}.txt &
#done

