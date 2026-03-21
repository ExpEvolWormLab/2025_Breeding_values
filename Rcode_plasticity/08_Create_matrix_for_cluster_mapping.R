rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
library(parallel)
library(bigsnpr)
vect_P_traits=c("beta")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

## Load Cemee data
#load('RData/Genotype_info_Cemee.RData')
#snp_positions=as.data.frame(snp_positions)
#snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
#snp_positions$POS=as.numeric(as.character(snp_positions$POS))

## Re-order
#genotype_info=genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
#snp_positions=snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

genotypes_gemma <- read.table(file="data_BLUE_plasticity//genotypes_for_gemma.txt")

for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,1)){
print(threshold_clust)
## Load cluster data
all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged=all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]


## Subset
genotypes_gemma_sign=genotypes_gemma[genotypes_gemma$V1%in%all_cluster_merged$rs,]
### We are going to create a new matrix for gemma mapping - with a line per cluster
## Order groups/clusters
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
all_cluster_merged$group = as.numeric(factor(as.character(all_cluster_merged$group),levels=unique(all_cluster_merged$group)))

table(genotypes_gemma_sign$V1==all_cluster_merged$rs)
geno_matrix_groups = matrix(data=NA,nrow=max(all_cluster_merged$group),ncol=ncol(genotypes_gemma_sign)-3)

#for(i in 1:max(all_cluster_merged$group)){
#  if(i==1000*round(i/1000)){
#    print(i)
#  }
#  if(sum(all_cluster_merged$group==i)==1){
#    geno_matrix_groups[i,]=as.numeric(genotypes_gemma_sign[all_cluster_merged$group==i,4:ncol(geno_group)])
#  }else{
#    geno_group <- genotypes_gemma_sign[all_cluster_merged$group==i,]
#    ## Instead, just take the mean
#    geno_matrix_groups[i,] = colMeans(genotypes_gemma_sign[,4:ncol(geno_group)],na.rm=TRUE)
#  }
#}

genotypes_gemma_sign <- as.data.table(genotypes_gemma_sign)
genotypes_gemma_sign[, group := all_cluster_merged$group]
res <- genotypes_gemma_sign[, lapply(.SD, mean, na.rm = TRUE), by = group, .SDcols = 4:(ncol(genotypes_gemma_sign)-1)]
# Convert to matrix (removing "group" column)
geno_matrix_groups <- as.matrix(res[order(group), -"group"])
rm(res);gc()
## Unknown
geno_matrix_groups[is.na(geno_matrix_groups)] <- 3
#table(geno_matrix_groups==3) # Rare
table(rowSums(geno_matrix_groups==3)>50) # No one has more than 50 NAs ?

### Create CHR/POS info for the clusters
df_clusters = data.frame(group=names(tapply((all_cluster_merged$POS),all_cluster_merged$group,mean)),
                         CHROM= tapply(as.character(all_cluster_merged$CHROM),all_cluster_merged$group,unique),
                         POS= floor(tapply((all_cluster_merged$POS),all_cluster_merged$group,mean)))
## Ensure unique clusterID by adding a gp info
df_clusters$clusterID = paste(df_clusters$group,df_clusters$CHROM,df_clusters$POS,sep='_')
head(df_clusters)

### BIMBAM FORMAT
bimbam <- data.frame(id=df_clusters$clusterID,major="A",minor="a")
bimbam <- cbind(bimbam,geno_matrix_groups)

write.table(bimbam,paste0("data_BLUE_clusters_plasticity/genotypes_clusters_R",round(10*threshold_clust),".txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(df_clusters,paste0("data_BLUE_clusters_plasticity/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),row.names=FALSE,quote=FALSE)

}

### Run Gemma
###
# Use the latest installed gemma version
#conda activate gemma-env 

### Only takes one column at a time
#for j in 2 4 6 8 9 10
#do
#for i in 1 2
#do
#cut -f $i data_BLUE_plasticity/gemma_pheno.txt > data_BLUE_plasticity/split/gemma_pheno_T${i}.txt
#nohup gemma -g data_BLUE_clusters_plasticity/genotypes_clusters_R${j}.txt -k data_BLUE_plasticity/kinship.txt -lmm 4 -p data_BLUE_plasticity/split/gemma_pheno_T${i}.txt -o gwas_results_clusters_BLUE_plasticity_R${j}_T${i} > gemma_output_BLUE_plasticity_clusters_txt/output_gemma_clusters_R${j}_t${i}.txt &
#done
#wait
#done
