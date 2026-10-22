rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
library(parallel)
library(bigsnpr)
vect_P_traits=c("beta")

## Load Cemee data
#load('RData/Genotype_info_Cemee.RData')
#snp_positions=as.data.frame(snp_positions)
#snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
#snp_positions$POS=as.numeric(as.character(snp_positions$POS))

## Re-order
#genotype_info=genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
#snp_positions=snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

# Old genotypes
#genotypes_gemma2 <- read.table(file="data_BLUE/genotypes_for_gemma.txt")
#env_Grid=new.env()
#load('RData/Preparatory_GLMM.RData',envir=env_Grid)
#line_labels=rownames(env_Grid$genotype_matrix)
#pos_labels=colnames(env_Grid$genotype_matrix)
#genotypes_gemma=cbind(pos_labels,rep("A",length(pos_labels)),rep("a",length(pos_labels)),t(env_Grid$genotype_matrix))
#genotypes_gemma=data.frame(genotypes_gemma)
#genotypes_gemma$V2=genotypes_gemma2$V2
#genotypes_gemma$V3=genotypes_gemma2$V3
#head(genotypes_gemma[,1:5])
#write.table(genotypes_gemma,file="data_BLUE_clusters_CHR//genotypes_for_gemma.txt",sep='\t',row.names=FALSE,quote=FALSE)
#rm(env_Grid,genotypes_gemma2);gc()
genotypes_gemma <- read.table(file="data_BLUE_clusters_CHR/genotypes_for_gemma.txt",h=TRUE)

for(threshold_clust in c(0.2,0.4,0.6,0.8)){
print(threshold_clust)
## Load cluster data
all_cluster_merged=read.table(paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged=all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

## Subset
genotypes_gemma_sign=genotypes_gemma[genotypes_gemma$pos_labels%in%all_cluster_merged$rs,]
### We are going to create a new matrix for gemma mapping - with a line per cluster
## Order groups/clusters
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
all_cluster_merged$group = as.numeric(factor(as.character(all_cluster_merged$group),levels=unique(all_cluster_merged$group)))

table(genotypes_gemma_sign$pos_labels==all_cluster_merged$rs)
geno_matrix_groups = matrix(data=NA,nrow=max(all_cluster_merged$group),ncol=ncol(genotypes_gemma_sign)-3)

genotypes_gemma_sign <- data.frame(genotypes_gemma_sign)
# Sets NA for non 0/2
for(i in 4:ncol(genotypes_gemma_sign)){
  genotypes_gemma_sign[,i][!as.numeric(genotypes_gemma_sign[,i])%in%c(0,2)]=NA
}

### For the groups with multiple SNPs, given sign of the beta
group_to_change=tapply(all_cluster_merged$beta_T1,all_cluster_merged$group,function(x){length(unique(sign(x)))})
group_to_change=names(group_to_change)[as.numeric(group_to_change)==2]
rs_to_change=NULL
## Within these groups, we need to pick the ones that are different
for(i in group_to_change){
  sign_beta = sign(mean(subset(all_cluster_merged,group==i)$beta_T1))
  rs_to_change=c(rs_to_change,subset(all_cluster_merged,group==i & sign(beta_T1)!=sign_beta)$rs)
}
## Change polarization when needed
for(i in rs_to_change){
  idx.pol=genotypes_gemma_sign$pos_labels==i
  tmp.A=genotypes_gemma_sign$V2[idx.pol]
  genotypes_gemma_sign$V2[idx.pol]=genotypes_gemma_sign$V3[idx.pol]
  genotypes_gemma_sign$V3[idx.pol]=tmp.A
  genotypes_gemma_sign[idx.pol,4:ncol(genotypes_gemma_sign)]=2-genotypes_gemma_sign[idx.pol,4:ncol(genotypes_gemma_sign)]
}

## Save the cluster polarization
write.table(genotypes_gemma_sign[,1:3],paste0("data_clustPolarisation_NaCl/Beta_genotypes_clusters_pol_R",round(10*threshold_clust),".txt"),quote=FALSE,row.names=FALSE)

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

### BIMBAM FORMAT - OLD
#bimbam <- data.frame(id=df_clusters$clusterID,major="A",minor="a")
#bimbam <- cbind(bimbam,geno_matrix_groups)
#write.table(bimbam,paste0("data_BLUE_GridLMM_NaCl/genotypes_clusters_R",round(10*threshold_clust),".txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(df_clusters,paste0("data_BLUE_GridLMM_NaCl/Beta_data_clusters_BLUE_R",round(10*threshold_clust),".txt"),row.names=FALSE,quote=FALSE)

rownames(geno_matrix_groups)=df_clusters$clusterID
write.table(geno_matrix_groups,paste0("data_BLUE_GridLMM_NaCl/Beta_genotypes_matrix_clusters_R",round(10*threshold_clust),".txt"),quote=FALSE)
}
