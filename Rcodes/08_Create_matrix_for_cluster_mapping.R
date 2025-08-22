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
all_cluster_merged = read.table(file=paste0('txt/all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
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

### We will load the former bed files from plink, update the genotype matrix and
### save a new one

bed_original <- snp_readBed("pink_binary/plink_binary.bed")
Cemee_bed_file <- snp_attach(bed_original)

n_ind <- 1001
n_snp <- nrow(geno_matrix_groups)

## 3. Create an FBM.code256 object on disk
## This will produce two files: geno.bk (data) and geno.rds (descriptor)
G <- FBM.code256(
  nrow = n_ind, 
  ncol = n_snp,
  code = bigsnpr::CODE_012,
  backingfile = "geno_matrix_groups"
)

## 4. Fill the FBM with your genotype codes
G[,] <- t(geno_matrix_groups)

## 5. Create corresponding .fam and .bim data frames
fam <- data.frame(
  family.ID = Cemee_bed_file$fam$family.ID,
  sample.ID = Cemee_bed_file$fam$sample.ID,
  paternal.ID=0,
  maternal.ID=0,
  sex = 0,
  affection = -9
)

bim <- data.frame(
  chromosome = df_clusters$CHROM,
  marker.ID = df_clusters$clusterID,
  genetic.dist = 0,
  physical.pos = df_clusters$POS,
  allele1 = "A",
  allele2 = "G"
)


# 6. Combine into a bigSNP object
obj.bigSNP <- list(genotypes = G, fam = fam, map = bim)
class(obj.bigSNP) <- "bigSNP"

# 7. Save as .rds (descriptor file)
saveRDS(obj.bigSNP, "geno.rds")
snp_writeBed(obj.bigSNP, bedfile = "plink_binary_clusters.bed")


### Run Gemma
###
# Use the latest installed gemma version
#conda activate gemma-env 

### Only takes one column at a time
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
#do
#cut -f $i phenotypic_data/gemma_pheno.txt > phenotypic_data/gemma_pheno_T${i}.txt
#nohup gemma -bfile pink_binary_clusters/plink_binary_clusters -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p phenotypic_data/gemma_pheno_T${i}.txt -o gwas_results_clusters_T${i} > gemma_output_txt/output_gemma_clusters_t${i}.txt &
#done

