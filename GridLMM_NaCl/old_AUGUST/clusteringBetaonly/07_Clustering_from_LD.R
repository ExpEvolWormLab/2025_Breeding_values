rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
library(parallel)
library(qqman)
vect_P_traits=c("beta")

trait_ID=0

env_label='NaCl'
### Load all results and concatenate them
  for(traitJ in vect_P_traits){
    print(traitJ)
    trait_ID=trait_ID+1
    # concatenate all the CHR
    all_GWAS.temp=NULL
    for(i in c("I","II","III","IV","V","X")){
      all_GWAS.temp=rbind(all_GWAS.temp,read.table(paste0("GridLMM_results_NaCl/gwas_results_",traitJ,"_CHR_",i,"_",env_label,".txt"),h=TRUE)[,c("X_ID","beta.5","p_value_REML")])
      
    }
    
    names(all_GWAS.temp)=c("rs",paste0("beta_T",trait_ID),paste0("pvalue_T",trait_ID))

    all_GWAS.temp$chr = factor(tstrsplit(all_GWAS.temp$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
    all_GWAS.temp$ps = as.numeric(tstrsplit(all_GWAS.temp$rs,"_")[[2]])
    all_GWAS.temp=all_GWAS.temp[order(all_GWAS.temp$chr,all_GWAS.temp$ps),]
    
    data_manhattan = data.frame(CHR=as.numeric(all_GWAS.temp$chr),BP=all_GWAS.temp$ps,SNP=c(1:nrow(all_GWAS.temp)),
                                P=-log10(p.adjust(all_GWAS.temp[,paste0("pvalue_T",trait_ID)],method="BH")))
    
    data_manhattan=subset(data_manhattan,CHR!=7)
    
    if(trait_ID==1){
      all_GWAS=all_GWAS.temp
    }else{
      all_GWAS=merge(all_GWAS,all_GWAS.temp)
    }
    
  }

dim(all_GWAS)
save(list=c("all_GWAS"),file='GridLMM_results_NaCl/Beta_Gwas.RData')

all_GWAS$is_sign_T1=p.adjust(all_GWAS$pvalue_T1,method="BH")<0.2
table(all_GWAS$is_sign_T1) # 8k SNPs

###############################################
#### Create the clusters and save - Merged ####
###############################################

## Re-order
all_GWAS=all_GWAS[order(all_GWAS$chr,all_GWAS$ps),]
all_GWAS$CHROM=all_GWAS$chr

#####################################
## We need the matrix of genotypes ##
#####################################
load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))
# Remove MtDNA
snp_positions=subset(snp_positions,ID%in%all_GWAS$rs)
snp.id.sign <- all_GWAS$rs[all_GWAS$is_sign_T1]

## Re-order
genotype_info=genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
snp_positions=snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

table(snp.id.sign%in%snp_positions$ID) # all present
g_sign=genotype_info[snp_positions$ID%in%snp.id.sign,]
snp_positions_sign=subset(snp_positions,ID%in%snp.id.sign)

snps_pos <- data.frame(CHR=snp_positions_sign$CHROM,pos=snp_positions_sign$POS)
snps_pos$CHR=as.factor(snps_pos$CHR)
snps_pos$pos=as.numeric(snps_pos$pos)

## Finally, we need to polarize the GWAS effects according to the genotype info
table(all_GWAS$rs==snp_positions$ID)

all_GWAS_sign <- subset(all_GWAS,is_sign_T1)
table(snp_positions_sign$ID==all_GWAS_sign$rs) # In the same order #
dim(snp_positions_sign)

### Save the txt file with the significant positions

write.table(snp_positions_sign,file="data_BLUE_GridLMM_NaCl/sign_positions_Beta.txt",quote=FALSE,col.names = FALSE,row.names=FALSE,sep='\t')
system("bcftools view -T data_BLUE_GridLMM_NaCl/sign_positions_Beta.txt -O z /Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz >  data_BLUE_GridLMM_NaCl/Cemee2025.sign.GWAS_BLUE_BETA.filt.vcf.gz")
system("vcftools --gzvcf data_BLUE_GridLMM_NaCl/Cemee2025.sign.GWAS_BLUE_BETA.filt.vcf.gz --hap-r2 --ld-window 400 --out data_BLUE_GridLMM_NaCl/Beta_ld_window_100")

### Now load the linkage blocks from plink

ld_snps <- read.table("data_BLUE_GridLMM_NaCl/Beta_ld_window_100.hap.ld",h=TRUE)
dim(ld_snps)

### Create clusters from these
for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,0.99)){


########################################################################
########## RUN ALL CODE IN A FUNCTION WITH MULTIPLE THRESHOLD ##########
########################################################################

# Create an undirected graph
ld_snps$ID1 <- paste(ld_snps$CHR,ld_snps$POS1,sep="_")
ld_snps$ID2 <- paste(ld_snps$CHR,ld_snps$POS2,sep="_")
edges_clust <- subset(ld_snps, R.2 > threshold_clust)

# Function to get groups per chromosome
get_groups_per_chr <- function(data_chr) {
  g <- graph_from_data_frame(data_chr[, c("ID1", "ID2")], directed = FALSE)
  comp <- components(g)
  data.frame(
    ID1 = names(comp$membership),
    cluster = comp$membership,
    stringsAsFactors = FALSE
  )
}

# Apply per chromosome and combine results
clusters_by_chr <- edges_clust %>%
  group_by(CHR) %>% get_groups_per_chr
row.names(clusters_by_chr)=NULL
clusters_by_chr$CHROM=tstrsplit(clusters_by_chr$ID1,"_")[[1]]
clusters_by_chr$POS=tstrsplit(clusters_by_chr$ID1,"_")[[2]]

### Add the cluster info to the GWAS results
all_GWAS_sign$POS=all_GWAS_sign$ps
dim(all_GWAS_sign)
all_GWAS_sign_with_clust <- (merge(all_GWAS_sign,clusters_by_chr,all.x=TRUE))
dim(all_GWAS_sign_with_clust)
length(unique(all_GWAS_sign_with_clust$cluster)) 
sum(is.na(all_GWAS_sign_with_clust$cluster))

clust_length <- tapply(all_GWAS_sign_with_clust$ps,all_GWAS_sign_with_clust$clust,function(x){
  return(max(x)-min(x))
})
max(clust_length)/10^6 # 1.76 Mb
hist(log10(clust_length),main="",xlab=c("Log10 cluster length"))

clust_mean <- tapply(all_GWAS_sign_with_clust$ps,all_GWAS_sign_with_clust$clust,function(x){
  mean(x)
})

max_clust_nb=max(all_GWAS_sign_with_clust$cluster,na.rm=TRUE)
all_cluster_merged=all_GWAS_sign_with_clust
all_cluster_merged$group=all_cluster_merged$cluster

plot(density(as.numeric(table(all_cluster_merged$group))),log="x")
summary(as.numeric(table(all_cluster_merged$group)))

### Add the unique SNPs
all_cluster_merged[is.na(all_cluster_merged$cluster),]$group=all_cluster_merged[is.na(all_cluster_merged$cluster),]$rs
table(is.na(all_cluster_merged$group))
## Number of clusters initially ?  # 639 for R=0.2
length(unique(all_GWAS_sign_with_clust$cluster)) + sum(is.na(all_GWAS_sign_with_clust$cluster))
## Number of clusters in the end ?
length(unique(all_cluster_merged$group)) # 733

all_cluster_merged=subset(all_cluster_merged,select= -ID1)
### Replace group ID by a number, first sort all_cluster_merged
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
all_cluster_merged$group=as.numeric(factor(as.character(all_cluster_merged$group),levels=unique(all_cluster_merged$group)))
#barplot(log(table(table(all_cluster_merged$group))))

for(idx_chr in sort(unique(all_cluster_merged$CHROM))){

  print(idx_chr)
  # No singletons
  df.CHR=subset(all_cluster_merged,CHROM==idx_chr & !is.na(cluster))
  if(nrow(df.CHR)>0){
  df.CHR$group=as.numeric(as.factor(as.character(df.CHR$group)))
  # Singletons
  df.CHR.single=subset(all_cluster_merged,CHROM==idx_chr & is.na(cluster))
  df.CHR.single$group=as.numeric(as.factor(as.character(df.CHR.single$group)))
  
  pdf(paste0("plots_BLUEs_GridLMM_NaCl/Beta_cluster_CHR_R",round(10*threshold_clust),"_",idx_chr,"_Merged-noPruning.pdf"),h=10,w=20)
  
  all_col_CHR <- sample(rainbow(max(df.CHR$cluster)))
  all_pch_CHR <- sample(1:25,max(df.CHR$cluster),replace=TRUE)

  kept_chr = all_GWAS$CHROM==idx_chr
  i_env="NaCl"
  i=1
  i_list_env=0
  
      plot(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="n",ylab="SNP effect",xlab="CHR pos (bp)",main=paste0("CHR ",idx_chr,' - Trait: ',vect_P_traits[i]," - ",i_env))
      abline(h=0,lty=2,col="orange",lwd=2)
      vcol=rep(adjustcolor('grey',alpha=.3),nrow(all_GWAS))
      points(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="p",col=vcol[kept_chr])
      
      # Plot singletons
      df.CHR.single.trait <- subset(df.CHR.single,df.CHR.single[,paste0("is_sign_T",(i+i_list_env))])
      points(df.CHR.single.trait$ps,df.CHR.single.trait[,paste0("beta_T",(i+i_list_env))],bg="black",pch=21)
      
      ## Plot clusters (no singletons)
      df.CHR.trait <- subset(df.CHR,df.CHR[,paste0("is_sign_T",(i+i_list_env))])
      points(df.CHR.trait$ps,df.CHR.trait[,paste0("beta_T",(i+i_list_env))],col=all_col_CHR[df.CHR.trait$cluster],pch=all_pch_CHR[df.CHR.trait$cluster])
      
    }
  
  dev.off()
}

#### Output the clusters
write.table(all_cluster_merged,file=paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',row.names=FALSE,quote=FALSE)
}
### Export statistics on the clustering ?
