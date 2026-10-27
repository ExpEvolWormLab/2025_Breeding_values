rm(list=ls());gc()
library(qqman)
library(data.table)
library(dplyr)
library(igraph)
library(parallel)

load('RData_GridLMM_NaCl/Compute_S_parallel.RData')
### Manhattan plot ###


### We are also going to cluster the SNPs that are significant for selection

results$P_corrected=p.adjust(results$P,method='BH')
table(results$P_corrected<0.1) # 68k SNPs
table(results$P_corrected<0.01) # 28k SNPs

data_manhattan = data.frame(CHR=tstrsplit(results$rs,"[.]")[[1]],
                            BP=as.numeric(tstrsplit(results$rs,"[.]")[[2]]),
                            P=as.numeric(results$P_corrected))

data_manhattan$CHR=as.numeric(factor(data_manhattan$CHR,levels=c("I","II","III","IV","V","X")))
data_manhattan=na.omit(data_manhattan)
data_manhattan$SNP=1:nrow(data_manhattan)

manhattan(subset(data_manhattan,P<0.1))

## We will set the cutoff at q<0.00001
#log10(0.00001)

data_manhattan$ID=paste(c("I","II","III","IV","V","X")[data_manhattan$CHR],data_manhattan$BP,sep='_')
sign_positions = subset(data_manhattan, P<0.00001)
dim(sign_positions) # 3.3k SNPs

#####################################
## We need the matrix of genotypes ##
#####################################

load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))
# Remove MtDNA
snp_positions=subset(snp_positions,ID%in%sign_positions$ID)
snp.id.sign = sign_positions$ID

## Re-order
genotype_info=genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
snp_positions=snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

table(snp.id.sign%in%snp_positions$ID) # all present
g_sign=genotype_info[snp_positions$ID%in%snp.id.sign,]
snp_positions_sign=subset(snp_positions,ID%in%snp.id.sign)

snps_pos <- data.frame(CHR=snp_positions_sign$CHROM,pos=snp_positions_sign$POS)
snps_pos$CHR=as.factor(snps_pos$CHR)
snps_pos$pos=as.numeric(snps_pos$pos)

### Save the txt file with the significant positions

write.table(snp_positions_sign,file="data_BLUE_GridLMM_NaCl/sign_positions_Selection.txt",quote=FALSE,col.names = FALSE,row.names=FALSE,sep='\t')
system("bcftools view -T data_BLUE_GridLMM_NaCl/sign_positions_Selection.txt -O z /Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz >  data_BLUE_GridLMM_NaCl/Cemee2025.sign.GWAS_BLUE.filt_SELECTION.vcf.gz")
system("vcftools --gzvcf data_BLUE_GridLMM_NaCl/Cemee2025.sign.GWAS_BLUE.filt_SELECTION.vcf.gz --hap-r2 --ld-window 100 --out data_BLUE_GridLMM_NaCl/ld_window_100_SELECTION")

### Now load the linkage blocks from plink

ld_snps <- read.table("data_BLUE_GridLMM_NaCl/ld_window_100_SELECTION.hap.ld",h=TRUE)
dim(ld_snps)

results$CHROM=(tstrsplit(results$rs,"[.]")[[1]])
results$POS=as.numeric(tstrsplit(results$rs,"[.]")[[2]])
results$ID=paste(results$CHROM,results$POS,sep='_')
results_sign = subset(results,ID%in%subset(data_manhattan, P<0.00001)$ID)
### Create clusters from these
for(threshold_clust in c(0.2,0.4,0.5,0.6,0.8,0.9,0.99)){
  
  
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
  dim(results_sign)
  results_with_clust <- (merge(results_sign,clusters_by_chr,all.x=TRUE))
  dim(results_with_clust)
  length(unique(results_with_clust$cluster)) 
  sum(is.na(results_with_clust$cluster))
  
  clust_length <- tapply(results_with_clust$POS,results_with_clust$clust,function(x){
    return(max(x)-min(x))
  })
  max(clust_length)/10^6 # 2.5 Mb
  hist(log10(clust_length),main="",xlab=c("Log10 cluster length"))
  
  clust_mean <- tapply(results_with_clust$POS,results_with_clust$clust,function(x){
    mean(x)
  })
  
  ###############################################################################################
  ### Split the clusters if they are two far apparts ?
  ### Split the clusters if they affect the same trait in opposite direction
  ###  Yes, as they are polarised ok ?
  ### (traits in each env. are considered different traits)
  ###############################################################################################
  
  max_clust_nb=max(results_with_clust$cluster,na.rm=TRUE)
  results_with_clust$group=results_with_clust$cluster
  all_cluster_merged=results_with_clust
  dim(results_with_clust)
  
  plot(density(as.numeric(table(all_cluster_merged$group))),log="x")
  summary(as.numeric(table(all_cluster_merged$group)))
  
  ### Add the unique SNPs
  all_cluster_merged[is.na(all_cluster_merged$cluster),]$group=all_cluster_merged[is.na(all_cluster_merged$cluster),]$rs
  table(is.na(all_cluster_merged$group))
  ## Number of clusters initially ?  # 167 for R=0.2
  length(unique(results_with_clust$cluster))+ sum(is.na(results_with_clust$cluster))
  ## Number of clusters in the end ?
  length(unique(all_cluster_merged$group)) # 166
  
  all_cluster_merged=subset(all_cluster_merged,select= -ID1)
  ### Replace group ID by a number, first sort all_cluster_merged
  all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
  all_cluster_merged$group=as.numeric(factor(as.character(all_cluster_merged$group),levels=unique(all_cluster_merged$group)))
  
  pdf(paste0("plots_BLUEs_GridLMM_NaCl/Cluster_SELECTION_allCHR_R",round(10*threshold_clust),"_","Merged-noPruning.pdf"),h=10,w=20)
  par(mfrow=c(2,3))
  
  for(idx_chr in sort(unique(all_cluster_merged$CHROM))){
    
    print(idx_chr)
    # No singletons
    df.CHR=subset(all_cluster_merged,CHROM==idx_chr & !is.na(cluster))
    if(nrow(df.CHR)>0){
      df.CHR$group=as.numeric(as.factor(as.character(df.CHR$group)))
      # Singletons
      df.CHR.single=subset(all_cluster_merged,CHROM==idx_chr & is.na(cluster))
      df.CHR.single$group=as.numeric(as.factor(as.character(df.CHR.single$group)))
      
      
      all_col_CHR <- sample(rainbow(max(df.CHR$cluster)))
      all_pch_CHR <- sample(1:25,max(df.CHR$cluster),replace=TRUE)
      
      kept_chr = results$CHROM==idx_chr
      i_env="NaCl"
      i_list_env <- 0
      
        plot(results$POS[kept_chr],-log10(results$P_corrected[kept_chr]),pch=16,cex=.8,type="n",ylab="SNP effect",xlab="CHR pos (bp)",main=paste0("CHR ",idx_chr))
        #abline(h=0,lty=2,col="orange",lwd=2)
        vcol=rep(adjustcolor('grey',alpha=.3),nrow(results))
        points(results$POS[kept_chr],-log10(results$P_corrected[kept_chr]),pch=16,cex=.8,type="p",col=vcol[kept_chr])
        
        # Plot singletons
        points(df.CHR.single$POS,-log10(df.CHR.single$P_corrected),bg="black",pch=21)
        
        ## Plot clusters (no singletons)
        points(df.CHR$POS,-log10(df.CHR$P_corrected),col=all_col_CHR[df.CHR$cluster],pch=all_pch_CHR[df.CHR$cluster])
        
      }
    }
  dev.off()
  #### Output the clusters
  write.table(all_cluster_merged,file=paste0("txt_GridLMM_NaCl/all_sign_SELECTION_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',row.names=FALSE,quote=FALSE)
}
### Export statistics on the clustering ?












