rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
library(parallel)
library(qqman)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

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
    
    pdf(paste0("plots_GridLMM_NaCl/Manhattan_GWAS_Trait_",traitJ,"_",env_label,".pdf"),w=10)
    manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
              highlight=data_manhattan$SNP[data_manhattan$is_sign],cex=.6,suggestiveline=-log10(0.2),genomewideline= -log10(0.3))
              abline(h=-log10(0.05),col="red")
    dev.off()
    
    if(trait_ID==1){
      all_GWAS=all_GWAS.temp
    }else{
      all_GWAS=merge(all_GWAS,all_GWAS.temp)
    }
    
  }

dim(all_GWAS)
save(list=c("all_GWAS"),file='GridLMM_results_NaCl/All_Gwas.RData')

for(i in 1:8){
  all_GWAS[,paste0("is_sign_T",i)] = p.adjust(all_GWAS[,paste0("pvalue_T",i)],method="BH")<0.2
  print(table(all_GWAS[,paste0("is_sign_T",i)]))
}

sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:8)]
table(rowSums(sign_all_df)>=1) #  33 k SNPS at least sign for one trait
table(rowSums(sign_all_df)) # up to 6 out of 8 traits

###############################################
#### Create the clusters and save - Merged ####
###############################################

## Re-order
all_GWAS=all_GWAS[order(all_GWAS$chr,all_GWAS$ps),]
all_GWAS$CHROM=all_GWAS$chr

sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:8)]
table(rowSums(sign_all_df)>=1) # ~ 33 k SNPS at least sign for one trait
all_GWAS$sign_overall <- rowSums(sign_all_df)>=1

#####################################
## We need the matrix of genotypes ##
#####################################
load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))
# Remove MtDNA
snp_positions=subset(snp_positions,ID%in%all_GWAS$rs)
snp.id.sign <- all_GWAS$rs[rowSums(sign_all_df)>=1]

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


### OLD FROM GEMMA
#table(all_GWAS$allele0==snp_positions$REF | all_GWAS$allele1==snp_positions$REF)
#table(all_GWAS$allele0==snp_positions$ALT | all_GWAS$allele1==snp_positions$ALT)

### Revert the sign of the betas when the polarization is changed
#vect_to_revert <- all_GWAS$allele0!=snp_positions$REF

#all_GWAS$allele0[vect_to_revert]=snp_positions$REF[vect_to_revert]
#all_GWAS$allele1[vect_to_revert]=snp_positions$ALT[vect_to_revert]
#all_GWAS$af.NaCl[vect_to_revert]=1-all_GWAS$af.NaCl[vect_to_revert]
#all_GWAS$af.NGM[vect_to_revert]=1-all_GWAS$af.NGM[vect_to_revert]

#for(i_col in paste0("beta_T",1:8)){
#  all_GWAS[,i_col] = -all_GWAS[,i_col]
#}

all_GWAS_sign <- subset(all_GWAS,sign_overall)
table(snp_positions_sign$ID==all_GWAS_sign$rs) # In the same order #

dim(snp_positions_sign)

### Save the txt file with the significant positions

write.table(snp_positions_sign,file="data_BLUE_GridLMM_NaCl/sign_positions.txt",quote=FALSE,col.names = FALSE,row.names=FALSE,sep='\t')
system("bcftools view -T data_BLUE_GridLMM_NaCl/sign_positions.txt -O z /Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz >  data_BLUE_GridLMM_NaCl/Cemee2025.sign.GWAS_BLUE.filt.vcf.gz")
system("vcftools --gzvcf data_BLUE_GridLMM_NaCl/Cemee2025.sign.GWAS_BLUE.filt.vcf.gz --hap-r2 --ld-window 100 --out data_BLUE_GridLMM_NaCl/ld_window_100")

### Now load the linkage blocks from plink

ld_snps <- read.table("data_BLUE_GridLMM_NaCl/ld_window_100.hap.ld",h=TRUE)
dim(ld_snps)

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
all_GWAS_sign$POS=all_GWAS_sign$ps
dim(all_GWAS_sign)
all_GWAS_sign_with_clust <- (merge(all_GWAS_sign,clusters_by_chr,all.x=TRUE))
dim(all_GWAS_sign_with_clust)
length(unique(all_GWAS_sign_with_clust$cluster)) 
sum(is.na(all_GWAS_sign_with_clust$cluster))

clust_length <- tapply(all_GWAS_sign_with_clust$ps,all_GWAS_sign_with_clust$clust,function(x){
  return(max(x)-min(x))
})
max(clust_length)/10^6 # 1.9 Mb
hist(log10(clust_length),main="",xlab=c("Log10 cluster length"))

clust_mean <- tapply(all_GWAS_sign_with_clust$ps,all_GWAS_sign_with_clust$clust,function(x){
  mean(x)
})

###############################################################################################
### Split the clusters if they are two far apparts ?
### Split the clusters if they affect the same trait in opposite direction
###  Yes, as they are polarised ok ?
### (traits in each env. are considered different traits)
###############################################################################################

max_clust_nb=max(all_GWAS_sign_with_clust$cluster,na.rm=TRUE)

split_clusters <- function(sign_all){
    ## Get their effects for each trait
    sign_direction=sign_all[,paste0('beta_T',1:8)]
    ## Becomes null if not significant
    sign_direction[!sign_all[,paste0('is_sign_T',1:8)]]=0
    names(sign_direction)=vect_P_traits
    
    for(i_trait in vect_P_traits){
      sign_direction[,i_trait]=sign(sign_direction[,i_trait])
       }
    
    sign_direction=cbind(sign_all$rs,sign_direction)
    names(sign_direction)[1]="rs"
    
    # test compatible lines
    compatible_lines <- function(r1, r2) {
      all(r1 == r2 | r1==0 | r2==0)
    }
    # identify groups
    current_group <- 1
    sign_direction$group=NA
    ## Reorder to put the max. non-zero first
    sign_direction$nonzero=rowSums(sign_direction[,vect_P_traits]!=0)
    sign_direction=sign_direction[order(sign_direction$nonzero,decreasing=TRUE),]
    
    ### Simplify the matrix in unique blocks
    sign_direction_simple <- unique(sign_direction[,2:ncol(sign_direction)])
  
    for (i_snp in 1:nrow(sign_direction_simple)){
      if(is.na(sign_direction_simple$group[i_snp])){
        sign_direction_simple$group[i_snp] <- current_group
        if(i_snp<nrow(sign_direction_simple)){
          for(j_snp in ((i_snp+1):nrow(sign_direction_simple))){
            if(is.na(sign_direction_simple$group[j_snp]) && compatible_lines(sign_direction_simple[i_snp,vect_P_traits], sign_direction_simple[j_snp, vect_P_traits])){
              ## Make sure it is compatible with all previous current group as well ##
              all_compatible=TRUE
              for(old_group in c(1:nrow(sign_direction_simple))[sign_direction_simple$group==current_group & !is.na(sign_direction_simple$group)]){
                if(!compatible_lines(sign_direction_simple[old_group,vect_P_traits], sign_direction_simple[j_snp, vect_P_traits]))
                  all_compatible=FALSE 
              }
              if(all_compatible){
                sign_direction_simple$group[j_snp] <- current_group
              }
            }
          }
        }
        current_group <-   current_group + 1
      }
    }
    sign_direction=(merge(subset(sign_direction,select= -group),sign_direction_simple))
    return(sign_direction)
}

clusters_split <- mclapply(split(all_GWAS_sign_with_clust, all_GWAS_sign_with_clust$cluster), split_clusters, mc.cores = detectCores()-2)

## Give a unique group name to each group within the list
clusters_split_gp_renamed <- lapply(clusters_split,function(x){
  x$group=paste(x$rs[1],x$group,sep="_")
  return(x)
})

all_cluster <- do.call(rbind, clusters_split_gp_renamed)
all_cluster=subset(all_cluster,select= -nonzero)

all_cluster_merged=merge(all_GWAS_sign_with_clust,all_cluster[,c("rs","group")],all.x=TRUE)
dim(all_cluster)
dim(all_cluster_merged)
dim(all_GWAS_sign_with_clust)

plot(density(as.numeric(table(all_cluster_merged$group))),log="x")
summary(as.numeric(table(all_cluster_merged$group)))

### Add the unique SNPs
all_cluster_merged[is.na(all_cluster_merged$cluster),]$group=all_cluster_merged[is.na(all_cluster_merged$cluster),]$rs
table(is.na(all_cluster_merged$group))
## Number of clusters initially ?  # 639 for R=0.2
length(unique(all_GWAS_sign_with_clust$cluster))+ sum(is.na(all_GWAS_sign_with_clust$cluster))
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
  
  pdf(paste0("plots_BLUEs_GridLMM_NaCl/Cluster_CHR_R",round(10*threshold_clust),"_",idx_chr,"_Merged-noPruning.pdf"),h=10,w=20)
  par(mfrow=c(2,8))
  
  all_col_CHR <- sample(rainbow(max(df.CHR$cluster)))
  all_pch_CHR <- sample(1:25,max(df.CHR$cluster),replace=TRUE)

  kept_chr = all_GWAS$CHROM==idx_chr
  i_env="NaCl"
  i_list_env <- 0
  
    for(i in 1:8){
      
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
}
#### Output the clusters
write.table(all_cluster_merged,file=paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',row.names=FALSE,quote=FALSE)
}
### Export statistics on the clustering ?
