rm(list=ls())
gc()
library(data.table)
library(cluster)
library(dplyr)
library(igraph)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))


load("RData/GWAS_results_concatenated.RData")
###############################################
#### Create the clusters and save - Merged ####
###############################################

all_GWAS$CHROM=factor(tstrsplit(all_GWAS$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
## Re-order
all_GWAS=all_GWAS[order(all_GWAS$CHROM,all_GWAS$ps),]

sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:14)]
table(rowSums(sign_all_df)>=1) # ~108k SNPS at least sign for one trait

all_GWAS$sign_overall <- rowSums(sign_all_df)>=1


#####################################
## We need the matrix of genotypes ##
#####################################
load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))

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

table(all_GWAS$allele0==snp_positions$REF | all_GWAS$allele1==snp_positions$REF)
table(all_GWAS$allele0==snp_positions$ALT | all_GWAS$allele1==snp_positions$ALT)

### Revert the sign of the betas when the polarization is changed
vect_to_revert <- all_GWAS$allele0!=snp_positions$REF

all_GWAS$allele0[vect_to_revert]=snp_positions$REF[vect_to_revert]
all_GWAS$allele1[vect_to_revert]=snp_positions$ALT[vect_to_revert]
all_GWAS$af.NaCl[vect_to_revert]=1-all_GWAS$af.NaCl[vect_to_revert]
all_GWAS$af.NGM[vect_to_revert]=1-all_GWAS$af.NGM[vect_to_revert]

for(i_col in paste0("beta_T",1:14)){
  all_GWAS[,i_col] = -all_GWAS[,i_col]
}

### Save the txt file with the significant positions
write.table(snps_pos,file="data/sign_positions.txt",quote=FALSE,col.names = FALSE,row.names=FALSE,sep='\t')
#bcftools view -T sign_positions.txt -O z Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz > Final_VCF/Cemee2025.sign.GWAS.filt.vcf.gz

all_GWAS_sign <- subset(all_GWAS,sign_overall)

#barplot(table(rowSums(all_GWAS_sign[,paste0("is_sign_T",1:14)])))
#par(mfrow=c(1,2))
#rowS_NGM <- table(rowSums(all_GWAS_sign[,paste0("is_sign_T",1:7)]))
#rowS_NaCl <-table(rowSums(all_GWAS_sign[,paste0("is_sign_T",8:14)]))

#plot(rowS_NGM,type='b',col='cornflowerblue',ylim=c(50,60000),log="y")
#points(rowS_NaCl,type='b',col='darkorange')
#mean((rowSums(all_GWAS_sign[,paste0("is_sign_T",1:7)]))[(rowSums(all_GWAS_sign[,paste0("is_sign_T",1:7)]))>0])
#mean((rowSums(all_GWAS_sign[,paste0("is_sign_T",8:14)]))[(rowSums(all_GWAS_sign[,paste0("is_sign_T",8:14)]))>0])

table(snp_positions_sign$ID==all_GWAS_sign$rs) # In the same order #

### Now load the linkage blocks from plink

ld_snps <- read.table("linkage_diseq_filt/linkage_diseq_filt/ld_window_100s.hap.ld",h=TRUE)
### Create clusters from these
threshold_clust <- 0.8

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
length(unique(all_GWAS_sign_with_clust$cluster)) # 5993 clusters / 4691 clusters / 4358 clusters
sum(is.na(all_GWAS_sign_with_clust$cluster)) # ??  / 3759  / 3444 unique SNPs with no clusters

clust_length <- tapply(all_GWAS_sign_with_clust$ps,all_GWAS_sign_with_clust$clust,function(x){
  I(max(x)-min(x))
})
summary(clust_length)
clust_mean <- tapply(all_GWAS_sign_with_clust$ps,all_GWAS_sign_with_clust$clust,function(x){
  mean(x)
})
#plot(clust_mean,clust_length)
#################################################################
### Split the clusters if they are two far apparts ?
### Split the clusters if they affect the same trait in opposite direction
### (traits in each env. are considered different traits)
#################################################################

  #temp.snp.pos <- snps_pos
  #temp.snp.pos$idx = 1:nrow(temp.snp.pos)
  max_clust_nb=max(all_GWAS_sign_with_clust$cluster,na.rm=TRUE)

# test compatible lines
compatible_lines <- function(r1, r2) {
  all(r1 == r2 | r1==0 | r2==0)
}

  for(p in 1:max_clust_nb){
    print(p)
    # Get the position of the SNPs in the cluster
    sign_all=subset(all_GWAS_sign_with_clust,cluster==p)
    ## Get their effects for each trait
    sign_direction=sign_all[,paste0('beta_T',1:14)]
    ## Becomes null if not significant
    sign_direction[!sign_all[,paste0('is_sign_T',1:14)]]=0
    names(sign_direction)=vect_P_traits_bothEnv
    
    for(i_trait in vect_P_traits_bothEnv){
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
    sign_direction$nonzero=rowSums(sign_direction[,vect_P_traits_bothEnv]!=0)
    sign_direction=sign_direction[order(sign_direction$nonzero,decreasing=TRUE),]
    for (i_snp in 1:nrow(sign_direction)){
      if(is.na(sign_direction$group[i_snp])){
        sign_direction$group[i_snp] <- current_group
        if(i_snp<nrow(sign_direction)){
          for(j_snp in ((i_snp+1):nrow(sign_direction))){
            if(is.na(sign_direction$group[j_snp]) && compatible_lines(sign_direction[i_snp,vect_P_traits_bothEnv], sign_direction[j_snp, vect_P_traits_bothEnv])){
              ## Make sure it is compatible with all previous current group as well ##
              all_compatible=TRUE
              for(old_group in c(1:nrow(sign_direction))[sign_direction$group==current_group & !is.na(sign_direction$group)]){
                if(!compatible_lines(sign_direction[old_group,vect_P_traits_bothEnv], sign_direction[j_snp, vect_P_traits_bothEnv]))
                  all_compatible=FALSE 
              }
              if(all_compatible){
                sign_direction$group[j_snp] <- current_group
              }
            }
          }
        }
        current_group <-   current_group + 1
      }
    }
    current_group=current_group-1
    if(current_group>1){
      for(i_cg in 2:current_group){
        max_clust_nb=max_clust_nb+1
        all_GWAS_sign[all_GWAS_sign$rs%in%subset(sign_direction,group==i_cg)$rs,]$clust=max_clust_nb
      }
    }
  }  
   max_clust_nb # 27 -> 40  
  
  # for(p in unique(cluster_CHR)){
  #   # Get the position of the SNPs in the cluster
  #   temp.snp.pos.2=subset(temp.snp.pos,CHR==idx_chr & pos%in%subset(df.CHR,clust==p)$pos)
  #   if(nrow(temp.snp.pos.2)>1){
  #   # Distance between the SNPs
  #   for(pp in 2:nrow(temp.snp.pos.2)){
  #     dist.temp <- (temp.snp.pos.2$idx[pp]-temp.snp.pos.2$idx[pp-1])
  #     if(dist.temp>250){
  #       ## Replace all subsequent clust values
  #       print("create new cluster")
  #       max_clust_nb=max_clust_nb+1
  #       df.CHR$clust[df.CHR$pos%in% temp.snp.pos.2[pp:nrow(temp.snp.pos.2),]$pos]=max_clust_nb
  #     }
  #   }
  # }
  # }
  # max_clust_nb # 50
  #################################################################
  pdf(paste0("plots/Cluster_CHR_",idx_chr,"Merged-noPruning.pdf"),h=10,w=20)
  par(mfrow=c(2,6))
  
  all_col_CHR <- sample(rainbow(max(df.CHR$clust)))
  all_pch_CHR <- sample(1:25,max(df.CHR$clust),replace=TRUE)
  kept_chr=snps_pos$CHR==idx_chr
  for(i_env in c("NGM","NaCl")){
    if(i_env=='NGM'){
      i_list_env <- list_NGM_pruned3
    }else{
      i_list_env <- list_NaCl_pruned3
    }
    for(i in 1:6){
      
      plot(snps_pos$pos[kept_chr],i_list_env[[i]]$temp1[kept_chr],pch=16,cex=.8,type="n",ylab="SNP effect",xlab="CHR pos (bp)",main=paste0("CHR ",idx_chr,' - Trait: ',vect_P_traits[i]," - ",i_env))
      abline(h=0,lty=2,col="orange",lwd=2)
      vcol=rep(adjustcolor('grey',alpha=.3),nrow(snps_pos))
      points(snps_pos$pos[kept_chr],i_list_env[[i]]$temp1[kept_chr],pch=16,cex=.8,type="p",col=vcol[kept_chr])
      
      ## merge with cluster ?
      df=merge(df.CHR[,c("CHR","pos","clust")],cbind(snps_pos,data.frame(temp1=i_list_env[[i]]$temp1,sign=i_list_env[[i]]$sign)))
      df=subset(df,sign)
      
      df$clust2=as.numeric(as.factor(df$clust))
      points(df$pos,df$temp1,col=all_col_CHR[df$clust],pch=all_pch_CHR[df$clust])
      
    }
  }
  dev.off()
  write.table(df.CHR[,1:3],file=paste0('txt/all_clusters_CHR',idx_chr,'_Merged-no_Pruning.txt'),sep='\t',row.names=FALSE,quote=FALSE)
  
  df.with.effects=df.CHR
  df.with.signs=df.CHR
  for(i in 1:6){
    df.with.effects=merge(df.with.effects,cbind(snps_pos,data.frame(temp1=list_NGM_pruned3[[i]]$temp1)))
    names(df.with.effects)[ncol(df.with.effects)]=vect_P_traits_bothEnv[i]
    
    df.with.signs=merge(df.with.signs,cbind(snps_pos,data.frame(temp1=list_NGM_pruned3[[i]]$sign)))
    names(df.with.signs)[ncol(df.with.signs)]=paste0(vect_P_traits_bothEnv[i],"_sign")
    
  }
  for(i in 1:6){
    df.with.effects=merge(df.with.effects,cbind(snps_pos,data.frame(temp1=list_NaCl_pruned3[[i]]$temp1)))
    names(df.with.effects)[ncol(df.with.effects)]=vect_P_traits_bothEnv[i+6]
    
    df.with.signs=merge(df.with.signs,cbind(snps_pos,data.frame(temp1=list_NaCl_pruned3[[i]]$sign)))
    names(df.with.signs)[ncol(df.with.signs)]=paste0(vect_P_traits_bothEnv[i+6],"_sign")
    
  }
  
  write.table(df.with.effects,file=paste0('txt/all_clusters_CHR',idx_chr,'_Merged_with_effects_noPruning.txt'),sep='\t',row.names=FALSE,quote=FALSE)
  write.table(df.with.signs,file=paste0('txt/all_clusters_CHR',idx_chr,'_Merged_with_sign_noPruning.txt'),sep='\t',row.names=FALSE,quote=FALSE)
} # End of the clustering