rm(list=ls());gc()
library(data.table)
library(lme4)

threshold_clust =0.2
load(paste0("RData_GridLMM_NaCl/Beta_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
load(file="GridLMM_results_NaCl/Beta_Gwas.RData")

vect_P_traits=paste0("beta","_NaCl")

all_GWAS_clusters$sign_overall <- sum(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits)])
dim(all_GWAS_clusters)
table(all_GWAS_clusters$is_sign_beta_NaCl)
table(all_GWAS_clusters$sign_overall)

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
  #read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

df_clusters=read.table(paste0("data_BLUE_GridLMM_NaCl/Beta_data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters) # 219 clusters for R2

### Add the tip/arm/center tag
rkmap <- read.table("Caendr_filtering/rockman_kruglya_map_liftOver.txt",h=TRUE)
all_GWAS$CHROM=tstrsplit(all_GWAS$rs,"_")[[1]]
table(all_GWAS$CHROM)
all_GWAS$chr_type=NA
for(i in 1:nrow(rkmap)){
  if(sum(all_GWAS$CHROM==rkmap$CHROM[i] & all_GWAS$ps>=rkmap$Start[i] & all_GWAS$ps<=rkmap$End[i])>0){
  all_GWAS$chr_type[all_GWAS$CHROM==rkmap$CHROM[i] & all_GWAS$ps>=rkmap$Start[i] & all_GWAS$ps<=rkmap$End[i]]=rkmap$type[i]
  }
}
table(all_GWAS$chr_type)
table(is.na(all_GWAS$chr_type))
summary(subset(all_GWAS,is.na(chr_type) & CHROM=='V')[,c("CHROM","ps")]$ps) # Some tips on CHR V
all_GWAS$chr_type[all_GWAS$CHROM=="V" & all_GWAS$ps>=20341612]="right_tip"
table(subset(all_GWAS,is.na(chr_type))$CHROM) # Only MtDNA left

## Now also add the group ID to the GWAS
dim(all_GWAS)
dim(all_cluster_merged)
all_GWAS_sign <- (merge(all_GWAS,all_cluster_merged[,c("rs","group","is_sign_T1","POS")]))
dim(all_GWAS_sign)
names(all_cluster_merged)[!names(all_cluster_merged)%in%names(all_GWAS_sign)]

### How much is in Caendr  - Founders
#library(vcfR)
#founders_Caendr <- read.vcfR("Caendr_filtering/WI.20250625.hard-filter.isotype.CemeeFounders_poly.vcf.gz")
#caendr_pos <- as.data.frame(founders_Caendr@fix[,1:2])

#sign_for_Caendr <- all_GWAS_sign
#names(sign_for_Caendr)[2:5]=c("CHROM","POS","ALT","REF")
#dim(merge(sign_for_Caendr,data.frame(founders_Caendr@fix)))
#dim(merge(sign_for_Caendr,data.frame(founders_Caendr@fix[,1:2])))

gp_nbSNPs=tapply(all_GWAS_sign$group,all_GWAS_sign$group,length)
gp_min=tapply(all_GWAS_sign$ps,all_GWAS_sign$group,min)
gp_max=tapply(all_GWAS_sign$ps,all_GWAS_sign$group,max)
gp_type=tapply(all_GWAS_sign$chr_type,all_GWAS_sign$group,unique)
table(as.numeric(lapply(gp_type,length)))

#gp_type[as.numeric(lapply(gp_type,length))==2] # We decide center>arm>tip
for(i in c(1:length(gp_type))[as.numeric(lapply(gp_type,length))==2]){
  if(gp_type[[i]][1]=="left_tip"){
    gp_type[[i]]='left_arm'
  }else if(gp_type[[i]][1]=="left_arm"){
    gp_type[[i]]='center'
  }else{
    gp_type[[i]]=gp_type[[i]][[2]]
  }
}
table(as.numeric(lapply(gp_type,length))) # OK
gp_type=unlist(gp_type)

gp_length=data.frame(group=names(gp_min),clust_length=as.numeric(gp_max-gp_min),nb_SNPs=as.numeric(gp_nbSNPs),chr_type=gp_type)
plot(gp_length$nb_SNPs,gp_length$clust_length,pch=16,cex=.8,col=c("red","cornflowerblue","lightgreen","lightgreen","cornflowerblue")[as.numeric(as.factor(gp_length$chr_type))])

## ClusterID & group are the same thing
all_GWAS_clusters=merge(all_GWAS_clusters,gp_length)
head(all_GWAS_clusters)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]

hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+clust_length+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic
