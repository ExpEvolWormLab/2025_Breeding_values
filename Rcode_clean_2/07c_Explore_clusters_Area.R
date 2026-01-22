rm(list=ls());gc()
library(gplots)
library(limma) 

threshold_clust=0.6
all_cluster_merged=read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)


i_Trait = 7
i_chr = "II"

pdf("Test_CHRII_Area_NGM_R06.pdf",w=12)
all_cluster_merged_CHRI = subset(all_cluster_merged,CHROM==i_chr)
plot(all_cluster_merged_CHRI$POS,all_cluster_merged_CHRI[,paste0("beta_T",i_Trait)],pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3)
all_cluster_merged_CHRI_nosign = all_cluster_merged_CHRI[!all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
all_cluster_merged_CHRI_sign = all_cluster_merged_CHRI[all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
points(all_cluster_merged_CHRI_nosign$POS,all_cluster_merged_CHRI_nosign[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor("grey",alpha=.7))
gp_length_1 = (names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)==1])
gp_length_sup2 = names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)>1]

all_cluster_merged_CHRI_gp1 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_1)
all_cluster_merged_CHRI_gp2 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_sup2)
points(all_cluster_merged_CHRI_gp1$POS,all_cluster_merged_CHRI_gp1[,paste0("beta_T",i_Trait)],pch=16,col="black")
vect_col_GP <- sample(rainbow(length(gp_length_sup2)))[as.numeric(as.factor(as.character(all_cluster_merged_CHRI_gp2$group)))]

points(all_cluster_merged_CHRI_gp2$POS,all_cluster_merged_CHRI_gp2[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor(vect_col_GP,alpha=0.8))
dev.off()




i_Trait = 14
i_chr = "II"

pdf("Test_CHRII_Area_NaCl_R06.pdf",w=12)
all_cluster_merged_CHRI = subset(all_cluster_merged,CHROM==i_chr)
plot(all_cluster_merged_CHRI$POS,all_cluster_merged_CHRI[,paste0("beta_T",i_Trait)],pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3)
all_cluster_merged_CHRI_nosign = all_cluster_merged_CHRI[!all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
all_cluster_merged_CHRI_sign = all_cluster_merged_CHRI[all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
points(all_cluster_merged_CHRI_nosign$POS,all_cluster_merged_CHRI_nosign[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor("grey",alpha=.7))
gp_length_1 = (names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)==1])
gp_length_sup2 = names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)>1]

all_cluster_merged_CHRI_gp1 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_1)
all_cluster_merged_CHRI_gp2 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_sup2)
points(all_cluster_merged_CHRI_gp1$POS,all_cluster_merged_CHRI_gp1[,paste0("beta_T",i_Trait)],pch=16,col="black")
vect_col_GP <- sample(rainbow(length(gp_length_sup2)))[as.numeric(as.factor(as.character(all_cluster_merged_CHRI_gp2$group)))]

points(all_cluster_merged_CHRI_gp2$POS,all_cluster_merged_CHRI_gp2[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor(vect_col_GP,alpha=0.8))
dev.off()


#### Zoom IN
pdf("Zoom_CHRII_Area.pdf",w=12,h=9)
par(mfrow=c(2,1))
for(i_Trait in c(7,14)){
all_cluster_merged_CHRI = subset(all_cluster_merged,CHROM==i_chr & POS>13000000 & POS<14000000)
plot(all_cluster_merged_CHRI$POS,all_cluster_merged_CHRI[,paste0("beta_T",i_Trait)],pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3)
all_cluster_merged_CHRI_nosign = all_cluster_merged_CHRI[!all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
all_cluster_merged_CHRI_sign = all_cluster_merged_CHRI[all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
points(all_cluster_merged_CHRI_nosign$POS,all_cluster_merged_CHRI_nosign[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor("grey",alpha=.7))
gp_length_1 = (names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)==1])
gp_length_sup2 = names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)>1]

all_cluster_merged_CHRI_gp1 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_1)
all_cluster_merged_CHRI_gp2 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_sup2)
points(all_cluster_merged_CHRI_gp1$POS,all_cluster_merged_CHRI_gp1[,paste0("beta_T",i_Trait)],pch=16,col="black")
vect_col_GP <- sample(rainbow(length(gp_length_sup2)))[as.numeric(as.factor(as.character(all_cluster_merged_CHRI_gp2$group)))]
points(all_cluster_merged_CHRI_gp2$POS,all_cluster_merged_CHRI_gp2[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor(vect_col_GP,alpha=0.8))

}
dev.off()


plot(all_cluster_merged_CHRI$beta_T7,all_cluster_merged_CHRI$beta_T14,pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3);abline(v=0,lty=2,col="orange",lwd=3)
points(all_cluster_merged_CHRI_nosign$beta_T7,all_cluster_merged_CHRI_nosign$beta_T14,pch=16,col=adjustcolor('grey',alpha=0.7))
points(all_cluster_merged_CHRI_gp1$beta_T7,all_cluster_merged_CHRI_gp1$beta_T14,pch=16,col="black")
points(all_cluster_merged_CHRI_gp2$beta_T7,all_cluster_merged_CHRI_gp2$beta_T14,pch=16,col=adjustcolor(vect_col_GP,alpha=0.8))

#### Three overlapping QTLs in NGML - body size
pdf("Overlapping_QTLS_Area_NaCl.pdf")
i_Trait=14
all_cluster_merged_CHRI = subset(all_cluster_merged,CHROM==i_chr & POS>13500000 & POS<13600000)
plot(all_cluster_merged_CHRI$POS,all_cluster_merged_CHRI[,paste0("beta_T",i_Trait)],pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3)
all_cluster_merged_CHRI_nosign = all_cluster_merged_CHRI[!all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
all_cluster_merged_CHRI_sign = all_cluster_merged_CHRI[all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
points(all_cluster_merged_CHRI_nosign$POS,all_cluster_merged_CHRI_nosign[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor("grey",alpha=.7))
gp_length_1 = (names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)==1])
gp_length_sup2 = names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)>1]

all_cluster_merged_CHRI_gp1 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_1)
all_cluster_merged_CHRI_gp2 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_sup2)
points(all_cluster_merged_CHRI_gp1$POS,all_cluster_merged_CHRI_gp1[,paste0("beta_T",i_Trait)],pch=16,col="black")
vect_col_GP <- sample(rainbow(length(gp_length_sup2)))[as.numeric(as.factor(as.character(all_cluster_merged_CHRI_gp2$group)))]
points(all_cluster_merged_CHRI_gp2$POS,all_cluster_merged_CHRI_gp2[,paste0("beta_T",i_Trait)],pch=16,col=vect_col_GP)
dev.off()

pdf("Overlapping_QTLS_Area_NGM.pdf")
i_Trait=7
all_cluster_merged_CHRI = subset(all_cluster_merged,CHROM==i_chr & POS>13500000 & POS<13600000)
plot(all_cluster_merged_CHRI$POS,all_cluster_merged_CHRI[,paste0("beta_T",i_Trait)],pch=16,type='n')
abline(h=0,lty=2,col="orange",lwd=3)
all_cluster_merged_CHRI_nosign = all_cluster_merged_CHRI[!all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
all_cluster_merged_CHRI_sign = all_cluster_merged_CHRI[all_cluster_merged_CHRI[,paste0("is_sign_T",i_Trait)],]
points(all_cluster_merged_CHRI_nosign$POS,all_cluster_merged_CHRI_nosign[,paste0("beta_T",i_Trait)],pch=16,col=adjustcolor("grey",alpha=.7))
gp_length_1 = (names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)==1])
gp_length_sup2 = names(table(all_cluster_merged_CHRI_sign$group))[table(all_cluster_merged_CHRI_sign$group)>1]

all_cluster_merged_CHRI_gp1 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_1)
all_cluster_merged_CHRI_gp2 = subset(all_cluster_merged_CHRI_sign,group%in%gp_length_sup2)
points(all_cluster_merged_CHRI_gp1$POS,all_cluster_merged_CHRI_gp1[,paste0("beta_T",i_Trait)],pch=16,col="black")
vect_col_GP <- sample(rainbow(length(gp_length_sup2)))[as.numeric(as.factor(as.character(all_cluster_merged_CHRI_gp2$group)))]
points(all_cluster_merged_CHRI_gp2$POS,all_cluster_merged_CHRI_gp2[,paste0("beta_T",i_Trait)],pch=16,col=vect_col_GP)
dev.off()

all_cluster_merged_CHRI_gp2$col_clust = vect_col_GP

### Frequency in the panel of the different clusters + mean phenotypic differences between groups ?
## Load phenotypes
phenotypes = read.table("data_BLUE/gemma_pheno.txt",h=FALSE)
geno = read.table("data_BLUE/genotypes_for_gemma.txt")

### Keep only relevant positions
geno_sub = subset(geno,V1%in%all_cluster_merged_CHRI$rs)
freq_group=NULL
## For each main cluster ?
 for(i_group in unique(all_cluster_merged_CHRI_gp2$group)){
 geno_sub_group = subset(geno_sub,V1%in%subset(all_cluster_merged_CHRI_gp2,group==i_group)$rs)
 #dim(geno_sub_group) # 14
 if(is.null(freq_group)){
   freq_group <- data.frame(group=i_group,rs=geno_sub_group$V1,
                            freq=apply(geno_sub_group[,4:ncol(geno_sub_group)],1,mean,na.rm=TRUE),
                            col=unique(subset(all_cluster_merged_CHRI_gp2,group==i_group)$col_clust))
 }else{
 freq_group <- rbind(freq_group,data.frame(group=i_group,rs=geno_sub_group$V1,freq=apply(geno_sub_group[,4:ncol(geno_sub_group)],1,mean,na.rm=TRUE),col=unique(subset(all_cluster_merged_CHRI_gp2,group==i_group)$col_clust)))
 }
 } 
freq_group$group=as.factor(freq_group$group)
pdf("freq_overlapping_QTLs.pdf")
boxplot(freq_group$freq~freq_group$group,col=unique(freq_group[,c("group","col")])$col)
dev.off()

## Overlap in the lines for the clusters ?
i_group=1522
geno_sub_group = subset(geno_sub,V1%in%subset(all_cluster_merged_CHRI_gp2,group==i_group)$rs)
all_clust_1 <- colSums(geno_sub_group[,4:ncol(geno_sub_group)])/2==nrow(geno_sub_group)

i_group=1524
geno_sub_group = subset(geno_sub,V1%in%subset(all_cluster_merged_CHRI_gp2,group==i_group)$rs)
all_clust_2 <- colSums(geno_sub_group[,4:ncol(geno_sub_group)])/2==nrow(geno_sub_group)

i_group=1527
geno_sub_group = subset(geno_sub,V1%in%subset(all_cluster_merged_CHRI_gp2,group==i_group)$rs)
all_clust_3 <- colSums(geno_sub_group[,4:ncol(geno_sub_group)])/2==nrow(geno_sub_group)

i_group=1529
geno_sub_group = subset(geno_sub,V1%in%subset(all_cluster_merged_CHRI_gp2,group==i_group)$rs)
all_clust_4 <- colSums(geno_sub_group[,4:ncol(geno_sub_group)])/2==nrow(geno_sub_group)

### Venn Diagram

vennDiagram(vennCounts(cbind(all_clust_1,all_clust_2,all_clust_3)),main='All')

## Load the clusters results
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters)

all_GWAS_clusters_signZOOM = subset(all_GWAS_clusters,group%in%all_cluster_merged_CHRI_gp2$group)
all_GWAS_clusters_signZOOM$is_sign_T7
all_GWAS_clusters_signZOOM$is_sign_T14

dim(all_GWAS_clusters_signZOOM)


### Mean phenotypes in CeMEE

i_group=1524
geno_sub_group = subset(geno_sub,V1%in%subset(all_cluster_merged_CHRI_gp2,group==i_group)$rs)
boxplot(phenotypes[geno_sub_group[1,4:ncol(geno_sub_group)]==0,7],phenotypes[geno_sub_group[1,4:ncol(geno_sub_group)]==2,7])




