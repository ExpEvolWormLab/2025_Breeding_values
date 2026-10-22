rm(list=ls());gc()
library(qqman)
library(data.table)
library(dplyr)
library(igraph)
library(parallel)
threshold_clust=.2

selection_clusters = read.table(paste0("txt_GridLMM_NaCl/all_sign_SELECTION_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
dim(selection_clusters) # 3k selected SNPs

diff_signs = tapply(selection_clusters$s,selection_clusters$group,function(x){length(unique(sign(x)))})
groups_diff_signs = names(diff_signs)[as.numeric(diff_signs)==2]

load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
GWAS_results = read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
table(GWAS_results$sign_overall) # 33k sign SNPs

## Overlap between the two ?

dim(subset(selection_clusters,ID%in%GWAS_results$rs)) # 895 out of 3k / More than just by chance ?

length(unique(subset(selection_clusters,ID%in%GWAS_results$rs)$group)) # 30 clusters out of 109
#SF
length(unique(subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T1)$rs)$group)) # 11 clusters
#SB
length(unique(subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T2)$rs)$group)) # 14 clusters
#Area
length(unique(subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T7)$rs)$group)) # 7 clusters

#Beta
length(unique(subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T8)$rs)$group)) # 14 clusters

all_sign_clusters = data.frame(group=unique(subset(selection_clusters,ID%in%GWAS_results$rs)$group))
all_sign_clusters$isSF = all_sign_clusters$group%in%subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T1)$rs)$group
all_sign_clusters$isFS = all_sign_clusters$group%in%subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T2)$rs)$group
all_sign_clusters$isArea = all_sign_clusters$group%in%subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T7)$rs)$group
all_sign_clusters$isBeta = all_sign_clusters$group%in%subset(selection_clusters,ID%in%subset(GWAS_results,is_sign_T8)$rs)$group

head(all_sign_clusters)

## Ad all other clusters and do a Venn Diagram ?
all_non_sign_clusters=data.frame(group=unique(subset(selection_clusters,!group%in%all_sign_clusters$group)$group),
                                 isSF=FALSE,isFS=FALSE,isArea=FALSE,isBeta=FALSE)
all_clusters=rbind(all_sign_clusters,all_non_sign_clusters)
library(gplots);library(limma)
pdf("VennDiagram_Gwas_Selection.pdf",h=3,w=5)

vennDiagram(vennCounts(all_clusters[,2:4]),names=c("SF","FS","Area"))
## Compare it to the distribution of shared GWAS clusters
vennDiagram(vennCounts(subset(all_GWAS_clusters,is_sign_SF_NaCl | is_sign_SB_NaCl | is_sign_area.F_NaCl)[,c("is_sign_SF_NaCl","is_sign_SB_NaCl", "is_sign_area.F_NaCl")]),names=c("SF","FS","Area"))
dev.off()

par(mfrow=c(1,1))
vennDiagram(vennCounts(all_clusters[,2:5]),names=c("SF","FS","Area","Beta"))

table(rowSums(all_sign_clusters[,2:5])) # 1 cluster affect all 4 traits
table(rowSums(all_sign_clusters[,2:4])) # only 1 cluster affect all 3; only 4 clusters are missed by ignoring beta
subset(all_sign_clusters,rowSums(all_sign_clusters[,2:4])==3) # group 151

dim(subset(selection_clusters,group==151)) # 94 SNPs

dim(subset(GWAS_results,rs%in%subset(selection_clusters,group==151)$ID)) # 39 SNPs significant in this group for GWAS

## Manhattan plot with highlighted significant SNPs for GWAS
temp.env=new.env()
load('RData_GridLMM_NaCl/Compute_S_parallel.RData',envir=temp.env)
results_clust = temp.env$results
rm(temp.env);gc()

results_clust$P_corrected=p.adjust(results_clust$P,method='BH')
results_clust$CHR=tstrsplit(results_clust$rs,"[.]")[[1]]
results_clust$BP=as.numeric(tstrsplit(results_clust$rs,"[.]")[[2]])
table(results_clust$P_corrected<0.1) # 68k SNPs
table(results_clust$P_corrected<0.01) # 28k SNPs

data_manhattan = data.frame(CHR=results_clust$CHR,
                            BP=results_clust$BP,
                            P=as.numeric(results_clust$P_corrected))
data_manhattan$CHR=as.numeric(factor(data_manhattan$CHR,levels=c("I","II","III","IV","V","X")))
data_manhattan=na.omit(data_manhattan)
data_manhattan$SNP=paste(c("I","II","III","IV","V","X")[data_manhattan$CHR],data_manhattan$BP,sep='_')

### All manhattan
manhattan(subset(data_manhattan,P<0.1))

selected_at_least_one = unique(subset(all_sign_clusters,rowSums(all_sign_clusters[,2:4])>=1)$group)
selected_IDs = subset(GWAS_results,is_sign_T1 | is_sign_T2 | is_sign_T7)$rs

## Highlight significant for 1 of the three selected traits 
manhattan(subset(data_manhattan,P<0.1),highlight = selected_IDs)

### Do a plot for each chromosome ?
library(Polychrome)
library(colorspace)
cols <- glasbey.colors(30)
# Convert to HCL and keep reasonably dark colors
hcl <- hex2RGB(cols) |> as("polarLUV")
cols <- cols[hcl@coords[, "L"] < 75]


chr_idx = "IV"
pdf("Manhattan_plots_selection_with_GWAS.pdf",h=5,w=8)
for(chr_idx in c("I","II","III","IV","V","X")){
temp.man = subset(results_clust,CHR==chr_idx)
temp.man.group = subset(selection_clusters,CHROM==chr_idx)
temp.man.group$gp2=as.numeric(as.factor(temp.man.group$group))
#library(Polychrome)
#all_col_CHR <- sample(glasbey.colors(20))#sample(rainbow(length(unique(temp.man.group$cluster))))


all_col_CHR <- sample(glasbey.colors(length(unique(temp.man.group$group))))

plot(temp.man$BP,-log10(temp.man$P_corrected),pch=16,cex=.8,col=adjustcolor("gray",alpha=.2),bty="n",las=1,
     xlab=c("Chromosome Position"),ylab=c("-log10 Q value"),main=paste0("Chromosome :  ",chr_idx))
abline(h=5,lwd=1.5,lty=2)
## Add significant cluster

temp.man.group.non.GWAS = subset(temp.man.group,!ID%in%subset(GWAS_results,is_sign_T1 | is_sign_T2 | is_sign_T7)$rs)
points(temp.man.group.non.GWAS$POS,-log10(temp.man.group.non.GWAS$P_corrected),pch=16,col=adjustcolor(all_col_CHR[temp.man.group.non.GWAS$gp2],alpha=.4))


### And the ones with significant GWAS hits
temp.man.group.GWAS = subset(temp.man.group,ID%in%subset(GWAS_results,is_sign_T1 | is_sign_T2 | is_sign_T7)$rs)
dim(temp.man.group.GWAS)
length(unique(temp.man.group.GWAS$group)) # 2 clusters
points(temp.man.group.GWAS$POS,-log10(temp.man.group.GWAS$P_corrected),pch=8,col=all_col_CHR[temp.man.group.GWAS$gp2],cex=1.2)

}
dev.off()


## How many per CHR ?
for(chr_idx in c("I","II","III","IV","V","X")){
  temp.man = subset(results_clust,CHR==chr_idx)
  temp.man.group = subset(selection_clusters,CHROM==chr_idx)
  temp.man.group.GWAS = subset(temp.man.group,ID%in%subset(GWAS_results,is_sign_T1 | is_sign_T2 | is_sign_T7)$rs)
  print(chr_idx)
  print(length(unique(temp.man.group$group)))
  print(length(unique(temp.man.group.GWAS$group)))
}

dim(all_GWAS_clusters)
head(all_GWAS_clusters)
dim(GWAS_results)

### We can export a list of GWAS clusters that are significant - here meaning at least one SNP selected
dim(subset(GWAS_results,rs%in%selection_clusters$ID)) # 895 SNPs
length(unique(subset(GWAS_results,rs%in%selection_clusters$ID)$group)) # 40 GWAS clusters
write.table(unique(subset(GWAS_results,rs%in%selection_clusters$ID)$group),file='Selected_GWAS_clusters.txt')




