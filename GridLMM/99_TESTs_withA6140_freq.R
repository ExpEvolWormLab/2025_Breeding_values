rm(list=ls());gc()
library(gplots)
library(limma) 

threshold_clust=0.8
all_cluster_merged = read.table(paste0("txt_BLUEs_GridLMM_CHR//all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))

#######################################
### Load the A6140 allele frequencies ?
#######################################

a6140_freq = read.table("data/A6140_frequencies.txt",sep='\t',h=TRUE)[,c(1:4,6,7)]
names(a6140_freq)[5:6]=c("a6140_freq","a6140_cov")
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_cluster_merged)
dim(all_GWAS_clusters)

## Polarization ##
temp.env=new.env()
load('RData/Preparatory_GLMM.RData',envir=temp.env)
geno_matrix = temp.env$genotype_matrix
rm(temp.env);gc()
A6160_Cemee_freq = data.frame(rs=colnames(geno_matrix),freq_Cemee=as.numeric(colMeans(geno_matrix[substring(rownames(geno_matrix),1,5)=="A6140",])))

### Some are not the same alleles ? ->> get the SNPs
temp.env=new.env()
load('RData/Genotype_info_Cemee.RData',envir=temp.env)
snp_positions=as.data.frame(temp.env$snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))
names(snp_positions)[3]="rs"
rm(temp.env);gc()
A6160_Cemee_freq=merge(A6160_Cemee_freq,snp_positions[,c("rs","REF","ALT")])

head(A6160_Cemee_freq)
  
all_cluster_merged_with_freq <- merge(all_cluster_merged,a6140_freq)
all_cluster_merged_with_freq=merge(all_cluster_merged_with_freq,A6160_Cemee_freq)

table(all_cluster_merged_with_freq$REF==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$ALT==all_cluster_merged_with_freq$minor)
table(all_cluster_merged_with_freq$ALT==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$REF==all_cluster_merged_with_freq$minor)

to_remove = !(all_cluster_merged_with_freq$REF==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$REF==all_cluster_merged_with_freq$minor) |
  ! (all_cluster_merged_with_freq$ALT==all_cluster_merged_with_freq$major | all_cluster_merged_with_freq$ALT==all_cluster_merged_with_freq$minor)

all_cluster_merged_with_freq=subset(all_cluster_merged_with_freq,!to_remove)
dim(all_cluster_merged_with_freq) # 73.250

## Now re-polarize ##

all_cluster_merged_with_freq[all_cluster_merged_with_freq$REF==all_cluster_merged_with_freq$minor,]$a6140_freq = 
  1- (all_cluster_merged_with_freq[all_cluster_merged_with_freq$REF==all_cluster_merged_with_freq$minor,]$a6140_freq)

plot(all_cluster_merged_with_freq$a6140_freq,all_cluster_merged_with_freq$freq_Cemee/2,xlab="Freq in A6140",ylab="Freq in CeMEE",
     pch=16,col=adjustcolor("grey",alpha=.3),type="n")

disjoint=abs(all_cluster_merged_with_freq$a6140_freq-all_cluster_merged_with_freq$freq_Cemee/2)>.15
points(all_cluster_merged_with_freq$a6140_freq[disjoint],I(all_cluster_merged_with_freq$freq_Cemee/2)[disjoint],xlab="Freq in A6140",ylab="Freq in CeMEE",
     pch=16,col=adjustcolor("red",alpha=.3))
points(all_cluster_merged_with_freq$a6140_freq[!disjoint],I(all_cluster_merged_with_freq$freq_Cemee/2)[!disjoint],xlab="Freq in A6140",ylab="Freq in CeMEE",
       pch=16,col=adjustcolor("grey",alpha=.3))
table(disjoint)


## Group frequencies
group_Freq <- tapply(all_cluster_merged_with_freq$a6140_freq,all_cluster_merged_with_freq$group,mean,na.rm=TRUE)
group_Freq=data.frame(group=names(group_Freq),a6140_freq=1-as.numeric(group_Freq))
dim(merge(all_GWAS_clusters,group_Freq)) # 736 with frequencies

group_Freq_noD <- tapply(all_cluster_merged_with_freq$a6140_freq[!disjoint],all_cluster_merged_with_freq$group[!disjoint],mean,na.rm=TRUE)
group_Freq_noD=data.frame(group=names(group_Freq_noD),a6140_freq=1-as.numeric(group_Freq_noD))
dim(merge(all_GWAS_clusters,group_Freq_noD)) # 676 out of 736 with frequencies

all_GWAS_clusters_with_freq = merge(all_GWAS_clusters,group_Freq)
all_GWAS_clusters_with_freq_noD = merge(all_GWAS_clusters,group_Freq_noD)

## Load the frequencies of the clusters as well
genotyped_matrix_clust <- read.table(paste0("data_BLUE_clusters_CHR/genotypes_matrix_clusters_R",round(10*threshold_clust),".txt"))
A6160_Cemee_freq_clusters = data.frame(rs=rownames(genotyped_matrix_clust),freq_Cemee=2-as.numeric(rowMeans(genotyped_matrix_clust[,substring(colnames(genotyped_matrix_clust),1,5)=="A6140"])))

all_GWAS_clusters_with_freq=merge(all_GWAS_clusters_with_freq,A6160_Cemee_freq_clusters)
all_GWAS_clusters_with_freq_noD=merge(all_GWAS_clusters_with_freq_noD,A6160_Cemee_freq_clusters)

plot(all_GWAS_clusters_with_freq$a6140_freq,all_GWAS_clusters_with_freq$freq_Cemee/2,pch=16)
abline(a=0,b=1)
plot(all_GWAS_clusters_with_freq_noD$a6140_freq,all_GWAS_clusters_with_freq_noD$freq_Cemee/2,pch=16)
abline(a=0,b=1)

#plot(density(subset(all_cluster_merged,is_sign_area.F_NGM)$beta_area.F_NGM))
#lines(density(subset(all_cluster_merged,is_sign_area.F_NaCl)$beta_area.F_NaCl),col="orange")

#boxplot(c(abs(subset(all_cluster_merged,is_sign_area.F_NGM)$beta_area.F_NGM),
#        abs(subset(all_cluster_merged,is_sign_area.F_NaCl)$beta_area.F_NaCl))~rep(c("NGM","NaCl"),c(sum(all_cluster_merged$is_sign_area.F_NGM),sum(all_cluster_merged$is_sign_area.F_NaCl))))

#summary(abs(subset(all_cluster_merged,is_sign_area.F_NGM)$beta_area.F_NGM))
#summary(abs(subset(all_cluster_merged,is_sign_area.F_NaCl)$beta_area.F_NaCl))

### Should be done at the cluster level ###

plot(density(subset(all_GWAS_clusters,is_sign_area.F_NGM)$beta_area.F_NGM))
lines(density(subset(all_GWAS_clusters,is_sign_area.F_NaCl)$beta_area.F_NaCl),col="orange")

boxplot(c(abs(subset(all_GWAS_clusters,is_sign_area.F_NGM)$beta_area.F_NGM),
          abs(subset(all_GWAS_clusters,is_sign_area.F_NaCl)$beta_area.F_NaCl))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_area.F_NGM),sum(all_GWAS_clusters$is_sign_area.F_NaCl))))

#boxplot(c((subset(all_GWAS_clusters,is_sign_area.F_NGM & beta_area.F_NGM<0)$beta_area.F_NGM),
#          (subset(all_GWAS_clusters,is_sign_area.F_NaCl & beta_area.F_NaCl<0)$beta_area.F_NaCl))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_area.F_NGM & all_GWAS_clusters$beta_area.F_NGM<0),sum(all_GWAS_clusters$is_sign_area.F_NaCl & all_GWAS_clusters$beta_area.F_NaCl<0))))

#boxplot(c((subset(all_GWAS_clusters,is_sign_area.F_NGM & beta_area.F_NGM>0)$beta_area.F_NGM),
#          (subset(all_GWAS_clusters,is_sign_area.F_NaCl & beta_area.F_NaCl>0)$beta_area.F_NaCl))~rep(c("NGM","NaCl"),c(sum(all_GWAS_clusters$is_sign_area.F_NGM & all_GWAS_clusters$beta_area.F_NGM>0),sum(all_GWAS_clusters$is_sign_area.F_NaCl & all_GWAS_clusters$beta_area.F_NaCl>0))))

summary(abs(subset(all_GWAS_clusters,is_sign_area.F_NGM)$beta_area.F_NGM))
summary(abs(subset(all_GWAS_clusters,is_sign_area.F_NaCl)$beta_area.F_NaCl))


### Allele frequencies ?
all_GWAS_clusters_with_freq_area.F_NGM=subset(all_GWAS_clusters_with_freq,is_sign_area.F_NGM)
all_GWAS_clusters_with_freq_area.F_NGM$pq.A6140 = all_GWAS_clusters_with_freq_area.F_NGM$a6140_freq * (1-all_GWAS_clusters_with_freq_area.F_NGM$a6140_freq)
plot(all_GWAS_clusters_with_freq_area.F_NGM$pq.A6140,abs(all_GWAS_clusters_with_freq_area.F_NGM$beta_area.F_NGM),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_area.F_NGM$beta_area.F_NGM)+2])
plot(all_GWAS_clusters_with_freq_area.F_NGM$pq.A6140,1/(all_GWAS_clusters_with_freq_area.F_NGM$beta_area.F_NGM^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_area.F_NGM$beta_area.F_NGM)+2])


all_GWAS_clusters_with_freq_area.F_NaCl=subset(all_GWAS_clusters_with_freq,is_sign_area.F_NaCl)
all_GWAS_clusters_with_freq_area.F_NaCl$pq.A6140 = all_GWAS_clusters_with_freq_area.F_NaCl$a6140_freq * (1-all_GWAS_clusters_with_freq_area.F_NaCl$a6140_freq)

plot(all_GWAS_clusters_with_freq_area.F_NaCl$pq.A6140,abs(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl)+2])




all_GWAS_clusters_with_freq_area.F_NaCl=subset(all_GWAS_clusters_with_freq,is_sign_area.F_NaCl)
plot(all_GWAS_clusters_with_freq_area.F_NaCl$a6140_freq,all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl,pch=16)
all_GWAS_clusters_with_freq_area.F_NaCl$pq.A6140 = all_GWAS_clusters_with_freq_area.F_NaCl$a6140_freq * (1-all_GWAS_clusters_with_freq_area.F_NaCl$a6140_freq)
plot(all_GWAS_clusters_with_freq_area.F_NaCl$pq.A6140,abs(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl)+2])
plot(all_GWAS_clusters_with_freq_area.F_NaCl$pq.A6140,1/(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl)+2])

test.mod <- lm(I(1/all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl^2) ~ pq.A6140 + as.factor(sign(all_GWAS_clusters_with_freq_area.F_NaCl$beta_area.F_NaCl)),data=all_GWAS_clusters_with_freq_area.F_NaCl)
summary(test.mod)

### All SNPs ?
all_GWAS_clusters_with_freq$pq.A6140=all_GWAS_clusters_with_freq$a6140_freq * (1-all_GWAS_clusters_with_freq$a6140_freq)
all_GWAS_clusters_with_freq_noD$pq.A6140=all_GWAS_clusters_with_freq_noD$a6140_freq * (1-all_GWAS_clusters_with_freq_noD$a6140_freq)
all_GWAS_clusters_with_freq_noD$pq.Cemee=all_GWAS_clusters_with_freq_noD$freq_Cemee/2 * (1-all_GWAS_clusters_with_freq_noD$freq_Cemee/2)
plot(all_GWAS_clusters_with_freq$pq.A6140,1/(all_GWAS_clusters_with_freq$beta_area.F_NGM^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_with_freq$beta_area.F_NGM)+2],ylim=c(0.1,200000000),log="y")

plot(all_GWAS_clusters_with_freq_noD$pq.A6140,abs(all_GWAS_clusters_with_freq_noD$beta_area.F_NGM),pch=16,cex=1)
plot(all_GWAS_clusters_with_freq_noD$a6140_freq,abs(all_GWAS_clusters_with_freq_noD$beta_area.F_NGM),pch=16,cex=1)

summary(lm((all_GWAS_clusters_with_freq_noD$beta_area.F_NGM)^2~all_GWAS_clusters_with_freq_noD$pq.A6140))

### Now plot both T7/T14
all_GWAS_clusters_Area=subset(all_GWAS_clusters_with_freq,is_sign_area.F_NGM | is_sign_area.F_NaCl)
all_GWAS_clusters_Area$pq.A6140 = all_GWAS_clusters_Area$a6140_freq * (1-all_GWAS_clusters_Area$a6140_freq)

all_GWAS_clusters_Area$sign_effect_retained = all_GWAS_clusters_Area$beta_area.F_NGM
all_GWAS_clusters_Area$sign_effect_retained[all_GWAS_clusters_Area$is_sign_area.F_NaCl & !all_GWAS_clusters_Area$is_sign_area.F_NGM] = all_GWAS_clusters_Area$beta_area.F_NaCl[all_GWAS_clusters_Area$is_sign_area.F_NaCl & !all_GWAS_clusters_Area$is_sign_area.F_NGM]

plot(all_GWAS_clusters_Area$pq.A6140,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_area.F_NGM)+2])
plot(all_GWAS_clusters_Area$pq.A6140,all_GWAS_clusters_Area$sign_effect_retained,pch=16,col=c("blue","red","red")[1+all_GWAS_clusters_Area$is_sign_area.F_NGM])
plot(all_GWAS_clusters_Area$pq.A6140,1/(all_GWAS_clusters_Area$sign_effect_retained^2),pch=16,col=c("blue","black","red")[sign(all_GWAS_clusters_Area$beta_area.F_NGM)+2])




