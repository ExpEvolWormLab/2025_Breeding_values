rm(list=ls())
library(data.table)
library(poolSeq)
library(qqman)

### COMPARE THE RESULTS WITH OTHER A6140 - /mnt/data3/mallard/2025_Extra_sequencing/populations   ######
### I ended up creating a new sync file with some focal populations of interest ########################

### What should we do ? Two different directions : analyize the results (CMH ?) and do the overlap ?
### Or ask whether the significant SNPs for size are enriched for selected variants ? ###

### Alternative : Map the rotated phenotypes and ensure whether they are increasing ?


### Extract the corresponding sync file
a=read.table('population_sync/bam_to_sync.txt')$V1
a=(tstrsplit(a,'_'))[[1]]
a=(tstrsplit(a,'.bam'))[[1]]
a[1:3]=c("A00_early","A00","A00_late")

### I will now create the list of columns I want to keep in R and run a terminal command from R
### we want all the CA populations and the first three columns
### CA columns  index export ##

### Create a table with the corresponding populations

populations=data.frame(population=a,ancestor = "A00")
populations$generations = c(0,0,0,0,36,36,36,36)
populations$replicate = c(1:8)
rm(a)

### Load allele frequencies from PoolSeq data ###
### Load the sync file
GA.sync=read.sync('population_sync/A00_GA_poly_Cemee.sync',gen=populations$generations,repl=populations$replicate)

### First extract allele frequencies and coverage
all.af = af(GA.sync,gen=populations$generations,repl=populations$replicate)
all.cov = coverage(GA.sync,gen=populations$generations,repl=populations$replicate)
### Extract counts ##
## Reference allele counts
all.counts = round(all.af*all.cov)
## alternative allele counts
all.counts.alt = round((1-all.af)*all.cov)

### Position and chromosome information
snps_positions.sync=data.frame(ID=row.names(all.counts))
snps_positions.sync$CHROM=tstrsplit(snps_positions.sync$ID,"[.]")[[1]]
snps_positions.sync$POS=tstrsplit(snps_positions.sync$ID,"[.]")[[2]]
head(snps_positions.sync)
dim(snps_positions.sync) # 307,600 SNPS

## PCA between samples ?

pca1 <- prcomp(t(na.omit(all.af)),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population,col=rep(c("blue","red"),each=4)) ### Ok, we should only use the two others then

### CMH / Manhattan plot ?

vect_A00 <- sample(c(2,4),4,replace=TRUE)
vect_GA36 <- 5:8
GA36_cmh <- cmh.test(t(all.counts[,vect_A00]),t(all.counts.alt[,vect_A00]),t(all.counts[,vect_GA36]),t(all.counts.alt[,vect_GA36]),min.cov=30,max.cov=300)

CMH_GA36 = na.omit(data.frame(CHR=as.numeric(factor(snps_positions.sync$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync$POS),SNP=c(1:nrow(snps_positions.sync)),P=-log10(GA36_cmh)))
manhattan(subset(CMH_GA36,P>2), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE)

### Top peaks ?
int_SNPs <- CMH_GA36$CHR==1
plot(CMH_GA36$BP[int_SNPs],CMH_GA36$P[int_SNPs],pch=16)
subset(CMH_GA36,CHR==1 & P>15)

int_freq <- all.af[snps_positions.sync$CHROM=='I' & as.numeric(snps_positions.sync$POS)%in%subset(CMH_GA36,CHR==1 & P>15)$BP,]

par(mfrow=c(3,3))
for(i in 1:nrow(int_freq)){
  plot(int_freq[i,c(2,4:8)]~populations$generations[c(2,4:8)],pch=16,ylim=c(0,1))
}

### Allele frequency changes in GA populations

vect_medians_cov <- rowMeans(all.cov[,c(2,4:8)])
kept_cov <- vect_medians_cov>=50 & vect_medians_cov<=250

GA_af_changes <- all.af[kept_cov,c(5:8)] - rowMeans(all.af[kept_cov,c(2,4)])
GA_af_changes=na.omit(GA_af_changes)
dim(GA_af_changes) # 275,117 SNPs

plot(GA_af_changes[,1],GA_af_changes[,2],pch=16,col=adjustcolor("grey",alpha=.6))
heatmap(cor(na.omit(GA_af_changes)))

snps_positions.sync$rs=paste(snps_positions.sync$CHROM,snps_positions.sync$POS,sep='_')


### We want to estimate whether some founder haplotype are rising - how to link them with the GWAS haplotypes ?

### Load cluster data
all_cluster_merged = read.table(file=paste0('txt_BLUEs/all_sign_GWAS_with_clusters_R02.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged=all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
dim(all_cluster_merged)
length(unique(all_cluster_merged$group))

### 932k clusters

# Convert data.frames to data.tables for speed
all_cluster_merged <- as.data.table(all_cluster_merged)
snps_positions_sync <- as.data.table(snps_positions.sync)

# Convert to data.table if not already
all_cluster_merged <- as.data.table(all_cluster_merged)
snps_positions_sync.kept=data.frame(rs=row.names(GA_af_changes))
dim(snps_positions_sync.kept)
snps_positions_sync.kept$CHROM = tstrsplit(snps_positions_sync.kept$rs,'[.]')[[1]]
snps_positions_sync.kept$POS = as.numeric(tstrsplit(snps_positions_sync.kept$rs,'[.]')[[2]])
snps_positions_sync.kept$rs = paste(snps_positions_sync.kept$CHROM,snps_positions_sync.kept$POS,sep="_")
snps_positions_sync.kept <- as.data.table(snps_positions_sync.kept)

GA_af_changes=data.frame(GA_af_changes)
names(GA_af_changes)=paste("GA",1:4,"36",sep='')
GA_af_changes$rs=snps_positions_sync.kept$rs
GA_af_changes <- as.data.table(GA_af_changes) # Make sure all.af has column 'rs'

# Optionally add a key on rs to all_af for faster lookup
setkey(GA_af_changes, rs)

# For each group in all_cluster_merged, get all rs for that group
group_rs <- split(all_cluster_merged$rs, all_cluster_merged$group)

# Function to compute mean af per SNP group
calc_group_af <- function(rs_vector) {
  # Select matching rows from all_af
  matched <- GA_af_changes[rs %in% rs_vector]
  if (nrow(matched) > 1) {
    return(colMeans(matched[, !"rs", with = FALSE]))
  } else if (nrow(matched) == 1) {
    return(as.numeric(matched[, !"rs", with = FALSE]))
  } else {
    return(rep(NA, ncol(GA_af_changes) - 1))
  }
}

# Apply across groups (no error)
temp_GA_af <- lapply(group_rs, calc_group_af)
clusters_af = do.call(rbind,temp_GA_af)
colnames(clusters_af)=paste("GA",1:4,"36",sep='')

clusters_af=data.frame(clusters_af)
clusters_af$group=as.numeric(names(group_rs))

heatmap(cor(na.omit(clusters_af[,1:4])))

## Load GWAS results
load("RData_BLUEs//GWAS_cluster_results_BLUE_concatenated_R02.RData")
load(file="RData_BLUEs/GWAS_results_concatenated.RData")
GWAS_with_clust <- read.table("txt_BLUEs/all_sign_GWAS_with_clusters.txt",sep="\t",header=TRUE)

gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)])
all_GWAS_clusters$sign_NGM <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:7)])
all_GWAS_clusters$sign_NaCl <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",8:14)])

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters_R02.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

df_clusters=read.table("data_BLUE/data_clusters_BLUE_R02.txt",h=TRUE)
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)

dim(all_GWAS_clusters)
dim(merge(clusters_af,all_GWAS_clusters))

all_GWAS_clusters=merge(all_GWAS_clusters,clusters_af)
head(all_GWAS_clusters)

plot(all_GWAS_clusters$beta_T14,all_GWAS_clusters$GA136,pch=16,type="n")
abline(h=0,lty=2,col="darkorange",lwd=2)
points(all_GWAS_clusters$beta_T14[!all_GWAS_clusters$is_sign_T14],all_GWAS_clusters$GA136[!all_GWAS_clusters$is_sign_T14],pch=16)
points(all_GWAS_clusters$beta_T14[all_GWAS_clusters$is_sign_T14],all_GWAS_clusters$GA136[all_GWAS_clusters$is_sign_T14],pch=16,col="red")


plot(all_GWAS_clusters$beta_T7,all_GWAS_clusters$GA136,pch=16,type="n")
abline(h=0,lty=2,col="darkorange",lwd=2)
points(all_GWAS_clusters$beta_T7[!all_GWAS_clusters$is_sign_T7],all_GWAS_clusters$GA136[!all_GWAS_clusters$is_sign_T7],pch=16)
points(all_GWAS_clusters$beta_T7[all_GWAS_clusters$is_sign_T7],all_GWAS_clusters$GA136[all_GWAS_clusters$is_sign_T7],pch=16,col="red")


test=subset(all_GWAS_clusters,is_sign_T14)
test$GA_mean = (test$GA136 + test$GA236 + test$GA436)/3
plot((test$GA_mean),test$beta_T14)
abline(h=0)
abline(v=0)

test=subset(all_GWAS_clusters,is_sign_T7)
test$GA_mean = (test$GA136 + test$GA236 + test$GA436)/3
plot((test$GA_mean),test$beta_T7)
abline(h=0)
abline(v=0)

sum((test$GA136)*test$beta_T14,na.rm=TRUE)
sum((test$GA236)*test$beta_T14,na.rm=TRUE)
sum((test$GA336)*test$beta_T14,na.rm=TRUE)
sum((test$GA436)*test$beta_T14,na.rm=TRUE)

sum((test$GA136)*test$beta_T7,na.rm=TRUE)
sum((test$GA236)*test$beta_T7,na.rm=TRUE)
sum((test$GA336)*test$beta_T7,na.rm=TRUE)
sum((test$GA436)*test$beta_T7,na.rm=TRUE)


pdf("Trait_change_prediction_GA1_R02_sign.pdf")
bb=list()
for(k in 1:14){
  #k=14
  test=subset(all_GWAS_clusters,all_GWAS_clusters[,paste0("is_sign_T",k)])
  #test=all_GWAS_clusters
  
  bb_temp=NULL
  for(i in 1:10000){
    bb_temp=c(bb_temp,sum(sample(test$GA136)*test[,paste0("beta_T",k)],na.rm=TRUE))
  }
  bb[[k]]=bb_temp
  hist(bb_temp,main=paste("Trait Nb: ",k))
  abline(v=sum(test$GA136*test[,paste0("beta_T",k)],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:14){
  print(k)
  print(sum(test$GA136*test[,paste0("beta_T",k)],na.rm=TRUE)<sort(bb[[k]])[length(bb[[k]])*0.05] | sum(test$GA136*test[,paste0("beta_T",k)],na.rm=TRUE)>sort(bb[[k]])[length(bb[[k]])*0.95])
}


pdf("Trait_change_prediction_GA2_R02_sign.pdf")
cc=list()
for(k in 1:14){
  #test=all_GWAS_clusters
  test=subset(all_GWAS_clusters,all_GWAS_clusters[,paste0("is_sign_T",k)])
  cc_temp=NULL
  for(i in 1:10000){
    cc_temp=c(cc_temp,sum(sample(test$GA236)*test[,paste0("beta_T",k)],na.rm=TRUE))
  }
  cc[[k]]=cc_temp
  hist(cc_temp,main=paste("Trait Nb: ",k))
  abline(v=sum(test$GA236*test[,paste0("beta_T",k)],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:14){
  print(k)
  print(sum(test$GA236*test[,paste0("beta_T",k)],na.rm=TRUE)<sort(cc[[k]])[length(cc[[k]])*0.05] | sum(test$GA236*test[,paste0("beta_T",k)],na.rm=TRUE)>sort(cc[[k]])[length(cc[[k]])*0.95])
}

########

pdf("Trait_change_prediction_GA3_R02_sign.pdf")
dd=list()
for(k in 1:14){
  #test=all_GWAS_clusters
  test=subset(all_GWAS_clusters,all_GWAS_clusters[,paste0("is_sign_T",k)])
  dd_temp=NULL
  for(i in 1:10000){
    dd_temp=c(dd_temp,sum(sample(test$GA336)*test[,paste0("beta_T",k)],na.rm=TRUE))
  }
  dd[[k]]=dd_temp
  hist(dd_temp,main=paste("Trait Nb: ",k))
  abline(v=sum(test$GA336*test[,paste0("beta_T",k)],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:14){
  print(k)
  print(sum(test$GA336*test[,paste0("beta_T",k)],na.rm=TRUE)<sort(dd[[k]])[length(dd[[k]])*0.05] | sum(test$GA336*test[,paste0("beta_T",k)],na.rm=TRUE)>sort(dd[[k]])[length(dd[[k]])*0.95])
}

pdf("Trait_change_prediction_GA4_R02_sign.pdf")
ee=list()
for(k in 1:14){
  #test=all_GWAS_clusters
  test=subset(all_GWAS_clusters,all_GWAS_clusters[,paste0("is_sign_T",k)])
  ee_temp=NULL
  for(i in 1:10000){
    ee_temp=c(ee_temp,sum(sample(test$GA436)*test[,paste0("beta_T",k)],na.rm=TRUE))
  }
  ee[[k]]=ee_temp
  hist(ee_temp,main=paste("Trait Nb: ",k))
  abline(v=sum(test$GA436*test[,paste0("beta_T",k)],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:14){
  print(k)
  print(sum(test$GA436*test[,paste0("beta_T",k)],na.rm=TRUE)<sort(ee[[k]])[length(ee[[k]])*0.05] | sum(test$GA436*test[,paste0("beta_T",k)],na.rm=TRUE)>sort(ee[[k]])[length(ee[[k]])*0.95])
}

### What explains whether a SNP/cluster increase ? 