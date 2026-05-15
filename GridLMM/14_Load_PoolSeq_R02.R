rm(list=ls())
library(data.table)
library(poolSeq)
library(qqman)
library(parallel)

### Extract the corresponding sync file
a=read.table('Poolseq_GA_2026/all_bam_list.txt')$V1
a=(tstrsplit(a,'_'))[[1]]
a=(tstrsplit(a,'.bam'))[[1]]
a[1:3]=c("A00_early","A00","A00_late")

## Removed 2 populations
a=a[c(2,4:15)]

### I will now create the list of columns I want to keep in R and run a terminal command from R
### we want all the CA populations and the first three columns
### CA columns  index export ##

### Create a table with the corresponding populations

populations=data.frame(population=a,ancestor = "A00")
populations$generations = c(0,0,10,50,36,10,50,36,10,36,10,50,36)
populations$replicate = c(1:2,3,3,3,4,4,4,5,5,6,6,6)
rm(a)

### Load allele frequencies from PoolSeq data ###
### Load the sync file
GA.sync=read.sync('Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync',gen=populations$generations,repl=populations$replicate)

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
dim(snps_positions.sync) # 693,377 SNPS

## Re-order populations
populations=populations[order(populations$replicate,populations$generations),]
populations$median_cov=colMedians(all.cov)

hist(colMedians(all.cov[,c(2,4:ncol(all.cov))]))
sort(colMedians(all.cov))

plot(density(all.cov[,3]),xlim=c(0,200),ylim=c(0,0.1))
for(i in 2:13){
  lines(density(all.cov[,i]),xlim=c(0,200))
}


## PCA between samples ?

pca1 <- prcomp(t(na.omit(all.af)),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population,col=rep(c("blue","red"),c(2,11))) ### Ok, we should only use the two others then

# Remove A6140
pca1 <- prcomp(t(na.omit(all.af[,c(1,3:13)])),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population[c(1,3:13)],col=rep(c("blue","red"),c(1,11))) ### Ok, we should only use the two others then


pca1 <- prcomp(t(na.omit(all.af[,c(1,3:13)])),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population[c(1,3:13)],col=rep(c("blue","red","gold","green","magenta"),c(1,3,3,2,3))) ### Ok, we should only use the two others then


### Strange variation in coverage ?
rel_cov=all.cov
for(i in 1:ncol(all.cov)){
  rel_cov[,i]=rel_cov[,i]/mean(rel_cov[,i],na.rm=TRUE)
}

snps_positions.sync$minCov <- rowMins(rel_cov)
snps_positions.sync$maxCov <- rowMaxs(rel_cov)
snps_positions.sync$medianCov <- rowMedians(rel_cov)

#plot(snps_positions.sync$minCov,snps_positions.sync$maxCov,pch=16,col=adjustcolor("grey",alpha=.5))

### We can only keep alleles with normal coverage variation ? 
#kept.alleles = maxCov>20 & maxCov<80 & minCov>1.5
#table(kept.alleles)
#table(snps_positions.sync$CHROM[kept.alleles]) # MtDNA

## No MtDNA
#temp=subset(snps_positions.sync,CHROM!="MtDNA")
#plot(temp$minCov,temp$maxCov,pch=16,col=adjustcolor("grey",alpha=.5))
#hist(log10(snps_positions.sync$medianCov))

kept.alleles = snps_positions.sync$maxCov<=4 & snps_positions.sync$medianCov<5 & snps_positions.sync$medianCov>.15
table(kept.alleles) # removed 4k SNPS including 2k MtDNA

#hist(log10(snps_positions.sync$medianCov[kept.alleles]))

all.af=all.af[rownames(all.af)%in%snps_positions.sync$ID[kept.alleles],]
all.cov=all.cov[rownames(all.cov)%in%snps_positions.sync$ID[kept.alleles],]
all.counts=all.counts[rownames(all.counts)%in%snps_positions.sync$ID[kept.alleles],]
all.counts.alt=all.counts.alt[rownames(all.counts.alt)%in%snps_positions.sync$ID[kept.alleles],]

snps_positions.sync=snps_positions.sync[kept.alleles,]

### Redo the PCA ?

pca1 <- prcomp(t(na.omit(all.af)),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population,col=rep(c("blue","red"),c(2,11))) ### Ok, we should only use the two others then

pca1 <- prcomp(t(na.omit(all.af[,c(1,3:13)])),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population[c(1,3:13)],col=rep(c("blue","red","gold","green","magenta"),c(1,3,3,2,3))) ### Ok, we should only use the two others then


## Keep only the Caendr SNPs ?
#library(vcfR)
#founders_Caendr <- read.vcfR("Caendr_filtering/WI.20250625.hard-filter.isotype.CemeeFounders_poly.vcf.gz")
#caendr_pos <- as.data.frame(founders_Caendr@fix[,1:2])
#write.table(caendr_pos,file='data/caendr_SNPs.txt',quote=FALSE,sep='\t',row.names=FALSE)
caendr_pos=read.table('data/caendr_SNPs.txt',sep='\t',h=TRUE)

all.af.caendr=all.af[rownames(all.af)%in%paste(caendr_pos$CHROM,caendr_pos$POS,sep='.'),]
all.cov.caendr=all.cov[rownames(all.cov)%in%paste(caendr_pos$CHROM,caendr_pos$POS,sep='.'),]
all.counts.caendr=all.counts[rownames(all.counts)%in%paste(caendr_pos$CHROM,caendr_pos$POS,sep='.'),]
all.counts.alt.caendr=all.counts.alt[rownames(all.counts.alt)%in%paste(caendr_pos$CHROM,caendr_pos$POS,sep='.'),]
snps_positions.sync.caendr=subset(snps_positions.sync,ID%in%rownames(all.counts.alt.caendr))
dim(snps_positions.sync.caendr)

dim(all.af.caendr) # 273k

pca1 <- prcomp(t(na.omit(all.af.caendr)),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population,col=rep(c("blue","red"),c(2,11))) ### Ok, we should only use the two others then

pca1 <- prcomp(t(na.omit(all.af.caendr[,c(1,3:13)])),center=TRUE,scale=FALSE)
plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
text(pca1$x[,1],pca1$x[,2],labels=populations$population[c(1,3:13)],col=rep(c("blue","red","gold","green","magenta"),c(1,3,3,2,3))) ### Ok, we should only use the two others then


### We will try the GLMS
library(lme4)
safe_glmer <- function(formula, data) {
  fit <- try(
    glmer(formula, data = data, family = "binomial"),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(c(NA, NA))
  
  coef <- summary(fit)$coef
  coef[2, c(1, 4)]
}

populations$replicate_fact=as.factor(populations$replicate)

vect_All=rep(TRUE,13)
vect_GA1=populations$replicate%in%c(1,2,3)
vect_GA2=populations$replicate%in%c(1,2,4)
vect_GA4=populations$replicate%in%c(1,2,6)

populations_GA1=subset(populations,replicate%in%c(1,2,3))
populations_GA2=subset(populations,replicate%in%c(1,2,4))
populations_GA4=subset(populations,replicate%in%c(1,2,6))

# res <- mclapply(1:nrow(all.counts.caendr), function(i) {
#   
#   f_All <- cbind(all.counts.caendr[i, vect_All],
#                  all.counts.alt.caendr[i, vect_All]) ~ generations + (1|replicate_fact)
#   
#   f_GA1 <- cbind(all.counts.caendr[i, vect_GA1],
#                  all.counts.alt.caendr[i, vect_GA1]) ~ generations + (1|replicate_fact)
#   
#   f_GA2 <- cbind(all.counts.caendr[i, vect_GA2],
#                  all.counts.alt.caendr[i, vect_GA2]) ~ generations + (1|replicate_fact)
# 
#   f_GA4 <- cbind(all.counts.caendr[i, vect_GA4],
#                  all.counts.alt.caendr[i, vect_GA4]) ~ generations + (1|replicate_fact)
#   
#   list(
#     All = safe_glmer(f_All, populations),
#     GA1 = safe_glmer(f_GA1, populations_GA1),
#     GA2 = safe_glmer(f_GA2, populations_GA2),
#     GA4 = safe_glmer(f_GA4, populations_GA4)
#   )
#   
# }, mc.cores = detectCores() - 2)
# 
# save(list=ls(),file='SaveGLM.RData')
load('SaveGLM.RData')

snps_positions.sync.caendr=subset(snps_positions.sync,ID%in%rownames(all.counts.alt.caendr))

results_All <- do.call(rbind, lapply(res, `[[`, "All"))
results_GA1 <- do.call(rbind, lapply(res, `[[`, "GA1"))
results_GA2 <- do.call(rbind, lapply(res, `[[`, "GA2"))
results_GA4 <- do.call(rbind, lapply(res, `[[`, "GA4"))

snps_positions.sync.caendr$GA_ALL=results_All[,2]
snps_positions.sync.caendr$GA1=results_GA1[,2]
snps_positions.sync.caendr$GA2=results_GA2[,2]
snps_positions.sync.caendr$GA4=results_GA4[,2]

GLM_All = (data.frame(CHR=as.numeric(factor(snps_positions.sync.caendr$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync.caendr$POS),SNP=c(1:nrow(snps_positions.sync.caendr)),P=-log10(snps_positions.sync.caendr$GA_ALL)))
GLM_GA1 = (data.frame(CHR=as.numeric(factor(snps_positions.sync.caendr$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync.caendr$POS),SNP=c(1:nrow(snps_positions.sync.caendr)),P=-log10(snps_positions.sync.caendr$GA1)))
GLM_GA2 = (data.frame(CHR=as.numeric(factor(snps_positions.sync.caendr$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync.caendr$POS),SNP=c(1:nrow(snps_positions.sync.caendr)),P=-log10(snps_positions.sync.caendr$GA2)))
GLM_GA4 = (data.frame(CHR=as.numeric(factor(snps_positions.sync.caendr$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync.caendr$POS),SNP=c(1:nrow(snps_positions.sync.caendr)),P=-log10(snps_positions.sync.caendr$GA4)))

manhattan(subset(GLM_All,P>2 & is.finite(P)), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE)
manhattan(subset(GLM_GA1,P>2 & is.finite(P)), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE)

#subset(GLM_All,P>100)
#int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='X' & as.numeric(snps_positions.sync.caendr$POS)%in%subset(GLM_All,CHR==6 & P>200 & is.finite(P))$BP,]

int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='II' 
                          & as.numeric(snps_positions.sync.caendr$POS)==13649880,]
plot(int_freq~populations$generations,pch=16,ylim=c(0,1))

int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='IV' 
                          & as.numeric(snps_positions.sync.caendr$POS)==2129332,]
plot(int_freq~populations$generations,pch=16,ylim=c(0,1))

#######################################
### Keep only > 1% all the time ? #####
#######################################

base_freq = rowMins(all.af.caendr)
base_freq = rowMins(all.af.caendr)
kept_base_freq = base_freq>0.01

manhattan(subset(GLM_All,P>2 & is.finite(P) & kept_base_freq), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE)
manhattan(subset(GLM_GA1,P>2 & is.finite(P) & kept_base_freq), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE)


subset(GLM_GA1,P>30 & is.finite(P) & kept_base_freq)
int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='V' 
                          & as.numeric(snps_positions.sync.caendr$POS)==15966464,]
plot(int_freq~populations$generations,pch=16,ylim=c(0,1))

int_POS <- subset(GLM_GA1,P>10 & is.finite(P) & kept_base_freq & CHR==5 & BP>15000000 & BP<17000000)$BP

int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='V' 
                          & as.numeric(snps_positions.sync.caendr$POS)%in%int_POS,]

plot(colMeans(int_freq)~populations$generations,type="n",ylim=c(0,1))
for(i in 1:nrow(int_freq)){
  points(int_freq[i,1:2]~c(0,3),pch=16,type='b')  
  for(j in 3:6){
  points(int_freq[i,populations$replicate==j]~populations$generations[populations$replicate==j],type='b',col=c("black","black","red","blue","darkgreen","orange")[j])
  }
}
temp.glm=subset(GLM_All,P>2 & is.finite(P) & kept_base_freq & CHR==1 & BP<3000000 & BP>2500000)
plot(temp.glm$BP,temp.glm$P,pch=16 )


int_POS <- subset(GLM_All,P>9 & is.finite(P) & kept_base_freq & CHR==1 & BP<3000000 & BP>2500000)$BP
int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='I' 
                          & as.numeric(snps_positions.sync.caendr$POS)%in%int_POS,]

plot(colMeans(int_freq)~populations$generations,type="n",ylim=c(0,1))
for(i in 1:nrow(int_freq)){
  points(int_freq[i,1:2]~c(0,3),pch=16,type='b')  
  for(j in 3:6){
    points(int_freq[i,populations$replicate==j]~populations$generations[populations$replicate==j],type='b',col=c("black","black","red","blue","darkgreen","orange")[j])
  }
}

temp.glm=subset(GLM_All,P>2 & is.finite(P) & kept_base_freq & CHR==3 & BP<8000000 & BP>4000000)
plot(temp.glm$BP,temp.glm$P,pch=16)

int_POS <- subset(GLM_All,P>10 & is.finite(P) & kept_base_freq & CHR==3  & BP<8000000 & BP>4000000)$BP
int_freq <- all.af.caendr[snps_positions.sync.caendr$CHROM=='III' 
                          & as.numeric(snps_positions.sync.caendr$POS)%in%int_POS,]

## Increasing means
diff=rowMeans(int_freq[,populations$generations==50])-rowMeans(int_freq[,populations$generations==10])
int_freq[diff<0,]=1-int_freq[diff<0,]

plot(colMeans(int_freq)~populations$generations,type="n",ylim=c(0,1))
for(i in 1:nrow(int_freq)){
  points(int_freq[i,1:2]~c(0,3),pch=16,type='b')  
  for(j in 3:6){
    points(int_freq[i,populations$replicate==j]~populations$generations[populations$replicate==j],type='b',col=c("black","black","red","blue","darkgreen","orange")[j])
  }
}


### We want to estimate whether some founder haplotype are rising - how to link them with the GWAS haplotypes ?
### Load cluster data
threshold_clust=0.2
all_cluster_merged = read.table(paste0("txt_BLUEs_GridLMM_CHR//all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged=all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
dim(all_cluster_merged)
length(unique(all_cluster_merged$group))
### 4.1k clusters

### TEST
all_cluster_merged$CHR=as.numeric(factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA")))
all_cluster_merged$BP=all_cluster_merged$POS
dim(merge(all_cluster_merged,GLM_All))
clust_with_Man <- merge(all_cluster_merged,GLM_All)

## Are they in clusters ?
#pos_selected <- GLM_All[nb_repl_selected>=1,c("CHR","BP")]
#pos_selected$is_selected = TRUE
#pos_selected$nb_clust = nb_repl_selected[nb_repl_selected>=1]

is_selected=p.adjust(10^(-GLM_All$P),method="BH")
qqnorm(-log10(is_selected[!is.na(is_selected) & is.finite(is_selected) & is_selected>0]))
table(is_selected)

pos_selected=snps_positions.sync.caendr
pos_selected$rs=paste(pos_selected$CHROM,pos_selected$POS,sep='_')
pos_selected$is_selected=(is_selected<0.01 & !is.na(is_selected) & is.finite(is_selected) & is_selected>0)

pos_selected$nb_clust=1


all_cluster_merged$pleiotropy <- rowSums(all_cluster_merged[,paste0("is_sign_T",1:16)],na.rm=TRUE)
all_cluster_merged$pleiotropy_NaCl <- rowSums(all_cluster_merged[,paste0("is_sign_T",9:16)],na.rm=TRUE)
all_cluster_merged$pleiotropy_NGM <- rowSums(all_cluster_merged[,paste0("is_sign_T",1:8)],na.rm=TRUE)
all_cluster_merged_selected <- merge(all_cluster_merged,pos_selected)

all_cluster_merged_selected_uniq_clust <- unique(all_cluster_merged_selected[,c("is_selected","nb_clust","group","pleiotropy","pleiotropy_NaCl","pleiotropy_NGM")])
dim(all_cluster_merged_selected_uniq_clust)
#boxplot(all_cluster_merged_selected_uniq_clust$pleiotropy~all_cluster_merged_selected_uniq_clust$nb_clust)
#boxplot(all_cluster_merged_selected_uniq_clust$pleiotropy_NaCl~all_cluster_merged_selected_uniq_clust$nb_clust)
#boxplot(all_cluster_merged_selected_uniq_clust$pleiotropy_NGM~all_cluster_merged_selected_uniq_clust$nb_clust)

#tapply(all_cluster_merged_selected_uniq_clust$pleiotropy,all_cluster_merged_selected_uniq_clust$nb_clust,mean)
#tapply(all_cluster_merged_selected_uniq_clust$pleiotropy_NaCl,all_cluster_merged_selected_uniq_clust$nb_clust,mean)
#tapply(all_cluster_merged_selected_uniq_clust$pleiotropy_NGM,all_cluster_merged_selected_uniq_clust$nb_clust,mean)

mean(all_cluster_merged_selected_uniq_clust$pleiotropy)
mean(all_cluster_merged_selected_uniq_clust$pleiotropy_NaCl)
mean(all_cluster_merged_selected_uniq_clust$pleiotropy_NGM)

### Top peaks ?

### Allele frequency changes in GA populations
vect_medians_cov <- rowMeans(all.cov.caendr)
kept_cov <- vect_medians_cov>=50 & vect_medians_cov<=250

idx_GA1_4_F10 = c(3,6,11)
idx_GA1_4_F50 = c(4,7,12)

GA_af_changes <- all.af.caendr[kept_cov,idx_GA1_4_F50] - (all.af.caendr[kept_cov,idx_GA1_4_F10])
GA_af_changes=na.omit(GA_af_changes)
dim(GA_af_changes) # 275,117 SNPs

plot(GA_af_changes[,1],GA_af_changes[,2],pch=16,col=adjustcolor("grey",alpha=.6))
heatmap(cor(na.omit(GA_af_changes)))

snps_positions.sync.caendr$rs=paste(snps_positions.sync.caendr$CHROM,snps_positions.sync.caendr$POS,sep='_')

kCHR=6; cCHR="X"
plot(subset(GLM_All,CHR==kCHR & kept_base_freq)$BP,subset(GLM_All,CHR==kCHR & kept_base_freq)$P,pch=21,bg=adjustcolor("grey",alpha=.5))
kept_pos=subset(GLM_All,CHR==kCHR & kept_base_freq)$BP
points(subset(clust_with_Man,CHROM==cCHR & BP%in%kept_pos)$BP,subset(clust_with_Man,CHROM==cCHR & BP%in%kept_pos)$P,pch=16,col=adjustcolor("red",alpha=.3))

kCHR=6; cCHR="X"
plot(subset(GLM_All,CHR==kCHR & kept_base_freq)$BP,subset(GLM_All,CHR==kCHR & kept_base_freq)$P,pch=21,bg=adjustcolor("grey",alpha=.5))
kept_pos=subset(GLM_All,CHR==kCHR & kept_base_freq)$BP
clust_with_Man_BetaNaCl=subset(clust_with_Man,is_sign_T16)
points(subset(clust_with_Man_BetaNaCl,CHROM==cCHR & BP%in%kept_pos)$BP,subset(clust_with_Man_BetaNaCl,CHROM==cCHR & BP%in%kept_pos)$P,pch=16,col=adjustcolor("red",alpha=.3))

# Convert data.frames to data.tables for speed
all_cluster_merged <- as.data.table(all_cluster_merged)
snps_positions.sync.caendr <- as.data.table(snps_positions.sync.caendr)

# Convert to data.table if not already
snps_positions_sync.kept=data.frame(rs=row.names(GA_af_changes))
dim(snps_positions_sync.kept)
snps_positions_sync.kept$CHROM = tstrsplit(snps_positions_sync.kept$rs,'[.]')[[1]]
snps_positions_sync.kept$POS = as.numeric(tstrsplit(snps_positions_sync.kept$rs,'[.]')[[2]])
snps_positions_sync.kept$rs = paste(snps_positions_sync.kept$CHROM,snps_positions_sync.kept$POS,sep="_")
snps_positions_sync.kept <- as.data.table(snps_positions_sync.kept)

GA_af_changes=data.frame(GA_af_changes)
names(GA_af_changes)=paste("GA",c(1,2,4),"36",sep='')
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
colnames(clusters_af)=paste("GA",c(1,2,4),"50",sep='')

clusters_af=data.frame(clusters_af)
clusters_af$group=as.numeric(names(group_rs))

heatmap(cor(na.omit(clusters_af[,1:3])))

## Load GWAS results

load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
GWAS_with_clust <- read.table("txt_BLUEs/all_sign_GWAS_with_clusters_R8.txt",sep="\t",header=TRUE)

#gemma_pheno=read.table('data_BLUE/gemma_pheno.txt',sep='\t',h=FALSE)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv)])
all_GWAS_clusters$sign_NGM <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv[1:8])])
all_GWAS_clusters$sign_NaCl <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv[9:16])])

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(paste0("txt_BLUEs_GridLMM_CHR//all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters)
dim(merge(clusters_af,all_GWAS_clusters))

all_GWAS_clusters=merge(all_GWAS_clusters,clusters_af)
head(all_GWAS_clusters)

plot(all_GWAS_clusters$beta_area.F_NaCl,all_GWAS_clusters$GA150,pch=16,type="n")
abline(h=0,lty=2,col="darkorange",lwd=2)
points(all_GWAS_clusters$beta_area.F_NaCl[!all_GWAS_clusters$is_sign_area.F_NaCl],all_GWAS_clusters$GA150[!all_GWAS_clusters$is_sign_area.F_NaCl],pch=16)
points(all_GWAS_clusters$beta_area.F_NaCl[all_GWAS_clusters$is_sign_area.F_NaCl],all_GWAS_clusters$GA150[all_GWAS_clusters$is_sign_area.F_NaCl],pch=16,col="red")


plot(all_GWAS_clusters$beta_area.F_NGM,all_GWAS_clusters$GA150,pch=16,type="n")
abline(h=0,lty=2,col="darkorange",lwd=2)
points(all_GWAS_clusters$beta_area.F_NGM[!all_GWAS_clusters$is_sign_area.F_NGM],all_GWAS_clusters$GA150[!all_GWAS_clusters$is_sign_area.F_NGM],pch=16)
points(all_GWAS_clusters$beta_area.F_NGM[all_GWAS_clusters$is_sign_area.F_NGM],all_GWAS_clusters$GA150[all_GWAS_clusters$is_sign_area.F_NGM],pch=16,col="red")


test=subset(all_GWAS_clusters,is_sign_area.F_NaCl)
test$GA_mean = (test$GA150 + test$GA250 + test$GA450)/3
plot((test$GA_mean),test$beta_area.F_NaCl)
abline(h=0)
abline(v=0)

test=subset(all_GWAS_clusters,is_sign_area.F_NGM)
test$GA_mean = (test$GA150 + test$GA250 + test$GA450)/3
plot((test$GA_mean),test$beta_area.F_NGM)
abline(h=0)
abline(v=0)

sum((test$GA150)*test$beta_area.F_NaCl,na.rm=TRUE)
sum((test$GA250)*test$beta_area.F_NaCl,na.rm=TRUE)
sum((test$GA450)*test$beta_area.F_NaCl,na.rm=TRUE)

sum((test$GA150)*test$beta_area.F_NGM,na.rm=TRUE)
sum((test$GA250)*test$beta_area.F_NGM,na.rm=TRUE)
sum((test$GA450)*test$beta_area.F_NGM,na.rm=TRUE)

aa=0
for(i in 1:10000){
aa=aa+sum(rnorm(length(test$beta_area.F_NaCl),0,0.09)*test$beta_area.F_NaCl,na.rm=TRUE)
}
aa/10000

pdf("Trait_change_prediction_GA1_R02.pdf")
bb=list()
for(k in 1:16){
  print(k)
test=all_GWAS_clusters

bb_temp=NULL
for(i in 1:10000){
  bb_temp=c(bb_temp,sum(sample(test$GA150)*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE))
}
bb[[k]]=bb_temp
hist(bb_temp,main=paste("Trait Nb: ",k))
abline(v=sum(test$GA150*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:16){
  print(k)
  print(sum(test$GA150*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE)<sort(bb[[k]])[length(bb[[k]])*0.05] | sum(test$GA150*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE)>sort(bb[[k]])[length(bb[[k]])*0.95])
}


pdf("Trait_change_prediction_GA2_R02.pdf")
cc=list()
for(k in 1:16){
  test=all_GWAS_clusters
  cc_temp=NULL
  for(i in 1:10000){
    cc_temp=c(cc_temp,sum(sample(test$GA250)*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE))
  }
  cc[[k]]=cc_temp
  hist(cc_temp,main=paste("Trait : ",vect_P_traits_bothEnv[k]))
  abline(v=sum(test$GA250*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:16){
  print(k)
  print(sum(test$GA250*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE)<sort(cc[[k]])[length(cc[[k]])*0.05] | sum(test$GA250*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE)>sort(cc[[k]])[length(cc[[k]])*0.95])
}

########


pdf("Trait_change_prediction_GA4_R02.pdf")
ee=list()
for(k in 1:16){
  test=all_GWAS_clusters
  ee_temp=NULL
  for(i in 1:10000){
    ee_temp=c(ee_temp,sum(sample(test$GA450)*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE))
  }
  ee[[k]]=ee_temp
  hist(ee_temp,main=paste("Trait : ",vect_P_traits_bothEnv[k]))
  abline(v=sum(test$GA450*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE),col='red')
}
dev.off()


for(k in 1:16){
  print(k)
  print(sum(test$GA450*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE)<sort(ee[[k]])[length(ee[[k]])*0.05] | sum(test$GA450*test[,paste0("beta_",vect_P_traits_bothEnv[k])],na.rm=TRUE)>sort(ee[[k]])[length(ee[[k]])*0.95])
}


### Save temp
save(list=ls(),file='temp_R02.RData')

### What explains whether a SNP/cluster increase ? effect (not really) ? pleiotropy ? position ? / neighbouring SNPs ?




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
all_GWAS_sign <- (merge(all_GWAS,all_cluster_merged[,c("rs","group","sign_overall","POS")]))
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
## Only 28 with 2 types
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
head(gp_length)

plot(gp_length$nb_SNPs,gp_length$clust_length,pch=16,cex=.8,col=c("red","cornflowerblue","lightgreen","lightgreen","cornflowerblue")[as.numeric(as.factor(gp_length$chr_type))])

## ClusterID & group are the same thing
all_GWAS_clusters$group=all_GWAS_clusters$clusterID
all_GWAS_clusters=merge(all_GWAS_clusters,gp_length)
head(all_GWAS_clusters)
#plot(all_GWAS_clusters$ps,all_GWAS_clusters$clust_length/1000000,pch=16)

### We will now replace minor af
temp.af.NGM=1-all_GWAS_clusters$af.NGM
temp.af.NaCl=1-all_GWAS_clusters$af.NaCl
#temp.af.NGM[is.na(temp.af.NGM)]=1;temp.af.NaCl[is.na(temp.af.NaCl)]=1
all_GWAS_clusters$af.NGM[!is.na(temp.af.NGM) & all_GWAS_clusters$af.NGM>temp.af.NGM]=temp.af.NGM[!is.na(temp.af.NGM) & all_GWAS_clusters$af.NGM>temp.af.NGM]
all_GWAS_clusters$af.NaCl[!is.na(temp.af.NaCl) & all_GWAS_clusters$af.NaCl>temp.af.NaCl]=temp.af.NaCl[!is.na(temp.af.NaCl) & all_GWAS_clusters$af.NaCl>temp.af.NaCl]

plot((all_GWAS_clusters$af.NGM),(all_GWAS_clusters$af.NaCl),pch=16)


boxplot(all_GWAS_clusters$af.NGM~all_GWAS_clusters$sign_NGM)
boxplot(all_GWAS_clusters$af.NaCl~all_GWAS_clusters$sign_NaCl)

tapply(all_GWAS_clusters$af.NGM,all_GWAS_clusters$sign_NGM,mean)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]

summary(glm(sign_NGM~af.NGM,data=all_GWAS_clusters,family="quasipoisson"))

af.per.pleio <- rbind(
 data.frame(sign=all_GWAS_clusters$sign_NGM,env_label='NGM',af=all_GWAS_clusters$af.NGM,chr_type2=all_GWAS_clusters$chr_type2),
 data.frame(sign=all_GWAS_clusters$sign_NaCl,env_label='NaCl',af=all_GWAS_clusters$af.NaCl,chr_type2=all_GWAS_clusters$chr_type2))

af.per.pleio$sign=as.factor(af.per.pleio$sign)
af.per.pleio$sign_num=as.numeric(as.character(as.factor((af.per.pleio$sign))))

af.per.pleio=na.omit(af.per.pleio)
## only take minor alleles
af.per.pleio$af[af.per.pleio$af>0.5] = 1-af.per.pleio$af[af.per.pleio$af>0.5]

af.per.pleio$af.logit <- log(af.per.pleio$af/(1-af.per.pleio$af))

model.af <- (lm(af.logit~env_label*sign+chr_type2*env_label,data=subset(af.per.pleio,sign!=0)))
summary(model.af)
model.af2 <- (lm(af.logit~env_label*sign+chr_type2,data=subset(af.per.pleio,sign!=0)))
anova(model.af,model.af2)
summary(model.af2)
model.af3 <- (lm(af.logit~env_label+sign+chr_type2*env_label,data=subset(af.per.pleio,sign!=0)))
anova(model.af,model.af3)
model.af4 <- (lm(af.logit~chr_type2*env_label,data=subset(af.per.pleio,sign!=0)))
anova(model.af3,model.af4)
model.af5 <- (lm(af.logit~chr_type2+env_label,data=subset(af.per.pleio,sign!=0)))
anova(model.af5,model.af3)

summary(model.af3)

model.af <- (lm(af.logit~env_label*sign_num,data=subset(af.per.pleio,sign!=0)))
model.af2 <- (lm(af.logit~env_label+sign_num,data=subset(af.per.pleio,sign!=0)))
model.af3 <- (lm(af.logit~sign_num,data=subset(af.per.pleio,sign!=0)))
summary(model.af)
anova(model.af,model.af2)
anova(model.af2,model.af3)
anova(model.af2,model.af4)
hist(af.per.pleio$af*(1-af.per.pleio$af),n=40)

af.per.pleio.NGM <- subset(af.per.pleio,env_label=="NGM")
af.per.pleio.NaCl <- subset(af.per.pleio,env_label=="NaCl")
sd_NGM <- as.numeric(tapply(af.per.pleio.NGM$af,af.per.pleio.NGM$sign,sd))
mean_NGM <- as.numeric(tapply(af.per.pleio.NGM$af,af.per.pleio.NGM$sign,mean))

sd_NaCl <- as.numeric(tapply(af.per.pleio.NaCl$af,af.per.pleio.NaCl$sign,sd))
mean_NaCl <- as.numeric(tapply(af.per.pleio.NaCl$af,af.per.pleio.NaCl$sign,mean))

plot(mean_NGM~c(c(0:7)+0.1),pch=16,ylim=c(0,1),
     las=1,bty="n",ylab="Mean allele frequencies",xlab="Pleiotropy")
#arrows(c(c(0:7)+0.1),mean_NGM-2*sd_NGM,c(c(0:7)+0.1),mean_NGM+2*sd_NGM,code=3,angle=90)

points(mean_NaCl~c(c(0:7)-.1),pch=16,col="orange",bty="n",las=1)
#arrows(c(c(0:7)-0.1),mean_NaCl-2*sd_NaCl,c(c(0:7)-0.1),mean_NaCl+2*sd_NaCl,code=3,angle=90)


hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+clust_length+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic

boxplot(all_GWAS_clusters$clust_length+1~all_GWAS_clusters$sign_overall,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NGM,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NaCl,log="y")

### Add the pleiotropy
dim(merge(all_GWAS_sign,unique(all_GWAS_clusters[,c("group","sign_overall","sign_NaCl","sign_NGM")])))

#pheno_NGM=read.table("Breeding_values_NGM.txt",sep='\t',h=TRUE)
#pheno_NaCl=read.table("Breeding_values_NaCl.txt",sep='\t',h=TRUE)
gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

# Clusters
sign_all_df_clust <- all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)]
table(rowSums(sign_all_df_clust)>=1) # ~5411 Clusters / 1360 are no longer significant
all_GWAS_clusters$sign_overall <- rowSums(sign_all_df_clust)>=1

# SNPs
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:14)]
table(rowSums(sign_all_df)>=1) # 107k SNPs
all_GWAS$sign_overall <- rowSums(sign_all_df)>=1

rowS_NGM_clust <- table(rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:7)]))
rowS_NaCl_clust <-table(rowSums(all_GWAS_clusters[,paste0("is_sign_T",8:14)]))

rowS_NGM <- table(rowSums(all_GWAS[,paste0("is_sign_T",1:7)]))
rowS_NaCl <-table(rowSums(all_GWAS[,paste0("is_sign_T",8:14)]))

### Create a G-matrix from the SNPs ?
var_traits <- NULL
for(i in 1:14){
  if(i %in% 1:7){
    var_traits = c(var_traits, 
                 sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)]],na.rm=TRUE))
  }else{
    var_traits = c(var_traits, 
                   sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af.NaCl*(1-all_GWAS_clusters$af.NaCl))[all_GWAS_clusters[,paste0("is_sign_T",i)]],na.rm=TRUE))
  }
}
plot(var_traits[1:7],var_traits[8:14],pch=16,asp=1)
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)


### If we split center/arm ?
### Create a G-matrix from the SNPs ?
var_traits_center <- NULL; var_traits_arm <- NULL

for(i in 1:14){
  if(i %in% 1:7){
    var_traits_center = c(var_traits_center, 
                   sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE))
    var_traits_arm = c(var_traits_arm, 
                          sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE))
  }else{
    var_traits_center = c(var_traits_center, 
                   sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af.NaCl*(1-all_GWAS_clusters$af.NaCl))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE))
    var_traits_arm = c(var_traits_arm, 
                          sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af.NaCl*(1-all_GWAS_clusters$af.NaCl))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE))
  }
}

plot(var_traits_center,var_traits_arm,pch=16,col=rep(c("blue","red"),each=7))
abline(a=0,b=1); abline(h=0,lty=2); abline(v=0,lty=2) # Much more variance in the arms

####### Covariance between the traits ?

G_mat = matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_arm =matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_center =matrix(rep(NA,196),nrow=14,ncol=14)

for(i in 1:14){
  for(j in 1:14){
    G_mat[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)]],na.rm=TRUE)
    G_mat_center[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE)
    G_mat_arm[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE)
  
  }
}
table(all_GWAS_clusters$chr_type2)
6324/1181
boxplot(all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM)~all_GWAS_clusters$chr_type2)
tapply(all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM),all_GWAS_clusters$chr_type2,mean) # similar frequencies
tapply(abs(all_GWAS_clusters$beta_T7),all_GWAS_clusters$chr_type2,mean) # strong
boxplot(abs(all_GWAS_clusters$beta_T2)~all_GWAS_clusters$chr_type2) # strong

all_GWAS_clusters$varianceT1 <- all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM)*(all_GWAS_clusters$beta_T1)^2
all_GWAS_clusters$varianceT7 <- all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM)*(all_GWAS_clusters$beta_T7)^2

summary(glm((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,family="Gamma"))
summary(glm((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,family="inverse.gaussian"))
boxplot((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,log="y")

a=subset(all_GWAS_clusters,!is.na(all_GWAS_clusters$chr_type2) & varianceT7<0.00051)
plot(density(a$varianceT7[a$chr_type2=='center']),ylim=c(0,40000))
lines(density(a$varianceT7[a$chr_type2=='arm']),col='red')

plot(G_mat_center,G_mat_arm,pch=16)
abline(a=0,b=1); abline(h=0,lty=2); abline(v=0,lty=2) # Much more (co)-variance in the arms

rbind(var_traits[1:7],var_traits[8:14])

plot(cor(all_GWAS_clusters[,paste0("beta_T",1:7)]),
     cor(na.omit(all_GWAS_clusters[,paste0("beta_T",8:14)])))
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)

## Export all volcano plots
pdf("Volcano_plots_clusters_BLUEs.pdf")
vect_col=c(adjustcolor("grey",alpha=0.3),adjustcolor("firebrick3",alpha=1))
for(i in 1:14){
  plot(-log10(all_GWAS_clusters[,paste0("p_wald_T",i)])~all_GWAS_clusters[,paste0("beta_T",i)],
       pch=16,xlab='Beta coefficient',ylab="-log(p)",
       col=vect_col[1+as.numeric(all_GWAS_clusters[,paste0("is_sign_T",i)])],
       main=paste0('Trait: ',vect_P_traits_bothEnv[i]))
}
dev.off()

i=7
plot(-log10(all_GWAS_clusters[,paste0("p_wald_T",i)])~all_GWAS_clusters[,paste0("beta_T",i)],
     pch=16,xlab='Beta coefficient',ylab="-log(p)",
     col=vect_col[1+as.numeric(all_GWAS_clusters[,paste0("is_sign_T",i)])],
     main=paste0('Trait: ',vect_P_traits_bothEnv[i]))
subset(all_GWAS,p_wald_T7<0.0000001)[,1:5]
-log10(0.00000001)

plot((all_GWAS_clusters$beta_T7),(all_GWAS_clusters$beta_T14),pch=16)
plot((all_GWAS_clusters$beta_T7),(all_GWAS_clusters$beta_T14),pch=16,xlim=c(-.1,.1),ylim=c(-.1,.1))


### Heatmap plots between beta/sign ?
heatmap(cor(na.omit(all_GWAS_clusters[,paste0("beta_T",1:14)])),
        labRow=vect_P_traits_bothEnv,
        labCol=vect_P_traits_bothEnv)

## With the non-significant
plot(rowS_NGM/sum(rowS_NGM),type='b',col='cornflowerblue',ylim=c(.0001,1),log="y")
points(rowS_NaCl/sum(rowS_NaCl),type='b',col='darkorange')
points(rowS_NaCl_clust/sum(rowS_NaCl_clust),type='b',col='darkorange',lty=2)
points(rowS_NGM_clust/sum(rowS_NGM_clust),type='b',col='cornflowerblue',lty=2)

## Without the non-significant
pdf('plots_BLUEs/Pleiotropie_NGM_NaCl_BLUEs.pdf')
rowS_NGM2=c(rowS_NGM[2:7],0);rowS_NGM_clust2=c(rowS_NGM_clust[2:7],0)
rowS_NaCl2=rowS_NaCl[2:8];rowS_NaCl_clust2=c(rowS_NaCl_clust[2:7],0)
plot(rowS_NGM2/sum(rowS_NGM2),type='b',col='cornflowerblue',ylim=c(.0001,1),log="y",
     ylab="Number of significant SNPs/Clusters",xlab='Number of traits')
points(rowS_NaCl2/sum(rowS_NaCl2),type='b',col='darkorange')
points(rowS_NaCl_clust2/sum(rowS_NaCl_clust2),type='b',col='darkorange',lty=2)
points(rowS_NGM_clust2/sum(rowS_NGM_clust2),type='b',col='cornflowerblue',lty=2)
dev.off()

pdf('plots_BLUEs/Pleiotropie_NGM_NaCl_absolute.pdf')
plot(rowS_NGM2,type='b',col='cornflowerblue',ylim=c(1,45000),log="y",
     ylab="Number of significant SNPs/Clusters",xlab='Number of traits')
points(rowS_NaCl2,type='b',col='darkorange')
points(rowS_NaCl_clust2,type='b',col='darkorange',lty=2)
points(rowS_NGM_clust2,type='b',col='cornflowerblue',lty=2)
dev.off()

(rowS_NaCl2/sum(rowS_NaCl2))/(rowS_NGM2/sum(rowS_NGM2))
# 12 times more SNPs with all 7 traits in NaCl
sum(rowS_NGM2)/sum(rowS_NaCl2) # Same number of clusters


library(gplots)
library(limma) 
vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_T1),as.numeric(all_GWAS_clusters$is_sign_T8))),
  names=c("SF NGM","SF NaCl"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_T7),as.numeric(all_GWAS_clusters$is_sign_T14))),
  names=c("SF NGM","SF NaCl"),main='All')

vennDiagram(vennCounts(
  all_GWAS_clusters[,paste0("is_sign_T",c(14,11))]),
  names=c("test","Area NGM"),main='All')

### Could there be SNPs in the same genes but not the same variants 
### affecting the traits ? -> would be nice to give examples ????

#### Now the Neff to estimate the effect ~pleiotropy ?
## For each cluster, compute a level of pleiotropy

all_GWAS_clusters_sign <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:14)])>0,])
all_GWAS_clusters_sign_NGM <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:7)])>0,])
all_GWAS_clusters_sign_NaCl <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_T",8:14)])>0,])

nb_sign_clust_NGM <- nrow(all_GWAS_clusters_sign)
nb_sign_clust_NaCl <- nrow(all_GWAS_clusters_sign)
nb_sign_clust <- nrow(all_GWAS_clusters_sign)

pleiot.idx=data.frame(NGM=rep(0,nb_sign_clust),NaCl=rep(0,nb_sign_clust))
pleiot.idx$NGM.cor=0;pleiot.idx$NaCl.cor=0

pleiot.idx_mixed=data.frame(N=rep(0,nb_sign_clust),Neff=rep(0,nb_sign_clust))

#Gmat_NGM = cov(na.omit(all_GWAS_clusters[,paste0("beta_T",1:7)]))
#Gmat_NaCl = cov(na.omit(all_GWAS_clusters[,paste0("beta_T",8:14)]))
#G_mat_all=cov(na.omit(all_GWAS_clusters[,paste0("beta_T",1:14)]))

Gmat_NGM = cov(na.omit(gemma_pheno[,1:7]))
Gmat_NaCl = cov(na.omit(gemma_pheno[,8:14]))
G_mat_all= cov(na.omit(gemma_pheno[,1:14]))


#head(gemma_pheno)

#plot(cov(na.omit(all_GWAS_clusters[,paste0("beta_T",1:7)])),cov(na.omit(gemma_pheno[,1:7])))
#abline(a=0,b=1)
#plot(G_mat[1:7,1:7],cov(na.omit(gemma_pheno[,1:7])))
#abline(a=0,b=1)

for(i in 1:nrow(pleiot.idx_mixed)){
  ## NGM first
  pleiot.idx$NGM[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:7)],na.rm=TRUE)
  if(pleiot.idx$NGM[i]==1){
    pleiot.idx$NGM.cor[i]=1 
  }else if(pleiot.idx$NGM[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:7)])
    std_devs.NGM.temp <- sqrt(diag(Gmat_NGM[idx.sign,idx.sign]))
    cor_matrix.NGM.temp <- Gmat_NGM[idx.sign,idx.sign] / (std_devs.NGM.temp %o% std_devs.NGM.temp)
    pleiot.idx$NGM.cor[i]=pleiot.idx$NGM[i] - var(eigen(cor_matrix.NGM.temp)$values)
  }
  
  ## NaCl second
  pleiot.idx$NaCl[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_T",8:14)],na.rm=TRUE)
  if(pleiot.idx$NaCl[i]==1){
    pleiot.idx$NaCl.cor[i]=1 
  }else if(pleiot.idx$NaCl[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_T",8:14)])
    std_devs.NaCl.temp <- sqrt(diag(Gmat_NaCl[idx.sign,idx.sign]))
    cor_matrix.NaCl.temp <- Gmat_NaCl[idx.sign,idx.sign] / (std_devs.NaCl.temp %o% std_devs.NaCl.temp)
    pleiot.idx$NaCl.cor[i]=pleiot.idx$NaCl[i] - var(eigen(cor_matrix.NaCl.temp)$values)
  }
  
  ### Mixed
  pleiot.idx_mixed$N[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:14)],na.rm=TRUE)
  
  if(pleiot.idx_mixed$N[i]==1){
    pleiot.idx_mixed$Neff[i]=1  
  }else if(pleiot.idx_mixed$N[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:14)])
    idx.sign[is.na(idx.sign)]=FALSE
    std_devs.all.temp <- sqrt(diag(G_mat_all[idx.sign,idx.sign]))
    cor_matrix.all.temp <- G_mat_all[idx.sign,idx.sign] / (std_devs.all.temp %o% std_devs.all.temp)
    pleiot.idx_mixed$Neff[i]=pleiot.idx_mixed$N[i] - var(eigen(cor_matrix.all.temp)$values)
  }
}

par(mfrow=c(1,3))
hist(pleiot.idx$NGM.cor,main='NGM',xlim=c(0,6),ylim=c(0,3000))
hist(pleiot.idx$NaCl.cor,main='NaCl',xlim=c(0,6),ylim=c(0,3000))
hist(pleiot.idx_mixed$Neff,main='ALl 12',xlim=c(0,12),ylim=c(0,5000))

##### Now compute the summed effect of the top SNP of each signicant cluster ? ####

### We need to compute the standardized allelic effect ###
all_effects <- all_GWAS_clusters_sign[,paste0('beta_T',1:14)]
all_effects_std <- NULL
for(i in 1:14){
  all_effects_std[i] =
    sd(all_GWAS_clusters_sign[all_GWAS_clusters_sign[,paste0("is_sign_T",i)],paste0("beta_T",i)])
}

all_GWAS_clusters_sign_std <- all_GWAS_clusters_sign
for(i in 1:14){
  all_GWAS_clusters_sign_std[,paste0('beta_T',i)]=all_GWAS_clusters_sign_std[,paste0('beta_T',i)]/all_effects_std[i]
}

 euclidian_summed <- NULL;euclidian_NGM <- NULL;euclidian_NaCl <- NULL
 manh_summed <- NULL; manh_NGM <- NULL; manh_NaCl <- NULL
 for(i in 1:nrow(all_GWAS_clusters_sign_std)){
   euclidian_summed <- c(euclidian_summed,
                         sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:14)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:14)])]^2)))
   euclidian_NGM[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:7)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:7)])]^2))
   euclidian_NaCl[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',8:14)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',8:14)])]^2))
   
   manh_summed <- c(manh_summed,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:14)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:14)])])))
   manh_NGM=c(manh_NGM,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:7)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:7)])])))
   manh_NaCl=c(manh_NaCl,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',8:14)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',8:14)])])))
   
 }

 
 ### With nls - replace all
 pdf("plots_BLUEs/Scaling_pleiotropy_Neff_BLUE.pdf",h=10,w=12)
 
 is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
 is_center <- all_GWAS_clusters_sign$chr_type2=="center"
 
 temp=pleiot.idx_mixed$Neff
 coefs_euc <- coef(nls(euclidian_summed~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$Neff[is_arm]
 coefs_euc_arm <- coef(nls(euclidian_summed[is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$Neff[is_center]
 coefs_euc_center <- coef(nls(euclidian_summed[is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 temp.NGM <- euclidian_NGM>0
 temp=pleiot.idx$NGM.cor[temp.NGM]
 coefs_euc.NGM <- coef(nls(euclidian_NGM[temp.NGM]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_arm]
 ## a=0 here
 coefs_euc.NGM.arm <- c(0,coef(nls(euclidian_NGM[temp.NGM & is_arm]~b*temp^c,start=list(b=0.8,c=.5))))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_center]
 coefs_euc.NGM.center <- coef(nls(euclidian_NGM[temp.NGM & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 temp.NaCl <- euclidian_NaCl>0
 temp=pleiot.idx$NaCl.cor[temp.NaCl]
 coefs_euc.NaCl <- coef(nls(euclidian_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_euc.NaCl.arm <- coef(nls(euclidian_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_euc.NaCl.center <- coef(nls(euclidian_summed[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 #par(mfrow=c(2,3))
 plot(euclidian_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,11),ylim=c(0,14),type='n')
 points(euclidian_summed[is_arm]~pleiot.idx_mixed$Neff[is_arm],pch=16,col='green')
 points(euclidian_summed[is_center]~pleiot.idx_mixed$Neff[is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc[1]+(coefs_euc[2])*(seq(0,12,0.1)^coefs_euc[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc_arm[1]+(coefs_euc_arm[2])*(seq(0,12,0.1)^coefs_euc_arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc_center[1]+(coefs_euc_center[2])*(seq(0,12,0.1)^coefs_euc_center[3])),col="lightblue")
 
 abline(lm(euclidian_summed~pleiot.idx_mixed$Neff))
 text(2,12,paste0("c=",round(coefs_euc[3],digits=3)))
 text(2,11,paste0("c=",round(coefs_euc_arm[3],digits=3)),col="green")
 text(2,10,paste0("c=",round(coefs_euc_center[3],digits=3)),col='lightblue')
 
 plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NGM',xlim=c(0,7),ylim=c(0,10),
      type='n')
 points(euclidian_NGM[temp.NGM & is_arm]~pleiot.idx$NGM.cor[temp.NGM & is_arm],pch=16,col='green')
 points(euclidian_NGM[temp.NGM & is_center]~pleiot.idx$NGM.cor[temp.NGM & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NGM[1]+(coefs_euc.NGM[2])*(seq(0,12,0.1)^coefs_euc.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NGM.arm[1]+(coefs_euc.NGM.arm[2])*(seq(0,12,0.1)^coefs_euc.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NGM.center[1]+(coefs_euc.NGM.center[2])*(seq(0,12,0.1)^coefs_euc.NGM.center[3])),col="lightblue")
 
 abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]))
 text(1.5,7.5,paste0("c=",round(coefs_euc.NGM[3],digits=3)))
 text(1.5,6.5,paste0("c=",round(coefs_euc.NGM.arm[3],digits=3)),col="green")
 text(1.5,5.5,paste0("c=",round(coefs_euc.NGM.center[3],digits=3)),col="lightblue")
 
 plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,10),type='n')
 points(euclidian_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=16,col='green')
 points(euclidian_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NaCl[1]+(coefs_euc.NaCl[2])*(seq(0,12,0.1)^coefs_euc.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.arm[1]+(coefs_euc.NaCl.arm[2])*(seq(0,12,0.1)^coefs_euc.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.center[1]+(coefs_euc.NaCl.center[2])*(seq(0,12,0.1)^coefs_euc.NaCl.center[3])),col="lightblue")
 
 abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))
 text(1,8.5,paste0("c=",round(coefs_euc.NaCl[3],digits=3)))
 text(1,7.5,paste0("c=",round(coefs_euc.NaCl.arm[3],digits=3)),col="green")
 text(1,6.5,paste0("c=",round(coefs_euc.NaCl.center[3],digits=3)),col='lightblue')
 
 temp=pleiot.idx_mixed$Neff
 coefs_man <- coef(nls(manh_summed~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx_mixed$Neff[is_arm]
 coefs_man.arm <- coef(nls(manh_summed[is_arm]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx_mixed$Neff[is_center]
 coefs_man.center <- coef(nls(manh_summed[is_center]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
 temp.NGM <- manh_NGM>0
 temp=pleiot.idx$NGM.cor[temp.NGM]
 coefs_man.NGM <- coef(nls(manh_NGM[temp.NGM]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_arm]
 coefs_man.NGM.arm <- coef(nls(manh_NGM[temp.NGM & is_arm]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_center]
 coefs_man.NGM.center <- coef(nls(manh_NGM[temp.NGM & is_center]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
 temp.NaCl <- manh_NaCl>0
 temp=pleiot.idx$NaCl.cor[temp.NaCl]
 coefs_man.NaCl <- coef(nls(manh_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 
 plot(manh_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TM ',xlab="Neff",xlim=c(0,11),ylim=c(0,40),type='n')
 points(manh_summed[is_arm]~pleiot.idx_mixed$Neff[is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_summed[is_center]~pleiot.idx_mixed$Neff[is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man[1]+(coefs_man[2])*(seq(0,12,0.1)^coefs_man[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.arm[1]+(coefs_man.arm[2])*(seq(0,12,0.1)^coefs_man.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.center[1]+(coefs_man.center[2])*(seq(0,12,0.1)^coefs_man.center[3])),col="lightblue")
 
 abline(lm(manh_summed~pleiot.idx_mixed$Neff))
 text(2,30,paste0("c=",round(coefs_man[3],digits=3)))
 text(2,28,paste0("c=",round(coefs_man.arm[3],digits=3)),col='green')
 text(2,26,paste0("c=",round(coefs_man.center[3],digits=3)),col='lightblue')
 
 plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(1,7),ylim=c(0,25),type='n')
 points(manh_NGM[temp.NGM & is_arm]~pleiot.idx$NGM.cor[temp.NGM & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NGM[temp.NGM & is_center]~pleiot.idx$NGM.cor[temp.NGM & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NGM[1]+(coefs_man.NGM[2])*(seq(0,12,0.1)^coefs_man.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NGM.arm[1]+(coefs_man.NGM.arm[2])*(seq(0,12,0.1)^coefs_man.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NGM.center[1]+(coefs_man.NGM.center[2])*(seq(0,12,0.1)^coefs_man.NGM.center[3])),col="lightblue")
 
 abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]))
 text(1.5,20,paste0("c=",round(coefs_man.NGM[3],digits=3)))
 text(1.5,18,paste0("c=",round(coefs_man.NGM.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NGM.center[3],digits=3)),col='lightblue')
 
 plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,7),ylim=c(0,25),type='n')
 points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(0,12,0.1)^coefs_man.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(0,12,0.1)^coefs_man.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(0,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue")
 
 abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))
 text(1.5,20,paste0("c=",round(coefs_man.NaCl[3],digits=3)))
 text(1.5,18,paste0("c=",round(coefs_man.NaCl.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NaCl.center[3],digits=3)),col='lightblue')
 
 
 ### With N instead of Neff - do a linear model ?
 
 vect_chr_type <- all_GWAS_clusters_sign_std$chr_type2
 model_pleiotropie <- lm(log(manh_summed)~log(pleiot.idx_mixed$Neff)*vect_chr_type)
 model_pleiotropie2 <- lm(log(manh_summed)~log(pleiot.idx_mixed$Neff)+vect_chr_type)
 anova(model_pleiotropie,model_pleiotropie2)
 
 
 dev.off()
 
 # 
# 
# 
# coefs_euc <- as.numeric((lm(log(euclidian_summed)~log(pleiot.idx_mixed$Neff)))$coef)
# temp.NGM <- euclidian_NGM>0
# coefs_euc.NGM <- as.numeric((lm(log(euclidian_NGM[temp.NGM])~log(pleiot.idx$NGM.cor[temp.NGM])))$coef)
# temp.NaCl <- euclidian_NaCl>0
# coefs_euc.NaCl <- as.numeric((lm(log(euclidian_NaCl[temp.NaCl])~log(pleiot.idx$NaCl.cor[temp.NaCl])))$coef)
# 
# par(mfrow=c(2,3))
# plot(euclidian_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,11),ylim=c(0,14))
# lines(seq(0,12,0.1),I(exp(coefs_euc[1])*(seq(0,12,0.1)^coefs_euc[2])),col="red")
# abline(lm(euclidian_summed~pleiot.idx_mixed$Neff-1))
# text(2,12,paste0("b=",round(coefs_euc[2],digits=3)))
# 
# plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NGM',xlim=c(0,5.5),ylim=c(0,10))
# lines(seq(0,12,0.1),I(exp(coefs_euc.NGM[1])*(seq(0,12,0.1)^coefs_euc.NGM[2])),col="red")
# abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]-1))
# text(1.5,7.5,paste0("b=",round(coefs_euc.NGM[2],digits=3)))
# 
# plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,10))
# lines(seq(0,12,0.1),I(exp(coefs_euc.NaCl[1])*(seq(0,12,0.1)^coefs_euc.NaCl[2])),col="red")
# abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]-1))
# text(1,8.5,paste0("b=",round(coefs_euc.NaCl[2],digits=3)))
# 
# # mod1 <- lm(log(manh_summed)~log(pleiot.idx_mixed$Neff))
# # mod2 <- lm(log(manh_summed)~offset(log(pleiot.idx_mixed$Neff)))
# # summary(mod1)
# # summary(mod2)
# # anova(mod1,mod2)
# # AIC(mod1)
# # AIC(mod2)
# 
# coefs_man <- as.numeric((lm(log(manh_summed)~log(pleiot.idx_mixed$Neff)))$coef)
# temp.NGM <- manh_NGM>0
# coefs_man.NGM <- as.numeric((lm(log(manh_NGM[temp.NGM])~log(pleiot.idx$NGM.cor[temp.NGM])))$coef)
# temp.NaCl <- manh_NaCl>0
# coefs_man.NaCl <- as.numeric((lm(log(manh_NaCl[temp.NaCl])~log(pleiot.idx$NaCl.cor[temp.NaCl])))$coef)
# 
# plot(manh_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TM ',xlab="Neff",xlim=c(0,11),ylim=c(0,40))
# lines(seq(0,12,0.1),I(exp(coefs_man[1])*(seq(0,12,0.1)^coefs_man[2])),col="red")
# abline(lm(manh_summed~pleiot.idx_mixed$Neff-1))
# text(2,30,paste0("b=",round(coefs_man[2],digits=3)))
# 
# plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(0,5.5),ylim=c(0,25))
# lines(seq(0,12,0.1),I(exp(coefs_man.NGM[1])*(seq(0,12,0.1)^coefs_man.NGM[2])),col="red")
# abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]-1))
# text(1.5,20,paste0("b=",round(coefs_man.NGM[2],digits=3)))
# 
# plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,25))
# lines(seq(0,12,0.1),I(exp(coefs_man.NaCl[1])*(seq(0,12,0.1)^coefs_man.NaCl[2])),col="red")
# abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]-1))
# text(1,20,paste0("b=",round(coefs_man.NaCl[2],digits=3)))
# 
# dev.off()
# 
# 
# 
# #################
# 
# ### The same with N ? ##
# pdf("plots/Scaling_pleiotropy_N.pdf",h=14,w=12)
# 
# coefs_euc <- as.numeric((lm(log(euclidian_summed)~log(pleiot.idx_mixed$N)))$coef)
# temp.NGM <- euclidian_NGM>0
# coefs_euc.NGM <- as.numeric((lm(log(euclidian_NGM[temp.NGM])~log(pleiot.idx$NGM[temp.NGM])))$coef)
# temp.NaCl <- euclidian_NaCl>0
# coefs_euc.NaCl <- as.numeric((lm(log(euclidian_NaCl[temp.NaCl])~log(pleiot.idx$NaCl.[temp.NaCl])))$coef)
# 
# par(mfrow=c(2,3))
# plot(euclidian_summed~pleiot.idx_mixed$N,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,13),ylim=c(0,14))
# lines(seq(0,12,0.1),I(exp(coefs_euc[1])*(seq(0,12,0.1)^coefs_euc[2])),col="red")
# abline(lm(euclidian_summed~pleiot.idx_mixed$N-1))
# text(2,12,paste0("b=",round(coefs_euc[2],digits=3)))
# 
# plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NGM',xlim=c(0,6.5),ylim=c(0,10))
# lines(seq(0,12,0.1),I(exp(coefs_euc.NGM[1])*(seq(0,12,0.1)^coefs_euc.NGM[2])),col="red")
# abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM]-1))
# text(1.5,7.5,paste0("b=",round(coefs_euc.NGM[2],digits=3)))
# 
# plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,6.5),ylim=c(0,10))
# lines(seq(0,12,0.1),I(exp(coefs_euc.NaCl[1])*(seq(0,12,0.1)^coefs_euc.NaCl[2])),col="red")
# abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]-1))
# text(1,8.5,paste0("b=",round(coefs_euc.NaCl[2],digits=3)))
# 
# coefs_man <- as.numeric((lm(log(manh_summed)~log(pleiot.idx_mixed$N)))$coef)
# temp.NGM <- manh_NGM>0
# coefs_man.NGM <- as.numeric((lm(log(manh_NGM[temp.NGM])~log(pleiot.idx$NGM[temp.NGM])))$coef)
# temp.NaCl <- manh_NaCl>0
# coefs_man.NaCl <- as.numeric((lm(log(manh_NaCl[temp.NaCl])~log(pleiot.idx$NaCl[temp.NaCl])))$coef)
# 
# plot(manh_summed~pleiot.idx_mixed$N,pch=16,ylab='TM ',xlab="Neff",xlim=c(0,12),ylim=c(0,40))
# lines(seq(0,12,0.1),I(exp(coefs_man[1])*(seq(0,12,0.1)^coefs_man[2])),col="red")
# abline(lm(manh_summed~pleiot.idx_mixed$N-1))
# text(2,30,paste0("b=",round(coefs_man[2],digits=3)))
# 
# plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(0,6.5),ylim=c(0,25))
# lines(seq(0,12,0.1),I(exp(coefs_man.NGM[1])*(seq(0,12,0.1)^coefs_man.NGM[2])),col="red")
# abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.[temp.NGM]-1))
# text(1.5,20,paste0("b=",round(coefs_man.NGM[2],digits=3)))
# 
# plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,6.5),ylim=c(0,25))
# lines(seq(0,12,0.1),I(exp(coefs_man.NaCl[1])*(seq(0,12,0.1)^coefs_man.NaCl[2])),col="red")
# abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]-1))
# text(1,20,paste0("b=",round(coefs_man.NaCl[2],digits=3)))
# 
# dev.off()


########################################################################
 ########################################################################
 ########################################################################
 
 ### Create a G-matrix from the SNPs - depending on pleiotropy ?
 sign_total <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)])
G_mat = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center =matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_P1 = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm_P1 =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center_P1 =matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_P2 = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm_P2 =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center_P2 =matrix(rep(NA,196),nrow=14,ncol=14) 
G_mat_P3 = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm_P3 =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center_P3 =matrix(rep(NA,196),nrow=14,ncol=14) 
vect_P1 <- sign_total==1 & all_GWAS_clusters$clust_length<10000
vect_P2 <- sign_total%in%c(2,3) & all_GWAS_clusters$clust_length<10000
vect_P3 <- sign_total>3 & all_GWAS_clusters$clust_length<10000


 for(i in 1:14){
   for(j in 1:14){
     G_mat[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE)
     G_mat_arm[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE)
     

     G_mat_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_T",i)]*all_GWAS_clusters[vect_P1,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P1]*(1-all_GWAS_clusters$af.NGM[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P1,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_T",i)]*all_GWAS_clusters[vect_P1,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P1]*(1-all_GWAS_clusters$af.NGM[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P1,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P1]=="center"],na.rm=TRUE)
     G_mat_arm_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_T",i)]*all_GWAS_clusters[vect_P1,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P1]*(1-all_GWAS_clusters$af.NGM[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P1,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P1]=="arm"],na.rm=TRUE)
     
     G_mat_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_T",i)]*all_GWAS_clusters[vect_P2,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P2]*(1-all_GWAS_clusters$af.NGM[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P2,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_T",i)]*all_GWAS_clusters[vect_P2,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P2]*(1-all_GWAS_clusters$af.NGM[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P2,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P2]=="center"],na.rm=TRUE)
     G_mat_arm_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_T",i)]*all_GWAS_clusters[vect_P2,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P2]*(1-all_GWAS_clusters$af.NGM[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P2,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P2]=="arm"],na.rm=TRUE)
     
     G_mat_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_T",i)]*all_GWAS_clusters[vect_P3,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P3]*(1-all_GWAS_clusters$af.NGM[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P3,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_T",i)]*all_GWAS_clusters[vect_P3,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P3]*(1-all_GWAS_clusters$af.NGM[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P3,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P3]=="center"],na.rm=TRUE)
     G_mat_arm_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_T",i)]*all_GWAS_clusters[vect_P3,paste0("beta_T",j)]*all_GWAS_clusters$af.NGM[vect_P3]*(1-all_GWAS_clusters$af.NGM[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P3,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P3]=="arm"],na.rm=TRUE)
     
        }
 }
 
 
plot(scale(eigen(G_mat)$vectors[,1]),scale(eigen(G_mat_P1)$vectors[,1]))

 plot(G_mat[1:7,1:7],G_mat_P1[1:7,1:7])
 
 allPcor=NULL
 
 for(i in 1:13){
   for(j in (i+1):14){
     allPcor=rbind(allPcor,
     data.frame(T1=i,T2=j,Pcor=G_mat_all[i,j],nbClust=sum(rowSums(all_GWAS_clusters_sign[,paste0("is_sign_T",c(i,j))])==2),
                nbClust_arm=sum(rowSums(all_GWAS_clusters_sign[all_GWAS_clusters_sign$chr_type2=="arm",paste0("is_sign_T",c(i,j))])==2),
                nbClust_center=sum(rowSums(all_GWAS_clusters_sign[all_GWAS_clusters_sign$chr_type2=="center",paste0("is_sign_T",c(i,j))])==2))
     )
   }
 }
 allPcor$same_env = 0
 allPcor$same_env[allPcor$T1%in%(1:7) & allPcor$T2%in%(1:7)]=1
 allPcor$same_env[allPcor$T1%in%(8:14) & allPcor$T2%in%(1:14)]=2
 
 allPcor$same_T = 0
 allPcor$same_T[allPcor$T1==(allPcor$T2-7)]=1
 
 plot(abs(allPcor$Pcor),allPcor$nbClust,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 plot(abs(allPcor$Pcor),allPcor$nbClust_arm,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 plot(abs(allPcor$Pcor),allPcor$nbClust_center,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 
 model_sameT <- lm(nbClust~abs(Pcor),data=subset(allPcor,same_T==1))
 summary(model_sameT)
 
 allPcor$same_env=as.factor(allPcor$same_env)
 allPcor$same_T=as.factor(allPcor$same_T)
 
 allPcor$same_env_NGM <- allPcor$same_env==1
 allPcor$same_env_NaCl <- allPcor$same_env==2
 
 model <- lm(allPcor$nbClust~abs(allPcor$Pcor)*(allPcor$same_T+allPcor$same_env_NGM+allPcor$same_env_NaCl))
 summary(model)
 drop1(model)
 
 model <- lm(allPcor$nbClust~abs(allPcor$Pcor)*as.factor(allPcor$same_T)+allPcor$same_env_NGM+allPcor$same_env_NaCl)
 summary(model)
 drop1(model,test='F')
 
 model2 <- lm(allPcor$nbClust~abs(allPcor$Pcor)*as.factor(allPcor$same_T)+allPcor$same_env_NaCl)
 summary(model2)
 drop1(model2,test='F')
 anova(model,model2,test="F")
 
 model3 <- lm(allPcor$nbClust~abs(allPcor$Pcor)*as.factor(allPcor$same_T))
 summary(model3)
 drop1(model3,test='F')
 
 
allPcor$same_env=as.numeric(allPcor$same_env)-1
 plot(abs(allPcor$Pcor),allPcor$nbClust,pch=c(16,11)[as.numeric(allPcor$same_T)],
      col=c("black","blue","orange")[as.numeric(allPcor$same_env)],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 abline(a=coef(model3)[1]+coef(model3)[3],b=coef(model3)[2]+coef(model3)[4])
 abline(a=coef(model3)[1],b=coef(model3)[2],col="blue")
 abline(a=(coef(model2)[1]+coef(model2)[4]),b=coef(model2)[2],col="orange")
 
 model2_arm <- lm(allPcor$nbClust~abs(allPcor$Pcor)+as.factor(allPcor$same_env))
 drop1(model2_arm,test='F')
 model2_center <- lm(allPcor$nbClust~abs(allPcor$Pcor)+as.factor(allPcor$same_env))
 drop1(model2_center,test='F')
 
 
 
 