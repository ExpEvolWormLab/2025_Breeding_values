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

## Output a file for the haplotype pipeline
pools = data.frame(sample_id=populations$population,sync_file="GA_populations_polymorphisms_v2.sync",
                   sync_column=4:(nrow(populations)+3),generation=populations$generations,replicate=populations$replicate)
write.table(pools,file="pools.tsv",sep='\t',row.names=FALSE,quote=FALSE)

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

#plot(density(all.cov[,3]),xlim=c(0,200),ylim=c(0,0.1))
#for(i in 2:13){
#  lines(density(all.cov[,i]),xlim=c(0,200))
#}


## PCA between samples ?

#pca1 <- prcomp(t(na.omit(all.af)),center=TRUE,scale=FALSE)
#plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
#text(pca1$x[,1],pca1$x[,2],labels=populations$population,col=rep(c("blue","red"),c(2,11))) ### Ok, we should only use the two others then

# Remove A6140
#pca1 <- prcomp(t(na.omit(all.af[,c(1,3:13)])),center=TRUE,scale=FALSE)
#plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
#text(pca1$x[,1],pca1$x[,2],labels=populations$population[c(1,3:13)],col=rep(c("blue","red"),c(1,11))) ### Ok, we should only use the two others then

#pca1 <- prcomp(t(na.omit(all.af[,c(1,3:13)])),center=TRUE,scale=FALSE)
#plot(pca1$x[,1],pca1$x[,2],col=rep(c("blue","red"),each=4),pch=16,type='n')
#text(pca1$x[,1],pca1$x[,2],labels=populations$population[c(1,3:13)],col=rep(c("blue","red","gold","green","magenta"),c(1,3,3,2,3))) ### Ok, we should only use the two others then


### Strange variation in coverage ?
rel_cov=all.cov
for(i in 1:ncol(all.cov)){
  rel_cov[,i]=rel_cov[,i]/mean(rel_cov[,i],na.rm=TRUE)
}

snps_positions.sync$minCov <- rowMins(rel_cov)
snps_positions.sync$maxCov <- rowMaxs(rel_cov)
snps_positions.sync$medianCov <- rowMedians(rel_cov)

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
snps_positions.sync.caendr=subset(snps_positions.sync,ID%in%rownames(all.counts.alt.caendr))
load('SaveGLM.RData')

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

pdf("GLM_plot.pdf",h=5)
manhattan(subset(GLM_All,P>2 & is.finite(P) & kept_base_freq), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,suggestiveline = FALSE,genomewideline = FALSE)
dev.off()

### Superpose the clusters ?
load('GridLMM_results_NaCl/All_Gwas.RData')
for(i in 1:8){
  all_GWAS[,paste0("is_sign_T",i)] = p.adjust(all_GWAS[,paste0("pvalue_T",i)],method="BH")<0.2
  print(table(all_GWAS[,paste0("is_sign_T",i)]))
}
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:8)]
all_GWAS$sign_overall <- rowSums(sign_all_df)>=1


sign_clusters = subset(all_GWAS,is_sign_T7)
sign_clusters=sign_clusters[,c("chr","ps")]
names(sign_clusters)=c("CHR","BP")
sign_clusters$CHR=as.numeric(factor(sign_clusters$CHR,levels=c("I","II","III","IV","V","X","MtDNA")))

GLM_All$clust=FALSE
GLM_All$clust[paste(GLM_All$CHR,GLM_All$BP)%in%paste(sign_clusters$CHR,sign_clusters$BP)]=TRUE

pdf("GLM_plot_withClust_AreaNaCl.pdf",h=5)
manhattan(subset(GLM_All,P>1 & is.finite(P) & kept_base_freq),
          chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,suggestiveline = FALSE,genomewideline = FALSE,
          highlight = subset(GLM_All,P>1 & is.finite(P) & kept_base_freq & clust)$SNP)
dev.off()


sign_clusters = subset(all_GWAS,is_sign_T8)
sign_clusters=sign_clusters[,c("chr","ps")]
names(sign_clusters)=c("CHR","BP")
sign_clusters$CHR=as.numeric(factor(sign_clusters$CHR,levels=c("I","II","III","IV","V","X","MtDNA")))

GLM_All$clust=FALSE
GLM_All$clust[paste(GLM_All$CHR,GLM_All$BP)%in%paste(sign_clusters$CHR,sign_clusters$BP)]=TRUE

pdf("GLM_plot_withClust_BetaNaCl.pdf",h=5)
manhattan(subset(GLM_All,P>1 & is.finite(P) & kept_base_freq),
          chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,suggestiveline = FALSE,genomewideline = FALSE,
          highlight = subset(GLM_All,P>1 & is.finite(P) & kept_base_freq & clust)$SNP)
dev.off()




