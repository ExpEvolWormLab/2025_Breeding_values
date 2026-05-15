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
snps_positions.sync.caendr=snps_positions.sync[kept.alleles,]

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

res <- mclapply(1:nrow(all.counts.caendr), function(i) {
  
  f_All <- cbind(all.counts.caendr[i, vect_All],
                 all.counts.alt.caendr[i, vect_All]) ~ generations + (1|replicate_fact)
  
  f_GA1 <- cbind(all.counts.caendr[i, vect_GA1],
                 all.counts.alt.caendr[i, vect_GA1]) ~ generations + (1|replicate_fact)
  
  f_GA2 <- cbind(all.counts.caendr[i, vect_GA2],
                 all.counts.alt.caendr[i, vect_GA2]) ~ generations + (1|replicate_fact)
  
  f_GA4 <- cbind(all.counts.caendr[i, vect_GA4],
                 all.counts.alt.caendr[i, vect_GA4]) ~ generations + (1|replicate_fact)
  
  list(
    All = safe_glmer(f_All, populations),
    GA1 = safe_glmer(f_GA1, populations_GA1),
    GA2 = safe_glmer(f_GA2, populations_GA2),
    GA4 = safe_glmer(f_GA4, populations_GA4)
  )
  
}, mc.cores = detectCores() - 2)

save(list=ls(),file='SaveGLM.RData')
