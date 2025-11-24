rm(list=ls())
library(data.table)
library(poolSeq)
library(qqman)
library(gplots)
library(limma) 

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


### Now, we will run the FET and set-up a significance threshold
### chisq GA1 / Manhattan plot
all.counts.noNa <- na.omit(all.counts)
all.counts.alt.noNa <- na.omit(all.counts.noNa)

A00 <- floor(rowMeans(all.counts[,c(2,4)],na.rm=TRUE))
a00 <- floor(rowMeans(all.counts.alt[,c(2,4)],na.rm=TRUE))

A_GA136 <- all.counts[,5];a_GA136 <- all.counts.alt[,5]
A_GA236 <- all.counts[,6];a_GA236 <- all.counts.alt[,6]
A_GA336 <- all.counts[,7];a_GA336 <- all.counts.alt[,7]
A_GA436 <- all.counts[,8];a_GA436 <- all.counts.alt[,8]

GA136_cmh <- chi.sq.test(A00,a00,A_GA136,a_GA136 ,min.cov=30,max.cov=300)
GA236_cmh <- chi.sq.test(A00,a00,A_GA236,a_GA236 ,min.cov=30,max.cov=300)
GA336_cmh <- chi.sq.test(A00,a00,A_GA336,a_GA336 ,min.cov=30,max.cov=300)
GA436_cmh <- chi.sq.test(A00,a00,A_GA436,a_GA436 ,min.cov=30,max.cov=300)

CHI_GA136 = data.frame(CHR=as.numeric(factor(snps_positions.sync$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync$POS),SNP=c(1:nrow(snps_positions.sync)),P=-log10(GA136_cmh),Padj = -log10(p.adjust(GA136_cmh,method ='BH')))
CHI_GA236 = data.frame(CHR=as.numeric(factor(snps_positions.sync$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync$POS),SNP=c(1:nrow(snps_positions.sync)),P=-log10(GA236_cmh),Padj = -log10(p.adjust(GA236_cmh,method ='BH')))
CHI_GA336 = data.frame(CHR=as.numeric(factor(snps_positions.sync$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync$POS),SNP=c(1:nrow(snps_positions.sync)),P=-log10(GA336_cmh),Padj = -log10(p.adjust(GA336_cmh,method ='BH')))
CHI_GA436 = data.frame(CHR=as.numeric(factor(snps_positions.sync$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=as.numeric(snps_positions.sync$POS),SNP=c(1:nrow(snps_positions.sync)),P=-log10(GA436_cmh),Padj = -log10(p.adjust(GA436_cmh,method ='BH')))

manhattan(CHI_GA136, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,highlight = subset(CHI_GA136,P>1.3)$SNP)
manhattan(CHI_GA136, chr="CHR", bp="BP", snp="SNP", p="Padj" ,logp=FALSE,highlight = subset(CHI_GA136,Padj>1.3)$SNP)

table(CHI_GA136$P>1.3) # 29k SNPs
table(CHI_GA236$P>1.3) # 27k SNPs
table(CHI_GA336$P>1.3) # 37k SNPs
table(CHI_GA436$P>1.3) # 29k SNPs

#######
vennDiagram(vennCounts(
  cbind(as.numeric(CHI_GA136$P>1.3),
        as.numeric(CHI_GA236$P>1.3),
        as.numeric(CHI_GA336$P>1.3),
        as.numeric(CHI_GA436$P>1.3))),names=paste0("GA",1:4),main='All')





