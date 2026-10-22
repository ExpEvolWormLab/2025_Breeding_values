rm(list=ls())
gc()

load('RData/Genotype_info_Cemee.RData')

founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]

## Import the more complete Cemee Data
snp_positions=as.data.frame(snp_positions)
snp_positions$POS=as.numeric(snp_positions$POS)

snp_positions$CHROM_num=as.numeric(factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA")))
table(snp_positions$CHROM_num)

### Add the position info on the genotype matrix - simpler later
row.names(genotype_info)=snp_positions$ID

all_BP_peak <- c(6637695, 6322080,6826421)
chr_peak = 5
length(all_BP_peak)

##### From the trajectories I think we first keep these SNPS:

max(all_BP_peak)-min(all_BP_peak) #500 k
length(all_BP_peak) # 3 SNPs

#PART 1: investigate haplotypes
#1. select your region - let's take 50kb on each side
start <- min(all_BP_peak)-50000 #put here start bp, organizing as a list so later we can put all interesting regions and loop below analyis
end <- max(all_BP_peak)+50000
CHR <-  c('V') #put here chromosome names
i = 1 #the region we want to start with. now we just have 1. R counts from 1.  
#easier:
region <- c(start[i], end[i])

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR[i] & snp_positions$POS>=start[i] & snp_positions$POS<=end[i],]
snp_positions_region <- subset(snp_positions,CHROM==CHR[i] & POS>=start[i] & POS<=end[i])

dim(genotype_info_region);dim(snp_positions_region)
## For region 1, 1k SNPs!

#how many haplotypes -might want to remove founders before?
unique_haplotypes <- unique(t(genotype_info_region)) # will this work? a unique combination of 0 and 1?
## -> Yes !
dim(unique_haplotypes)
# Here 255 combinaison of these SNPs
heatmap(unique_haplotypes,Colv=NA,scale='none')

### Ok, we see some are the founders haplotypes (N2 of course but not only)
founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]
founders_peak <- founders_haplo[snp_positions$CHROM==CHR[i] & snp_positions$POS>=start[i] & snp_positions$POS<=end[i],]

dim(unique(founders_peak)) # 11 different SNPs only
dim(unique(genotype_info_region)) # 266 unique SNPs in Cemee (recombination! / or sequencing error, we must be carefull)

dim(unique(t(founders_peak))) # 11 different haplotypes, can we differentiate groups ?
heatmap(t(founders_peak),Colv=NA,scale='none')
# Four main groups of haplotype. Is one of them selected ? -> Yes see after, the first one with 5 strains
heatmap(t(founders_peak),Colv=NA,scale='none')

### OK unique to CB4856 ###
dim(subset(founders_peak,founders_peak[,5]==2 & rowSums(founders_peak)==2)) # 386 SNPs
### OK unique to MY16 - None ###
dim(subset(founders_peak,founders_peak[,11]==2 & rowSums(founders_peak)==2)) # 386 SNPs


founders_peak <- founders_haplo[snp_positions$CHROM==CHR[i] & snp_positions$POS%in%all_BP_peak,]
founders_peak

