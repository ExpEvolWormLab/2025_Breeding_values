rm(list=ls());gc()

# load the poolseq data now
poolseq_env=new.env()
load(file='SaveGLM.RData',envir=poolseq_env)
populations=poolseq_env$populations
all.af.caendr=poolseq_env$all.af
all.cov.caendr=poolseq_env$all.cov
all.counts.caendr=poolseq_env$all.counts
all.counts.alt.caendr=poolseq_env$all.counts.alt
snps_positions.sync.caendr=poolseq_env$snps_positions.sync

### All select SNPS in interesting blocks

selected = read.table(file='Area.F_clusters.txt',h=TRUE,sep='\t')
table(selected$CHROM)

### SNPs on CHROM II

CHR="II"
snp_block=subset(selected,CHROM==CHR)$POS # 2 blocks
snp_all=all.af.caendr[row.names(all.af.caendr)%in%paste0(CHR,".",snp_block),]
dim(snp_all) # 36 SNPs - 5 are missing from one block---> All 6 from the second group are here

plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  tmp=populations$replicate_fact%in%c(1,i)
  lines(snp_all[1,tmp]~populations$generations[tmp])
}

plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  for(p in 1:36){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}

load('RData/Genotype_info_Cemee.RData')

founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]

## Import the more complete Cemee Data
snp_positions=as.data.frame(snp_positions)
snp_positions$POS=as.numeric(snp_positions$POS)

snp_positions$CHROM_num=as.numeric(factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA")))
table(snp_positions$CHROM_num)

### Add the position info on the genotype matrix - simpler later
row.names(genotype_info)=snp_positions$ID
chr_peak = 2

## Second block now
all_BP_peak <- snp_block[36:41]
length(all_BP_peak)
genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')

## Extend the region / take all intermediate SNPs ?

max(all_BP_peak)-min(all_BP_peak) # 80kb
length(all_BP_peak) # 65 SNPs

#PART 1: investigate haplotypes
#1. select your region - let's take 10kb on each side
start <- min(all_BP_peak)-10000 #put here start bp, organizing as a list so later we can put all interesting regions and loop below analyis
end <- max(all_BP_peak)+10000
region <- c(start, end)

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 349 SNPs
## For region 2, 349 SNPs only in a small region 

library(data.table)
plot(as.numeric(tstrsplit(rownames(genotype_info_region),'_')[[2]]))
abline(h=min(all_BP_peak),col='red');abline(h=max(all_BP_peak),col='red')

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 349 SNPs left - only half


#how many haplotypes -might want to remove founders before?
unique_haplotypes <- unique(t(genotype_info_region)) # will this work? a unique combination of 0 and 1?
## -> Yes !
dim(unique_haplotypes)
# Here 433 combinaison of these SNPs
heatmap(unique_haplotypes,Colv=NA,scale='none')  ### Highly recombined !
# All
heatmap(t(genotype_info_region),Colv=NA,scale="none")

### Unique A6140
kept_CeMEE = substring(colnames(genotype_info_region),1,2)%in%c("A6")
unique_haplotypes <- unique(t(genotype_info_region[,kept_CeMEE])) # will this work? a unique combination of 0 and 1?
heatmap(unique_haplotypes,Colv=NA,scale='none')

#####?##################################################################################################################################
### Ok, we see some are the founders haplotypes (N2 of course but not only)
founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]
founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]

dim(unique(founders_peak)) # 39 different SNPs only
dim(unique(genotype_info_region)) # 260 unique SNPs in Cemee (recombination! / or sequencing error, we must be carefull)

dim(unique(t(founders_peak))) # 166 different haplotypes, can we differentiate groups ?
heatmap(t(founders_peak),Colv=NA,scale='none')
# At least 5 main groups of haplotype !
founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(founders_peak),Colv=NA,scale='none')

### Unique to JU319
dim(founders_peak[rowSums(founders_peak)<=2 & founders_peak[,7]==2,]) # 4 SNPs only are private - others are shared with AB1
SNPs_JU319 = rownames(founders_peak[rowSums(founders_peak)<=2 & founders_peak[,7]==2,])

founders_haplo_JU319 <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_JU319,'_')[[2]],]
heatmap(t(founders_haplo_JU319),Colv=NA,scale='none')

genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_JU319,'_')[[2]],]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')

snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("II.",tstrsplit(SNPs_JU319,"_")[[2]]))
plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  for(p in 1:4){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}
par(mfrow=c(1,3))
for(i in c(3,4,6)){
plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
  for(p in 1:nrow(snp_all_sub)){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}

mapp_env=new.env()
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData',envir=mapp_env)
## JU319 specific SNPs
g1 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%SNPs_JU319]
is_G1 = (rowSums(g1)==8)
table(is_G1) # 42
## all others
is_G2 = (rowSums(g1)==0)
table(is_G2) # 885

phen = mapp_env$all_BLUEs_genotyped_NaCl
phen=subset(phen,pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1 | is_G2])
length(unique(phen$pop_label)) # 927

phen$isG1 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1]
phen$isG2 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G2]

new_geno = mapp_env$genotype_matrix[,1:2]
colnames(new_geno)=c("G1","G2")
new_geno[,1]=is_G1;new_geno[,2]=(is_G2 | is_G1)

library(AGHmatrix)
GRM= Gmatrix(mapp_env$genotype_pruned_matrix, method="VanRaden") 

library(GridLMM)
gwas.NaCl.area = GridLMM_GWAS(formula = area.F ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                             test_formula = ~1, reduced_formula = ~0,
                             data = phen, X = new_geno,
                             X_ID = 'pop_label', relmat = list(pop_label = mapp_env$GRM_sub), mc.cores = 1 ,method = 'REML')

gwas.NaCl.beta = GridLMM_GWAS(formula = beta ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                              test_formula = ~1, reduced_formula = ~0,
                              data = phen, X = new_geno,
                              X_ID = 'pop_label', relmat = list(pop_label = mapp_env$GRM_sub), mc.cores = 1 ,method = 'REML')

gwas.NaCl.area$results
gwas.NaCl.beta$results


### Regress ?
library(lme4)
mod1 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1 ,data=phen)
drop1(mod1,test='Chisq')

mod1.area = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1 ,data=phen)
drop1(mod1.area,test='Chisq')

mod_pred1 = lmer(beta ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(beta))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2")]))
pred_beta=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$beta = pred_beta
## Remove weird ones
#new_df2=subset(new_df, !(isG2 & !isG1))
new_df2=new_df
new_df2$geno = 'N2etc'
new_df2$geno[new_df2$isG1] = 'JU319'
boxplot(new_df2$beta~new_df2$geno,ylab="Beta phenotype",xlab="")


##### And for body area  - or others?
mod_pred1 = lmer(area.F ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2")]))
pred_area=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$area = pred_area
## Remove weird ones
#new_df2=subset(new_df, !(isG2 & !isG1))
new_df2=new_df
new_df2$geno = 'N2etc'
new_df2$geno[new_df2$isG2 & !new_df2$isG1] = 'JU319'
boxplot(new_df2$area~new_df2$geno,ylab="Area phenotype",xlab="")










