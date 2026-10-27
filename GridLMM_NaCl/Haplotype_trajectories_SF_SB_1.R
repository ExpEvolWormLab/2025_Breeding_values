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
### SNPs on CHROM IV

CHR="IV"
snp_block=8892663
snp_all=all.af.caendr[row.names(all.af.caendr)%in%paste0(CHR,".",snp_block),]
length(snp_all) # OK

plot(snp_all~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  tmp=populations$replicate_fact%in%c(1,i)
  lines(snp_all[tmp]~populations$generations[tmp])
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

## First block only
all_BP_peak <- snp_block
length(all_BP_peak)
genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
#heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')
table(genotype_info_region_exact)
## Extend the region / take all intermediate SNPs ?

#max(all_BP_peak)-min(all_BP_peak) # Only 1 SNP
#length(all_BP_peak)

#PART 1: investigate haplotypes
#1. select your region - let's take 10kb on each side
start <- min(all_BP_peak)-100000 #put here start bp, organizing as a list so later we can put all interesting regions and loop below analyis
end <- max(all_BP_peak)+100000
region <- c(start, end)

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 156 SNPs

library(data.table)
plot(as.numeric(tstrsplit(rownames(genotype_info_region),'_')[[2]]))
abline(h=min(all_BP_peak),col='red');abline(h=max(all_BP_peak),col='red')

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 856 SNPs left - only half

#how many haplotypes -might want to remove founders before?
unique_haplotypes <- unique(t(genotype_info_region)) # will this work? a unique combination of 0 and 1?
## -> Yes !
dim(unique_haplotypes)
# Here 56 combinaison of these SNPs
heatmap(unique_haplotypes,Colv=NA,scale='none')
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

dim(unique(founders_peak)) # 6 different SNPs only
dim(unique(genotype_info_region)) # 72 unique SNPs in Cemee (recombination! / or sequencing error, we must be carefull)

dim(unique(t(founders_peak))) # 5 different haplotypes, can we differentiate groups ?
heatmap(t(founders_peak),Colv=NA,scale='none')

### IN which background do we found the mutation ?
ID_geno <- genotype_info_region[row.names(genotype_info_region)=="IV_8892663",]==2
geno_with_mut <- genotype_info_region[,ID_geno]
heatmap(t(geno_with_mut),Colv=NA,scale='none')
heatmap(t(cbind(founders_peak,geno_with_mut)[,1:50]),Colv=NA,Rowv=NA,scale='none')



snp_all_sub=snp_all
plot(snp_all_sub~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all[tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  
}
par(mfrow=c(1,1))
plot(snp_all_sub~populations$generations,type='n',ylim=c(0,1))
for(i in c(3,4,6)){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
}

### The two new mutations ?

CHR="IV"
snp_block=c(8892663,10499135)
snp_all=all.af.caendr[row.names(all.af.caendr)%in%paste0(CHR,".",snp_block),]
dim(snp_all) # OK
all_BP_peak=snp_block
genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')
table(genotype_info_region_exact)
## Extend the region / take all intermediate SNPs ?

start <- 8892663-500000 #put here start bp, organizing as a list so later we can put all interesting regions and loop below analyis
end <- 10499135+500000
region <- c(start, end)

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 1754 SNPs

plot(as.numeric(tstrsplit(rownames(genotype_info_region),'_')[[2]]))
abline(h=8892663,col='red');abline(h=10499135,col='red')

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 1754 SNPs left - only half

#how many haplotypes -might want to remove founders before?
unique_haplotypes <- unique(t(genotype_info_region)) # will this work? a unique combination of 0 and 1?
## -> Yes !
dim(unique_haplotypes)
# Here 565 combinaison of these SNPs
heatmap(unique_haplotypes,Colv=NA,scale='none')
# All
heatmap(t(genotype_info_region),Colv=NA,scale="none")





### STOPPED HERE



## They don't seem very correlated.. a group close to 0 at G10
SNPs_low_start = rownames(snp_all_sub[rowSums(snp_all_sub[,c(3,6)])==0,])
SNPs_high_start = rownames(snp_all_sub[rowSums(snp_all_sub[,c(3,6)])>0,])
all_SNPs =SNPs_high_start = rownames(snp_all_sub)

founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_low_start,"[.]")[[2]],]
heatmap(t(founders_peak),Colv=NA,scale='none')

founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_high_start,"[.]")[[2]],]
heatmap(t(founders_peak),Colv=NA,scale='none')


founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(all_SNPs,"[.]")[[2]],]
heatmap(t(founders_peak),Colv=NA,scale='none')

genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(all_SNPs,"[.]")[[2]],]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')
# almost no recombination - small region though

mapp_env=new.env()
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData',envir=mapp_env)
## MY16 specific SNPs
g1 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%paste0("II_",tstrsplit(all_SNPs,"[.]")[[2]])]
is_G1 = (rowSums(g1)==40)
table(is_G1) # 46
## all others
g2 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%paste0("II_",tstrsplit(SNPs_high_start,"[.]")[[2]])]
is_G2 = (rowSums(g2)==0)
table(is_G2) # 870

phen = mapp_env$all_BLUEs_genotyped_NaCl
phen=subset(phen,pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1 | is_G2])
length(unique(phen$pop_label)) # 916

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
#mod1 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1*isG2) ,data=phen)
mod2 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1, data=phen)
drop1(mod2,test='Chisq')
summary(mod2)

mod_pred1 = lmer(beta ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(beta))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2")]))
pred_beta=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$beta = pred_beta
## Remove weird ones
#new_df2=subset(new_df, !(isG2 & !isG1))
new_df2=new_df
new_df2$geno = 'non-MY16'
new_df2$geno[new_df2$isG1] = 'MY16'
boxplot(new_df2$beta~new_df2$geno,ylab="Beta phenotype",xlab="")


##### And for body area  - or others
mod2 = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1, data=phen)
drop1(mod2,test='Chisq')
summary(mod2)

mod_pred1 = lmer(area.F ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1")]))
pred_area=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$area = pred_area
## Remove weird ones

new_df2=new_df
new_df2$geno = 'non-MY16'
new_df2$geno[new_df2$isG1] = 'MY16'
boxplot(new_df2$area~new_df2$geno,ylab="Beta phenotype",xlab="")

vect_traits=c("SF","SB","FS","FB","BS","BF")

par(mfrow=c(2,3))
for(i in 1:6){
phen$y=phen[,vect_traits[i]]
mod_pred1 = lmer(y ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1")]))
new_df$pred_i=predict(mod_pred1,newdata=new_df,re.form=NA)
names(new_df)[ncol(new_df)]=vect_traits[i]
## Remove weird ones

new_df2=new_df
new_df2$geno = 'non-MY16'
new_df2$geno[new_df2$isG1] = 'MY16'
boxplot(new_df2[,vect_traits[i]]~new_df2$geno,ylab="Beta phenotype",xlab="")
}






