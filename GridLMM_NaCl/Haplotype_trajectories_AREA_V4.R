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

### SNPs on CHROM X

CHR="X"
snp_block=subset(selected,CHROM==CHR)$POS # 1 block
snp_all=all.af.caendr[row.names(all.af.caendr)%in%paste0(CHR,".",snp_block),]
dim(snp_all) # 278 SNPs - 2 are missing from one block

plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  tmp=populations$replicate_fact%in%c(1,i)
  lines(snp_all[1,tmp]~populations$generations[tmp])
}

plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in 3:6){
  for(p in 1:278){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}

## Clearly multiple groups of frequencies /// This block should be broken I assume

load('RData/Genotype_info_Cemee.RData')
founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]

## Import the more complete Cemee Data
snp_positions=as.data.frame(snp_positions)
snp_positions$POS=as.numeric(snp_positions$POS)

snp_positions$CHROM_num=as.numeric(factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA")))
table(snp_positions$CHROM_num)

### Add the position info on the genotype matrix - simpler later
row.names(genotype_info)=snp_positions$ID
chr_peak = 6

all_BP_peak <- snp_block
length(all_BP_peak) # 280 SNPs
genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')

## Extend the region / take all intermediate SNPs ?

max(all_BP_peak)-min(all_BP_peak) # 825kb
length(all_BP_peak) # 280 SNPs

#PART 1: investigate haplotypes
#1. select your region - let's take 10kb on each side
start <- min(all_BP_peak)-100000 #put here start bp, organizing as a list so later we can put all interesting regions and loop below analyis
end <- max(all_BP_peak)+100000
region <- c(start, end)

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 2695 SNPs


library(data.table)
plot(as.numeric(tstrsplit(rownames(genotype_info_region),'_')[[2]]))
abline(h=min(all_BP_peak),col='red');abline(h=max(all_BP_peak),col='red')


#how many haplotypes -might want to remove founders before?
unique_haplotypes <- unique(t(genotype_info_region)) # will this work? a unique combination of 0 and 1?
## -> Yes !
dim(unique_haplotypes)
# Here 890 combinaison of these SNPs
heatmap(unique_haplotypes,Colv=NA,scale='none')
# All
#heatmap(t(genotype_info_region),Colv=NA,scale="none")

### Unique A6140
kept_CeMEE = substring(colnames(genotype_info_region),1,2)%in%c("A6")
unique_haplotypes <- unique(t(genotype_info_region[,kept_CeMEE])) # will this work? a unique combination of 0 and 1?
heatmap(unique_haplotypes,Colv=NA,scale='none')

#####?##################################################################################################################################
### Ok, we see some are the founders haplotypes (N2 of course but not only)
founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]
founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]

dim(unique(founders_peak)) # 94 different SNPs only
dim(unique(genotype_info_region)) # 1325 unique SNPs in Cemee (recombination! / or sequencing error, we must be carefull)

dim(unique(t(founders_peak))) # 15 different haplotypes, can we differentiate groups ?
heatmap(t(founders_peak),Colv=NA,scale='none')
# At least 6 main groups of haplotype. 
founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(founders_peak),Colv=NA,scale='none') # Issue with polarization ?
#heatmap(unique(t(founders_peak)),Colv=NA,scale='none') # Issue with polarization ?
### How to proceed ?

## MY1 haplotype ?
dim(founders_peak[rowSums(founders_peak)==6 & founders_peak[,4]==2 & founders_peak[,6]==2 & founders_peak[,10]==2,]) # 29 SNPs are private to this first group
SNPs_MY1 = rownames(founders_peak[rowSums(founders_peak)==6 & founders_peak[,4]==2 & founders_peak[,6]==2 & founders_peak[,10]==2,])
## CB4856 haplotype ?
dim(founders_peak[rowSums(founders_peak)==2 & founders_peak[,5]==2,]) # 27 SNPs are private to CB4856
SNPs_CB4856 = rownames(founders_peak[rowSums(founders_peak)==2 & founders_peak[,5]==2,])
## JU319/MY16
dim(founders_peak[rowSums(founders_peak)==4 & founders_peak[,7]==2 & founders_peak[,11]==2 ,]) # 0 SNPs are private to this group
## PB306 / JU345
dim(founders_peak[rowSums(founders_peak)==2 & founders_peak[,8]==2,]) # 0 SNPs are private to JU345
dim(founders_peak[rowSums(founders_peak)==2 & founders_peak[,13]==2,]) # 0 SNPs are private to PB306
# AB1 & CB4852
dim(founders_peak[rowSums(founders_peak)==4 & founders_peak[,3]==2 & founders_peak[,3]==2,]) # 5 SNPs are private to this first group
SNPs_AB1 = rownames(founders_peak[rowSums(founders_peak)==4 & founders_peak[,3]==2 & founders_peak[,3]==2,])
### We have three genotypes now ... (four alleles + undetermined)

founders_haplo_MY1 <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_MY1,'_')[[2]],]
heatmap(t(founders_haplo_MY1),Colv=NA,scale='none')

founders_haplo_CB4856 <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_CB4856,'_')[[2]],]
heatmap(t(founders_haplo_CB4856),Colv=NA,scale='none')

founders_haplo_AB1 <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_AB1,'_')[[2]],]
heatmap(t(founders_haplo_AB1),Colv=NA,scale='none')

genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_MY1,'_')[[2]],]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')

all_snps = c(SNPs_MY1,SNPs_CB4856,SNPs_AB1)

snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("X.",tstrsplit(all_snps,"_")[[2]]))
par(mfrow=c(1,3))
for(i in c(3,4,6)){
plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
  for(p in 1:nrow(snp_all_sub)){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}

#### By founders ?
par(mfrow=c(1,3))
snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("X.",tstrsplit(SNPs_MY1,"_")[[2]]))
for(i in c(3,4,6)){
  plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
  for(p in 1:nrow(snp_all_sub)){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}
snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("X.",tstrsplit(SNPs_CB4856,"_")[[2]]))
for(i in c(3,4,6)){
  plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
  for(p in 1:nrow(snp_all_sub)){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}
snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("X.",tstrsplit(SNPs_AB1,"_")[[2]]))
for(i in c(3,4,6)){
  plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
  for(p in 1:nrow(snp_all_sub)){
    tmp=populations$replicate_fact%in%c(1,i)
    lines(snp_all_sub[p,tmp]~populations$generations[tmp],col=c("black","black","red","blue","orange","darkgreen")[i])
  }
}

mapp_env=new.env()
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData',envir=mapp_env)
## MY1 specific SNPs
g1 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%SNPs_MY1]
is_G1 = (rowSums(g1)==58)
table(is_G1) # 283
## CB4856
g2 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%SNPs_CB4856]
is_G2 = (rowSums(g2)==54)
table(is_G2) # 57
## AB1
g3 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%SNPs_AB1]
is_G3 = (rowSums(g3)==10)
table(is_G3) # 236

## N2 - none of the other ones ?
is_G4 = (rowSums(g1)==0 & rowSums(g2)==0 & rowSums(g3)==0)
table(is_G4) # 222

phen = mapp_env$all_BLUEs_genotyped_NaCl
phen=subset(phen,pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1 | is_G2 | is_G3 | is_G4])
length(unique(phen$pop_label)) # 737

phen$isG1 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1]
phen$isG2 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G2]
phen$isG3 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G3]
phen$isG4 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G4]

 library(gplots)
 library(limma) 
 vennDiagram(vennCounts(
   cbind(as.numeric(phen$isG1),
         as.numeric(phen$isG2),
         as.numeric(phen$isG3))),
   names=c("G1","G2","G3"),main='All')

new_geno = mapp_env$genotype_matrix[,1:3]
colnames(new_geno)=c("G1","G2","G3")
new_geno[,1]=is_G1;new_geno[,2]=(is_G2);new_geno[,3]=(is_G3)

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
mod1 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2+isG3)^2 ,data=phen)
drop1(mod1,test='Chisq')
mod2 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2+isG3) ,data=phen)
drop1(mod2,test='Chisq')
summary(mod2)

mod2a = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1, data=subset(phen,!isG2 & !isG3))
drop1(mod2a,test='Chisq')
mod2b = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG2, data=subset(phen,!isG1 & !isG3))
drop1(mod2b,test='Chisq')
mod2c = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG3, data=subset(phen,!isG2 & !isG1))
drop1(mod2c,test='Chisq')

mod3=lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1*isG3, data=subset(phen,!isG2))
drop1(mod3,test='Chisq')
mod3b=lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1+isG3, data=subset(phen,!isG2))
drop1(mod3b,test='Chisq')

mod3.area=lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2+isG3)^2 ,data=phen)
drop1(mod3.area,test='Chisq')

mod3.area_test=lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1+isG2*isG3 ,data=phen)
drop1(mod3.area_test,test='Chisq')
summary(mod3.area_test)

mod3.area2=lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2+isG3) ,data=phen)
drop1(mod3.area2,test='Chisq')

mod_pred1 = lmer(beta ~ (-1) + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(beta))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2","isG3","isG4")]))
pred_beta=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$beta = pred_beta
## Remove weird ones
#new_df2=subset(new_df, !(isG2 & !isG1))
new_df2=new_df
new_df2$geno = 'mixed'
new_df2$geno[new_df2$isG1 & !new_df2$isG2 & !new_df2$isG3] = 'MY1'
new_df2$geno[new_df2$isG1 & !new_df2$isG2 & new_df2$isG3] = 'MY1 & AB1'
new_df2$geno[new_df2$isG1 & new_df2$isG2 & !new_df2$isG3] = 'MY1 & CB4856'
new_df2$geno[!new_df2$isG1 & new_df2$isG2 & new_df2$isG3] = 'CB4856 & AB1'
new_df2$geno[new_df2$isG2 & !new_df2$isG1 & !new_df2$isG3] = 'CB4856'
new_df2$geno[new_df2$isG2 & new_df2$isG1 & new_df2$isG3] = 'MY1 & AB1 & CB4856'
new_df2$geno[new_df2$isG3 & !new_df2$isG1 & !new_df2$isG2] = 'AB1'
new_df2$geno[new_df2$isG4] = 'N2'

boxplot(new_df2$beta~new_df2$geno,ylab="Beta phenotype",xlab="")

##### And for body area  - or others
mod2 = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1, data=phen)
drop1(mod2,test='Chisq')
summary(mod2)

mod_pred1 = lmer(area.F ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2","isG3","isG4")]))
pred_area=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$area = pred_area
## Remove weird ones

new_df2=new_df
new_df2$geno = 'mixed'
new_df2$geno[new_df2$isG1 & !new_df2$isG2 & !new_df2$isG3] = 'MY1'
new_df2$geno[new_df2$isG1 & !new_df2$isG2 & new_df2$isG3] = 'MY1 & AB1'
new_df2$geno[new_df2$isG1 & new_df2$isG2 & !new_df2$isG3] = 'MY1 & CB4856'
new_df2$geno[!new_df2$isG1 & new_df2$isG2 & new_df2$isG3] = 'CB4856 & AB1'
new_df2$geno[new_df2$isG2 & !new_df2$isG1 & !new_df2$isG3] = 'CB4856'
new_df2$geno[new_df2$isG2 & new_df2$isG1 & new_df2$isG3] = 'MY1 & AB1 & CB4856'
new_df2$geno[new_df2$isG3 & !new_df2$isG1 & !new_df2$isG2] = 'AB1'
new_df2$geno[new_df2$isG4] = 'N2'

boxplot(new_df2$area~new_df2$geno,ylab="Area phenotype",xlab="")

vect_traits=c("SF","SB","FS","FB","BS","BF")

par(mfrow=c(2,3))
for(i in 1:6){
phen$y=phen[,vect_traits[i]]
mod2 = lmer(y~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + isG1, data=phen)
print(drop1(mod2,test='Chisq'))
mod_pred1 = lmer(y ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1")]))
new_df$pred_i=predict(mod_pred1,newdata=new_df,re.form=NA)
names(new_df)[ncol(new_df)]=vect_traits[i]
## Remove weird ones

new_df2=new_df
new_df2$geno = 'non-MY1'
new_df2$geno[new_df2$isG1] = 'MY1'
boxplot(new_df2[,vect_traits[i]]~new_df2$geno,ylab="TR phenotype",xlab="")
}






