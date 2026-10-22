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

selected = read.table(file='Selected_clusters.txt',h=TRUE,sep='\t')
table(selected$CHROM)

### SNPs on CHROM I

CHR="V"
snp_block=subset(selected,CHROM==CHR)$POS
snp_all=all.af.caendr[row.names(all.af.caendr)%in%paste0(CHR,".",snp_block),]
dim(snp_all) # 15 SNPs - all here

plot(snp_all[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in c(3,4,6)){
  for(p in 1:9){
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
chr_peak = 1
all_BP_peak <- snp_block

length(all_BP_peak)
genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')

## Extend the region / take all intermediate SNPs ?

max(all_BP_peak)-min(all_BP_peak) # 1,232 MB
length(all_BP_peak) # 15 SNPs only are rising?

#PART 1: investigate haplotypes
#1. select your region - let's take 50kb on each side
start <- min(all_BP_peak)-1000000 #put here start bp, organizing as a list so later we can put all interesting regions and loop below analyis
end <- max(all_BP_peak)+40000
region <- c(start, end)

genotype_info_region <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
snp_positions_region <- subset(snp_positions,CHROM==CHR & POS>=start & POS<=end)
dim(genotype_info_region);dim(snp_positions_region) # 1882 SNPs

library(data.table)
plot(as.numeric(tstrsplit(rownames(genotype_info_region),'_')[[2]]))
abline(h=min(all_BP_peak),col='red');abline(h=max(all_BP_peak),col='red')

#how many haplotypes -might want to remove founders before?
unique_haplotypes <- unique(t(genotype_info_region)) # will this work? a unique combination of 0 and 1?
dim(unique_haplotypes)
# Here 835 combinaison of these SNPs
heatmap(unique_haplotypes,Colv=NA,scale='none')

### Unique A6140
kept_CeMEE = substring(colnames(genotype_info_region),1,2)%in%c("A6")
unique_haplotypes <- unique(t(genotype_info_region[,kept_CeMEE])) # will this work? a unique combination of 0 and 1?
dim(unique_haplotypes)
heatmap(unique_haplotypes,Colv=NA,scale='none')

unique_haplotypes <- unique(t((genotype_info_region))) # will this work? a unique combination of 0 and 1?
dim(unique_haplotypes)
heatmap(unique_haplotypes,Colv=NA,scale='none')

row_selected = (rowSums(unique_haplotypes[,paste0("V_",c(6322080,6637695,6826421))])==6)
heatmap(unique_haplotypes[row_selected,],Colv=NA,scale='none')

row_selected = (rowSums(unique_haplotypes[,paste0("V_",c(6322080,6637695,6826421))])==6)| !substring(rownames(unique_haplotypes),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")
heatmap(unique_haplotypes[row_selected,],Colv=NA,scale='none')
table(row_selected)

row_selected = (rowSums(unique_haplotypes[,paste0("V_",c(6637695,6826421))])==4)# | !substring(rownames(unique_haplotypes),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")
heatmap(unique(unique_haplotypes[row_selected,]),Colv=NA,scale='none')

row_selected = rowSums((unique_haplotypes)[,paste0("V_",all_BP_peak[c(1:4,6,7,10:15)])])==0
heatmap(unique(unique_haplotypes[row_selected,]),Colv=NA,scale='none')

row_selected=row_selected & (rowSums(unique_haplotypes[,paste0("V_",c(6322080,6637695,6826421))])==6)
heatmap(unique(unique_haplotypes[row_selected,]),Colv=NA,scale='none')
table(row_selected)

row_selected = (rowSums(t(genotype_info_region_exact)[,paste0("V_",c(6322080,6637695,6826421))])==6)
heatmap(t(genotype_info_region_exact)[row_selected,],Colv=NA,scale='none')
table(row_selected) # 41 but only 19 have not the other ones

heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')

#####?##################################################################################################################################
### Ok, we see some are the founders haplotypes (N2 of course but not only)
founders_haplo <- genotype_info[,!substring(colnames(genotype_info),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")]
founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS>=start & snp_positions$POS<=end,]
dim(unique(founders_peak)) # 25 different SNPs only
heatmap(t(founders_peak),Colv=NA,scale='none')
# Three main groups of haplotype. Is one of them selected ?
founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%all_BP_peak,]
heatmap(t(founders_peak),Colv=NA,scale='none')
## One SNP is just missing !
row_selected = (rowSums(unique_haplotypes[,paste0("V_",c(6637695,6826421))])==4)# | !substring(rownames(unique_haplotypes),1,2)%in%c("A6","CA","GA","GT","GM","LR","SM","EE","GF")
heatmap((unique_haplotypes[row_selected,]),Colv=NA,scale='none')

### Unique to CB4507 ?
dim(founders_peak[rowSums(founders_peak)<=4 & founders_peak[,2]==2,]) # 18 SNPs only
SNPs_all_three = paste0("V_",c(6322080,6637695,6826421))
SNPs_other_missing = paste0("V_",all_BP_peak[c(1:4,6,7,10:15)])

snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("V.",tstrsplit(SNPs_all_three,"_")[[2]]))
plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in c(3,4,6)){
  for(p in 1:3){
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

### Other ones
snp_all_sub=subset(snp_all,row.names(snp_all)%in%paste0("V.",tstrsplit(SNPs_other_missing,"_")[[2]]))
plot(snp_all_sub[1,]~populations$generations,type='n',ylim=c(0,1))
for(i in c(3,4,6)){
  for(p in 1:length(SNPs_other_missing)){
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

#####

founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_all_three,"[_]")[[2]],]
heatmap(t(founders_peak),Colv=NA,scale='none')

founders_peak <- founders_haplo[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(SNPs_other_missing,"[_]")[[2]],]
heatmap(t(founders_peak),Colv=NA,scale='none')

### Evidence for new mutation ? weird but OK
row_selected = rowSums((unique_haplotypes)[,paste0("V_",all_BP_peak[c(1:4,6,7,10:15)])])==0 | rownames(unique_haplotypes)%in%c("N2","CB4507")
row_selected=row_selected & (rowSums(unique_haplotypes[,paste0("V_",c(6322080,6637695,6826421))])>=6)
heatmap((unique_haplotypes[row_selected,]),Colv=NA,scale='none')

row_selected = rowSums((unique_haplotypes)[,paste0("V_",all_BP_peak[c(1:4,6,7,10:15)])])==0 | rownames(unique_haplotypes)%in%c("N2","CB4507")
row_selected=row_selected & (rowSums(unique_haplotypes[,paste0("V_",c(6322080,6637695,6826421))])>=6)
heatmap((unique_haplotypes[row_selected,]),Colv=NA,scale='none')


#genotype_info_region_exact <- genotype_info[snp_positions$CHROM==CHR & snp_positions$POS%in%tstrsplit(c(SNPs_all_three,SNPs_other_missing),"[_]")[[2]],]
#dim(genotype_info_region_exact)
heatmap(t(genotype_info_region_exact),Colv=NA,scale='none')
row_selected = genotype_info_region_exact[paste0("V_",6322080),]==2
heatmap(t(genotype_info_region_exact[,row_selected]),Colv=NA,scale='none')
############# OK then we can re-genotype according to these 3  ?

mapp_env=new.env()
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData',envir=mapp_env)

## All three ?
g1 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%paste0("V_",tstrsplit(SNPs_all_three,"[_]")[[2]])]
is_G1 = (rowSums(g1)==6)
table(is_G1) # 40

## Only the new mutation ?
is_G1b = g1[,1]==2
table(is_G1b) # 64

## All three AND not the other ones
g2 = mapp_env$genotype_matrix[,colnames(mapp_env$genotype_matrix)%in%paste0("V_",tstrsplit(SNPs_other_missing,"[_]")[[2]])]
is_G2 = (rowSums(g2)==0)
table(is_G2) # 38

is_G3 = (rowSums(g2)==24)
sum(is_G3)+sum(is_G2 | is_G1 | is_G1b) # 749 out of 960

phen = mapp_env$all_BLUEs_genotyped_NaCl
phen=subset(phen,pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1 | is_G2 | is_G3 | is_G3b])
length(unique(phen$pop_label))

phen$isG1 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1]
phen$isG1b = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G1b]
phen$isG2 = phen$pop_label%in%rownames(mapp_env$genotype_matrix)[is_G2]

new_geno = mapp_env$genotype_matrix[,1:3]
colnames(new_geno)=c("G1","G1b","G2")
new_geno[,1]=is_G1;new_geno[,2]=(is_G1b);new_geno[,3]=(is_G2)

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

gwas.NaCl.SF = GridLMM_GWAS(formula = SF ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                            test_formula = ~1, reduced_formula = ~0,
                            data = phen, X = new_geno,
                            X_ID = 'pop_label', relmat = list(pop_label = mapp_env$GRM_sub), mc.cores = 1 ,method = 'REML')

gwas.NaCl.SB = GridLMM_GWAS(formula = SB ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                            test_formula = ~1, reduced_formula = ~0,
                            data = phen, X = new_geno,
                            X_ID = 'pop_label', relmat = list(pop_label = mapp_env$GRM_sub), mc.cores = 1 ,method = 'REML')


gwas.NaCl.area$results
gwas.NaCl.beta$results

### Regress ?
library(lme4)
mod1 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2)^2 ,data=phen)
drop1(mod1,test='Chisq')
mod2 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) +(isG1+isG2), data=phen)
anova(mod1,mod2) #sign !

##### If I do it with G1b
mod1b = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1b+isG2)^2 ,data=phen)
drop1(mod1b,test='Chisq')
mod2b = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) +(isG1b+isG2), data=phen)
anova(mod1b,mod2b) #sign !

# And for Area?
mod1.a = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1*isG2) ,data=phen)
mod2.a = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) +(isG1+isG2), data=phen)
anova(mod1.a,mod2.a) #sign !

mod1.ab = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1b*isG2) ,data=phen)
mod2.ab = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) +(isG1b+isG2), data=phen)
anova(mod1.ab,mod2.ab) #not sign again


mod3.a = lmer(area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) +(isG1), data=phen)
anova(mod2.a,mod3.a) # not sign 
mod4 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id), data=phen)
anova(mod3,mod4) # not sign ?
drop1(mod2,test='Chisq')

mod_pred1 = lmer(beta ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(beta))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG1b","isG2")]))
pred_beta=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$beta = pred_beta
## Remove weird ones
#new_df2=subset(new_df, !(isG2 & !isG1))
new_df2=new_df
new_df2$geno = '1.None (876)'
new_df2$geno[new_df2$isG1b] = '2.New mut.(24)'
new_df2$geno[new_df2$isG1] = '3.Three mut.(21)'
new_df2$geno[new_df2$isG1 & new_df2$isG2] = '4.Full cluster (19)'
boxplot(new_df2$beta~new_df2$geno,ylab="Beta phenotype",xlab="")
abline(h=mean(subset(new_df2,geno=='1.None (876)')$beta),col="red")

##### And for body area  - or others?
mod_pred1 = lmer(area.F ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2","isG1b")]))
pred_area=predict(mod_pred1,newdata=new_df,re.form=NA)
new_df$area = pred_area
## Remove weird ones
#new_df2=subset(new_df, !(isG2 & !isG1))
new_df2=new_df
new_df2$geno = '1.None (876)'
new_df2$geno[new_df2$isG1b] = '2.New mut.(24)'
new_df2$geno[new_df2$isG1] = '3.Three mut.(21)'
new_df2$geno[new_df2$isG1 & new_df2$isG2] = '4.Full cluster (19)'
boxplot(new_df2$area~new_df2$geno,ylab="Area phenotype",xlab="")
abline(h=mean(subset(new_df2,geno=='1.None (876)')$area),col="red")

vect_traits=c("SF","SB","FS","FB","BS","BF")
par(mfrow=c(2,3))
for(i in 1:6){
  phen$y=phen[,vect_traits[i]]
  mod2 = lmer(y~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id)+isG1+isG2, data=phen)
  print(drop1(mod2,test='Chisq'))
  mod_pred1 = lmer(y ~ -1 + logD + rel_humidity + temperature + (1|date_str) + pop_label + (1|session_id) ,data=phen)
  new_df=data.frame(pop_label=sort(unique(subset(phen,!is.na(area.F))$pop_label)),temperature=0,rel_humidity=0,logD=0)
  new_df=merge(new_df,unique(phen[,c("pop_label","isG1","isG2","isG1b")]))
  new_df$pred_i=predict(mod_pred1,newdata=new_df,re.form=NA)
  names(new_df)[ncol(new_df)]=vect_traits[i]
  ## Remove weird ones
  
  new_df2=new_df
  new_df2$geno = '1.None (876)'
  new_df2$geno[new_df2$isG1b] = '2.New mut.(24)'
  new_df2$geno[new_df2$isG1] = '3.Three mut.(21)'
  new_df2$geno[new_df2$isG1 & new_df2$isG2] = '4.Full cluster (19)' 
  boxplot(new_df2[,vect_traits[i]]~new_df2$geno,ylab="TR phenotype",xlab="")
  abline(h=mean(subset(new_df2,geno=='1.None (876)')[,vect_traits[i]]),col="red")
}

rm(mapp_env);rm(poolseq_env)
rm(all.af.caendr);rm(all.counts.alt.caendr);rm(all.cov.caendr);rm(all.counts.caendr)
rm(genotype_info);rm(genotype_info_region);rm(GRM)
rm(founders_haplo);rm(unique_haplotypes)
gc()


### STOPPED HERE


## Epistasis

#test.mod = GridLMM(formula = beta ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1*isG2),
#                             test_formula = ~1, reduced_formula = ~0,
#                             data = phen,
#                             X_ID = 'pop_label', relmat = list(pop_label = GRM), mc.cores = 1 ,method = 'REML')

#X = genotype_matrix_CHR[,is.NA==0]

#library(sommer)

#fit <- mmer(
#  fixed  = beta ~ 1 + phen$logD + phen$rel_humidity + phen$temperature + phen$isG1*phen$isG2,
#  random = ~ vs(pop_label, Gu = GRM) + session_id + date_str,
#  rcov   = ~ units,
#  data   = phen)

#fit.no.int <- mmer(
#  fixed = beta ~ logD + rel_humidity + temperature +
#    isG1 + isG2,
#  random = ~ vs(pop_label, Gu = GRM) + session_id + date_str,
#  rcov = ~ units,
#  data = phen
#)

#anova(fit,fit.no.int)

#mod1 = lmer(beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1*isG2) ,data=phen)

#prior <- list(G = list(G1 = list(V = 1, nu = 0.002),G2 = list(V = 1, nu = 0.002),G3 = list(V = 1, nu = 0.002)),
#  R = list(V = 1, nu = 0.002))

##TOO SLOW
#inv.mat = as(solve(GRM),"dgCMatrix")
#fit.mcmc <- MCMCglmm(beta ~ isG1 * isG2 + logD + rel_humidity+ temperature -1,
#  random = ~ pop_label + date_str + session_id,
#  ginverse = list(pop_label = inv.mat),
#  data = phen,prior = prior,nitt = 1000,burnin = 300,thin = 10,family="gaussian")


#all(levels(factor(phen$pop_label)) %in% rownames(GRM))
phen$pop_label <- factor(phen$pop_label)
library(lme4qtl)

## Beta - evidence of epistatic effects
fit.lmer <- relmatLmer(
  beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1*isG2) ,
  data=subset(phen,!is.na(beta)),relmat = list(pop_label = GRM))

fit.lmer2 <- relmatLmer(
  beta~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2) ,
  data=subset(phen,!is.na(beta)),relmat = list(pop_label = GRM))
anova(fit.lmer,fit.lmer2,test='Chisq')

## Area - evidence of epistatic effects
fit.lmer.area <- relmatLmer(
  area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1*isG2) ,
  data=subset(phen,!is.na(beta)),relmat = list(pop_label = GRM))

fit.lmer.area2 <- relmatLmer(
  area.F~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id) + (isG1+isG2) ,
  data=subset(phen,!is.na(beta)),relmat = list(pop_label = GRM))
anova(fit.lmer.area,fit.lmer.area2,test='Chisq')

pred_G = predict(fit.lmer.area,data.frame(logD=0,rel_humidity=0,temperature=0,isG1=c(FALSE,FALSE,TRUE,TRUE),isG2=c(FALSE,TRUE,TRUE,FALSE)),re.form=NA,se.fit=TRUE)

coef(fit.lmer.area)
save(list=ls(),file='Haplotype_trajV3.RData')
plot(1:4,pred_G$fit,pch=16,ylim=c(-.2,.5),las=1,bty="n",xlab="",ylab="Area")
arrows(1:4,pred_G$fit-2*pred_G$se.fit,1:4,pred_G$fit+2*pred_G$se.fit,code=3,angle=90)






