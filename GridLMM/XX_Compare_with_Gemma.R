setwd("~/Documents/GitHub/2025_Breeding_values/")
library(data.table)
library(qqman)
FS_NGM=read.table("GridLMM_results_mixed/gwas_results_area.F_NGM.txt",h=TRUE)

names(FS_NGM)[2]="rs"
traitJ = read.table(paste0("gemma_results_BLUE/gwas_results_BLUE_T7.assoc.txt"),h=TRUE)
head(traitJ)

## Merge
merged=merge(traitJ,FS_NGM[,c("rs","beta.1","beta.2","beta.5","p_value_REML")])

plot(merged$beta,merged$beta.1,pch=16)
plot(merged$beta,merged$beta.5,pch=16)

plot(-log10(p.adjust(merged$p_wald,method="BH")),-log10(p.adjust(merged$p_value_REML,method="BH")),pch=16)
plot(-log10((merged$p_wald)),-log10(p.adjust(merged$p_value_REML,method="BH")),pch=16)
abline(h=-log10(0.3),col="red")

plot(merged$beta.1,-log10(p.adjust(merged$p_value_REML,method="BH")),pch=16)
abline(h=-log10(0.2),col="red")
abline(h=-log10(0.3),col="red")

plot(merged$beta,merged$beta.2,pch=16)
cor(merged[,c("beta",paste0("beta.",1:5))])

### Manhattan plot
merged$chr = tstrsplit(merged$rs,"_")[[1]]
merged$ps = as.numeric(tstrsplit(merged$rs,"_")[[2]])

data_manhattan = data.frame(CHR=as.numeric(factor(merged$chr,levels=c("I","II","III","IV","V","X","MtDNA"))),
                            BP=merged$ps,SNP=c(1:nrow(merged)),
                            P=-log10(merged$p_value_REML))
data_manhattan=subset(data_manhattan,CHR!=7)

manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          ,cex=.6,suggestiveline=FALSE,genomewideline=FALSE)







SB_NGM=read.table("GridLMM_results_mixed/gwas_results_area.F_NGM.txt",h=TRUE)
names(SB_NGM)[2]="rs"
traitJ = read.table(paste0("gemma_results_BLUE/gwas_results_BLUE_T2.assoc.txt"),h=TRUE)
merged=merge(traitJ,SB_NGM[,c("rs","beta.1","p_value_REML")])
plot(merged$beta,merged$beta.1,pch=16)

plot(-log10(merged$p_wald),-log10(merged$p_value_REML),pch=16)

plot(-log10(merged$p_wald),-log10(p.adjust(merged$p_value_REML,method="BH")),pch=16)
abline(h=-log10(0.05))

test1=p.adjust(merged$p_value_REML,method="BH")<0.05
test2=(merged$p_wald)<0.05

table(test1,test2)



beta_NGM=read.table("GridLMM_results_mixed/gwas_results_beta_NGM.txt",h=TRUE)
beta_NGM=subset(beta_NGM,!substring(beta_NGM$X_ID,1,1)=="M")
## Concatenate the LOO
beta_LOO_NGM=NULL
for(i in c("I","II","III","IV","V","X")){
  beta_LOO_NGM=rbind(beta_LOO_NGM,read.table(paste0("GridLMM_results_mixed_perCHR/gwas_results_beta_CHR_",i,"_NGM.txt"),h=TRUE))
  
}
dim(beta_LOO_NGM)
dim(beta_NGM)
table(beta_LOO_NGM$X_ID==beta_NGM$X_ID)
plot(beta_LOO_NGM$beta.1,beta_NGM$beta.1)
plot(-log10(beta_LOO_NGM$p_value_REML),-log10(beta_NGM$p_value_REML))
abline(a=0,b=1)
table(beta_LOO_NGM$p_value_REML<0.05)
table(beta_NGM$p_value_REML<0.05)

table(p.adjust(beta_LOO_NGM$p_value_REML,method="BH")<0.05)
table(p.adjust(beta_NGM$p_value_REML,method="BH")<0.05)

table(p.adjust(beta_LOO_NGM$p_value_REML,method="BH")<0.2)
table(p.adjust(beta_NGM$p_value_REML,method="BH")<0.2)


plot(beta_LOO_NGM$beta.1,-log10(p.adjust(beta_LOO_NGM$p_value_REML,method="BH")))
abline(h=-log10(0.05))
abline(h=-log10(0.2))







