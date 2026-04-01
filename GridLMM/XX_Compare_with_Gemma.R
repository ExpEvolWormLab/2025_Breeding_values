library(data.table)
library(qqman)
FS_NGM=read.table("GridLMM_results_mixed/gwas_results_area.F_NGM.txt",h=TRUE)

names(FS_NGM)[2]="rs"
traitJ = read.table(paste0("gemma_results_BLUE/gwas_results_BLUE_T7.assoc.txt"),h=TRUE)
head(traitJ)

## Merge
merged=merge(traitJ,FS_NGM[,c("rs","beta.1","beta.2","p_value_REML")])

plot(merged$beta,merged$beta.1,pch=16)

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