
FS_NGM=read.table("GridLMM_results/gwas_results_FS_NGM.txt",h=TRUE)

names(FS_NGM)[2]="rs"
traitJ = read.table(paste0("gemma_results_BLUE/gwas_results_BLUE_T1.assoc.txt"),h=TRUE)
head(traitJ)

## Merge
merged=merge(traitJ,FS_NGM[,c("rs","beta.1")])

plot(merged$beta,merged$beta.1,pch=16)
cor(merged[,c("beta",paste0("beta.",1:5))])

SB_NGM=read.table("GridLMM_results/gwas_results_SB_NGM.txt",h=TRUE)
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