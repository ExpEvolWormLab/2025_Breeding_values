rm(list=ls())
gc()
library(qqman)


t1=read.table("gemma_results_BLUE/gwas_results_BLUE_T7.assoc.txt",h=TRUE)
t2=read.table("gemma_results_BLUE/gwas_results_BLUE_T14.assoc.txt",h=TRUE)
pvals_fdr <-t1$p_wald# p.adjust(t1$p_wald, method = "BH")
plot(-log10(pvals_fdr),pch=16,cex=.6)

plot((t1$beta),(-log10(pvals_fdr)),pch=16,col=adjustcolor(c("grey","red")[1+as.numeric(-log10(pvals_fdr)>2)],alpha=.6),cex=.6)

plot((1:length(t1$beta))[(-log10(pvals_fdr)<=2)],t1$beta[(-log10(pvals_fdr)<=2)],col=adjustcolor("grey",alpha=.6),pch=16,cex=.6,
     ylim=c(min(t1$beta),max(t1$beta)))
points((1:length(t1$beta))[(-log10(pvals_fdr)>2)],t1$beta[(-log10(pvals_fdr)>2)],col=adjustcolor("red",alpha=.6),pch=16,cex=.6)


data_manhattan = data.frame(CHR=as.numeric(as.factor(t1$chr)),BP=as.numeric(t1$ps),SNP=c(1:nrow(t1)),
                         P=-log10(t1$p_wald))
manhattan(subset(data_manhattan,P>.8), chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE)
subset(data_manhattan,P>6)

sign1=t1$p_wald<0.05
sign2=t2$p_wald<0.05
plot(t1$beta[sign1 & sign2],t2$beta[sign1 & sign2],pch=16,cex=.6)
plot(t1$beta,t2$beta[sign1 & sign2],pch=16,cex=.6)
abline(v=0)
abline(h=0)

table(sign1+sign2)


