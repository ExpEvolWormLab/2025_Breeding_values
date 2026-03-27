
all_GWAS.temp=read.table("GridLMM_results_mixed_scaled/gwas_results_SF_NGM.txt",h=TRUE)[,c("X_ID","beta.1","p_value_REML")]
names(all_GWAS.temp)[1]="rs"
all_GWAS.temp$chr = factor(tstrsplit(all_GWAS.temp$rs,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA"))
all_GWAS.temp$ps = as.numeric(tstrsplit(all_GWAS.temp$rs,"_")[[2]])
all_GWAS.temp=all_GWAS.temp[order(all_GWAS.temp$chr,all_GWAS.temp$ps),]


head(all_GWAS.temp)
table(p.adjust(all_GWAS.temp$p_value_REML,method='BH')<0.3)

data_manhattan = data.frame(CHR=as.numeric(all_GWAS.temp$chr),BP=all_GWAS.temp$ps,SNP=c(1:nrow(all_GWAS.temp)),
                            P=-log10(all_GWAS.temp$p_value_REML))

data_manhattan = data.frame(CHR=as.numeric(all_GWAS.temp$chr),BP=all_GWAS.temp$ps,SNP=c(1:nrow(all_GWAS.temp)),
                            P=-log10(p.adjust(all_GWAS.temp$p_value_REML,method="fdr")))
data_manhattan=subset(data_manhattan,CHR!=7)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,
          highlight=data_manhattan$SNP[data_manhattan$is_sign],cex=.6,suggestiveline=-log10(0.05),genomewideline=FALSE)
