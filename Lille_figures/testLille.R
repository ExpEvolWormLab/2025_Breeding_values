
i_list_env=0
pdf("Plot_Lille_01.pdf",h=5)

plot(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="n",ylab="SNP effect",xlab="Physical position (bp)",main=paste0("CHR ",idx_chr,' - Trait: ',vect_P_traits[i]," - ",i_env),
     las=1,bty="n")
abline(h=0,lty=2,col="orange",lwd=2)
vcol=rep(adjustcolor('grey',alpha=.3),nrow(all_GWAS))
points(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="p",col=vcol[kept_chr])

# Plot singletons
df.CHR.single.trait <- subset(df.CHR.single,df.CHR.single[,paste0("is_sign_T",(i+i_list_env))])
points(df.CHR.single.trait$ps,df.CHR.single.trait[,paste0("beta_T",(i+i_list_env))],bg="black",pch=21,cex=.8)

## Plot clusters (no singletons)
df.CHR.trait <- subset(df.CHR,df.CHR[,paste0("is_sign_T",(i+i_list_env))])
points(df.CHR.trait$ps,df.CHR.trait[,paste0("beta_T",(i+i_list_env))],col=all_col_CHR[df.CHR.trait$cluster],pch=all_pch_CHR[df.CHR.trait$cluster])

dev.off()

i_list_env=7

pdf("Plot_Lille_02.pdf",h=5)

plot(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="n",ylab="SNP effect",xlab="Physical position (bp)",main=paste0("CHR ",idx_chr,' - Trait: ',vect_P_traits[i]," - ",i_env),
     las=1,bty="n")
abline(h=0,lty=2,col="orange",lwd=2)
vcol=rep(adjustcolor('grey',alpha=.3),nrow(all_GWAS))
points(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="p",col=vcol[kept_chr])

# Plot singletons
df.CHR.single.trait <- subset(df.CHR.single,df.CHR.single[,paste0("is_sign_T",(i+i_list_env))])
points(df.CHR.single.trait$ps,df.CHR.single.trait[,paste0("beta_T",(i+i_list_env))],bg="black",pch=21,cex=.8)

## Plot clusters (no singletons)
df.CHR.trait <- subset(df.CHR,df.CHR[,paste0("is_sign_T",(i+i_list_env))])
points(df.CHR.trait$ps,df.CHR.trait[,paste0("beta_T",(i+i_list_env))],col=all_col_CHR[df.CHR.trait$cluster],pch=all_pch_CHR[df.CHR.trait$cluster])

dev.off()