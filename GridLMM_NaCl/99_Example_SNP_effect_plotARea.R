

idx_chr="I"

print(idx_chr)
# No singletons
df.CHR=subset(all_cluster_merged,CHROM==idx_chr & !is.na(cluster))

  df.CHR$group=as.numeric(as.factor(as.character(df.CHR$group)))
  # Singletons
  df.CHR.single=subset(all_cluster_merged,CHROM==idx_chr & is.na(cluster))
  df.CHR.single$group=as.numeric(as.factor(as.character(df.CHR.single$group)))
  
 # pdf(paste0("plots_BLUEs_GridLMM_NaCl/Cluster_CHR_R",round(10*threshold_clust),"_",idx_chr,"_Merged-noPruning.pdf"),h=10,w=20)
  
  all_col_CHR <- sample(rainbow(max(df.CHR$cluster)))
  all_pch_CHR <- sample(1:25,max(df.CHR$cluster),replace=TRUE)
  
  kept_chr = all_GWAS$CHROM==idx_chr
  i_env="NaCl"
  i_list_env <- 0
  
  i=7
    pdf("Cluster_CHRI_Area.pdf",h=5)
    plot(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="n",ylab="SNP effect",xlab="CHR pos (bp)",main=paste0("CHR ",idx_chr,' - Trait: ',vect_P_traits[i]," - ",i_env))
    abline(h=0,lty=2,col="orange",lwd=2)
    vcol=rep(adjustcolor('grey',alpha=.3),nrow(all_GWAS))
    points(all_GWAS$ps[kept_chr],all_GWAS[kept_chr,paste0("beta_T",(i+i_list_env))],pch=16,cex=.8,type="p",col=vcol[kept_chr])
    
    # Plot singletons
    df.CHR.single.trait <- subset(df.CHR.single,df.CHR.single[,paste0("is_sign_T",(i+i_list_env))])
    points(df.CHR.single.trait$ps,df.CHR.single.trait[,paste0("beta_T",(i+i_list_env))],bg="black",pch=21)
    
    ## Plot clusters (no singletons)
    df.CHR.trait <- subset(df.CHR,df.CHR[,paste0("is_sign_T",(i+i_list_env))])
    points(df.CHR.trait$ps,df.CHR.trait[,paste0("beta_T",(i+i_list_env))],col=all_col_CHR[df.CHR.trait$cluster],pch=all_pch_CHR[df.CHR.trait$cluster])
    

  
  dev.off()