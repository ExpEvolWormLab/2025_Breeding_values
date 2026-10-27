rm(list=ls());gc()
library(data.table)
library(lme4)

threshold_clust =0.5
load(paste0("RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
load(file="GridLMM_results_NaCl/All_Gwas.RData")

vect_P_traits=paste0(c("SF",'SB','FS','FB','BS','BF',"area.F","beta"),"_NaCl")

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits)])
dim(all_GWAS_clusters)
table(all_GWAS_clusters$is_sign_area.F_NaCl)
table(all_GWAS_clusters$sign_overall)

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
  #read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

df_clusters=read.table(paste0("data_BLUE_GridLMM_NaCl/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters) # 1058 clusters for R2
### Add the tip/arm/center tag
rkmap <- read.table("Caendr_filtering/rockman_kruglya_map_liftOver.txt",h=TRUE)
all_GWAS$CHROM=tstrsplit(all_GWAS$rs,"_")[[1]]
table(all_GWAS$CHROM)
all_GWAS$chr_type=NA
for(i in 1:nrow(rkmap)){
  if(sum(all_GWAS$CHROM==rkmap$CHROM[i] & all_GWAS$ps>=rkmap$Start[i] & all_GWAS$ps<=rkmap$End[i])>0){
  all_GWAS$chr_type[all_GWAS$CHROM==rkmap$CHROM[i] & all_GWAS$ps>=rkmap$Start[i] & all_GWAS$ps<=rkmap$End[i]]=rkmap$type[i]
  }
}
table(all_GWAS$chr_type)
table(is.na(all_GWAS$chr_type))
summary(subset(all_GWAS,is.na(chr_type) & CHROM=='V')[,c("CHROM","ps")]$ps) # Some tips on CHR V
all_GWAS$chr_type[all_GWAS$CHROM=="V" & all_GWAS$ps>=20341612]="right_tip"
table(subset(all_GWAS,is.na(chr_type))$CHROM) # Only MtDNA left

## Now also add the group ID to the GWAS
dim(all_GWAS)
dim(all_cluster_merged)
all_GWAS_sign <- (merge(all_GWAS,all_cluster_merged[,c("rs","group","sign_overall","POS")]))
dim(all_GWAS_sign)
names(all_cluster_merged)[!names(all_cluster_merged)%in%names(all_GWAS_sign)]

### How much is in Caendr  - Founders
#library(vcfR)
#founders_Caendr <- read.vcfR("Caendr_filtering/WI.20250625.hard-filter.isotype.CemeeFounders_poly.vcf.gz")
#caendr_pos <- as.data.frame(founders_Caendr@fix[,1:2])

#sign_for_Caendr <- all_GWAS_sign
#names(sign_for_Caendr)[2:5]=c("CHROM","POS","ALT","REF")
#dim(merge(sign_for_Caendr,data.frame(founders_Caendr@fix)))
#dim(merge(sign_for_Caendr,data.frame(founders_Caendr@fix[,1:2])))

gp_nbSNPs=tapply(all_GWAS_sign$group,all_GWAS_sign$group,length)
gp_min=tapply(all_GWAS_sign$ps,all_GWAS_sign$group,min)
gp_max=tapply(all_GWAS_sign$ps,all_GWAS_sign$group,max)
gp_type=tapply(all_GWAS_sign$chr_type,all_GWAS_sign$group,unique)
table(as.numeric(lapply(gp_type,length)))

#gp_type[as.numeric(lapply(gp_type,length))==2] # We decide center>arm>tip
for(i in c(1:length(gp_type))[as.numeric(lapply(gp_type,length))==2]){
  if(gp_type[[i]][1]=="left_tip"){
    gp_type[[i]]='left_arm'
  }else if(gp_type[[i]][1]=="left_arm"){
    gp_type[[i]]='center'
  }else{
    gp_type[[i]]=gp_type[[i]][[2]]
  }
}
table(as.numeric(lapply(gp_type,length))) # OK
gp_type=unlist(gp_type)

gp_length=data.frame(group=names(gp_min),clust_length=as.numeric(gp_max-gp_min),nb_SNPs=as.numeric(gp_nbSNPs),chr_type=gp_type)
head(gp_length)

pdf("Cluster_description_R5.pdf",h=5,w=5)
plot(gp_length$nb_SNPs,gp_length$clust_length,pch=16,cex=.8,col=c("red","cornflowerblue","lightgreen","lightgreen","cornflowerblue")[as.numeric(as.factor(gp_length$chr_type))],
     log="xy",xlab="Number of SNPs in QTL",ylab="QTL length (bp)")
legend(50,100,c("center","arm","tip"),pch=16,col=c("red","cornflowerblue","lightgreen"))
dev.off()
pdf("Cluster_nbSNPs_R5.pdf",h=3,w=4)
hist(log10(gp_length$nb_SNPs),xlab="Nb SNPs (log10)",main='')
dev.off()
pdf("Cluster_length_R5.pdf",h=3,w=4)
hist(log10(gp_length$clust_length),xlab="Cluster length (log10 bp)",main='')
dev.off()

## ClusterID & group are the same thing
all_GWAS_clusters=merge(all_GWAS_clusters,gp_length)
head(all_GWAS_clusters)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]

hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+clust_length+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic

boxplot(all_GWAS_clusters$clust_length+1~all_GWAS_clusters$sign_overall,log="y")

### Add the pleiotropy
#dim(merge(all_GWAS_sign,unique(all_GWAS_clusters[,c("group","sign_overall","sign_NaCl","sign_NGM")])))

#pheno_NGM=read.table("Breeding_values_NGM.txt",sep='\t',h=TRUE)
#pheno_NaCl=read.table("Breeding_values_NaCl.txt",sep='\t',h=TRUE)
#gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)

# Clusters
sign_all_df_clust <- all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits)]
table(rowSums(sign_all_df_clust)>=1) # ~737 Clusters / 10 are no longer significant
all_GWAS_clusters$sign_overall <- rowSums(sign_all_df_clust)>=1

for(i in 1:8){
  all_GWAS[,paste0("is_sign_",vect_P_traits[i])] = p.adjust(all_GWAS[,paste0("pvalue_T",i)],method="BH")<0.2
  print(table(all_GWAS[,paste0("is_sign_",vect_P_traits[i])]))
}

# SNPs
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_",vect_P_traits)]
table(rowSums(sign_all_df)>=1) # 33k SNPs

all_GWAS$sign_overall <- rowSums(sign_all_df)>=1

rowS_NaCl_clust <-table(rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits)]))
rowS_NaCl <-table(rowSums(all_GWAS[,paste0("is_sign_",vect_P_traits)]))

vect_col=c(adjustcolor("grey",alpha=0.3),adjustcolor("firebrick3",alpha=1))
i=7
plot(-log10(all_GWAS_clusters[,paste0("p_value_REML.cor_",vect_P_traits[i])])~all_GWAS_clusters[,paste0("beta_",vect_P_traits[i])],
     pch=16,xlab='Beta coefficient',ylab="-log(p)",
     col=vect_col[1+as.numeric(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits[i])])],
     main=paste0('Trait: ',vect_P_traits[i]))

#plot((all_GWAS_clusters$beta_SF_NaCl),(all_GWAS_clusters$beta_area.F_NaCl),pch=16,xlim=c(-.1,.1),ylim=c(-.1,.1))


### Heatmap plots between beta/sign ?
heatmap(cor(na.omit(all_GWAS_clusters[,paste0("beta_",vect_P_traits)])),
        labRow=vect_P_traits,
        labCol=vect_P_traits)

## With the non-significant
plot(rowS_NaCl/sum(rowS_NaCl),type='b',col='cornflowerblue',ylim=c(.00000001,1),log="y",xlim=c(0,8))
points(rowS_NaCl_clust/sum(rowS_NaCl_clust),type='b',col='darkorange',lty=2)


pdf("Lille_Pleio_R5_clusters_NaCl.pdf",w=4,h=4)
plot(rowS_NaCl,type='b',col=adjustcolor('darkorange',alpha=.4),ylim=c(1,1000000),log="y",xlim=c(0,8),las=1,bty="n",xlab="",xaxt="n",ylab="")
axis(1,0:8)
points(rowS_NaCl_clust,type='b',col='darkorange',lty=2)
legend(4,1000000,lty=c(1,2),lwd=3,col=c("darkorange","darkorange"),c('SNPs',"LD blocks"))
dev.off()

library(gplots)
library(limma) 
vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_SF_NaCl),as.numeric(all_GWAS_clusters$is_sign_area.F_NaCl))),
  names=c("SF NaCl","Area NaCl"),main='All')

load('RData/Gmatrices_P.RData')
heatmap(G_mat_all)

#### Now the Neff to estimate the effect ~pleiotropy ?
## For each cluster, compute a level of pleiotropy

all_GWAS_clusters_sign_NaCl <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits)])>0,])
nb_sign_clust_NaCl <- nrow(all_GWAS_clusters_sign_NaCl)
all_GWAS_clusters_sign=all_GWAS_clusters_sign_NaCl
pleiot.idx=data.frame(NaCl=rep(0,nb_sign_clust_NaCl))
pleiot.idx$NaCl.cor=0

for(i in 1:nrow(pleiot.idx)){
  
  ## NaCl second
  pleiot.idx$NaCl[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits)],na.rm=TRUE)
  if(pleiot.idx$NaCl[i]==1){
    pleiot.idx$NaCl.cor[i]=1 
  }else if(pleiot.idx$NaCl[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits)])
    std_devs.NaCl.temp <- sqrt(diag(Gmat_NaCl[idx.sign,idx.sign]))
    cor_matrix.NaCl.temp <- Gmat_NaCl[idx.sign,idx.sign] / (std_devs.NaCl.temp %o% std_devs.NaCl.temp)
    pleiot.idx$NaCl.cor[i]=pleiot.idx$NaCl[i] - var(eigen(cor_matrix.NaCl.temp)$values)
  }else{
    pleiot.idx$NaCl.cor[i]=0 
  }
  
}

##### Now compute the summed effect of the top SNP of each signicant cluster ####
### We need to compute the standardized allelic effect ###
all_effects <- all_GWAS_clusters_sign[,paste0('beta_',vect_P_traits)]
all_effects_std <- NULL
for(i in 1:8){
  all_effects_std[i] =
    sd(all_GWAS_clusters_sign[all_GWAS_clusters_sign[,paste0("is_sign_",vect_P_traits[i])],paste0("beta_",vect_P_traits[i])])
}

all_GWAS_clusters_sign_std <- all_GWAS_clusters_sign
for(i in 1:8){
  all_GWAS_clusters_sign_std[,paste0('beta_',vect_P_traits[i])]=all_GWAS_clusters_sign_std[,paste0('beta_',vect_P_traits[i])]/all_effects_std[i]
}

 euclidian_NaCl <- NULL
  manh_NaCl <- NULL
 for(i in 1:nrow(all_GWAS_clusters_sign_NaCl)){
      euclidian_NaCl[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits)])]^2))
      manh_NaCl=c(manh_NaCl,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits)])])))
   
 }

 ### With nls - replace all
 pdf("plots_BLUEs_GridLMM_NaCl/Scaling_pleiotropy_Neff_BLUE_Arm_Centers_R5.pdf",h=10,w=12)
 
 is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
 is_center <- all_GWAS_clusters_sign$chr_type2=="center"
 
 temp.NaCl <- euclidian_NaCl>0
 temp=pleiot.idx$NaCl.cor[temp.NaCl]
 coefs_euc.NaCl <- coef(nls(euclidian_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_euc.NaCl.arm <- coef(nls(euclidian_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_euc.NaCl.center <- coef(nls(euclidian_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,10),type='n')
 points(euclidian_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=16,col='green')
 points(euclidian_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NaCl[1]+(coefs_euc.NaCl[2])*(seq(0,12,0.1)^coefs_euc.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.arm[1]+(coefs_euc.NaCl.arm[2])*(seq(0,12,0.1)^coefs_euc.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.center[1]+(coefs_euc.NaCl.center[2])*(seq(0,12,0.1)^coefs_euc.NaCl.center[3])),col="lightblue")
 
 abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))
 text(1,8.5,paste0("c=",round(coefs_euc.NaCl[3],digits=3)))
 text(1,7.5,paste0("c=",round(coefs_euc.NaCl.arm[3],digits=3)),col="green")
 text(1,6.5,paste0("c=",round(coefs_euc.NaCl.center[3],digits=3)),col='lightblue')
 
 temp.NaCl <- manh_NaCl>0
 temp=pleiot.idx$NaCl.cor[temp.NaCl]
 coefs_man.NaCl <- coef(nls(manh_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))

 plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,7),ylim=c(0,25),type='n')
 points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(0,12,0.1)^coefs_man.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(0,12,0.1)^coefs_man.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(0,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue")
 
 abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))
 text(1.5,20,paste0("c=",round(coefs_man.NaCl[3],digits=3)))
 text(1.5,18,paste0("c=",round(coefs_man.NaCl.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NaCl.center[3],digits=3)),col='lightblue')
 
 
 dev.off()
 
 # pdf("Extra_synergistic_NaCl_R2.pdf",h=5,w=5)
 # temp.NaCl <- manh_NaCl>0
 # temp=pleiot.idx$NaCl.cor[temp.NaCl]
 # coefs_man.NaCl <- coef(nls(manh_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 # temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 # coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 # temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 # coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 # 
 # plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='Euclidian summed effect',xlab="Neff (number of trait affected)",main='NaCl',
 #      xlim=c(0,5),ylim=c(0,25),type='n',bty="n",las=1)
 # points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 # points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 # 
 # lines(seq(1,12,0.1),seq(1,12,0.1)*coef(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))[2],lwd=3) 
 # lines(seq(1,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(1,12,0.1)^coefs_man.NaCl[3])),col="red",lwd=3)
 # lines(seq(1,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(1,12,0.1)^coefs_man.NaCl.arm[3])),col="green",lwd=3)
 # lines(seq(1,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(1,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue",lwd=3)
 # 
 # legend(0,25,bty="n",lwd=4,c("All clusters","Chromosome arms","Chromosome centers","No synergistic pleiotropy (hypothetic)"),col=c("red","green","lightblue","black"))
 # 
 # 
 # dev.off()
 # 
 
 
 pdf("Extra_synergistic_NaCl_R2.pdf",h=5,w=5)
 
 temp=pleiot.idx$NaCl.cor
 coefs_man.NaCl <- coef(nls(manh_NaCl~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 
 plot(manh_NaCl~pleiot.idx$NaCl.cor,pch=16,ylab='Euclidian summed effect',xlab="Neff (number of trait affected)",main='NaCl',
      xlim=c(0,5),ylim=c(0,25),type='n',bty="n",las=1)
 points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(1,12,0.1),seq(1,12,0.1)*coef(lm(manh_NaCl~pleiot.idx$NaCl.cor))[2],lwd=3) 
 lines(seq(1,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(1,12,0.1)^coefs_man.NaCl[3])),col="red",lwd=3)
 lines(seq(1,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(1,12,0.1)^coefs_man.NaCl.arm[3])),col="green",lwd=3)
 lines(seq(1,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(1,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue",lwd=3)
 
 legend(0,25,bty="n",lwd=4,c("All clusters","Chromosome arms","Chromosome centers","No synergistic pleiotropy (hypothetic)"),col=c("red","green","lightblue","black"))
 
 
 dev.off()
 
 
 
# #################
# 
# ### The same with N ? ##
 ### With nls - replace all
 pdf("plots_BLUEs_GridLMM_NaCl/Scaling_pleiotropy_N_BLUE_Arm_Centers_R2.pdf",h=10,w=12)
 
 is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
 is_center <- all_GWAS_clusters_sign$chr_type2=="center"
 
 temp.NaCl <- euclidian_NaCl>0
 temp=pleiot.idx$NaCl[temp.NaCl]
 coefs_euc.NaCl <- coef(nls(euclidian_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_arm]
 coefs_euc.NaCl.arm <- coef(nls(euclidian_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_center]
 coefs_euc.NaCl.center <- coef(nls(euclidian_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="N",main='NaCl',xlim=c(0,5.5),ylim=c(0,10),type='n')
 points(euclidian_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl[temp.NaCl & is_arm],pch=16,col='green')
 points(euclidian_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl[temp.NaCl & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NaCl[1]+(coefs_euc.NaCl[2])*(seq(0,12,0.1)^coefs_euc.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.arm[1]+(coefs_euc.NaCl.arm[2])*(seq(0,12,0.1)^coefs_euc.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.center[1]+(coefs_euc.NaCl.center[2])*(seq(0,12,0.1)^coefs_euc.NaCl.center[3])),col="lightblue")
 
 abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]))
 text(1,8.5,paste0("c=",round(coefs_euc.NaCl[3],digits=3)),col="red")
 text(1,7.5,paste0("c=",round(coefs_euc.NaCl.arm[3],digits=3)),col="green")
 text(1,6.5,paste0("c=",round(coefs_euc.NaCl.center[3],digits=3)),col='lightblue')
 
 temp.NaCl <- manh_NaCl>0
 temp=pleiot.idx$NaCl[temp.NaCl]
 coefs_man.NaCl <- coef(nls(manh_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_arm]
 coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_center]
 coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 
 plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TM',xlab="N",main='NaCl',xlim=c(0,7),ylim=c(0,25),type='n')
 points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('green',alpha=.6))
 points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(0,12,0.1)^coefs_man.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(0,12,0.1)^coefs_man.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(0,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue")
 
 abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]))
 text(1.5,20,paste0("c=",round(coefs_man.NaCl[3],digits=3)),col="red")
 text(1.5,18,paste0("c=",round(coefs_man.NaCl.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NaCl.center[3],digits=3)),col='lightblue')
 
 dev.off()
 

########################################################################
 ########################################################################
 ########################################################################
 
 