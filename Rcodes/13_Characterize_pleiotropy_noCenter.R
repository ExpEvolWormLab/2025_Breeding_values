rm(list=ls())
library(data.table)
load("RData/GWAS_cluster_results_concatenated.RData")
load(file="RData/GWAS_results_concatenated.RData")
GWAS_with_clust <- read.table("txt/all_sign_GWAS_with_clusters.txt",sep="\t",header=TRUE)

gemma_pheno=read.table('data/gemma_pheno.txt',sep='\t',h=FALSE)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)])
all_GWAS_clusters$sign_NGM <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:7)])
all_GWAS_clusters$sign_NaCl <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",8:14)])

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(file=paste0('txt/all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]
all_cluster_merged$group = as.numeric(factor(as.character(all_cluster_merged$group),levels=unique(all_cluster_merged$group)))

df_clusters = data.frame(clusterID=1:max(all_cluster_merged$group),
                         CHROM= tapply(as.character(all_cluster_merged$CHROM),all_cluster_merged$group,unique),
                         POS= tapply((all_cluster_merged$POS),all_cluster_merged$group,mean),group=tapply(as.character(all_cluster_merged$group),all_cluster_merged$group,unique))
df_clusters$ps=floor(df_clusters$POS)
all_GWAS_clusters$CHROM=all_GWAS_clusters$chr
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)

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
## Only 28 with 2 types
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

plot(gp_length$nb_SNPs,gp_length$clust_length,pch=16,cex=.8,col=c("red","cornflowerblue","lightgreen","lightgreen","cornflowerblue")[as.numeric(as.factor(gp_length$chr_type))])
#levels(as.factor(gp_length$chr_type))
all_GWAS_clusters=merge(all_GWAS_clusters,gp_length)
#plot(all_GWAS_clusters$ps,all_GWAS_clusters$clust_length/1000000,pch=16)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]

hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+log(clust_length+1)+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic

boxplot(all_GWAS_clusters$clust_length+1~all_GWAS_clusters$sign_overall,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NGM,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NaCl,log="y")

### Add the pleiotropy
dim(merge(all_GWAS_sign,unique(all_GWAS_clusters[,c("group","sign_overall","sign_NaCl","sign_NGM")])))


### Same but we will remove center/long clusters ?
plot(density(log10(1+all_GWAS_clusters$clust_length[!is.na(all_GWAS_clusters$chr_type2) & all_GWAS_clusters$chr_type2=="center"])))
lines(density(log10(1+all_GWAS_clusters$clust_length[!is.na(all_GWAS_clusters$chr_type2) & all_GWAS_clusters$chr_type2=="tip"])),col="blue")
lines(density(log10(1+all_GWAS_clusters$clust_length[!is.na(all_GWAS_clusters$chr_type2) & all_GWAS_clusters$chr_type2=="arm"])),col="orange")


plot(density(log10(all_GWAS_clusters$nb_SNPs[!is.na(all_GWAS_clusters$chr_type2) & all_GWAS_clusters$chr_type2=="center"])))
lines(density(log10(all_GWAS_clusters$nb_SNPs[!is.na(all_GWAS_clusters$chr_type2) & all_GWAS_clusters$chr_type2=="tip"])),col="blue")
lines(density(log10(all_GWAS_clusters$nb_SNPs[!is.na(all_GWAS_clusters$chr_type2) & all_GWAS_clusters$chr_type2=="arm"])),col="orange")


gemma_pheno=read.table('data/gemma_pheno.txt',sep='\t',h=FALSE)

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

# SNPs
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:14)]
table(rowSums(sign_all_df)>=1) # 107k SNPs
all_GWAS$sign_overall <- rowSums(sign_all_df)

## We remove the long cluster and the center clusters
table(all_GWAS_clusters$clust_length>10000 | all_GWAS_clusters$chr_type2!="arm") # 2577 clusters removed
#table(all_GWAS_clusters$nb_SNPs>=2) # 4323
#table(all_GWAS_clusters$clust_length>10000) # 1734
#table(all_GWAS_clusters$chr_type2!="arm") # 568

all_GWAS_clusters_filtered <- subset(all_GWAS_clusters,chr_type2=="arm" & clust_length<10000)
dim(all_GWAS_clusters_filtered) # 5160
table(all_GWAS_clusters_filtered$sign_overall>0)  # 3967 significant

model1 <- glm(sign_overall~log(nb_SNPs)+log(clust_length+1),
              family="quasipoisson",data=subset(all_GWAS_clusters_filtered,sign_overall>0 ))
summary(model1) # Yes, the more SNPs, the more pleiotropic
### Still significant.. but less

rowS_NGM_clust <- table(rowSums(all_GWAS_clusters_filtered[,paste0("is_sign_T",1:7)]))
rowS_NaCl_clust <-table(rowSums(all_GWAS_clusters_filtered[,paste0("is_sign_T",8:14)]))

#rowS_NGM <- table(rowSums(all_GWAS[,paste0("is_sign_T",1:7)]))
#rowS_NaCl <-table(rowSums(all_GWAS[,paste0("is_sign_T",8:14)]))

plot(cor(all_GWAS_clusters_filtered[,paste0("beta_T",1:7)]),
     cor(na.omit(all_GWAS_clusters_filtered[,paste0("beta_T",8:14)])))
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)

## Export all volcano plots
pdf("Volcano_plots_clusters.pdf")
vect_col=c(adjustcolor("grey",alpha=0.3),adjustcolor("firebrick3",alpha=1))
for(i in 1:14){
  plot(-log10(all_GWAS_clusters_filtered[,paste0("p_wald_T",i)])~all_GWAS_clusters_filtered[,paste0("beta_T",i)],
       pch=16,xlab='Beta coefficient',ylab="-log(p)",
       col=vect_col[1+as.numeric(all_GWAS_clusters_filtered[,paste0("is_sign_T",i)])],
       main=paste0('Trait: ',vect_P_traits_bothEnv[i]))
}
dev.off()

plot((all_GWAS_clusters_filtered$beta_T7),(all_GWAS_clusters_filtered$beta_T14),pch=16)


### Heatmap plots between beta/sign ?
heatmap(cor(na.omit(all_GWAS_clusters_filtered[,paste0("beta_T",1:14)])),
        labRow=vect_P_traits_bothEnv,
        labCol=vect_P_traits_bothEnv)

## With the non-significant
#plot(rowS_NGM/sum(rowS_NGM),type='b',col='cornflowerblue',ylim=c(.0001,1),log="y")
#points(rowS_NaCl/sum(rowS_NaCl),type='b',col='darkorange')
plot(rowS_NaCl_clust/sum(rowS_NaCl_clust),type='b',col='darkorange',lty=2,ylim=c(.0001,1),log="y")
points(rowS_NGM_clust/sum(rowS_NGM_clust),type='b',col='cornflowerblue',lty=2)

## Without the non-significant
pdf('plots/Pleiotropie_NGM_NaCl_nolargeClust.pdf')
rowS_NGM_clust2=rowS_NGM_clust[2:8]
rowS_NaCl_clust2=rowS_NaCl_clust[2:8]
plot(rowS_NGM_clust2/sum(rowS_NGM_clust2),type='b',col='cornflowerblue',ylim=c(.0001,1),log="y",
     ylab="Number of significant SNPs/Clusters",xlab='Number of traits')
points(rowS_NaCl_clust2/sum(rowS_NaCl_clust2),type='b',col='darkorange',lty=2)
dev.off()

pdf('plots/Pleiotropie_NGM_NaCl_absolute.pdf')
plot(rowS_NGM_clust2,type='b',col='cornflowerblue',ylim=c(1,45000),log="y",
     ylab="Number of significant SNPs/Clusters",xlab='Number of traits')
points(rowS_NaCl_clust2,type='b',col='darkorange',lty=2)
dev.off()

library(gplots)
library(limma) 
vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters_filtered$is_sign_T1),as.numeric(all_GWAS_clusters_filtered$is_sign_T8))),
  names=c("SF NGM","SF NaCl"),main='All')

vennDiagram(vennCounts(
  all_GWAS_clusters_filtered[,paste0("is_sign_T",c(14,11))]),
  names=c("test","Area NGM"),main='All')

### Could there be SNPs in the same genes but not the same variants 
### affecting the traits ? -> would be nice to give examples ????

#### Now the Neff to estimate the effect ~pleiotropy ?
## For each cluster, compute a level of pleiotropy

all_GWAS_clusters_filtered_sign <- na.omit(all_GWAS_clusters_filtered[rowSums(all_GWAS_clusters_filtered[,paste0("is_sign_T",1:14)])>0,])
all_GWAS_clusters_filtered_sign_NGM <- na.omit(all_GWAS_clusters_filtered[rowSums(all_GWAS_clusters_filtered[,paste0("is_sign_T",1:7)])>0,])
all_GWAS_clusters_filtered_sign_NaCl <- na.omit(all_GWAS_clusters_filtered[rowSums(all_GWAS_clusters_filtered[,paste0("is_sign_T",8:14)])>0,])

nb_sign_clust_NGM <- nrow(all_GWAS_clusters_filtered_sign)
nb_sign_clust_NaCl <- nrow(all_GWAS_clusters_filtered_sign)
nb_sign_clust <- nrow(all_GWAS_clusters_filtered_sign)

pleiot.idx=data.frame(NGM=rep(0,nb_sign_clust),NaCl=rep(0,nb_sign_clust))
pleiot.idx$NGM.cor=0;pleiot.idx$NaCl.cor=0

pleiot.idx_mixed=data.frame(N=rep(0,nb_sign_clust),Neff=rep(0,nb_sign_clust))

Gmat_NGM = cov(na.omit(all_GWAS_clusters_filtered[,paste0("beta_T",1:7)]))
Gmat_NaCl = cov(na.omit(all_GWAS_clusters_filtered[,paste0("beta_T",8:14)]))
G_mat_all=cov(na.omit(all_GWAS_clusters_filtered[,paste0("beta_T",1:14)]))

for(i in 1:nrow(pleiot.idx_mixed)){
  ## NGM first
  pleiot.idx$NGM[i]=sum(all_GWAS_clusters_filtered_sign[i,paste0("is_sign_T",1:7)],na.rm=TRUE)
  if(pleiot.idx$NGM[i]==1){
    pleiot.idx$NGM.cor[i]=1 
  }else if(pleiot.idx$NGM[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_filtered_sign[i,paste0("is_sign_T",1:7)])
    std_devs.NGM.temp <- sqrt(diag(Gmat_NGM[idx.sign,idx.sign]))
    cor_matrix.NGM.temp <- Gmat_NGM[idx.sign,idx.sign] / (std_devs.NGM.temp %o% std_devs.NGM.temp)
    pleiot.idx$NGM.cor[i]=pleiot.idx$NGM[i] - var(eigen(cor_matrix.NGM.temp)$values)
  }
  
  ## NaCl second
  pleiot.idx$NaCl[i]=sum(all_GWAS_clusters_filtered_sign[i,paste0("is_sign_T",8:14)],na.rm=TRUE)
  if(pleiot.idx$NaCl[i]==1){
    pleiot.idx$NaCl.cor[i]=1 
  }else if(pleiot.idx$NaCl[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_filtered_sign[i,paste0("is_sign_T",8:14)])
    std_devs.NaCl.temp <- sqrt(diag(Gmat_NaCl[idx.sign,idx.sign]))
    cor_matrix.NaCl.temp <- Gmat_NaCl[idx.sign,idx.sign] / (std_devs.NaCl.temp %o% std_devs.NaCl.temp)
    pleiot.idx$NaCl.cor[i]=pleiot.idx$NaCl[i] - var(eigen(cor_matrix.NaCl.temp)$values)
  }
  
  ### Mixed
  pleiot.idx_mixed$N[i]=sum(all_GWAS_clusters_filtered_sign[i,paste0("is_sign_T",1:14)],na.rm=TRUE)
  
  if(pleiot.idx_mixed$N[i]==1){
    pleiot.idx_mixed$Neff[i]=1  
  }else if(pleiot.idx_mixed$N[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_filtered_sign[i,paste0("is_sign_T",1:14)])
    idx.sign[is.na(idx.sign)]=FALSE
    std_devs.all.temp <- sqrt(diag(G_mat_all[idx.sign,idx.sign]))
    cor_matrix.all.temp <- G_mat_all[idx.sign,idx.sign] / (std_devs.all.temp %o% std_devs.all.temp)
    pleiot.idx_mixed$Neff[i]=pleiot.idx_mixed$N[i] - var(eigen(cor_matrix.all.temp)$values)
  }
}

par(mfrow=c(1,3))
hist(pleiot.idx$NGM.cor,main='NGM',xlim=c(0,6),ylim=c(0,3000))
hist(pleiot.idx$NaCl.cor,main='NaCl',xlim=c(0,6),ylim=c(0,3000))
hist(pleiot.idx_mixed$Neff,main='ALl 12',xlim=c(0,12),ylim=c(0,5000))

##### Now compute the summed effect of the top SNP of each signicant cluster ? ####

### We need to compute the standardized allelic effect ###
all_effects <- all_GWAS_clusters_filtered_sign[,paste0('beta_T',1:14)]
all_effects_std <- NULL
for(i in 1:14){
  all_effects_std[i] =
    sd(all_GWAS_clusters_filtered_sign[all_GWAS_clusters_filtered_sign[,paste0("is_sign_T",i)],paste0("beta_T",i)])
}

all_GWAS_clusters_filtered_sign_std <- all_GWAS_clusters_filtered_sign
for(i in 1:14){
  all_GWAS_clusters_filtered_sign_std[,paste0('beta_T',i)]=all_GWAS_clusters_filtered_sign_std[,paste0('beta_T',i)]/all_effects_std[i]
}

euclidian_summed <- NULL;euclidian_NGM <- NULL;euclidian_NaCl <- NULL
manh_summed <- NULL; manh_NGM <- NULL; manh_NaCl <- NULL
for(i in 1:nrow(all_GWAS_clusters_filtered_sign_std)){
  euclidian_summed <- c(euclidian_summed,sqrt(sum(as.numeric(all_GWAS_clusters_filtered_sign_std[i,paste0('beta_T',1:14)])[as.logical(all_GWAS_clusters_filtered_sign_std[i,paste0('is_sign_T',1:14)])]^2)))
  euclidian_NGM[i]=sqrt(sum(as.numeric(all_GWAS_clusters_filtered_sign_std[i,paste0('beta_T',1:7)])[as.logical(all_GWAS_clusters_filtered_sign_std[i,paste0('is_sign_T',1:7)])]^2))
  euclidian_NaCl[i]=sqrt(sum(as.numeric(all_GWAS_clusters_filtered_sign_std[i,paste0('beta_T',8:14)])[as.logical(all_GWAS_clusters_filtered_sign_std[i,paste0('is_sign_T',8:14)])]^2))
  
  manh_summed <- c(manh_summed,(sum(abs(as.numeric(all_GWAS_clusters_filtered_sign_std[i,paste0('beta_T',1:14)]))[as.logical(all_GWAS_clusters_filtered_sign_std[i,paste0('is_sign_T',1:14)])])))
  manh_NGM=c(manh_NGM,(sum(abs(as.numeric(all_GWAS_clusters_filtered_sign_std[i,paste0('beta_T',1:7)]))[as.logical(all_GWAS_clusters_filtered_sign_std[i,paste0('is_sign_T',1:7)])])))
  manh_NaCl=c(manh_NaCl,(sum(abs(as.numeric(all_GWAS_clusters_filtered_sign_std[i,paste0('beta_T',8:14)]))[as.logical(all_GWAS_clusters_filtered_sign_std[i,paste0('is_sign_T',8:14)])])))
  
}

pdf("plots/Scaling_pleiotropy_Neff_no_large_clust.pdf",h=10,w=12)

coefs_euc <- as.numeric((lm(log(euclidian_summed)~log(pleiot.idx_mixed$Neff)))$coef)
temp.NGM <- euclidian_NGM>0
coefs_euc.NGM <- as.numeric((lm(log(euclidian_NGM[temp.NGM])~log(pleiot.idx$NGM.cor[temp.NGM])))$coef)
temp.NaCl <- euclidian_NaCl>0
coefs_euc.NaCl <- as.numeric((lm(log(euclidian_NaCl[temp.NaCl])~log(pleiot.idx$NaCl.cor[temp.NaCl])))$coef)

par(mfrow=c(2,3))
plot(euclidian_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,11),ylim=c(0,14))
lines(seq(0,12,0.1),I(exp(coefs_euc[1])*(seq(0,12,0.1)^coefs_euc[2])),col="red")
abline(lm(euclidian_summed~pleiot.idx_mixed$Neff-1))
text(2,12,paste0("b=",round(coefs_euc[2],digits=3)))

plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NGM',xlim=c(0,5.5),ylim=c(0,10))
lines(seq(0,12,0.1),I(exp(coefs_euc.NGM[1])*(seq(0,12,0.1)^coefs_euc.NGM[2])),col="red")
abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]-1))
text(1.5,7.5,paste0("b=",round(coefs_euc.NGM[2],digits=3)))

plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,10))
lines(seq(0,12,0.1),I(exp(coefs_euc.NaCl[1])*(seq(0,12,0.1)^coefs_euc.NaCl[2])),col="red")
abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]-1))
text(1,8.5,paste0("b=",round(coefs_euc.NaCl[2],digits=3)))

# mod1 <- lm(log(manh_summed)~log(pleiot.idx_mixed$Neff))
# mod2 <- lm(log(manh_summed)~offset(log(pleiot.idx_mixed$Neff)))
# summary(mod1)
# summary(mod2)
# anova(mod1,mod2)
# AIC(mod1)
# AIC(mod2)

coefs_man <- as.numeric((lm(log(manh_summed)~log(pleiot.idx_mixed$Neff)))$coef)
temp.NGM <- manh_NGM>0
coefs_man.NGM <- as.numeric((lm(log(manh_NGM[temp.NGM])~log(pleiot.idx$NGM.cor[temp.NGM])))$coef)
temp.NaCl <- manh_NaCl>0
coefs_man.NaCl <- as.numeric((lm(log(manh_NaCl[temp.NaCl])~log(pleiot.idx$NaCl.cor[temp.NaCl])))$coef)

plot(manh_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TM ',xlab="Neff",xlim=c(0,11),ylim=c(0,40))
lines(seq(0,12,0.1),I(exp(coefs_man[1])*(seq(0,12,0.1)^coefs_man[2])),col="red")
abline(lm(manh_summed~pleiot.idx_mixed$Neff-1))
text(2,30,paste0("b=",round(coefs_man[2],digits=3)))

plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(0,5.5),ylim=c(0,25))
lines(seq(0,12,0.1),I(exp(coefs_man.NGM[1])*(seq(0,12,0.1)^coefs_man.NGM[2])),col="red")
abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]-1))
text(1.5,20,paste0("b=",round(coefs_man.NGM[2],digits=3)))

plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,25))
lines(seq(0,12,0.1),I(exp(coefs_man.NaCl[1])*(seq(0,12,0.1)^coefs_man.NaCl[2])),col="red")
abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]-1))
text(1,20,paste0("b=",round(coefs_man.NaCl[2],digits=3)))

dev.off()



#################

### The same with N ? ##
pdf("plots/Scaling_pleiotropy_N.pdf",h=14,w=12)

coefs_euc <- as.numeric((lm(log(euclidian_summed)~log(pleiot.idx_mixed$N)))$coef)
temp.NGM <- euclidian_NGM>0
coefs_euc.NGM <- as.numeric((lm(log(euclidian_NGM[temp.NGM])~log(pleiot.idx$NGM[temp.NGM])))$coef)
temp.NaCl <- euclidian_NaCl>0
coefs_euc.NaCl <- as.numeric((lm(log(euclidian_NaCl[temp.NaCl])~log(pleiot.idx$NaCl.[temp.NaCl])))$coef)

par(mfrow=c(2,3))
plot(euclidian_summed~pleiot.idx_mixed$N,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,13),ylim=c(0,14))
lines(seq(0,12,0.1),I(exp(coefs_euc[1])*(seq(0,12,0.1)^coefs_euc[2])),col="red")
abline(lm(euclidian_summed~pleiot.idx_mixed$N-1))
text(2,12,paste0("b=",round(coefs_euc[2],digits=3)))

plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NGM',xlim=c(0,6.5),ylim=c(0,10))
lines(seq(0,12,0.1),I(exp(coefs_euc.NGM[1])*(seq(0,12,0.1)^coefs_euc.NGM[2])),col="red")
abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM]-1))
text(1.5,7.5,paste0("b=",round(coefs_euc.NGM[2],digits=3)))

plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,6.5),ylim=c(0,10))
lines(seq(0,12,0.1),I(exp(coefs_euc.NaCl[1])*(seq(0,12,0.1)^coefs_euc.NaCl[2])),col="red")
abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]-1))
text(1,8.5,paste0("b=",round(coefs_euc.NaCl[2],digits=3)))

coefs_man <- as.numeric((lm(log(manh_summed)~log(pleiot.idx_mixed$N)))$coef)
temp.NGM <- manh_NGM>0
coefs_man.NGM <- as.numeric((lm(log(manh_NGM[temp.NGM])~log(pleiot.idx$NGM[temp.NGM])))$coef)
temp.NaCl <- manh_NaCl>0
coefs_man.NaCl <- as.numeric((lm(log(manh_NaCl[temp.NaCl])~log(pleiot.idx$NaCl[temp.NaCl])))$coef)

plot(manh_summed~pleiot.idx_mixed$N,pch=16,ylab='TM ',xlab="Neff",xlim=c(0,12),ylim=c(0,40))
lines(seq(0,12,0.1),I(exp(coefs_man[1])*(seq(0,12,0.1)^coefs_man[2])),col="red")
abline(lm(manh_summed~pleiot.idx_mixed$N-1))
text(2,30,paste0("b=",round(coefs_man[2],digits=3)))

plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(0,6.5),ylim=c(0,25))
lines(seq(0,12,0.1),I(exp(coefs_man.NGM[1])*(seq(0,12,0.1)^coefs_man.NGM[2])),col="red")
abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.[temp.NGM]-1))
text(1.5,20,paste0("b=",round(coefs_man.NGM[2],digits=3)))

plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,6.5),ylim=c(0,25))
lines(seq(0,12,0.1),I(exp(coefs_man.NaCl[1])*(seq(0,12,0.1)^coefs_man.NaCl[2])),col="red")
abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]-1))
text(1,20,paste0("b=",round(coefs_man.NaCl[2],digits=3)))

dev.off()


###### Now we 







