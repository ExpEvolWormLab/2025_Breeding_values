rm(list=ls());gc()
library(data.table)

threshold_clust = 0.4
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
load(file="RData/GWAS_results_concatenated.RData")
gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)
#GWAS_with_clust <- read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)])
all_GWAS_clusters$sign_NGM <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:7)])
all_GWAS_clusters$sign_NaCl <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",8:14)])

all_GWAS_clusters$p_wald_T1_FDR = p.adjust(all_GWAS_clusters$p_wald_T1,method="BH")
summary(all_GWAS_clusters$p_wald_T1_FDR)
table(all_GWAS_clusters$p_wald_T1_FDR<0.05)

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(paste0("txt_BLUEs/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
  #read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

df_clusters=read.table(paste0("data_BLUE_clusters/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters)
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

## ClusterID & group are the same thing
all_GWAS_clusters=merge(all_GWAS_clusters,gp_length)
head(all_GWAS_clusters)
#plot(all_GWAS_clusters$ps,all_GWAS_clusters$clust_length/1000000,pch=16)

### We will now replace minor af
temp.af_T7=1-all_GWAS_clusters$af_T7
temp.af_T14=1-all_GWAS_clusters$af_T14
#temp.af_T7[is.na(temp.af_T7)]=1;temp.af_T14[is.na(temp.af_T14)]=1
all_GWAS_clusters$af_T7[!is.na(temp.af_T7) & all_GWAS_clusters$af_T7>temp.af_T7]=temp.af_T7[!is.na(temp.af_T7) & all_GWAS_clusters$af_T7>temp.af_T7]
all_GWAS_clusters$af_T14[!is.na(temp.af_T14) & all_GWAS_clusters$af_T14>temp.af_T14]=temp.af_T14[!is.na(temp.af_T14) & all_GWAS_clusters$af_T14>temp.af_T14]

#plot((all_GWAS_clusters$af_T3),(all_GWAS_clusters$af_T14),pch=16)

boxplot(all_GWAS_clusters$af_T7~all_GWAS_clusters$sign_NGM)
boxplot(all_GWAS_clusters$af_T14~all_GWAS_clusters$sign_NaCl)

tapply(all_GWAS_clusters$af_T7,all_GWAS_clusters$sign_NGM,mean)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]

summary(glm(sign_NGM~af_T7,data=all_GWAS_clusters,family="quasipoisson"))

af.per.pleio <- rbind(
 data.frame(sign=all_GWAS_clusters$sign_NGM,env_label='NGM',af=all_GWAS_clusters$af_T7,chr_type2=all_GWAS_clusters$chr_type2),
 data.frame(sign=all_GWAS_clusters$sign_NaCl,env_label='NaCl',af=all_GWAS_clusters$af_T14,chr_type2=all_GWAS_clusters$chr_type2))

af.per.pleio$sign=as.factor(af.per.pleio$sign)
af.per.pleio$sign_num=as.numeric(as.character(as.factor((af.per.pleio$sign))))

af.per.pleio=na.omit(af.per.pleio)
## only take minor alleles
af.per.pleio$af[af.per.pleio$af>0.5] = 1-af.per.pleio$af[af.per.pleio$af>0.5]

af.per.pleio$af.logit <- log(af.per.pleio$af/(1-af.per.pleio$af))

model.af <- (lm(af.logit~env_label*sign+chr_type2*env_label,data=subset(af.per.pleio,sign!=0)))
summary(model.af)
model.af2 <- (lm(af.logit~env_label*sign+chr_type2,data=subset(af.per.pleio,sign!=0)))
anova(model.af,model.af2)
summary(model.af2)
model.af3 <- (lm(af.logit~env_label+sign+chr_type2*env_label,data=subset(af.per.pleio,sign!=0)))
anova(model.af,model.af3)
model.af4 <- (lm(af.logit~chr_type2*env_label,data=subset(af.per.pleio,sign!=0)))
anova(model.af3,model.af2)
model.af5 <- (lm(af.logit~chr_type2+env_label,data=subset(af.per.pleio,sign!=0)))
anova(model.af5,model.af3)

summary(model.af3)

model.af <- (lm(af.logit~env_label*sign_num,data=subset(af.per.pleio,sign!=0)))
model.af2 <- (lm(af.logit~env_label+sign_num,data=subset(af.per.pleio,sign!=0)))
model.af3 <- (lm(af.logit~sign_num,data=subset(af.per.pleio,sign!=0)))
summary(model.af)
anova(model.af,model.af2)
anova(model.af2,model.af3)
anova(model.af2,model.af4)
hist(af.per.pleio$af*(1-af.per.pleio$af),n=40)

af.per.pleio.NGM <- subset(af.per.pleio,env_label=="NGM")
af.per.pleio.NaCl <- subset(af.per.pleio,env_label=="NaCl")
sd_NGM <- as.numeric(tapply(af.per.pleio.NGM$af,af.per.pleio.NGM$sign,sd))
mean_NGM <- as.numeric(tapply(af.per.pleio.NGM$af,af.per.pleio.NGM$sign,mean))

sd_NaCl <- as.numeric(tapply(af.per.pleio.NaCl$af,af.per.pleio.NaCl$sign,sd))
mean_NaCl <- as.numeric(tapply(af.per.pleio.NaCl$af,af.per.pleio.NaCl$sign,mean))

plot(mean_NGM~c(c(0:6)+0.1),pch=16,ylim=c(0,1),
     las=1,bty="n",ylab="Mean allele frequencies",xlab="Pleiotropy")
#arrows(c(c(0:7)+0.1),mean_NGM-2*sd_NGM,c(c(0:7)+0.1),mean_NGM+2*sd_NGM,code=3,angle=90)

points(mean_NaCl~c(c(0:6)-.1),pch=16,col="orange",bty="n",las=1)
#arrows(c(c(0:7)-0.1),mean_NaCl-2*sd_NaCl,c(c(0:7)-0.1),mean_NaCl+2*sd_NaCl,code=3,angle=90)

hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+clust_length+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic

boxplot(all_GWAS_clusters$clust_length+1~all_GWAS_clusters$sign_overall,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NGM,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NaCl,log="y")

### Add the pleiotropy
dim(merge(all_GWAS_sign,unique(all_GWAS_clusters[,c("group","sign_overall","sign_NaCl","sign_NGM")])))

#pheno_NGM=read.table("Breeding_values_NGM.txt",sep='\t',h=TRUE)
#pheno_NaCl=read.table("Breeding_values_NaCl.txt",sep='\t',h=TRUE)
gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

# Clusters
sign_all_df_clust <- all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)]
table(rowSums(sign_all_df_clust)>=1) # ~5411 Clusters / 1360 are no longer significant
all_GWAS_clusters$sign_overall <- rowSums(sign_all_df_clust)>=1

# SNPs
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:14)]
table(rowSums(sign_all_df)>=1) # 107k SNPs
all_GWAS$sign_overall <- rowSums(sign_all_df)>=1

rowS_NGM_clust <- table(rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:7)]))
rowS_NaCl_clust <-table(rowSums(all_GWAS_clusters[,paste0("is_sign_T",8:14)]))

rowS_NGM <- table(rowSums(all_GWAS[,paste0("is_sign_T",1:7)]))
rowS_NaCl <-table(rowSums(all_GWAS[,paste0("is_sign_T",8:14)]))

### Create a G-matrix from the SNPs ?
var_traits <- NULL
for(i in 1:14){
  if(i %in% 1:7){
    var_traits = c(var_traits, 
                 sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)]],na.rm=TRUE))
  }else{
    var_traits = c(var_traits, 
                   sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af_T14*(1-all_GWAS_clusters$af_T14))[all_GWAS_clusters[,paste0("is_sign_T",i)]],na.rm=TRUE))
  }
}
plot(var_traits[1:7],var_traits[8:14],pch=16,asp=1)
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)


### If we split center/arm ?
### Create a G-matrix from the SNPs ?
var_traits_center <- NULL; var_traits_arm <- NULL

for(i in 1:14){
  if(i %in% 1:7){
    var_traits_center = c(var_traits_center, 
                   sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE))
    var_traits_arm = c(var_traits_arm, 
                          sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE))
  }else{
    var_traits_center = c(var_traits_center, 
                   sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af_T14*(1-all_GWAS_clusters$af_T14))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE))
    var_traits_arm = c(var_traits_arm, 
                          sum((2*(all_GWAS_clusters[,paste0("beta_T",i)])^2*all_GWAS_clusters$af_T14*(1-all_GWAS_clusters$af_T14))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE))
  }
}

plot(var_traits_center,var_traits_arm,pch=16,col=rep(c("blue","red"),each=7))
abline(a=0,b=1); abline(h=0,lty=2); abline(v=0,lty=2) # Much more variance in the arms

####### Covariance between the traits ?

G_mat = matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_arm =matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_center =matrix(rep(NA,196),nrow=14,ncol=14)

for(i in 1:14){
  for(j in 1:14){
    G_mat[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)]],na.rm=TRUE)
    G_mat_center[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE)
    G_mat_arm[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE)
  
  }
}
table(all_GWAS_clusters$chr_type2)

boxplot(all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7)~all_GWAS_clusters$chr_type2)
tapply(all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7),all_GWAS_clusters$chr_type2,mean) # similar frequencies
tapply(abs(all_GWAS_clusters$beta_T7),all_GWAS_clusters$chr_type2,mean) # strong
boxplot(abs(all_GWAS_clusters$beta_T2)~all_GWAS_clusters$chr_type2) # strong

all_GWAS_clusters$varianceT1 <- all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7)*(all_GWAS_clusters$beta_T1)^2
all_GWAS_clusters$varianceT7 <- all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7)*(all_GWAS_clusters$beta_T7)^2

summary(glm((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,family="Gamma"))
summary(glm((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,family="inverse.gaussian"))
boxplot((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,log="y")

a=subset(all_GWAS_clusters,!is.na(all_GWAS_clusters$chr_type2) & varianceT7<0.00051)
plot(density(a$varianceT7[a$chr_type2=='center']),ylim=c(0,40000))
lines(density(a$varianceT7[a$chr_type2=='arm']),col='red')

plot(G_mat_center,G_mat_arm,pch=16)
abline(a=0,b=1); abline(h=0,lty=2); abline(v=0,lty=2) # Much more (co)-variance in the arms

rbind(var_traits[1:7],var_traits[8:14])

plot(cor(all_GWAS_clusters[,paste0("beta_T",1:7)]),
     cor(na.omit(all_GWAS_clusters[,paste0("beta_T",8:14)])))
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)

vect_col=c(adjustcolor("grey",alpha=0.3),adjustcolor("firebrick3",alpha=1))
i=7
plot(-log10(all_GWAS_clusters[,paste0("p_wald_T",i)])~all_GWAS_clusters[,paste0("beta_T",i)],
     pch=16,xlab='Beta coefficient',ylab="-log(p)",
     col=vect_col[1+as.numeric(all_GWAS_clusters[,paste0("is_sign_T",i)])],
     main=paste0('Trait: ',vect_P_traits_bothEnv[i]))

plot((all_GWAS_clusters$beta_T7),(all_GWAS_clusters$beta_T14),pch=16)
plot((all_GWAS_clusters$beta_T7),(all_GWAS_clusters$beta_T14),pch=16,xlim=c(-.1,.1),ylim=c(-.1,.1))


### Heatmap plots between beta/sign ?
heatmap(cor(na.omit(all_GWAS_clusters[,paste0("beta_T",1:14)])),
        labRow=vect_P_traits_bothEnv,
        labCol=vect_P_traits_bothEnv)

## With the non-significant
plot(rowS_NGM/sum(rowS_NGM),type='b',col='cornflowerblue',ylim=c(.00001,1),log="y",xlim=c(0,7))
points(rowS_NaCl/sum(rowS_NaCl),type='b',col='darkorange')
points(rowS_NaCl_clust/sum(rowS_NaCl_clust),type='b',col='darkorange',lty=2)
points(rowS_NGM_clust/sum(rowS_NGM_clust),type='b',col='cornflowerblue',lty=2)

#(rowS_NaCl2/sum(rowS_NaCl2))/(rowS_NGM2/sum(rowS_NGM2))
# 12 times more SNPs with all 7 traits in NaCl
#sum(rowS_NGM2)/sum(rowS_NaCl2) # Same number of clusters

pdf("Lille_07.pdf",w=4,h=4)
plot(rowS_NGM,type='b',col='cornflowerblue',ylim=c(1,1000000),log="y",xlim=c(0,7),las=1,bty="n",xlab="",xaxt="n",ylab="")
axis(1,0:7)
points(rowS_NaCl,type='b',col='darkorange')
points(rowS_NaCl_clust,type='b',col='darkorange',lty=2)
points(rowS_NGM_clust,type='b',col='cornflowerblue',lty=2)
legend(4,1000000,pch=16,col=c("cornflowerblue","darkorange"),c('Low Salt',"High Salt"))
dev.off()

library(gplots)
library(limma) 
vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_T1),as.numeric(all_GWAS_clusters$is_sign_T8))),
  names=c("SF NGM","SF NaCl"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_T7),as.numeric(all_GWAS_clusters$is_sign_T14))),
  names=c("Area low Salt","Area High Salt"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_T7),as.numeric(all_GWAS_clusters$is_sign_T8))),
  names=c("Area low Salt","SF High Salt"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_T7),as.numeric(all_GWAS_clusters$is_sign_T9))),
  names=c("Area low Salt","SB High Salt"),main='All')

vennDiagram(vennCounts(
  all_GWAS_clusters[,paste0("is_sign_T",c(14,11))]),
  names=c("test","Area NGM"),main='All')

### Could there be SNPs in the same genes but not the same variants 
### affecting the traits ? -> would be nice to give examples ????

#### Now the Neff to estimate the effect ~pleiotropy ?
## For each cluster, compute a level of pleiotropy

all_GWAS_clusters_sign <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:14)])>0,])
all_GWAS_clusters_sign_NGM <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:7)])>0,])
all_GWAS_clusters_sign_NaCl <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_T",8:14)])>0,])

nb_sign_clust_NGM <- nrow(all_GWAS_clusters_sign)
nb_sign_clust_NaCl <- nrow(all_GWAS_clusters_sign)
nb_sign_clust <- nrow(all_GWAS_clusters_sign)

pleiot.idx=data.frame(NGM=rep(0,nb_sign_clust),NaCl=rep(0,nb_sign_clust))
pleiot.idx$NGM.cor=0;pleiot.idx$NaCl.cor=0

pleiot.idx_mixed=data.frame(N=rep(0,nb_sign_clust),Neff=rep(0,nb_sign_clust))

#Gmat_NGM = cov(na.omit(all_GWAS_clusters[,paste0("beta_T",1:7)]))
#Gmat_NaCl = cov(na.omit(all_GWAS_clusters[,paste0("beta_T",8:14)]))
#G_mat_all=cov(na.omit(all_GWAS_clusters[,paste0("beta_T",1:14)]))

Gmat_NGM = cov(na.omit(gemma_pheno[,1:7]))
Gmat_NaCl = cov(na.omit(gemma_pheno[,8:14]))
G_mat_all= cov(na.omit(gemma_pheno[,1:14]))


#head(gemma_pheno)

#plot(cov(na.omit(all_GWAS_clusters[,paste0("beta_T",1:7)])),cov(na.omit(gemma_pheno[,1:7])))
#abline(a=0,b=1)
#plot(G_mat[1:7,1:7],cov(na.omit(gemma_pheno[,1:7])))
#abline(a=0,b=1)

for(i in 1:nrow(pleiot.idx_mixed)){
  ## NGM first
  pleiot.idx$NGM[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:7)],na.rm=TRUE)
  if(pleiot.idx$NGM[i]==1){
    pleiot.idx$NGM.cor[i]=1 
  }else if(pleiot.idx$NGM[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:7)])
    std_devs.NGM.temp <- sqrt(diag(Gmat_NGM[idx.sign,idx.sign]))
    cor_matrix.NGM.temp <- Gmat_NGM[idx.sign,idx.sign] / (std_devs.NGM.temp %o% std_devs.NGM.temp)
    pleiot.idx$NGM.cor[i]=pleiot.idx$NGM[i] - var(eigen(cor_matrix.NGM.temp)$values)
  }
  
  ## NaCl second
  pleiot.idx$NaCl[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_T",8:14)],na.rm=TRUE)
  if(pleiot.idx$NaCl[i]==1){
    pleiot.idx$NaCl.cor[i]=1 
  }else if(pleiot.idx$NaCl[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_T",8:14)])
    std_devs.NaCl.temp <- sqrt(diag(Gmat_NaCl[idx.sign,idx.sign]))
    cor_matrix.NaCl.temp <- Gmat_NaCl[idx.sign,idx.sign] / (std_devs.NaCl.temp %o% std_devs.NaCl.temp)
    pleiot.idx$NaCl.cor[i]=pleiot.idx$NaCl[i] - var(eigen(cor_matrix.NaCl.temp)$values)
  }
  
  ### Mixed
  pleiot.idx_mixed$N[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:14)],na.rm=TRUE)
  
  if(pleiot.idx_mixed$N[i]==1){
    pleiot.idx_mixed$Neff[i]=1  
  }else if(pleiot.idx_mixed$N[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_T",1:14)])
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
all_effects <- all_GWAS_clusters_sign[,paste0('beta_T',1:14)]
all_effects_std <- NULL
for(i in 1:14){
  all_effects_std[i] =
    sd(all_GWAS_clusters_sign[all_GWAS_clusters_sign[,paste0("is_sign_T",i)],paste0("beta_T",i)])
}

all_GWAS_clusters_sign_std <- all_GWAS_clusters_sign
for(i in 1:14){
  all_GWAS_clusters_sign_std[,paste0('beta_T',i)]=all_GWAS_clusters_sign_std[,paste0('beta_T',i)]/all_effects_std[i]
}

 euclidian_summed <- NULL;euclidian_NGM <- NULL;euclidian_NaCl <- NULL
 manh_summed <- NULL; manh_NGM <- NULL; manh_NaCl <- NULL
 for(i in 1:nrow(all_GWAS_clusters_sign_std)){
   euclidian_summed <- c(euclidian_summed,
                         sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:14)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:14)])]^2)))
   euclidian_NGM[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:7)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:7)])]^2))
   euclidian_NaCl[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',8:14)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',8:14)])]^2))
   
   manh_summed <- c(manh_summed,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:14)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:14)])])))
   manh_NGM=c(manh_NGM,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',1:7)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',1:7)])])))
   manh_NaCl=c(manh_NaCl,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_T',8:14)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_T',8:14)])])))
   
 }

 
 ### With nls - replace all
 pdf("plots_BLUEs/Scaling_pleiotropy_Neff_BLUE_Arm_Centers.pdf",h=10,w=12)
 
 is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
 is_center <- all_GWAS_clusters_sign$chr_type2=="center"
 
 temp=pleiot.idx_mixed$Neff
 coefs_euc <- coef(nls(euclidian_summed~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$Neff[is_arm]
 coefs_euc_arm <- coef(nls(euclidian_summed[is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$Neff[is_center]
 coefs_euc_center <- coef(nls(euclidian_summed[is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 temp.NGM <- euclidian_NGM>0
 temp=pleiot.idx$NGM.cor[temp.NGM]
 coefs_euc.NGM <- coef(nls(euclidian_NGM[temp.NGM]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_arm]
 ## a=0 here
 coefs_euc.NGM.arm <- c(0,coef(nls(euclidian_NGM[temp.NGM & is_arm]~b*temp^c,start=list(b=0.8,c=.5))))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_center]
 coefs_euc.NGM.center <- coef(nls(euclidian_NGM[temp.NGM & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 temp.NaCl <- euclidian_NaCl>0
 temp=pleiot.idx$NaCl.cor[temp.NaCl]
 coefs_euc.NaCl <- coef(nls(euclidian_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_euc.NaCl.arm <- coef(nls(euclidian_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0.0,b=0.2,c=1.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_euc.NaCl.center <- coef(nls(euclidian_summed[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 #par(mfrow=c(2,3))
 plot(euclidian_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,11),ylim=c(0,14),type='n')
 points(euclidian_summed[is_arm]~pleiot.idx_mixed$Neff[is_arm],pch=16,col='green')
 points(euclidian_summed[is_center]~pleiot.idx_mixed$Neff[is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc[1]+(coefs_euc[2])*(seq(0,12,0.1)^coefs_euc[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc_arm[1]+(coefs_euc_arm[2])*(seq(0,12,0.1)^coefs_euc_arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc_center[1]+(coefs_euc_center[2])*(seq(0,12,0.1)^coefs_euc_center[3])),col="lightblue")
 
 abline(lm(euclidian_summed~pleiot.idx_mixed$Neff))
 text(2,12,paste0("c=",round(coefs_euc[3],digits=3)))
 text(2,11,paste0("c=",round(coefs_euc_arm[3],digits=3)),col="green")
 text(2,10,paste0("c=",round(coefs_euc_center[3],digits=3)),col='lightblue')
 
 plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NGM',xlim=c(0,7),ylim=c(0,10),
      type='n')
 points(euclidian_NGM[temp.NGM & is_arm]~pleiot.idx$NGM.cor[temp.NGM & is_arm],pch=16,col='green')
 points(euclidian_NGM[temp.NGM & is_center]~pleiot.idx$NGM.cor[temp.NGM & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NGM[1]+(coefs_euc.NGM[2])*(seq(0,12,0.1)^coefs_euc.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NGM.arm[1]+(coefs_euc.NGM.arm[2])*(seq(0,12,0.1)^coefs_euc.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NGM.center[1]+(coefs_euc.NGM.center[2])*(seq(0,12,0.1)^coefs_euc.NGM.center[3])),col="lightblue")
 
 abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]))
 text(1.5,7.5,paste0("c=",round(coefs_euc.NGM[3],digits=3)))
 text(1.5,6.5,paste0("c=",round(coefs_euc.NGM.arm[3],digits=3)),col="green")
 text(1.5,5.5,paste0("c=",round(coefs_euc.NGM.center[3],digits=3)),col="lightblue")
 
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
 
 temp=pleiot.idx_mixed$Neff
 coefs_man <- coef(nls(manh_summed~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx_mixed$Neff[is_arm]
 coefs_man.arm <- coef(nls(manh_summed[is_arm]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx_mixed$Neff[is_center]
 coefs_man.center <- coef(nls(manh_summed[is_center]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
 temp.NGM <- manh_NGM>0
 temp=pleiot.idx$NGM.cor[temp.NGM]
 coefs_man.NGM <- coef(nls(manh_NGM[temp.NGM]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_arm]
 coefs_man.NGM.arm <- coef(nls(manh_NGM[temp.NGM & is_arm]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx$NGM.cor[temp.NGM & is_center]
 coefs_man.NGM.center <- coef(nls(manh_NGM[temp.NGM & is_center]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
 temp.NaCl <- manh_NaCl>0
 temp=pleiot.idx$NaCl.cor[temp.NaCl]
 coefs_man.NaCl <- coef(nls(manh_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_arm]
 coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 
 plot(manh_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TM ',xlab="Neff",xlim=c(0,11),ylim=c(0,40),type='n')
 points(manh_summed[is_arm]~pleiot.idx_mixed$Neff[is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_summed[is_center]~pleiot.idx_mixed$Neff[is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man[1]+(coefs_man[2])*(seq(0,12,0.1)^coefs_man[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.arm[1]+(coefs_man.arm[2])*(seq(0,12,0.1)^coefs_man.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.center[1]+(coefs_man.center[2])*(seq(0,12,0.1)^coefs_man.center[3])),col="lightblue")
 
 abline(lm(manh_summed~pleiot.idx_mixed$Neff))
 text(2,30,paste0("c=",round(coefs_man[3],digits=3)))
 text(2,28,paste0("c=",round(coefs_man.arm[3],digits=3)),col='green')
 text(2,26,paste0("c=",round(coefs_man.center[3],digits=3)),col='lightblue')
 
 plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(1,7),ylim=c(0,25),type='n')
 points(manh_NGM[temp.NGM & is_arm]~pleiot.idx$NGM.cor[temp.NGM & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NGM[temp.NGM & is_center]~pleiot.idx$NGM.cor[temp.NGM & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NGM[1]+(coefs_man.NGM[2])*(seq(0,12,0.1)^coefs_man.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NGM.arm[1]+(coefs_man.NGM.arm[2])*(seq(0,12,0.1)^coefs_man.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NGM.center[1]+(coefs_man.NGM.center[2])*(seq(0,12,0.1)^coefs_man.NGM.center[3])),col="lightblue")
 
 abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]))
 text(1.5,20,paste0("c=",round(coefs_man.NGM[3],digits=3)))
 text(1.5,18,paste0("c=",round(coefs_man.NGM.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NGM.center[3],digits=3)),col='lightblue')
 
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
 
 
 ### With N instead of Neff - do a linear model ?
 
 vect_chr_type <- all_GWAS_clusters_sign_std$chr_type2
 model_pleiotropie <- lm(log(manh_summed)~log(pleiot.idx_mixed$Neff)*vect_chr_type)
 model_pleiotropie2 <- lm(log(manh_summed)~log(pleiot.idx_mixed$Neff)+vect_chr_type)
 anova(model_pleiotropie,model_pleiotropie2)
 
 
 dev.off()
 

# 
# #################
# 
# ### The same with N ? ##
 ### With nls - replace all
 pdf("plots_BLUEs/Scaling_pleiotropy_N_BLUE_Arm_Centers.pdf",h=10,w=12)
 
 is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
 is_center <- all_GWAS_clusters_sign$chr_type2=="center"
 
 temp=pleiot.idx_mixed$N
 coefs_euc <- coef(nls(euclidian_summed~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$N[is_arm]
 coefs_euc_arm <- coef(nls(euclidian_summed[is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$N[is_center]
 coefs_euc_center <- coef(nls(euclidian_summed[is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 temp.NGM <- euclidian_NGM>0
 temp=pleiot.idx$NGM[temp.NGM]
 coefs_euc.NGM <- coef(nls(euclidian_NGM[temp.NGM]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NGM[temp.NGM & is_arm]
 ## a=0 here
 coefs_euc.NGM.arm <- c(0,coef(nls(euclidian_NGM[temp.NGM & is_arm]~b*temp^c,start=list(b=0.8,c=.5))))
 temp=pleiot.idx$NGM[temp.NGM & is_center]
 coefs_euc.NGM.center <- coef(nls(euclidian_NGM[temp.NGM & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 temp.NaCl <- euclidian_NaCl>0
 temp=pleiot.idx$NaCl[temp.NaCl]
 coefs_euc.NaCl <- coef(nls(euclidian_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_arm]
 coefs_euc.NaCl.arm <- coef(nls(euclidian_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_center]
 coefs_euc.NaCl.center <- coef(nls(euclidian_summed[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 #par(mfrow=c(2,3))
 plot(euclidian_summed~pleiot.idx_mixed$N,pch=16,ylab='TE (euc. distance)',xlab="N",xlim=c(0,11),ylim=c(0,14),type='n')
 points(euclidian_summed[is_arm]~pleiot.idx_mixed$N[is_arm],pch=16,col='green')
 points(euclidian_summed[is_center]~pleiot.idx_mixed$N[is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc[1]+(coefs_euc[2])*(seq(0,12,0.1)^coefs_euc[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc_arm[1]+(coefs_euc_arm[2])*(seq(0,12,0.1)^coefs_euc_arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc_center[1]+(coefs_euc_center[2])*(seq(0,12,0.1)^coefs_euc_center[3])),col="lightblue")
 
 abline(lm(euclidian_summed~pleiot.idx_mixed$N))
 text(2,12,paste0("c=",round(coefs_euc[3],digits=3)),col='red')
 text(2,11,paste0("c=",round(coefs_euc_arm[3],digits=3)),col="green")
 text(2,10,paste0("c=",round(coefs_euc_center[3],digits=3)),col='lightblue')
 
 plot(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM],pch=16,ylab='TE (euc. distance)',xlab="N",main='NGM',xlim=c(0,7),ylim=c(0,10),
      type='n')
 points(euclidian_NGM[temp.NGM & is_arm]~pleiot.idx$NGM[temp.NGM & is_arm],pch=16,col='green')
 points(euclidian_NGM[temp.NGM & is_center]~pleiot.idx$NGM[temp.NGM & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NGM[1]+(coefs_euc.NGM[2])*(seq(0,12,0.1)^coefs_euc.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NGM.arm[1]+(coefs_euc.NGM.arm[2])*(seq(0,12,0.1)^coefs_euc.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NGM.center[1]+(coefs_euc.NGM.center[2])*(seq(0,12,0.1)^coefs_euc.NGM.center[3])),col="lightblue")
 
 abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM]))
 text(1.5,7.5,paste0("c=",round(coefs_euc.NGM[3],digits=3)),col='red')
 text(1.5,6.5,paste0("c=",round(coefs_euc.NGM.arm[3],digits=3)),col="green")
 text(1.5,5.5,paste0("c=",round(coefs_euc.NGM.center[3],digits=3)),col="lightblue")
 
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
 
 temp=pleiot.idx_mixed$N
 coefs_man <- coef(nls(manh_summed~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx_mixed$N[is_arm]
 coefs_man.arm <- coef(nls(manh_summed[is_arm]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx_mixed$N[is_center]
 coefs_man.center <- coef(nls(manh_summed[is_center]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
 temp.NGM <- manh_NGM>0
 temp=pleiot.idx$NGM[temp.NGM]
 coefs_man.NGM <- coef(nls(manh_NGM[temp.NGM]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx$NGM[temp.NGM & is_arm]
 coefs_man.NGM.arm <- coef(nls(manh_NGM[temp.NGM & is_arm]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 temp=pleiot.idx$NGM[temp.NGM & is_center]
 coefs_man.NGM.center <- coef(nls(manh_NGM[temp.NGM & is_center]~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
 temp.NaCl <- manh_NaCl>0
 temp=pleiot.idx$NaCl[temp.NaCl]
 coefs_man.NaCl <- coef(nls(manh_NaCl[temp.NaCl]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_arm]
 coefs_man.NaCl.arm <- coef(nls(manh_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 temp=pleiot.idx$NaCl[temp.NaCl & is_center]
 coefs_man.NaCl.center <- coef(nls(manh_NaCl[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=1,c=1)))
 
 plot(manh_summed~pleiot.idx_mixed$N,pch=16,ylab='TM ',xlab="N",xlim=c(0,11),ylim=c(0,40),type='n')
 points(manh_summed[is_arm]~pleiot.idx_mixed$N[is_arm],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('green',alpha=.6))
 points(manh_summed[is_center]~pleiot.idx_mixed$N[is_center],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man[1]+(coefs_man[2])*(seq(0,12,0.1)^coefs_man[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.arm[1]+(coefs_man.arm[2])*(seq(0,12,0.1)^coefs_man.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.center[1]+(coefs_man.center[2])*(seq(0,12,0.1)^coefs_man.center[3])),col="lightblue")
 
 abline(lm(manh_summed~pleiot.idx_mixed$N))
 text(2,30,paste0("c=",round(coefs_man[3],digits=3)),col="red")
 text(2,28,paste0("c=",round(coefs_man.arm[3],digits=3)),col='green')
 text(2,26,paste0("c=",round(coefs_man.center[3],digits=3)),col='lightblue')
 
 plot(manh_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM],pch=16,ylab='TM',xlab="N",main='NGM',xlim=c(1,7),ylim=c(0,25),type='n')
 points(manh_NGM[temp.NGM & is_arm]~pleiot.idx$NGM[temp.NGM & is_arm],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('green',alpha=.6))
 points(manh_NGM[temp.NGM & is_center]~pleiot.idx$NGM[temp.NGM & is_center],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NGM[1]+(coefs_man.NGM[2])*(seq(0,12,0.1)^coefs_man.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NGM.arm[1]+(coefs_man.NGM.arm[2])*(seq(0,12,0.1)^coefs_man.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NGM.center[1]+(coefs_man.NGM.center[2])*(seq(0,12,0.1)^coefs_man.NGM.center[3])),col="lightblue")
 
 abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM]))
 text(1.5,20,paste0("c=",round(coefs_man.NGM[3],digits=3)),col="red")
 text(1.5,18,paste0("c=",round(coefs_man.NGM.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NGM.center[3],digits=3)),col='lightblue')
 
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
 
 ### Create a G-matrix from the SNPs - depending on pleiotropy ?
 sign_total <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_T",1:14)])
G_mat = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center =matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_P1 = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm_P1 =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center_P1 =matrix(rep(NA,196),nrow=14,ncol=14)
G_mat_P2 = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm_P2 =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center_P2 =matrix(rep(NA,196),nrow=14,ncol=14) 
G_mat_P3 = matrix(rep(NA,196),nrow=14,ncol=14); G_mat_arm_P3 =matrix(rep(NA,196),nrow=14,ncol=14); G_mat_center_P3 =matrix(rep(NA,196),nrow=14,ncol=14) 
vect_P1 <- sign_total==1 & all_GWAS_clusters$clust_length<10000
vect_P2 <- sign_total%in%c(2,3) & all_GWAS_clusters$clust_length<10000
vect_P3 <- sign_total>3 & all_GWAS_clusters$clust_length<10000


 for(i in 1:14){
   for(j in 1:14){
     G_mat[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE)
     G_mat_arm[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_T",i)]*all_GWAS_clusters[,paste0("beta_T",j)]*all_GWAS_clusters$af_T7*(1-all_GWAS_clusters$af_T7))[all_GWAS_clusters[,paste0("is_sign_T",i)] & all_GWAS_clusters[,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE)
     

     G_mat_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_T",i)]*all_GWAS_clusters[vect_P1,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P1]*(1-all_GWAS_clusters$af_T7[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P1,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_T",i)]*all_GWAS_clusters[vect_P1,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P1]*(1-all_GWAS_clusters$af_T7[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P1,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P1]=="center"],na.rm=TRUE)
     G_mat_arm_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_T",i)]*all_GWAS_clusters[vect_P1,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P1]*(1-all_GWAS_clusters$af_T7[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P1,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P1]=="arm"],na.rm=TRUE)
     
     G_mat_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_T",i)]*all_GWAS_clusters[vect_P2,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P2]*(1-all_GWAS_clusters$af_T7[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P2,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_T",i)]*all_GWAS_clusters[vect_P2,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P2]*(1-all_GWAS_clusters$af_T7[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P2,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P2]=="center"],na.rm=TRUE)
     G_mat_arm_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_T",i)]*all_GWAS_clusters[vect_P2,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P2]*(1-all_GWAS_clusters$af_T7[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P2,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P2]=="arm"],na.rm=TRUE)
     
     G_mat_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_T",i)]*all_GWAS_clusters[vect_P3,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P3]*(1-all_GWAS_clusters$af_T7[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P3,paste0("is_sign_T",j)]],na.rm=TRUE)
     G_mat_center_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_T",i)]*all_GWAS_clusters[vect_P3,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P3]*(1-all_GWAS_clusters$af_T7[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P3,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P3]=="center"],na.rm=TRUE)
     G_mat_arm_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_T",i)]*all_GWAS_clusters[vect_P3,paste0("beta_T",j)]*all_GWAS_clusters$af_T7[vect_P3]*(1-all_GWAS_clusters$af_T7[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_T",i)] & all_GWAS_clusters[vect_P3,paste0("is_sign_T",j)] & all_GWAS_clusters$chr_type2[vect_P3]=="arm"],na.rm=TRUE)
     
        }
 }
 
 
plot(scale(eigen(G_mat)$vectors[,1]),scale(eigen(G_mat_P1)$vectors[,1]))

 plot(G_mat[1:7,1:7],G_mat_P1[1:7,1:7])
 
 allPcor=NULL
 
 for(i in 1:13){
   for(j in (i+1):14){
     allPcor=rbind(allPcor,
     data.frame(T1=i,T2=j,Pcor=G_mat_all[i,j],nbClust=sum(rowSums(all_GWAS_clusters_sign[,paste0("is_sign_T",c(i,j))])==2),
                nbClust_arm=sum(rowSums(all_GWAS_clusters_sign[all_GWAS_clusters_sign$chr_type2=="arm",paste0("is_sign_T",c(i,j))])==2),
                nbClust_center=sum(rowSums(all_GWAS_clusters_sign[all_GWAS_clusters_sign$chr_type2=="center",paste0("is_sign_T",c(i,j))])==2))
     )
   }
 }
 allPcor$same_env = 0
 allPcor$same_env[allPcor$T1%in%(1:7) & allPcor$T2%in%(1:7)]=1
 allPcor$same_env[allPcor$T1%in%(8:14) & allPcor$T2%in%(1:14)]=2
 
 allPcor$same_T = 0
 allPcor$same_T[allPcor$T1==(allPcor$T2-7)]=1
 
 
 par(mfrow=c(1,3))
 plot(abs(allPcor$Pcor),allPcor$nbClust,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 plot(abs(allPcor$Pcor),allPcor$nbClust_arm,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 plot(abs(allPcor$Pcor),allPcor$nbClust_center,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 
 pdf("Lille_05.pdf",w=4,h=4)
 vect_kept=allPcor$same_T==0 & allPcor$same_env>0
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.4),ylim=c(0,310),las=1,bty="n")
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2)
 legend(0,300,c("Low Salt","High Salt"),pch=16,col=c("blue","orange"))
 dev.off()
 
 pdf("Lille_06.pdf",w=4,h=4)
 vect_kept=allPcor$same_T==0
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.4),ylim=c(0,310),las=1,bty="n")
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2)
 legend(0,320,c("Low Salt","High Salt","Across environment"),pch=16,col=c("blue","orange","black"),cex=.75)
 dev.off()
 
 
 pdf("Lille_07.pdf",w=4,h=4)
 vect_kept=allPcor$same_T==0
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.4),ylim=c(0,310),las=1,bty="n")
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2)
 
 vect_kept=allPcor$same_T==1
 points(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col="red",xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.4),ylim=c(0,310))
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2,col="red")

 
 legend(0,320,c("Low Salt","High Salt","Across environment","Same trait"),pch=c(16,16,16,11),col=c("blue","orange","black","red"),cex=.6,ncol=2)
 dev.off()
 
 
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.4),ylim=c(0,310))
 
 
 
 model_sameT <- lm(nbClust~abs(Pcor),data=subset(allPcor,same_T==1))
 summary(model_sameT)
 
 allPcor$same_env=as.factor(allPcor$same_env)
 allPcor$same_T=as.factor(allPcor$same_T)
 
 allPcor$same_env_NGM <- allPcor$same_env==1
 allPcor$same_env_NaCl <- allPcor$same_env==2
 
 model1 <- lm(allPcor$nbClust~abs(allPcor$Pcor)*(allPcor$same_T+allPcor$same_env_NGM+allPcor$same_env_NaCl))
 summary(model1)
 drop1(model1,test='F')
 
 model1.1 <- lm(allPcor$nbClust~abs(allPcor$Pcor)*allPcor$same_T-allPcor$same_T)
 summary(model1.1)
 drop1(model1,test='F')
 
 
 
 
 model1 <- lm(allPcor$nbClust~abs(allPcor$Pcor)*(allPcor$same_T+allPcor$same_env_NGM+allPcor$same_env_NaCl)-abs(allPcor$Pcor):allPcor$same_T )
 summary(model1)
 drop1(model1,test='F')
 anova(model1,model2)
 
 model <- lm(allPcor$nbClust~abs(allPcor$Pcor)*(allPcor$same_T+allPcor$same_env_NGM+allPcor$same_env_NaCl)-abs(allPcor$Pcor):allPcor$same_T-allPcor$same_T  )
 summary(model)
 drop1(model,test='F')
 
 
 
 model2 <- lm(allPcor$nbClust~abs(allPcor$Pcor)*(allPcor$same_env_NGM+allPcor$same_env_NaCl))
 summary(model2)
 drop1(model2,test='F')
 anova(model1,model2,test="F")
 
 model3=model2
 coef(model3)
 
allPcor$same_env=as.numeric(allPcor$same_env)-1
 plot(abs(allPcor$Pcor),allPcor$nbClust,pch=c(16,11)[as.numeric(allPcor$same_T)],
      col=c("black","blue","orange")[2+as.numeric(allPcor$same_env)],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 abline(a=coef(model3)[1]+coef(model3)[3],b=coef(model3)[2]+coef(model3)[4])
 abline(a=coef(model3)[1],b=coef(model3)[2],col="blue")
 abline(a=(coef(model2)[1]+coef(model2)[4]),b=coef(model2)[2],col="orange")
 
 
 model2_arm <- lm(allPcor$nbClust_arm~abs(allPcor$Pcor)*(allPcor$same_env_NGM+allPcor$same_env_NaCl)) 
 drop1(model2_arm,test='F')
 model2_center <- lm(allPcor$nbClust_center~abs(allPcor$Pcor)*(allPcor$same_env_NGM+allPcor$same_env_NaCl)) 
 drop1(model2_center,test='F')
 
 
 
 