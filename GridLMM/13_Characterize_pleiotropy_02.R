rm(list=ls());gc()
library(data.table)
library(lme4)

threshold_clust = 0.2
load(paste0("RData/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
load(file="GridLMM_results_CHR/All_Gwas.RData")

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv)])
all_GWAS_clusters$sign_NGM <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv[1:8])])
all_GWAS_clusters$sign_NaCl <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv[9:16])])

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(paste0("txt_BLUEs_GridLMM_CHR//all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
  #read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

df_clusters=read.table(paste0("data_BLUE_clusters_CHR/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
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
## Only 37 with 2 types
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
#temp.af.NGM=1-all_GWAS_clusters$af.NGM
#temp.af.NaCl=1-all_GWAS_clusters$af.NaCl

#all_GWAS_clusters$af.NGM[!is.na(temp.af.NGM) & all_GWAS_clusters$af.NGM>temp.af.NGM]=temp.af.NGM[!is.na(temp.af.NGM) & all_GWAS_clusters$af.NGM>temp.af.NGM]
#all_GWAS_clusters$af.NaCl[!is.na(temp.af.NaCl) & all_GWAS_clusters$af.NaCl>temp.af.NaCl]=temp.af.NaCl[!is.na(temp.af.NaCl) & all_GWAS_clusters$af.NaCl>temp.af.NaCl]
#plot((all_GWAS_clusters$af.NGM),(all_GWAS_clusters$af.NaCl),pch=16)

#boxplot(all_GWAS_clusters$af.NGM~all_GWAS_clusters$sign_NGM)
#boxplot(all_GWAS_clusters$af.NaCl~all_GWAS_clusters$sign_NaCl)

#tapply(all_GWAS_clusters$af.NGM,all_GWAS_clusters$sign_NGM,mean)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]


hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+clust_length+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic

boxplot(all_GWAS_clusters$clust_length+1~all_GWAS_clusters$sign_overall,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NGM,log="y")
boxplot(all_GWAS_clusters$nb_SNPs~all_GWAS_clusters$sign_NaCl,log="y")

### Add the pleiotropy
#dim(merge(all_GWAS_sign,unique(all_GWAS_clusters[,c("group","sign_overall","sign_NaCl","sign_NGM")])))

#pheno_NGM=read.table("Breeding_values_NGM.txt",sep='\t',h=TRUE)
#pheno_NaCl=read.table("Breeding_values_NaCl.txt",sep='\t',h=TRUE)
#gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)

# Clusters
sign_all_df_clust <- all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv)]
table(rowSums(sign_all_df_clust)>=1) # ~737 Clusters / 10 are no longer significant
all_GWAS_clusters$sign_overall <- rowSums(sign_all_df_clust)>=1

for(i in 1:16){
  all_GWAS[,paste0("is_sign_",vect_P_traits_bothEnv[i])] = p.adjust(all_GWAS[,paste0("pvalue_T",i)],method="BH")<0.2
  print(table(all_GWAS[,paste0("is_sign_",vect_P_traits_bothEnv[i])]))
}

# SNPs
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_T",1:16)]
table(rowSums(sign_all_df)>=1) # 72k SNPs

all_GWAS$sign_overall <- rowSums(sign_all_df)>=1

rowS_NGM_clust <- table(rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[1:8])]))
rowS_NaCl_clust <-table(rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[9:16])]))

rowS_NGM <- table(rowSums(all_GWAS[,paste0("is_sign_",vect_P_traits_bothEnv[1:8])]))
rowS_NaCl <-table(rowSums(all_GWAS[,paste0("is_sign_",vect_P_traits_bothEnv[9:16])]))

### Create a G-matrix from the SNPs ?
#var_traits <- NULL
#for(i in 1:16){
#  if(i %in% 1:8){
#    var_traits = c(var_traits, 
#                 sum((2*(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])])^2*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])]],na.rm=TRUE))
#  }else{
#    var_traits = c(var_traits, 
#                   sum((2*(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])])^2*all_GWAS_clusters$af.NaCl*(1-all_GWAS_clusters$af.NaCl))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i]))]],na.rm=TRUE))
#  }
#}
#plot(var_traits[1:7],var_traits[8:14],pch=16,asp=1)
#abline(a=0,b=1);abline(h=0,lty=2);abline(v=0,lty=2)


### If we split center/arm ?
### Create a G-matrix from the SNPs ?
#var_traits_center <- NULL; var_traits_arm <- NULL

#for(i in 1:14){
#  if(i %in% 1:7){
#    var_traits_center = c(var_traits_center, 
#                   sum((2*(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])])^2*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE))
#    var_traits_arm = c(var_traits_arm, 
#                          sum((2*(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])])^2*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE))
#  }else{
#    var_traits_center = c(var_traits_center, 
#                   sum((2*(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])])^2*all_GWAS_clusters$af.NaCl*(1-all_GWAS_clusters$af.NaCl))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE))
#    var_traits_arm = c(var_traits_arm, 
#                          sum((2*(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])])^2*all_GWAS_clusters$af.NaCl*(1-all_GWAS_clusters$af.NaCl))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE))
#  }
#}

#plot(var_traits_center,var_traits_arm,pch=16,col=rep(c("blue","red"),each=7))
#abline(a=0,b=1); abline(h=0,lty=2); abline(v=0,lty=2) # Much more variance in the arms

####### Covariance between the traits ?

#G_mat = matrix(rep(NA,196),nrow=14,ncol=14)
#G_mat_arm =matrix(rep(NA,196),nrow=14,ncol=14)
#G_mat_center =matrix(rep(NA,196),nrow=14,ncol=14)

#for(i in 1:16){
#  for(j in 1:16){
#    G_mat[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[j])]],na.rm=TRUE)
#    G_mat_center[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE)
#    G_mat_arm[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE)
#  
#  }
#}
#table(all_GWAS_clusters$chr_type2)

#boxplot(all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM)~all_GWAS_clusters$chr_type2)
#tapply(all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM),all_GWAS_clusters$chr_type2,mean) # similar frequencies
#tapply(abs(all_GWAS_clusters$beta_T7),all_GWAS_clusters$chr_type2,mean) # strong
#boxplot(abs(all_GWAS_clusters$beta_T2)~all_GWAS_clusters$chr_type2) # strong

#all_GWAS_clusters$varianceT1 <- all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM)*(all_GWAS_clusters$beta_T1)^2
#all_GWAS_clusters$varianceT7 <- all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM)*(all_GWAS_clusters$beta_T7)^2
#summary(glm((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,family="Gamma"))
#summary(glm((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,family="inverse.gaussian"))
#boxplot((all_GWAS_clusters$varianceT7)~all_GWAS_clusters$chr_type2,log="y")

#a=subset(all_GWAS_clusters,!is.na(all_GWAS_clusters$chr_type2) & varianceT7<0.00051)
#plot(density(a$varianceT7[a$chr_type2=='center']),ylim=c(0,40000))
#lines(density(a$varianceT7[a$chr_type2=='arm']),col='red')

#plot(G_mat_center,G_mat_arm,pch=16)
#abline(a=0,b=1); abline(h=0,lty=2); abline(v=0,lty=2) # Much more (co)-variance in the arms
#rbind(var_traits[1:7],var_traits[8:14])

plot(cor(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[1:8])]),
     cor(na.omit(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[9:16])])))
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)

vect_col=c(adjustcolor("grey",alpha=0.3),adjustcolor("firebrick3",alpha=1))
i=7
plot(-log10(all_GWAS_clusters[,paste0("p_value_REML.cor_",vect_P_traits_bothEnv[i])])~all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])],
     pch=16,xlab='Beta coefficient',ylab="-log(p)",
     col=vect_col[1+as.numeric(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])])],
     main=paste0('Trait: ',vect_P_traits_bothEnv[i]))

plot((all_GWAS_clusters$beta_SF_NGM),(all_GWAS_clusters$beta_area.F_NaCl),pch=16,xlim=c(-.1,.1),ylim=c(-.1,.1))
plot((all_GWAS_clusters$beta_area.F_NGM),(all_GWAS_clusters$beta_area.F_NaCl),pch=16,xlim=c(-.1,.1),ylim=c(-.1,.1))

### Heatmap plots between beta/sign ?
heatmap(cor(na.omit(all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv)])),
        labRow=vect_P_traits_bothEnv,
        labCol=vect_P_traits_bothEnv)

## With the non-significant
plot(rowS_NGM/sum(rowS_NGM),type='b',col='cornflowerblue',ylim=c(.00000001,1),log="y",xlim=c(0,8))
points(rowS_NaCl/sum(rowS_NaCl),type='b',col='darkorange')
points(rowS_NaCl_clust/sum(rowS_NaCl_clust),type='b',col='darkorange',lty=2)
points(rowS_NGM_clust/sum(rowS_NGM_clust),type='b',col='cornflowerblue',lty=2)

#(rowS_NaCl2/sum(rowS_NaCl2))/(rowS_NGM2/sum(rowS_NGM2))
# 12 times more SNPs with all 7 traits in NaCl
#sum(rowS_NGM2)/sum(rowS_NaCl2) # Same number of clusters

pdf("Lille_07_R2.pdf",w=4,h=4)
plot(rowS_NGM,type='b',col='cornflowerblue',ylim=c(1,1000000),log="y",xlim=c(0,8),las=1,bty="n",xlab="",xaxt="n",ylab="")
axis(1,0:8)
points(rowS_NaCl,type='b',col='darkorange')
points(rowS_NaCl_clust,type='b',col='darkorange',lty=2)
points(rowS_NGM_clust,type='b',col='cornflowerblue',lty=2)
legend(4,1000000,pch=16,col=c("cornflowerblue","darkorange"),c('Low Salt',"High Salt"))
dev.off()


library(gplots)
library(limma) 
vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_SF_NGM),as.numeric(all_GWAS_clusters$is_sign_SF_NaCl))),
  names=c("SF NGM","SF NaCl"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_area.F_NGM),as.numeric(all_GWAS_clusters$is_sign_area.F_NaCl))),
  names=c("Area low Salt","Area High Salt"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_area.F_NGM),as.numeric(all_GWAS_clusters$is_sign_SF_NaCl))),
  names=c("Area low Salt","SF High Salt"),main='All')

vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_area.F_NGM),as.numeric(all_GWAS_clusters$is_sign_SB_NaCl))),
  names=c("Area low Salt","SB High Salt"),main='All')


### Could there be SNPs in the same genes but not the same variants
### affecting the traits ? -> would be nice to give examples ????


### ### ### ###  Can be reloaded easily
#glmm.env=new.env()
#load('RData/Preparatory_GLMM.RData',envir=glmm.env)
# data_NGM=glmm.env$all_BLUEs_genotyped_NGM
# data_NaCl=glmm.env$all_BLUEs_genotyped_NaCl
# 
# newdata_NGM=data.frame(logD=0,rel_humidity=0, temperature=0,
#                    pop_label=sort(unique(c(data_NGM$pop_label))))
# newdata_NaCl=data.frame(logD=0,rel_humidity=0, temperature=0,
#                        pop_label=sort(unique(c(data_NaCl$pop_label))))
# 
# all_traits_corrected=data.frame(pop_label=sort(unique(c(newdata_NGM$pop_label,newdata_NaCl$pop_label))))
# 
# for(i in vect_P_traits){
# print(i)
# data_NGM$y=scale(data_NGM[,i])
# data_NaCl$y=scale(data_NaCl[,i])
# trait_lmer <- lmer(y~ logD + rel_humidity + temperature  + pop_label + (1|date_str) + (1|session_id), data=data_NGM)
# all_traits_corrected.temp <- data.frame(pop_label=newdata_NGM$pop_label,trait=predict(trait_lmer,newdata_NGM,re.form=NA))
# names(all_traits_corrected.temp)[1]='pop_label'
# names(all_traits_corrected.temp)[2]=paste0(i,"_NGM")
# all_traits_corrected=merge(all_traits_corrected,all_traits_corrected.temp,all=TRUE)
# 
# trait_lmer <- lmer(y~ logD + rel_humidity + temperature  + pop_label + (1|date_str) + (1|session_id), data=data_NaCl)
# newdata_NaCl_temp=subset(newdata_NaCl,pop_label%in%subset(data_NaCl,!is.na(data_NaCl$y))$pop_label)
# all_traits_corrected.temp <- data.frame(pop_label=newdata_NaCl_temp$pop_label,trait=predict(trait_lmer,newdata_NaCl_temp,re.form=NA))
# names(all_traits_corrected.temp)[1]='pop_label'
# names(all_traits_corrected.temp)[2]=paste0(i,"_NaCl")
# all_traits_corrected=merge(all_traits_corrected,all_traits_corrected.temp,all=TRUE)
# }
# dim(all_traits_corrected)
# 
# Gmat_NGM = cov(na.omit(all_traits_corrected[,vect_P_traits_bothEnv[1:8]]))
# Gmat_NaCl = cov(na.omit(all_traits_corrected[,vect_P_traits_bothEnv[9:16]]))
# G_mat_all= cov(na.omit(all_traits_corrected[,vect_P_traits_bothEnv]))
# 
# save(list=c("all_traits_corrected","Gmat_NGM","Gmat_NaCl","G_mat_all"),file='RData/Gmatrices_P.RData')

load('RData/Gmatrices_P.RData')
heatmap(G_mat_all)

#### Now the Neff to estimate the effect ~pleiotropy ?
## For each cluster, compute a level of pleiotropy

all_GWAS_clusters_sign <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv)])>0,])
all_GWAS_clusters_sign_NGM <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[1:8])])>0,])
all_GWAS_clusters_sign_NaCl <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[9:16])])>0,])

nb_sign_clust_NGM <- nrow(all_GWAS_clusters_sign)
nb_sign_clust_NaCl <- nrow(all_GWAS_clusters_sign)
nb_sign_clust <- nrow(all_GWAS_clusters_sign)

pleiot.idx=data.frame(NGM=rep(0,nb_sign_clust),NaCl=rep(0,nb_sign_clust))
pleiot.idx$NGM.cor=0;pleiot.idx$NaCl.cor=0

pleiot.idx_mixed=data.frame(N=rep(0,nb_sign_clust),Neff=rep(0,nb_sign_clust))

for(i in 1:nrow(pleiot.idx_mixed)){

  ## NGM first
  pleiot.idx$NGM[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits_bothEnv[1:8])],na.rm=TRUE)
  
  if(pleiot.idx$NGM[i]==1){
    pleiot.idx$NGM.cor[i]=1 
  }else if(pleiot.idx$NGM[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits_bothEnv[1:8])])
    std_devs.NGM.temp <- sqrt(diag(Gmat_NGM[idx.sign,idx.sign]))
    cor_matrix.NGM.temp <- Gmat_NGM[idx.sign,idx.sign] / (std_devs.NGM.temp %o% std_devs.NGM.temp)
    pleiot.idx$NGM.cor[i]=pleiot.idx$NGM[i] - var(eigen(cor_matrix.NGM.temp)$values)
  }else{
    pleiot.idx$NGM.cor[i]=0 
  }
  
  ## NaCl second
  pleiot.idx$NaCl[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits_bothEnv[9:16])],na.rm=TRUE)
  if(pleiot.idx$NaCl[i]==1){
    pleiot.idx$NaCl.cor[i]=1 
  }else if(pleiot.idx$NaCl[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits_bothEnv[9:16])])
    std_devs.NaCl.temp <- sqrt(diag(Gmat_NaCl[idx.sign,idx.sign]))
    cor_matrix.NaCl.temp <- Gmat_NaCl[idx.sign,idx.sign] / (std_devs.NaCl.temp %o% std_devs.NaCl.temp)
    pleiot.idx$NaCl.cor[i]=pleiot.idx$NaCl[i] - var(eigen(cor_matrix.NaCl.temp)$values)
  }else{
    pleiot.idx$NaCl.cor[i]=0 
  }
  
  ### Mixed
  pleiot.idx_mixed$N[i]=sum(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits_bothEnv)],na.rm=TRUE)
  
  if(pleiot.idx_mixed$N[i]==1){
    pleiot.idx_mixed$Neff[i]=1  
  }else if(pleiot.idx_mixed$N[i]>1){
    idx.sign <- as.logical(all_GWAS_clusters_sign[i,paste0("is_sign_",vect_P_traits_bothEnv)])
    idx.sign[is.na(idx.sign)]=FALSE
    std_devs.all.temp <- sqrt(diag(G_mat_all[idx.sign,idx.sign]))
    cor_matrix.all.temp <- G_mat_all[idx.sign,idx.sign] / (std_devs.all.temp %o% std_devs.all.temp)
    pleiot.idx_mixed$Neff[i]=pleiot.idx_mixed$N[i] - var(eigen(cor_matrix.all.temp)$values)
  }
}

#par(mfrow=c(1,3))
#hist(pleiot.idx$NGM.cor,main='NGM',xlim=c(0,6),ylim=c(0,3000))
#hist(pleiot.idx$NaCl.cor,main='NaCl',xlim=c(0,6),ylim=c(0,3000))
#hist(pleiot.idx_mixed$Neff,main='ALl 12',xlim=c(0,12),ylim=c(0,5000))

##### Now compute the summed effect of the top SNP of each signicant cluster ? ####

### We need to compute the standardized allelic effect ###
all_effects <- all_GWAS_clusters_sign[,paste0('beta_',vect_P_traits_bothEnv)]
all_effects_std <- NULL
for(i in 1:16){
  all_effects_std[i] =
    sd(all_GWAS_clusters_sign[all_GWAS_clusters_sign[,paste0("is_sign_",vect_P_traits_bothEnv[i])],paste0("beta_",vect_P_traits_bothEnv[i])])
}

all_GWAS_clusters_sign_std <- all_GWAS_clusters_sign
for(i in 1:16){
  all_GWAS_clusters_sign_std[,paste0('beta_',vect_P_traits_bothEnv[i])]=all_GWAS_clusters_sign_std[,paste0('beta_',vect_P_traits_bothEnv[i])]/all_effects_std[i]
}

 euclidian_summed <- NULL;euclidian_NGM <- NULL;euclidian_NaCl <- NULL
 manh_summed <- NULL; manh_NGM <- NULL; manh_NaCl <- NULL
 for(i in 1:nrow(all_GWAS_clusters_sign_std)){
   euclidian_summed <- c(euclidian_summed,
                         sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits_bothEnv)])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits_bothEnv)])]^2)))
   euclidian_NGM[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits_bothEnv[1:8])])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits_bothEnv[1:8])])]^2))
   euclidian_NaCl[i]=sqrt(sum(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits_bothEnv[9:16])])[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits_bothEnv[9:16])])]^2))
   
   manh_summed <- c(manh_summed,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits_bothEnv)]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits_bothEnv)])])))
   manh_NGM=c(manh_NGM,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits_bothEnv[1:8])]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits_bothEnv[1:8])])])))
   manh_NaCl=c(manh_NaCl,(sum(abs(as.numeric(all_GWAS_clusters_sign_std[i,paste0('beta_',vect_P_traits_bothEnv[9:16])]))[as.logical(all_GWAS_clusters_sign_std[i,paste0('is_sign_',vect_P_traits_bothEnv[9:16])])])))
   
 }

 ### With nls - replace all
 pdf("plots_BLUEs/Scaling_pleiotropy_Neff_BLUE_Arm_Centers_R2.pdf",h=10,w=12)
 
 is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
 is_center <- all_GWAS_clusters_sign$chr_type2=="center"
 
 temp=pleiot.idx_mixed$Neff[pleiot.idx_mixed$Neff>0]
 temp2 = euclidian_summed[pleiot.idx_mixed$Neff>0]
 
 coefs_euc <- coef(nls(temp2~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx_mixed$Neff[is_arm & pleiot.idx_mixed$Neff>0]
 temp2=euclidian_summed[is_arm & pleiot.idx_mixed$Neff>0]
 coefs_euc_arm <- coef(nls(temp2~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))

 temp=pleiot.idx_mixed$Neff[is_center & pleiot.idx_mixed$Neff>0]
 temp2=euclidian_summed[is_center & pleiot.idx_mixed$Neff>0]
 coefs_euc_center <- coef(nls(temp2~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
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
 coefs_euc.NaCl.arm <- coef(nls(euclidian_NaCl[temp.NaCl & is_arm]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 temp=pleiot.idx$NaCl.cor[temp.NaCl & is_center]
 coefs_euc.NaCl.center <- coef(nls(euclidian_summed[temp.NaCl & is_center]~a+b*temp^c,start=list(a=0,b=0.8,c=.5)))
 
 #par(mfrow=c(2,3))
 plot(euclidian_summed~pleiot.idx_mixed$Neff,pch=16,ylab='TE (euc. distance)',xlab="Neff",xlim=c(0,11),ylim=c(0,14),type='n')
 points(euclidian_summed[is_arm]~pleiot.idx_mixed$Neff[is_arm],pch=16,col='green')
 points(euclidian_summed[is_center]~pleiot.idx_mixed$Neff[is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc[1]+(coefs_euc[2])*(seq(0,12,0.1)^coefs_euc[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc_arm[1]+(coefs_euc_arm[2])*(seq(0,12,0.1)^coefs_euc_arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc_center[1]+(coefs_euc_center[2])*(seq(0,12,0.1)^coefs_euc_center[3])),col="lightblue")
 
 #abline(lm(euclidian_summed~pleiot.idx_mixed$Neff))
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
 
 #abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]))
 text(1.5,7.5,paste0("c=",round(coefs_euc.NGM[3],digits=3)))
 text(1.5,6.5,paste0("c=",round(coefs_euc.NGM.arm[3],digits=3)),col="green")
 text(1.5,5.5,paste0("c=",round(coefs_euc.NGM.center[3],digits=3)),col="lightblue")
 
 plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="Neff",main='NaCl',xlim=c(0,5.5),ylim=c(0,10),type='n')
 points(euclidian_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=16,col='green')
 points(euclidian_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NaCl[1]+(coefs_euc.NaCl[2])*(seq(0,12,0.1)^coefs_euc.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.arm[1]+(coefs_euc.NaCl.arm[2])*(seq(0,12,0.1)^coefs_euc.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.center[1]+(coefs_euc.NaCl.center[2])*(seq(0,12,0.1)^coefs_euc.NaCl.center[3])),col="lightblue")
 
 #abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))
 text(1,8.5,paste0("c=",round(coefs_euc.NaCl[3],digits=3)))
 text(1,7.5,paste0("c=",round(coefs_euc.NaCl.arm[3],digits=3)),col="green")
 text(1,6.5,paste0("c=",round(coefs_euc.NaCl.center[3],digits=3)),col='lightblue')
 
 temp=pleiot.idx_mixed$Neff[pleiot.idx_mixed$Neff>0]
 temp2=manh_summed[pleiot.idx_mixed$Neff>0]
 coefs_man <- coef(nls(temp2~a+b*temp^c,start=list(a=1,b=1,c=1)))

 temp=pleiot.idx_mixed$Neff[is_arm & pleiot.idx_mixed$Neff>0]
 temp2=manh_summed[pleiot.idx_mixed$Neff>0 & is_arm]
 coefs_man.arm <- coef(nls(temp2~a+b*temp^c,start=list(a=1,b=1,c=1)))

 temp=pleiot.idx_mixed$Neff[pleiot.idx_mixed$Neff>0 &is_center]
 temp2=manh_summed[pleiot.idx_mixed$Neff>0 & is_center]
 coefs_man.center <- coef(nls(temp2~a+b*temp^c,start=list(a=1,b=1,c=1)))
 
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
 
 #abline(lm(manh_summed~pleiot.idx_mixed$Neff))
 text(2,30,paste0("c=",round(coefs_man[3],digits=3)))
 text(2,28,paste0("c=",round(coefs_man.arm[3],digits=3)),col='green')
 text(2,26,paste0("c=",round(coefs_man.center[3],digits=3)),col='lightblue')
 
 plot(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM],pch=16,ylab='TM',xlab="Neff",main='NGM',xlim=c(1,7),ylim=c(0,25),type='n')
 points(manh_NGM[temp.NGM & is_arm]~pleiot.idx$NGM.cor[temp.NGM & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NGM[temp.NGM & is_center]~pleiot.idx$NGM.cor[temp.NGM & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NGM[1]+(coefs_man.NGM[2])*(seq(0,12,0.1)^coefs_man.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NGM.arm[1]+(coefs_man.NGM.arm[2])*(seq(0,12,0.1)^coefs_man.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NGM.center[1]+(coefs_man.NGM.center[2])*(seq(0,12,0.1)^coefs_man.NGM.center[3])),col="lightblue")
 
 #abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM.cor[temp.NGM]))
 text(1.5,20,paste0("c=",round(coefs_man.NGM[3],digits=3)))
 text(1.5,18,paste0("c=",round(coefs_man.NGM.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NGM.center[3],digits=3)),col='lightblue')
 
 plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl],pch=16,ylab='TM',xlab="Neff",main='NaCl',xlim=c(0,7),ylim=c(0,25),type='n')
 points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl.cor[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('green',alpha=.6))
 points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl.cor[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="Neff",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(0,12,0.1)^coefs_man.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(0,12,0.1)^coefs_man.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(0,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue")
 
 #abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl.cor[temp.NaCl]))
 text(1.5,20,paste0("c=",round(coefs_man.NaCl[3],digits=3)),col="red")
 text(1.5,18,paste0("c=",round(coefs_man.NaCl.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NaCl.center[3],digits=3)),col='lightblue')
 
 
 ### With N instead of Neff - do a linear model ?
 
 vect_chr_type <- all_GWAS_clusters_sign_std$chr_type2
 model_pleiotropie <- lm(log(manh_summed[pleiot.idx_mixed$Neff>0])~log(pleiot.idx_mixed$Neff[pleiot.idx_mixed$Neff>0])*vect_chr_type[pleiot.idx_mixed$Neff>0])
 model_pleiotropie2 <- lm(log(manh_summed[pleiot.idx_mixed$Neff>0])~log(pleiot.idx_mixed$Neff[pleiot.idx_mixed$Neff>0])+vect_chr_type[pleiot.idx_mixed$Neff>0])
 anova(model_pleiotropie,model_pleiotropie2)
 
 
 dev.off()
 

# 
# #################
# 
# ### The same with N ? ##
 ### With nls - replace all
 pdf("plots_BLUEs/Scaling_pleiotropy_N_BLUE_Arm_Centers_R2.pdf",h=10,w=12)
 
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
 
 #abline(lm(euclidian_summed~pleiot.idx_mixed$N))
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
 
 #abline(lm(euclidian_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM]))
 text(1.5,7.5,paste0("c=",round(coefs_euc.NGM[3],digits=3)),col='red')
 text(1.5,6.5,paste0("c=",round(coefs_euc.NGM.arm[3],digits=3)),col="green")
 text(1.5,5.5,paste0("c=",round(coefs_euc.NGM.center[3],digits=3)),col="lightblue")
 
 plot(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TE (euc. distance)',xlab="N",main='NaCl',xlim=c(0,5.5),ylim=c(0,10),type='n')
 points(euclidian_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl[temp.NaCl & is_arm],pch=16,col='green')
 points(euclidian_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl[temp.NaCl & is_center],pch=16,col='lightblue')
 
 lines(seq(0,12,0.1),I(coefs_euc.NaCl[1]+(coefs_euc.NaCl[2])*(seq(0,12,0.1)^coefs_euc.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.arm[1]+(coefs_euc.NaCl.arm[2])*(seq(0,12,0.1)^coefs_euc.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_euc.NaCl.center[1]+(coefs_euc.NaCl.center[2])*(seq(0,12,0.1)^coefs_euc.NaCl.center[3])),col="lightblue")
 
 #abline(lm(euclidian_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]))
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
 
 #abline(lm(manh_summed~pleiot.idx_mixed$N))
 text(2,30,paste0("c=",round(coefs_man[3],digits=3)),col="red")
 text(2,28,paste0("c=",round(coefs_man.arm[3],digits=3)),col='green')
 text(2,26,paste0("c=",round(coefs_man.center[3],digits=3)),col='lightblue')
 
 plot(manh_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM],pch=16,ylab='TM',xlab="N",main='NGM',xlim=c(1,7),ylim=c(0,25),type='n')
 points(manh_NGM[temp.NGM & is_arm]~pleiot.idx$NGM[temp.NGM & is_arm],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('green',alpha=.6))
 points(manh_NGM[temp.NGM & is_center]~pleiot.idx$NGM[temp.NGM & is_center],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NGM[1]+(coefs_man.NGM[2])*(seq(0,12,0.1)^coefs_man.NGM[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NGM.arm[1]+(coefs_man.NGM.arm[2])*(seq(0,12,0.1)^coefs_man.NGM.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NGM.center[1]+(coefs_man.NGM.center[2])*(seq(0,12,0.1)^coefs_man.NGM.center[3])),col="lightblue")
 
 #abline(lm(manh_NGM[temp.NGM]~pleiot.idx$NGM[temp.NGM]))
 text(1.5,20,paste0("c=",round(coefs_man.NGM[3],digits=3)),col="red")
 text(1.5,18,paste0("c=",round(coefs_man.NGM.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NGM.center[3],digits=3)),col='lightblue')
 
 plot(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl],pch=16,ylab='TM',xlab="N",main='NaCl',xlim=c(0,7),ylim=c(0,25),type='n')
 points(manh_NaCl[temp.NaCl & is_arm]~pleiot.idx$NaCl[temp.NaCl & is_arm],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('green',alpha=.6))
 points(manh_NaCl[temp.NaCl & is_center]~pleiot.idx$NaCl[temp.NaCl & is_center],pch=21,ylab='TM ',xlab="N",bg=adjustcolor('lightblue',alpha=.6))
 
 lines(seq(0,12,0.1),I(coefs_man.NaCl[1]+(coefs_man.NaCl[2])*(seq(0,12,0.1)^coefs_man.NaCl[3])),col="red")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.arm[1]+(coefs_man.NaCl.arm[2])*(seq(0,12,0.1)^coefs_man.NaCl.arm[3])),col="green")
 lines(seq(0,12,0.1),I(coefs_man.NaCl.center[1]+(coefs_man.NaCl.center[2])*(seq(0,12,0.1)^coefs_man.NaCl.center[3])),col="lightblue")
 
 #abline(lm(manh_NaCl[temp.NaCl]~pleiot.idx$NaCl[temp.NaCl]))
 text(1.5,20,paste0("c=",round(coefs_man.NaCl[3],digits=3)),col="red")
 text(1.5,18,paste0("c=",round(coefs_man.NaCl.arm[3],digits=3)),col='green')
 text(1.5,16,paste0("c=",round(coefs_man.NaCl.center[3],digits=3)),col='lightblue')
 

 dev.off()
 

########################################################################
 ########################################################################
 ########################################################################
 
 ### Create a G-matrix from the SNPs - depending on pleiotropy ?
 sign_total <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits_bothEnv)])
G_mat = matrix(rep(NA,256),nrow=16,ncol=16); G_mat_arm =matrix(rep(NA,256),nrow=16,ncol=16); G_mat_center =matrix(rep(NA,256),nrow=16,ncol=16)
G_mat_P1 = matrix(rep(NA,256),nrow=16,ncol=16); G_mat_arm_P1 =matrix(rep(NA,256),nrow=16,ncol=16); G_mat_center_P1 =matrix(rep(NA,256),nrow=16,ncol=16)
G_mat_P2 = matrix(rep(NA,256),nrow=16,ncol=16); G_mat_arm_P2 =matrix(rep(NA,256),nrow=16,ncol=16); G_mat_center_P2 =matrix(rep(NA,256),nrow=16,ncol=16) 
G_mat_P3 = matrix(rep(NA,256),nrow=16,ncol=16); G_mat_arm_P3 =matrix(rep(NA,256),nrow=16,ncol=16); G_mat_center_P3 =matrix(rep(NA,256),nrow=16,ncol=16) 
vect_P1 <- sign_total==1 & all_GWAS_clusters$clust_length<10000
vect_P2 <- sign_total%in%c(2,3) & all_GWAS_clusters$clust_length<10000
vect_P3 <- sign_total>3 & all_GWAS_clusters$clust_length<10000


 for(i in 1:16){
   for(j in 1:16){
     G_mat[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[j])]],na.rm=TRUE)
     G_mat_center[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2=="center"],na.rm=TRUE)
     G_mat_arm[i,j]=sum((2*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM*(1-all_GWAS_clusters$af.NGM))[all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2=="arm"],na.rm=TRUE)
     

     G_mat_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P1,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P1]*(1-all_GWAS_clusters$af.NGM[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P1,paste0("is_sign_",vect_P_traits_bothEnv[j])]],na.rm=TRUE)
     G_mat_center_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P1,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P1]*(1-all_GWAS_clusters$af.NGM[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P1,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2[vect_P1]=="center"],na.rm=TRUE)
     G_mat_arm_P1[i,j]=sum((2*all_GWAS_clusters[vect_P1,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P1,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P1]*(1-all_GWAS_clusters$af.NGM[vect_P1]))[all_GWAS_clusters[vect_P1,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P1,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2[vect_P1]=="arm"],na.rm=TRUE)
     
     G_mat_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P2,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P2]*(1-all_GWAS_clusters$af.NGM[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P2,paste0("is_sign_",vect_P_traits_bothEnv[j])]],na.rm=TRUE)
     G_mat_center_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P2,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P2]*(1-all_GWAS_clusters$af.NGM[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P2,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2[vect_P2]=="center"],na.rm=TRUE)
     G_mat_arm_P2[i,j]=sum((2*all_GWAS_clusters[vect_P2,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P2,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P2]*(1-all_GWAS_clusters$af.NGM[vect_P2]))[all_GWAS_clusters[vect_P2,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P2,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2[vect_P2]=="arm"],na.rm=TRUE)
     
     G_mat_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P3,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P3]*(1-all_GWAS_clusters$af.NGM[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P3,paste0("is_sign_",vect_P_traits_bothEnv[j])]],na.rm=TRUE)
     G_mat_center_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P3,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P3]*(1-all_GWAS_clusters$af.NGM[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P3,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2[vect_P3]=="center"],na.rm=TRUE)
     G_mat_arm_P3[i,j]=sum((2*all_GWAS_clusters[vect_P3,paste0("beta_",vect_P_traits_bothEnv[i])]*all_GWAS_clusters[vect_P3,paste0("beta_",vect_P_traits_bothEnv[j])]*all_GWAS_clusters$af.NGM[vect_P3]*(1-all_GWAS_clusters$af.NGM[vect_P3]))[all_GWAS_clusters[vect_P3,paste0("is_sign_",vect_P_traits_bothEnv[i])] & all_GWAS_clusters[vect_P3,paste0("is_sign_",vect_P_traits_bothEnv[j])] & all_GWAS_clusters$chr_type2[vect_P3]=="arm"],na.rm=TRUE)
     
        }
 }
 
 
plot(scale(eigen(G_mat)$vectors[,1]),scale(eigen(G_mat_P1)$vectors[,1]))

 plot(G_mat[1:8,1:8],G_mat_P1[1:8,1:8])
 
 allPcor=NULL
 
 for(i in 1:15){
   for(j in (i+1):16){
     allPcor=rbind(allPcor,
     data.frame(T1=i,T2=j,Pcor=G_mat_all[i,j],nbClust=sum(rowSums(all_GWAS_clusters_sign[,paste0("is_sign_",vect_P_traits_bothEnv[c(i,j)])])==2),
                nbClust_arm=sum(rowSums(all_GWAS_clusters_sign[all_GWAS_clusters_sign$chr_type2=="arm",paste0("is_sign_",vect_P_traits_bothEnv[c(i,j)])])==2),
                nbClust_center=sum(rowSums(all_GWAS_clusters_sign[all_GWAS_clusters_sign$chr_type2=="center",paste0("is_sign_",vect_P_traits_bothEnv[c(i,j)])])==2))
     )
   }
 }
 allPcor$same_env = 0
 allPcor$same_env[allPcor$T1%in%(1:8) & allPcor$T2%in%(1:8)]=1
 allPcor$same_env[allPcor$T1%in%(9:16) & allPcor$T2%in%(1:16)]=2
 
 allPcor$same_T = 0
 allPcor$same_T[allPcor$T1==(allPcor$T2-8)]=1
 
 
 par(mfrow=c(1,3))
 plot(abs(allPcor$Pcor),allPcor$nbClust,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 plot(abs(allPcor$Pcor),allPcor$nbClust_arm,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 plot(abs(allPcor$Pcor),allPcor$nbClust_center,pch=c(16,11)[1+allPcor$same_T],
      col=c("black","blue","orange")[1+allPcor$same_env],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters")
 
 
 pdf("Lille_05_R2.pdf",w=4,h=4)
 vect_kept=allPcor$same_T==0 & allPcor$same_env>0
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.45),ylim=c(0,200),las=1,bty="n")
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2)
 legend(0,300,c("Low Salt","High Salt"),pch=16,col=c("blue","orange"))
 dev.off()
 
 pdf("Lille_06_R2.pdf",w=4,h=4)
 vect_kept=allPcor$same_T==0
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.45),ylim=c(0,200),las=1,bty="n")
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2)
 legend(0.2,200,c("Low Salt","High Salt","Across environment"),pch=16,col=c("blue","orange","black"),cex=.75)
 dev.off()
 
 
 pdf("Lille_07b_R2.pdf",w=4,h=4)
 vect_kept=allPcor$same_T==0
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.6),ylim=c(0,200),las=1,bty="n")
 
 abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2)
 
 vect_kept=allPcor$same_T==1
 points(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col="red",xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.6),ylim=c(0,200))
 
 #abline(lm(allPcor$nbClust[vect_kept]~abs(allPcor$Pcor)[vect_kept]),lty=2,col="red")

 
 legend(0,320,c("Low Salt","High Salt","Across environment","Same trait"),pch=c(16,16,16,11),col=c("blue","orange","black","red"),cex=.6,ncol=2)
 dev.off()
 
 plot(abs(allPcor$Pcor)[vect_kept],allPcor$nbClust[vect_kept],pch=c(16,11)[1+allPcor$same_T[vect_kept]],
      col=c("black","blue","orange")[1+allPcor$same_env[vect_kept]],xlab="Phenotypic correlation (BV)",ylab="Nb of  shared pleiotropic clusters",
      xlim=c(0,.45),ylim=c(0,310))
 
 
 
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
 
 
 
 