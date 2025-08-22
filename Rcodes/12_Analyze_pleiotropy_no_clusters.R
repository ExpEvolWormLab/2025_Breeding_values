rm(list=ls());gc()
#load("RData/GWAS_cluster_results_concatenated.RData")
load(file="RData/GWAS_results_concatenated.RData")
gemma_pheno=read.table('data/gemma_pheno.txt',sep='\t',h=FALSE)

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area")
vect_P_traits_bothEnv <- c(paste0(vect_P_traits,"_NGM"),paste0(vect_P_traits,"_NaCl"))

### Could there be SNPs in the same genes but not the same variants 
### affecting the traits ? -> would be nice to give examples ????

#### Now the Neff to estimate the effect ~pleiotropy ?
## For each cluster, compute a level of pleiotropy

all_GWAS_sign <- na.omit(all_GWAS[rowSums(all_GWAS[,paste0("is_sign_T",1:14)])>0,])
all_GWAS_sign_NGM <- na.omit(all_GWAS[rowSums(all_GWAS[,paste0("is_sign_T",1:7)])>0,])
all_GWAS_sign_NaCl <- na.omit(all_GWAS[rowSums(all_GWAS[,paste0("is_sign_T",8:14)])>0,])

nb_sign_NGM <- nrow(all_GWAS_sign)
nb_sign_NaCl <- nrow(all_GWAS_sign)
nb_sign <- nrow(all_GWAS_sign)

pleiot.idx=data.frame(NGM=rep(0,nb_sign),NaCl=rep(0,nb_sign))
pleiot.idx$NGM.cor=0;pleiot.idx$NaCl.cor=0

pleiot.idx_mixed=data.frame(N=rep(0,nb_sign),Neff=rep(0,nb_sign))

Gmat_NGM = cov(na.omit(all_GWAS[,paste0("beta_T",1:7)]))
Gmat_NaCl = cov(na.omit(all_GWAS[,paste0("beta_T",8:14)]))
G_mat_all=cov(na.omit(all_GWAS[,paste0("beta_T",1:14)]))

for(i in 1:nrow(pleiot.idx_mixed)){
  ## NGM first
  pleiot.idx$NGM[i]=sum(all_GWAS_sign[i,paste0("is_sign_T",1:7)],na.rm=TRUE)
  if(pleiot.idx$NGM[i]==1){
    pleiot.idx$NGM.cor[i]=1 
  }else if(pleiot.idx$NGM[i]>1){
    idx.sign <- as.logical(all_GWAS_sign[i,paste0("is_sign_T",1:7)])
    std_devs.NGM.temp <- sqrt(diag(Gmat_NGM[idx.sign,idx.sign]))
    cor_matrix.NGM.temp <- Gmat_NGM[idx.sign,idx.sign] / (std_devs.NGM.temp %o% std_devs.NGM.temp)
    pleiot.idx$NGM.cor[i]=pleiot.idx$NGM[i] - var(eigen(cor_matrix.NGM.temp)$values)
  }
  
  ## NaCl second
  pleiot.idx$NaCl[i]=sum(all_GWAS_sign[i,paste0("is_sign_T",8:14)],na.rm=TRUE)
  if(pleiot.idx$NaCl[i]==1){
    pleiot.idx$NaCl.cor[i]=1 
  }else if(pleiot.idx$NaCl[i]>1){
    idx.sign <- as.logical(all_GWAS_sign[i,paste0("is_sign_T",8:14)])
    std_devs.NaCl.temp <- sqrt(diag(Gmat_NaCl[idx.sign,idx.sign]))
    cor_matrix.NaCl.temp <- Gmat_NaCl[idx.sign,idx.sign] / (std_devs.NaCl.temp %o% std_devs.NaCl.temp)
    pleiot.idx$NaCl.cor[i]=pleiot.idx$NaCl[i] - var(eigen(cor_matrix.NaCl.temp)$values)
  }
  
  ### Mixed
  pleiot.idx_mixed$N[i]=sum(all_GWAS_sign[i,paste0("is_sign_T",1:14)],na.rm=TRUE)
  
  if(pleiot.idx_mixed$N[i]==1){
    pleiot.idx_mixed$Neff[i]=1  
  }else if(pleiot.idx_mixed$N[i]>1){
    idx.sign <- as.logical(all_GWAS_sign[i,paste0("is_sign_T",1:14)])
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
all_effects <- all_GWAS_sign[,paste0('beta_T',1:14)]
all_effects_std <- NULL
for(i in 1:14){
  all_effects_std[i] =
    sd(all_GWAS_sign[all_GWAS_sign[,paste0("is_sign_T",i)],paste0("beta_T",i)])
}

all_GWAS_sign_std <- all_GWAS_sign
for(i in 1:14){
  all_GWAS_sign_std[,paste0('beta_T',i)]=all_GWAS_sign_std[,paste0('beta_T',i)]/all_effects_std[i]
}

euclidian_summed <- NULL;euclidian_NGM <- NULL;euclidian_NaCl <- NULL
manh_summed <- NULL; manh_NGM <- NULL; manh_NaCl <- NULL
for(i in 1:nrow(all_GWAS_sign_std)){
  euclidian_summed <- c(euclidian_summed,sqrt(sum(as.numeric(all_GWAS_sign_std[i,paste0('beta_T',1:14)])[as.logical(all_GWAS_sign_std[i,paste0('is_sign_T',1:14)])]^2)))
  euclidian_NGM[i]=sqrt(sum(as.numeric(all_GWAS_sign_std[i,paste0('beta_T',1:7)])[as.logical(all_GWAS_sign_std[i,paste0('is_sign_T',1:7)])]^2))
  euclidian_NaCl[i]=sqrt(sum(as.numeric(all_GWAS_sign_std[i,paste0('beta_T',8:14)])[as.logical(all_GWAS_sign_std[i,paste0('is_sign_T',8:14)])]^2))
  
  manh_summed <- c(manh_summed,(sum(abs(as.numeric(all_GWAS_sign_std[i,paste0('beta_T',1:14)]))[as.logical(all_GWAS_sign_std[i,paste0('is_sign_T',1:14)])])))
  manh_NGM=c(manh_NGM,(sum(abs(as.numeric(all_GWAS_sign_std[i,paste0('beta_T',1:7)]))[as.logical(all_GWAS_sign_std[i,paste0('is_sign_T',1:7)])])))
  manh_NaCl=c(manh_NaCl,(sum(abs(as.numeric(all_GWAS_sign_std[i,paste0('beta_T',8:14)]))[as.logical(all_GWAS_sign_std[i,paste0('is_sign_T',8:14)])])))
  
}

pdf("plots/Scaling_pleiotropy_Neff_no_clusters.pdf",h=10,w=12)

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
pdf("plots/Scaling_pleiotropy_N_no_clusters.pdf",h=14,w=12)

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









