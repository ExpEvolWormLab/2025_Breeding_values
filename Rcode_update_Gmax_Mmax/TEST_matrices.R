

load('RData/Gmat_sommer_NGM.RData')
load('RData/Gmat_sommer_NaCl.RData')

Gmat_NaCl <- fit.NaCl$sigma$`u:line`
Rmat_NaCl <- fit.NaCl$sigma$`u:units`

Gmat_NGM <- fit_NGM$sigma$`u:line`
Rmat_NGM <- fit_NGM$sigma$`u:units`

write.table(Gmat_NaCl,file='txt_BLUEs/gmat_NaCl.txt',col.names=FALSE,quote=FALSE)
write.table(Gmat_NGM,file='txt_BLUEs/gmat_NGM.txt',col.names=FALSE,quote=FALSE)

abs(eigen(Gmat_NaCl)$values)/sum(abs(eigen(Gmat_NaCl)$values))
abs(eigen(Gmat_NGM)$values)/sum(abs(eigen(Gmat_NGM)$values))

barplot(eigen(Gmat_NGM)$values)
barplot(eigen(Gmat_NaCl)$values)

eigen(Gmat_NGM)$vector[,1]
eigen(Gmat_NaCl)$vector[,2]

gamma_stab <- read.table("txt_BLUEs/gamma_NGM_stabilizing.txt")
beta_NaCl <- read.table("txt_BLUEs/betas_NaCl.txt")
#eigen(cov(phen_NaCl_sub[,2:8]))$vector

cor(eigen(Gmat_NGM)$vector[,1],eigen(Gmat_NaCl)$vector)
cor(eigen(Gmat_NaCl)$vector[,1],eigen(Gmat_NGM)$vector)

cor(gamma_stab$V1,eigen(Gmat_NGM)$vector)
cor(gamma_stab$V1,eigen(Gmat_NaCl)$vector)

cor(beta_NaCl$V2,eigen(Gmat_NGM)$vector)
cor(beta_NaCl$V2,eigen(Gmat_NaCl)$vector)

100*gamma_stab$V1%*%Gmat_NaCl%*%gamma_stab$V1/sum(abs(eigen(Gmat_NaCl)$values))
100*gamma_stab$V1%*%Gmat_NGM%*%gamma_stab$V1/sum(abs(eigen(Gmat_NGM)$values))

100*beta_NaCl$V2%*%Gmat_NaCl%*%beta_NaCl$V2/sum(abs(eigen(Gmat_NaCl)$values))
100*beta_NaCl$V2%*%Gmat_NGM%*%beta_NaCl$V2/sum(abs(eigen(Gmat_NGM)$values))


cor(beta,eigen(Gmat_NGM)$vector)
cor(beta,eigen(Gmat_NaCl)$vector)

cor(eigen(gamma)$vectors[,7],eigen(Gmat_NGM)$vector)
cor(eigen(gamma)$vectors[,7],eigen(Gmat_NaCl)$vector)



