rm(list=ls())
gc()
library(MCMCglmm)

load("NGM_MCMCglmm_10.RData")

### Export the breeding values for each trait ###
library(data.table)
col_subset <- colnames(model_MCMC.NGM[["Sol"]])[grepl('pop_label', colnames(model_MCMC.NGM[["Sol"]]))]
median_breeding <- (apply(model_MCMC.NGM[["Sol"]][,col_subset],2,median))

traitID <- substring(tstrsplit(names(median_breeding),'[.]')[[1]],6,20)
popLabel <- tstrsplit(names(median_breeding),'pop_label.')[[2]]
final_BV <- data.frame(traitID=traitID,pop_label=popLabel,BV=as.numeric(median_breeding))
write.table(final_BV,file="Breeding_values_NGM.txt",sep='\t',quote=FALSE,row.names=FALSE)

boxplot(final_BV$BV~final_BV$traitID)
boxplot(final_BV$BV[final_BV$traitID=="area"]~tstrsplit(final_BV$pop_label[final_BV$traitID=="area"],"L")[[1]])
boxplot(final_BV$BV[final_BV$traitID=="area"]~substring(final_BV$pop_label[final_BV$traitID=="area"],1,3))
boxplot(final_BV$BV[final_BV$traitID=="SB"]~substring(final_BV$pop_label[final_BV$traitID=="SB"],1,3))



### Export the breeding values for each trait - NaCl ###
load("NaCl_MCMCglmm_10.RData")
col_subset <- colnames(model_MCMC.NaCl[["Sol"]])[grepl('pop_label', colnames(model_MCMC.NaCl[["Sol"]]))]
median_breeding <- (apply(model_MCMC.NaCl[["Sol"]][,col_subset],2,median))

traitID <- substring(tstrsplit(names(median_breeding),'[.]')[[1]],6,20)
popLabel <- tstrsplit(names(median_breeding),'pop_label.')[[2]]
final_BV_NaCl <- data.frame(traitID=traitID,pop_label=popLabel,BV=as.numeric(median_breeding))
write.table(final_BV_NaCl,file="Breeding_values_NaCl.txt",sep='\t',quote=FALSE,row.names=FALSE)
boxplot(final_BV_NaCl$BV[final_BV_NaCl$traitID=="area"]~substring(final_BV_NaCl$pop_label[final_BV_NaCl$traitID=="area"],1,3))
