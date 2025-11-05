rm(list=ls());gc()
library(data.table)
library(MCMCglmm)
library(matrixStats)
##### Load the phenotypic data

phen_data <- read.table("data/all_phenotypes_from_movies.txt",sep='\t',h=TRUE)

# ## Create the session ID 
 phen_data$session_id <- substring(phen_data$date_str,1,4)
 phen_data$session_id[phen_data$session_id=="2014"] <- "2013"
 phen_data$session_id[substring(phen_data$block_id,1,1)=="F"] <- "F"
 phen_data$session_id[substring(phen_data$block_id,1,2)=="SB"] <- "SB"

 ## Set the env. covariates as numeric
 for(i in c("temperature","rel_humidity","logD")){
   phen_data[,i]=as.numeric(phen_data[,i])
 }

 phen_data=subset(phen_data,!is.na(phen_data$pop_label))  # Mostly MA lines and populations
 # Remove additonal populations
 phen_data=subset(phen_data,!phen_data$block_id%in%paste0("B",300:350))
 phen_data=subset(phen_data,!phen_data$pop_label=="A00")
 phen_data=subset(phen_data,!phen_data$pop_label=="CX12311L0")
 phen_data=subset(phen_data,!phen_data$pop_label=="LSJ2L0")
 
 ## Create the population ID
 phen_data$population = tstrsplit(phen_data$pop_label,"L")[[1]]
 phen_data$population[phen_data$population==""]=paste0("L",tstrsplit(phen_data$pop_label,"L")[[2]])[phen_data$population==""]
 phen_data$population[phen_data$population=="S"]="SL3"
 #phen_data$pop_label[phen_data$population=="LSJ2"]
 phen_data$population2=phen_data$population
 phen_data$population2=substring(phen_data$population,1,2)
 phen_data$population2[phen_data$population=="EEV1401"]="A6"

 phen_data$date_str = as.factor(as.character(phen_data$date_str))
 
 ### Remove NA values
 phen_data=subset(phen_data,!is.na(phen_data$temperature))
 phen_data=subset(phen_data,!is.na(phen_data$rel_humidity))
 phen_data=subset(phen_data,!is.na(phen_data$logD))
 
 for(i in c("temperature","rel_humidity","logD")){
   phen_data[,i]= (phen_data[,i]-mean(phen_data[,i],na.rm=TRUE))/sd(phen_data[,i])
 }
 
 vect_P_traits = c("SF", "SB", "FS", "FB", "BS","BF","area.F")

### Finally we scale the body area to match the sd of the transition rates
### currently in mm^2

mean(colSds(as.matrix(phen_data[phen_data$env_label=="NaCl",vect_P_traits[1:6]]))) # 0.65
mean(colSds(as.matrix(phen_data[phen_data$env_label=="NGM",vect_P_traits[1:6]]))) # 0.63
tapply(na.omit(phen_data)$area.F,na.omit(phen_data)$env_label,sd) # 0.014 / 0.019 in NGM - High Salt

phen_data_noNA <- na.omit(phen_data)

mean(c(mean(colSds(as.matrix(phen_data_noNA[phen_data_noNA$env_label=="NaCl",vect_P_traits[1:6]])))/tapply(phen_data_noNA$area.F,phen_data_noNA$env_label,sd)[1],
       mean(colSds(as.matrix(phen_data_noNA[phen_data_noNA$env_label=="NGM",vect_P_traits[1:6]])))/tapply(phen_data_noNA$area.F,phen_data_noNA$env_label,sd)[2]))

# Mean factor in each environment ~38
phen_data$area.F = phen_data$area.F * 38
kept_column_for_MCMCglmm <- c(vect_P_traits,
                              c("logD" , "rel_humidity" , "temperature" , "session_id",
                                "date_str","pop_label"))
data_for_GBLUP_NaCl <- na.omit(subset(phen_data,env_label=='NaCl')[,kept_column_for_MCMCglmm])
data_for_GBLUP_NGM <- na.omit(subset(phen_data,env_label=='NGM')[,kept_column_for_MCMCglmm])

dim(data_for_GBLUP_NaCl) # 2236
dim(data_for_GBLUP_NGM) # 3045

##############################@
nb_trait=7
phen.var = diag(nb_trait) * diag(var(subset(data_for_GBLUP_NGM, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var / 3, nu = nb_trait),
                           G2 = list(V = phen.var / 3, nu = nb_trait)), 
                  R = list(V = phen.var / 3, nu = nb_trait))


model_MCMC.NGM <- MCMCglmm(cbind(SF, SB, FS, FB, BS, BF,area.F) ~ (logD + rel_humidity + temperature + session_id) + pop_label + trait - 1, 
                           random = ~ us(trait):date_str, 
                           rcov = ~us(trait):units, 
                           family = rep("gaussian", nb_trait), 
                           prior = prior_mod, 
                           data = data_for_GBLUP_NGM, 
                           nitt = 1100000, 
                           burnin = 100000, 
                           thin = 1000, 
                           pr = TRUE)


save(list=ls(),file="NGM_MCMCglmm.RData")

# rownames(sol.conv2)

# #to check convergence

# 
# sol.conv <- heidel.diag(model_MCMC.NGM[["Sol"]][,col_subset])
# VCV.conv <- heidel.diag(model_MCMC.NGM[["VCV"]])
# 
# sol.conv2=sol.conv
# rownames(sol.conv2)=rep("",nrow(sol.conv2))
# # Let's look at autocorrelation and effective sample size
# autocorr.vect <- autocorr.diag(model_MCMC.NGM[["Sol"]][,col_subset])
# autocorr.plot(model_MCMC.NGM[["Sol"]])
# 
# 

