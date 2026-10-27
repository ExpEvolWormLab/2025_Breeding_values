rm(list=ls());gc()
library(vcfR)
library(ggplot2)
#install.packages("~/Downloads/rutilstimflutre-0.170.0.tar.gz")
library(rutilstimflutre)
library(ASRgenomics)
#devtools::install_github('deruncie/GridLMM')
library(GridLMM)
library(data.table)
library(matrixStats)
library(performance)
library(lme4)
library(AGHmatrix)
library(qqman)
##### Load the phenotypic data
#library(MegaLMM)

phen_data <- read.table("data/all_phenotypes_from_movies_fixed.txt",sep='\t',h=TRUE)

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
  phen_data[,i]= scale(phen_data[,i])
}

vect_P_traits = c("SF", "SB", "FS", "FB", "BS","BF","area.F")

# Mean factor in each environment ~38
kept_column <- c(vect_P_traits,
                 c("logD" , "rel_humidity" , "temperature" , "session_id",
                   "date_str","pop_label"))

data_for_BLUE_NaCl <- (subset(phen_data,env_label=='NaCl')[,kept_column])

#####################################################################################
#####################################################################################
#Load the fertility table from Noble 2017
load("data/fertility.rda")

mean_relative_fit <- subset(ecoefs,env=="NaCl")
mean_relative_fit$w <- exp(mean_relative_fit$fertility)/mean(exp(mean_relative_fit$fertility))

mean_relative_fit <- mean_relative_fit[,c("line","w")]
names(mean_relative_fit)[1] = 'pop_label'

dim(data_for_BLUE_NaCl)
data_for_BLUE_NaCl_with_w <- (merge(data_for_BLUE_NaCl, mean_relative_fit))
dim(data_for_BLUE_NaCl_with_w) # 506

nb_trait = 8
k=0

vect_P_traits <- c(vect_P_traits,"w")

## Scale
for(i in vect_P_traits[1:7]){
  data_for_BLUE_NaCl_with_w[,i]=scale(data_for_BLUE_NaCl_with_w[,i])
}

phen.var = diag(nb_trait) * diag(var(subset(data_for_BLUE_NaCl_with_w, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var/3, nu =1.002), G2 = list(V = phen.var/3, nu = 1.002)), 
                  R = list(V = phen.var/3, nu = 1.002))

model_MCMC <- MCMCglmm(cbind(c(SF, SB, FS, FB, BS, BF,area.F,w)) ~ temperature+rel_humidity+logD + session_id + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                       rcov = ~us(trait):units, prior = prior_mod,
                       family = rep("gaussian", nb_trait), data = data_for_BLUE_NaCl_with_w, verbose = TRUE,nitt=500000, burnin=50000,thin=100)


library(brms)


data_for_BLUE_NaCl_with_w$pop_label <- factor(data_for_BLUE_NaCl_with_w$pop_label)
data_for_BLUE_NaCl_with_w$date_str <- factor(data_for_BLUE_NaCl_with_w$date_str)


#model_brms <- brm(bf(mvbind(SF, SB, FS, FB, BS, BF, area.F, w) ~ 
#temperature +rel_humidity + logD + session_id + 0 + trait + (1 | p | pop_label) + (1 | q | date_str)) + 
#  set_rescor(TRUE),data = data_for_BLUE_NaCl_with_w,family = gaussian(),
#chains = 4,cores = 4,iter = 5000,warmup = 1000,seed = 123)

par(mfrow = c(7, 2), mar = c(2, 2, 1, 0))
plot(model_MCMC$Sol, auto.layout = F)
plot.acfs(model_MCMC$Sol)

post_dist = posterior.mode(model_MCMC$VCV)

VCV_with_w_NaCl <- list(Population = "A6_with_w", G1_mat = matrix(post_dist[1:nb_trait^2], 
                                                                                                        nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                            R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)



 s = VCV_with_w_NaCl$G1_mat[1:7,8]
 G = VCV_with_w_NaCl$G1_mat[1:7,1:7]
 beta_MCMC = solve(G)%*%s
 plot(beta_MCMC,beta_NaCl,pch=16)
 abline(a=0,b=1)
 beta_VCV=NULL
 for(i in 1:nrow(model_MCMC$VCV)){
     
       temp_Gzw = matrix(model_MCMC$VCV[i,1:64],8,8)
       beta_VCV = cbind(beta_VCV,solve(temp_Gzw[1:7,1:7])%*%temp_Gzw[1:7,8])
       
       }
HPDinterval(as.mcmc(t(beta_VCV)))

save(list=ls(),file='temp.beta.MCMC_test.RData')