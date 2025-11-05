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

library(lme4)

newdata=data.frame(logD=0,rel_humidity=0, temperature=0,
                   pop_label=sort(unique(data_for_GBLUP_NGM$pop_label)))
all_traits_corrected=NULL
for(i in vect_P_traits){
  print(i)
  data_for_GBLUP_NGM$trait=data_for_GBLUP_NGM[,i]
  trait_lmer <- lmer(trait~ logD + rel_humidity + temperature  + pop_label + (1|date_str) + (1|session_id), data=data_for_GBLUP_NGM)
  if(i=="SF"){
    all_traits_corrected <- data.frame(pop_label=newdata$pop_label,trait=predict(trait_lmer,newdata,re.form=NA))
    names(all_traits_corrected)[2]=i
  }else{
    all_traits_corrected <- cbind(all_traits_corrected,predict(trait_lmer,newdata,re.form=NA))
  }
  
  
}

names(all_traits_corrected)[2:8]=vect_P_traits
save(list=ls(),file="RData_BLUEs/NGM_LMER_BLUEs.RData")


