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

# Mean factor in each environment ~38
kept_column_for_MCMCglmm <- c(vect_P_traits,
                              c("logD" , "rel_humidity" , "temperature" , "session_id",
                                "date_str","pop_label"))

data_for_BLUE_NaCl <- na.omit(subset(phen_data,env_label=='NaCl')[,kept_column_for_MCMCglmm])
data_for_BLUE_NGM <- na.omit(subset(phen_data,env_label=='NGM')[,kept_column_for_MCMCglmm])

### Finally we scale all traits

for(i in vect_P_traits){
  data_for_BLUE_NaCl[,i] = scale(data_for_BLUE_NaCl[,i])
  data_for_BLUE_NGM[,i] = scale(data_for_BLUE_NGM[,i])
}

dim(data_for_BLUE_NaCl) # 2236
dim(data_for_BLUE_NGM) # 3045

###################################################################################################
######################################### Univariate models #######################################
###################################################################################################

library(lme4)

newdata=data.frame(logD=0,rel_humidity=0, temperature=0,
                   pop_label=sort(unique(data_for_BLUE_NaCl$pop_label)))

for(i in vect_P_traits){
  print(i)
    data_for_BLUE_NaCl$trait=data_for_BLUE_NaCl[,i]
    trait_lmer <- lmer(trait~ logD + rel_humidity + temperature  + pop_label + (1|date_str) + (1|session_id), data=data_for_BLUE_NaCl)
    if(i=="SF"){
      all_traits_corrected <- data.frame(pop_label=newdata$pop_label,trait=predict(trait_lmer,newdata,re.form=NA))
      names(all_traits_corrected)[2]=i
    }else{
      all_traits_corrected <- cbind(all_traits_corrected,predict(trait_lmer,newdata,re.form=NA))
    }
}
names(all_traits_corrected)[2:8]=vect_P_traits
head(all_traits_corrected)
save(list=ls(),file="RData_BLUEs/NaCl_LMER_BLUEs.RData")

newdata=data.frame(logD=0,rel_humidity=0, temperature=0,
                   pop_label=sort(unique(data_for_BLUE_NGM$pop_label)))

all_traits_corrected=NULL
for(i in vect_P_traits){
  print(i)
  data_for_BLUE_NGM$trait=data_for_BLUE_NGM[,i]
  trait_lmer <- lmer(trait~ logD + rel_humidity + temperature  + pop_label + (1|date_str) + (1|session_id), data=data_for_BLUE_NGM)
  if(i=="SF"){
    all_traits_corrected <- data.frame(pop_label=newdata$pop_label,trait=predict(trait_lmer,newdata,re.form=NA))
    names(all_traits_corrected)[2]=i
  }else{
    all_traits_corrected <- cbind(all_traits_corrected,predict(trait_lmer,newdata,re.form=NA))
  }
  
  
}

names(all_traits_corrected)[2:8]=vect_P_traits
save(list=ls(),file="RData_BLUEs/NGM_LMER_BLUEs.RData")
