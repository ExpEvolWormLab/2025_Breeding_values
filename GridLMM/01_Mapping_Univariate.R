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
# 
# phen_data <- read.table("data/all_phenotypes_from_movies_fixed.txt",sep='\t',h=TRUE)
# 
# # ## Create the session ID 
# phen_data$session_id <- substring(phen_data$date_str,1,4)
# phen_data$session_id[phen_data$session_id=="2014"] <- "2013"
# phen_data$session_id[substring(phen_data$block_id,1,1)=="F"] <- "F"
# phen_data$session_id[substring(phen_data$block_id,1,2)=="SB"] <- "SB"
# 
# ## Set the env. covariates as numeric
# for(i in c("temperature","rel_humidity","logD")){
#   phen_data[,i]=as.numeric(phen_data[,i])
# }
# 
# phen_data=subset(phen_data,!is.na(phen_data$pop_label))  # Mostly MA lines and populations
# # Remove additonal populations
# phen_data=subset(phen_data,!phen_data$block_id%in%paste0("B",300:350))
# phen_data=subset(phen_data,!phen_data$pop_label=="A00")
# phen_data=subset(phen_data,!phen_data$pop_label=="CX12311L0")
# phen_data=subset(phen_data,!phen_data$pop_label=="LSJ2L0")
# 
# ## Create the population ID
# phen_data$population = tstrsplit(phen_data$pop_label,"L")[[1]]
# phen_data$population[phen_data$population==""]=paste0("L",tstrsplit(phen_data$pop_label,"L")[[2]])[phen_data$population==""]
# phen_data$population[phen_data$population=="S"]="SL3"
# #phen_data$pop_label[phen_data$population=="LSJ2"]
# phen_data$population2=phen_data$population
# phen_data$population2=substring(phen_data$population,1,2)
# phen_data$population2[phen_data$population=="EEV1401"]="A6"
# 
# phen_data$date_str = as.factor(as.character(phen_data$date_str))
# 
# ### Remove NA values
# phen_data=subset(phen_data,!is.na(phen_data$temperature))
# phen_data=subset(phen_data,!is.na(phen_data$rel_humidity))
# phen_data=subset(phen_data,!is.na(phen_data$logD))
# 
# for(i in c("temperature","rel_humidity","logD")){
#   phen_data[,i]= scale(phen_data[,i])
# }
# 
# vect_P_traits = c("SF", "SB", "FS", "FB", "BS","BF","area.F")
# 
# # Mean factor in each environment ~38
# kept_column <- c(vect_P_traits,
#                  c("logD" , "rel_humidity" , "temperature" , "session_id",
#                    "date_str","pop_label"))
# 
# data_for_BLUE_NaCl <- (subset(phen_data,env_label=='NaCl')[,kept_column])
# data_for_BLUE_NGM <- (subset(phen_data,env_label=='NGM')[,kept_column])
# 
# ### Finally we scale all traits
# 
# dim(data_for_BLUE_NaCl) # 2335
# dim(data_for_BLUE_NGM) # 3046
# 
# data_for_BLUE_NaCl$env_label='NaCl'
# data_for_BLUE_NGM$env_label='NGM'
# all_BLUEs <- rbind(data_for_BLUE_NGM,data_for_BLUE_NaCl)
# dim(all_BLUEs) # 5381 movies
# 
# ## Load the isotype data
# removed_isotypes <- read.table("data/removed_isotypes_Cemee_2025.txt")
# corr_isotypes <- read.csv("data/correspondance_isotypes.csv")
# ### Remove the additional name info
# removed_isotypes$V1 <- gsub('CeMee', '', removed_isotypes$V1)
# removed_isotypes$V1 <- gsub('_sorted', '', removed_isotypes$V1)
# removed_isotypes=unique(removed_isotypes)
# corr_isotypes$Line <- gsub('CeMee', '', corr_isotypes$Line)
# corr_isotypes$Line <- gsub('_sorted', '', corr_isotypes$Line)
# corr_isotypes$Representant_line <- gsub('CeMee', '', corr_isotypes$Representant_line)
# corr_isotypes$Representant_line <- gsub('_sorted', '', corr_isotypes$Representant_line)
# corr_isotypes <- corr_isotypes[,c(1,3)]
# corr_isotypes=subset(corr_isotypes,Line!=Representant_line)
# table(removed_isotypes$V1%in%corr_isotypes$Representant_line)
# replaced_iso <- subset(corr_isotypes,Line%in%removed_isotypes$V1)
# 
# ### We will use the bimbam format which is simple for GEMMA
# ## First load the genotypes
# ## Load Cemee data
# 
# load('RData/Genotype_info_Cemee.RData')
# 
# snp_positions=as.data.frame(snp_positions)
# snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
# snp_positions$POS=as.numeric(as.character(snp_positions$POS))
# 
# ## Order the genotypes
# genotype_info <- genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
# snp_positions <- snp_positions[order(snp_positions$CHROM,snp_positions$POS),]
# 
# ### Simply replace the isotype lines
# non.isotype=NULL;pbm.line=NULL
# for(i in 1:nrow(all_BLUEs)){
#   temp.iso=subset(replaced_iso,Line==all_BLUEs$pop_label[i] | Representant_line==all_BLUEs$pop_label[i])
#   if(length(unique(temp.iso$Representant_line))==1){
#   all_BLUEs$pop_label[i]=unique(temp.iso$Representant_line)
#   }else if(nrow(temp.iso)==0){
#     non.isotype=c(non.isotype,all_BLUEs$pop_label[i])
#   }else{
#     pbm.line=c(pbm.line,all_BLUEs$pop_label[i])
#   }
# }
# 
# 
# ## Order the genotype - keep shared lines
# genotyped_lines <- colnames(genotype_info)
# 
# all_BLUEs_genotyped=subset(all_BLUEs,pop_label%in%genotyped_lines)
# dim(all_BLUEs_genotyped) # 4420 out of 5381
# phenotyped_lines=unique(all_BLUEs_genotyped$pop_label)
# genotype_info_kept <- genotype_info[,colnames(genotype_info)%in%phenotyped_lines]
# dim(genotype_info_kept) # 966 kept, OK
# 
# bimbam <- data.frame(id=snp_positions$ID,major=snp_positions$REF,minor=snp_positions$ALT)
# bimbam <- cbind(bimbam,genotype_info_kept)
# genotype_matrix=t(bimbam[,4:ncol(bimbam)])
# colnames(genotype_matrix)=bimbam$id
# 
# ## Load the kinship matrix
# library("vcfR")
# genotype_pruned <- read.vcfR('plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz')
# genotype_pruned_matrix <- apply(genotype_pruned@gt[,2:1002],1,function(x){
#   as.numeric(substring(x,1,1))
# })
# genotype_pruned_matrix=2*genotype_pruned_matrix
# strainID <- tstrsplit(colnames(genotype_pruned@gt)[2:1002],'_')[[1]]
# GRM <- Gmatrix((genotype_pruned_matrix), method="VanRaden") # 2112 / 2112
# dim(GRM)
# ## Re-order
# 
# GRM_sub = GRM[strainID %in% colnames(genotype_info_kept),strainID %in% colnames(genotype_info_kept)]
# strainID_sub = strainID[strainID %in% colnames(genotype_info_kept)]
# table(strainID_sub==rownames(genotype_matrix)) # Same order
# 
# ##############################
# 
# all_BLUEs_genotyped_NGM=subset(all_BLUEs_genotyped,env_label=="NGM")
# all_BLUEs_genotyped_NaCl=subset(all_BLUEs_genotyped,env_label=="NaCl")
# 
# colnames(GRM_sub)=rownames(genotype_matrix)
# rownames(GRM_sub)=rownames(genotype_matrix)
# 
# ## Remove SNPs with NAs
# is.NA=colSums(is.na(genotype_matrix))
# table(is.NA)


save(list=c("all_BLUEs_genotyped_NGM","all_BLUEs_genotyped_NaCl","genotype_matrix","is.NA","GRM_sub"),file='RData/GridLMM_mapping.RData')
rm(list=ls());gc()
load('RData/GridLMM_mapping.RData')
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F")

for(traitID in vect_P_traits){

all_BLUEs_genotyped_NGM$y=all_BLUEs_genotyped_NGM[,traitID]
all_BLUEs_genotyped_NaCl$y=all_BLUEs_genotyped_NaCl[,traitID]

all_BLUEs_genotyped_NGM$y= scale(all_BLUEs_genotyped_NGM$y)
all_BLUEs_genotyped_NaCl$y= scale(all_BLUEs_genotyped_NaCl$y)


gwas.NGM = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                         data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')
NGM_results <- (gwas.NGM$results)
write.table(NGM_results,paste0("GridLMM_results/gwas_results_",traitID,"_NGM.txt"),quote=FALSE,row.names=FALSE)

gwas.NaCl = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str)  + (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                        data = all_BLUEs_genotyped_NaCl, X = genotype_matrix[,is.NA==0], X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')
NaCl_results <- (gwas.NaCl$results)
write.table(NGM_results,paste0("GridLMM_results/gwas_results_",traitID,"_NaCl.txt"),quote=FALSE,row.names=FALSE)

}
