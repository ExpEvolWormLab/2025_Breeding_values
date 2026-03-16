rm(list=ls())
gc()
library(data.table)
library(AGHmatrix)
vect_P_traits= c("SF","SB","FS","FB","BS","BF","area")
# Load phenotypic data
BLUE1 <- read.table("RData_BLUEs_plasticity/BLUES_NGM.txt",sep='\t',h=TRUE)
BLUE2 <- read.table("RData_BLUEs_plasticity/BLUES_NaCl.txt",sep='\t',h=TRUE)

## Create a data frame

phen_NaCl = data.frame(line=BLUE2$pop_label,SF=NA,SB=NA,FS=NA,FB=NA,BS=NA,BF=NA,area=NA)
for(i in 1:nrow(phen_NaCl)){
  phen_NaCl[i,2:8]=as.numeric(subset(BLUE2,pop_label==phen_NaCl$line[i])$BLUE)
}
phen_NaCl=unique(phen_NaCl)
dim(phen_NaCl)

phen_NGM = data.frame(line=BLUE1$pop_label,SF=NA,SB=NA,FS=NA,FB=NA,BS=NA,BF=NA,area=NA)
for(i in 1:nrow(phen_NGM)){
  phen_NGM[i,2:8]=as.numeric(subset(BLUE1,pop_label==phen_NGM$line[i])$BLUE)
}
phen_NGM=unique(phen_NGM)
dim(phen_NGM)


#Load the fertility table from Noble 2017
#load("fertility.rda")
load("data/fertility.rda")
fertility_NaCl <- subset(ecoefs,env=="NaCl")
phen_NaCl_with_fitness = merge(phen_NaCl,fertility_NaCl[,1:2])
dim(phen_NaCl_with_fitness) # 251 measurements

phen_NaCl_with_fitness$rel_fit = exp(phen_NaCl_with_fitness$fertility)/mean(exp(phen_NaCl_with_fitness$fertility),na.rm=TRUE)

for(i in vect_P_traits){
  phen_NaCl_with_fitness[,i]=scale(phen_NaCl_with_fitness[,i])
}

model_beta <- lm(rel_fit ~ SF+SB+FS+FB+BS+BF+area, data = phen_NaCl_with_fitness)
beta_NaCl <- coef(model_beta)[vect_P_traits]

phen_NaCl$beta_trait = NULL
for(i in 1:nrow(phen_NaCl)){
  phen_NaCl$beta_trait[i]  = sum(phen_NaCl[i,vect_P_traits]*beta_NaCl)
}
hist(phen_NaCl$beta_trait)

phen_NGM$beta_trait = NULL
for(i in 1:nrow(phen_NGM)){
  phen_NGM$beta_trait[i]  = sum(phen_NGM[i,vect_P_traits]*beta_NaCl)
}
hist(phen_NGM$beta_trait)

names(phen_NGM)[2:9]=paste0(names(phen_NGM)[2:9],".NGM")
names(phen_NaCl)[2:9]=paste0(names(phen_NaCl)[2:9],".NaCl")

dim(phen_NGM);dim(phen_NaCl)
dim(merge(phen_NGM,phen_NaCl))

phen_all=merge(phen_NGM,phen_NaCl)
plot(phen_all$beta_trait.NGM,phen_all$beta_trait.NaCl)
mean(phen_all$beta_trait.NGM,na.rm=TRUE)
mean(phen_all$beta_trait.NaCl,na.rm=TRUE)

boxplot(c(phen_all$beta_trait.NGM,phen_all$beta_trait.NaCl)~rep(c("NGM","NaCl"),each=nrow(phen_all)))

# Load the isotype information
removed_isotypes <- read.table("data/removed_isotypes_Cemee_2025.txt")
corr_isotypes <- read.csv("data/correspondance_isotypes.csv")

### Remove the additional name info
removed_isotypes$V1 <- gsub('CeMee', '', removed_isotypes$V1)
removed_isotypes$V1 <- gsub('_sorted', '', removed_isotypes$V1)
removed_isotypes=unique(removed_isotypes)
corr_isotypes$Line <- gsub('CeMee', '', corr_isotypes$Line)
corr_isotypes$Line <- gsub('_sorted', '', corr_isotypes$Line)
corr_isotypes$Representant_line <- gsub('CeMee', '', corr_isotypes$Representant_line)
corr_isotypes$Representant_line <- gsub('_sorted', '', corr_isotypes$Representant_line)
corr_isotypes <- corr_isotypes[,c(1,3)]
corr_isotypes=subset(corr_isotypes,Line!=Representant_line)
table(removed_isotypes$V1%in%corr_isotypes$Representant_line)
replaced_iso <- subset(corr_isotypes,Line%in%removed_isotypes$V1)


#all_BLUEs$traitID <- factor(as.character(all_BLUEs$traitID),levels = 
#                              paste0(rep(c("SF","SB","FS",'FB',"BS","BF",'area'),2),rep(c("_NGM","_NaCl"),each=7)))

### We will use the bimbam format which is simple for GEMMA
## First load the genotypes
## Load Cemee data
load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))

## Order the genotypes
genotype_info <- genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
snp_positions <- snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

## Order the genotype - keep shared lines
genotyped_lines <- colnames(genotype_info)
all_lines <- sort(unique(phen_all$line))
phenotyped_lines = rep(TRUE,length(genotyped_lines))
gemma_pheno <- NULL


### Now we will create the phenotype file for GEMMA

k=0
for(i in genotyped_lines){
  k=k+1
  ## Is there multiple isotypes ?
  temp_iso <- subset(replaced_iso,i==Line | i==Representant_line)
  all_lines_ID <- unique(c(i,temp_iso$Line,temp_iso$Representant_line))
  # Phenotyped once ?
  if(sum(phen_all$line%in%all_lines_ID)==0){
    phenotyped_lines[k]=FALSE
  }else{
    temp_pheno <- subset(phen_all,line%in%all_lines_ID)
    ## Now extract all phenotypic data for these lines
    if(nrow(temp_pheno)==0){
      gemma_pheno=rbind(gemma_pheno,rep(NA,2))
    }else if(nrow(temp_pheno)==1){
      gemma_pheno <- rbind(gemma_pheno,c(temp_pheno[,"beta_trait.NGM"],temp_pheno[,"beta_trait.NaCl"]))
      
    }else{
      gemma_pheno <- rbind(gemma_pheno,c(mean(temp_pheno[,"beta_trait.NGM"]),mean(temp_pheno[,"beta_trait.NaCl"])))
    }
  }
}

table(phenotyped_lines) # 960 out of 966 TRUE
dim(gemma_pheno) # 960 OK 

genotype_info_kept <- genotype_info[,phenotyped_lines]

bimbam <- data.frame(id=snp_positions$ID,major=snp_positions$REF,minor=snp_positions$ALT)
bimbam <- cbind(bimbam,genotype_info_kept)
write.table(bimbam,"data_BLUE_plasticity/genotypes_for_gemma.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(gemma_pheno,file="data_BLUE_plasticity/gemma_pheno.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

gemma_pheno_withlines=data.frame(gemma_pheno)
gemma_pheno_withlines$strain=
  as.character(genotyped_lines[phenotyped_lines])

write.table(gemma_pheno_withlines,file="data_BLUE_plasticity/gemma_pheno_withLines.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

### Compute a GRM
##### Pruning
for(i in c(1,2,4,6,8,9,99)){
  system(paste0("mkdir plink_pruning/pruning0",i))
  system(paste0("/mnt/data2/desmarais/tool_plink/plink --vcf /projects/nematode_genetics/Cemee_data/vcf_2025/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.",i," --out plink_pruning/pruning0",i,"/CeMEE.pruned.0.",i," --allow-extra-chr --double-id"))
  system(paste0("/mnt/data2/desmarais/tool_plink/plink --vcf /projects/nematode_genetics/Cemee_data/vcf_2025/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract plink_pruning/pruning0",i,"/CeMEE.pruned.0.",i,".prune.in --recode vcf --out plink_pruning/pruning0",i,"/CeMEE.pruned.0.",i," --allow-extra-chr --double-id"))
  system(paste0("bgzip plink_pruning/pruning0",i,"/CeMEE.pruned.0.",i,".vcf"))
}

library("vcfR")
genotype_pruned <- read.vcfR('plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz')
genotype_pruned_matrix <- apply(genotype_pruned@gt[,2:1002],1,function(x){
  as.numeric(substring(x,1,1))
})
str(genotype_pruned_matrix)
genotype_pruned_matrix=2*genotype_pruned_matrix
strainID <- tstrsplit(colnames(genotype_pruned@gt)[2:1002],'_')[[1]]
GRM <- Gmatrix((genotype_pruned_matrix), method="VanRaden") # 2112 / 2112
dim(GRM)
## Re-order
GRM_sub = GRM[strainID %in% colnames(genotype_info_kept),strainID %in% colnames(genotype_info_kept)]
strainID_sub = strainID[strainID %in% colnames(genotype_info_kept)]
table(strainID_sub==colnames(genotype_info_kept)) # Same order

head(bimbam)
head(gemma_pheno)

write.table(GRM_sub,file="data_BLUE_plasticity/kinship.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

### Run Gemma
###
# Use the latest installed gemma version - use BIMBAM format
#conda activate gemma-env 

### Only takes one column at a time
#for i in 1 2
#do
#cut -f $i data_BLUE_plasticity/gemma_pheno.txt > data_BLUE_plasticity/split/gemma_pheno_T${i}.txt
#nohup gemma -g data_BLUE_plasticity/genotypes_for_gemma.txt -k data_BLUE_plasticity/kinship.txt -lmm 4 -p data_BLUE_plasticity/split/gemma_pheno_T${i}.txt -o gwas_results_BLUE_plasticity_T${i} > gemma_output_BLUE_plasticity_txt/output_gemma_t${i}.txt &
#done


mapping_T1 <- read.table("gemma_results_plasticity/gwas_results_BLUE_plasticity_T1.assoc.txt",h=TRUE)

mapping_T1$CHR=tstrsplit(mapping_T1$rs,"_")[[1]]
mapping_T1$BP=as.numeric(tstrsplit(mapping_T1$rs,"_")[[2]])

mapping_T2 <- read.table("gemma_results_plasticity/gwas_results_BLUE_plasticity_T2.assoc.txt",h=TRUE)
mapping_T2$CHR=tstrsplit(mapping_T2$rs,"_")[[1]]
mapping_T2$BP=as.numeric(tstrsplit(mapping_T2$rs,"_")[[2]])

data_manhattan = data.frame(CHR=as.numeric(factor(mapping_T1$CHR,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=mapping_T1$BP,SNP=c(1:nrow(mapping_T1)),
                            P=-log10(mapping_T1$p_wald))
data_manhattan=subset(data_manhattan,CHR!=7)
library(qqman)
manhattan(data_manhattan, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,cex=.6,suggestiveline=FALSE,genomewideline=FALSE)

##############################

data_manhattan2 = data.frame(CHR=as.numeric(factor(mapping_T2$CHR,levels=c("I","II","III","IV","V","X","MtDNA"))),BP=mapping_T2$BP,SNP=c(1:nrow(mapping_T2)),
                            P=-log10(mapping_T2$p_wald))
data_manhattan2=subset(data_manhattan2,CHR!=7)
library(qqman)
manhattan(data_manhattan2, chr="CHR", bp="BP", snp="SNP", p="P" ,logp=FALSE,cex=.6,suggestiveline=FALSE,genomewideline=FALSE)
subset(data_manhattan2,P>10)











