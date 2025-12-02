rm(list=ls())
gc()

load("data/fertility_with_n.rda")
head(ecoefs)

## load the phenotypes
BLUE1 <- read.table("BLUES_NGM.txt",sep='\t',h=TRUE)
BLUE2 <- read.table("BLUES_NaCl.txt",sep='\t',h=TRUE)

BLUE1$traitID=paste0(BLUE1$traitID,'_NGM')
BLUE2$traitID=paste0(BLUE2$traitID,'_NaCl')

all_BLUEs <- rbind(BLUE1,BLUE2)

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

gemma_pheno <- NULL

all_BLUEs$traitID <- factor(as.character(all_BLUEs$traitID),levels = 
                              paste0(rep(c("SF","SB","FS",'FB',"BS","BF",'area'),2),rep(c("_NGM","_NaCl"),each=7)))



## load the fam file
fam <- read.table("pink_binary/plink_binary.fam")
### Now we will create the phenotype file for GEMMA
for(i in fam$V1){
  #i="A6140L102"
  
  ## Is there multiple isotypes ?
  temp_iso <- subset(replaced_iso,i==Line | i==Representant_line)
  all_lines_ID <- unique(c(i,temp_iso$Line,temp_iso$Representant_line))
  
  ## Now extract all phenotypic data for these lines
  temp_pheno <- subset(all_BLUEs,pop_label%in%all_lines_ID)
  if(nrow(temp_pheno)==0){
    gemma_pheno=rbind(gemma_pheno,rep(NA,14))
  }else{
    gemma_pheno <- rbind(gemma_pheno,tapply(temp_pheno$BLUE,temp_pheno$traitID,mean,na.rm=TRUE))
  }
}
head(fam)
dim(gemma_pheno)
dim(fam)

### OK, now we need to recompute the beta vector and rotate the phenotypes ?

### Add the fertility data ?
head(ecoefs)
## Replace the isotypes ?
kept_lines <- unique(ecoefs$line[ecoefs$line%in%fam$V1])
length(kept_lines) # 204 out of 279
to_replace <- unique(ecoefs$line[!ecoefs$line%in%fam$V1])
#table(to_replace%in%c(corr_isotypes$Line,corr_isotypes$Representant_line)) # 40 can be rescued
table(to_replace%in%c(corr_isotypes$Line)) # 40 can be rescued

ecoefs$isotype=ecoefs$line
for(i in 1:nrow(ecoefs)){
  if(ecoefs$line[i]%in%corr_isotypes$Line){
    ecoefs$isotype[i] = subset(corr_isotypes,Line==ecoefs$line[i])$Representant_line
  }
}
length((ecoefs$isotype[ecoefs$isotype%in%fam$V1])) #215, ok remove duplicates

ecoefs$uniqID = paste(ecoefs$isotype,ecoefs$env,sep='_')
table(table(ecoefs$uniqID)) ## Up to 8, GM !

uniq_ecoefs <- NULL
for(i in unique(ecoefs$uniqID)){
  temp=subset(ecoefs,uniqID==i)
  if(length(uniq_ecoefs)==0){
    uniq_ecoefs=data.frame(uniqID=i,isotype=unique(temp$isotype),fertility=mean(temp$fertility),env=unique(temp$env))
  }else{
    uniq_ecoefs=rbind(uniq_ecoefs,data.frame(uniqID=i,isotype=unique(temp$isotype),fertility=mean(temp$fertility),env=unique(temp$env)))
  }
}

temp_pheno=as.data.frame(gemma_pheno)
### Merge with the traits
temp_pheno$isotype=as.character(fam$V1)

uniq_ecoefs=(merge(uniq_ecoefs,temp_pheno))

NaCl_regression = lm(fertility~SF_NaCl+SB_NaCl+FS_NaCl+FB_NaCl+BS_NaCl+BF_NaCl+area_NaCl,data=uniq_ecoefs)
vect_Beta_NaCl <- coef(NaCl_regression)[2:8]

### Now we compute the angle with the beta gradient coefficient ?
angle_theta <- function(x, y) {
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
  as.numeric(theta)
}

gemma_pheno_beta_angle  <- unlist(apply(gemma_pheno[,8:14],1,function(x){
  angle_theta(x,vect_Beta_NaCl)
}))


gemma_pheno_beta_angle[!is.na(gemma_pheno_beta_angle) & gemma_pheno_beta_angle>90]=180-gemma_pheno_beta_angle[!is.na(gemma_pheno_beta_angle) & gemma_pheno_beta_angle>90]
plot(gemma_pheno_beta_angle~gemma_pheno[,14])

gemma_pheno_fitness_effect <- rowSums(t(t(gemma_pheno[,8:14])*vect_Beta_NaCl))

plot(gemma_pheno_fitness_effect~gemma_pheno[,14])

### We will use the bimbam format which is simple for GEMMA
## First load the genotypes
## Load Cemee data
load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))
snp_positions <- snp_positions[order(snp_positions$CHROM,snp_positions$POS),]
head(snp_positions)

table(fam$V1==colnames(genotype_info)) # OK, all TRUE

bimbam <- data.frame(id=snp_positions$ID,major=snp_positions$REF,minor=snp_positions$ALT)
bimbam <- cbind(bimbam,genotype_info)
head(bimbam[,1:10])

write.table(bimbam,"data_BLUE/genotypes_for_gemma.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

fam$V6[rowSums(is.na(gemma_pheno))<14]=1
#Set to female
fam$V5=2
sort(subset(fam,V6<0)$V1)
write.table(fam,file="data_BLUE/plink_binary_corrected.fam",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(gemma_pheno,file="data_BLUE/gemma_pheno.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


### Create the kinship matrix using GEMMA
#mkdir kinship_matrix
## Rename the chromosome for GEMMA
#cut -f1 pink_binary_pruned_09/plink_binary.pruned09.bim | sort | uniq > Pruned/chrom-map_renamed.txt
#/mnt/data2/desmarais/tool_plink/plink --bfile pink_binary_pruned_09/plink_binary.pruned09 --update-chr Pruned/chrom-map_renamed.txt  --make-bed --out pink_binary_pruned_09/plink_binary.pruned09_renamed --allow-extra-chr

## Fix sex info 
#awk '{$5=2; $6=1; print}' pink_binary_pruned_09/plink_binary.pruned09.fam > pink_binary_pruned_09/plink_binary.pruned09_fixed.fam
#mv pink_binary_pruned_09/plink_binary.pruned09_fixed.fam pink_binary_pruned_09/plink_binary.pruned09.fam
#gemma -bfile pink_binary_pruned_09/plink_binary.pruned09 -gk 1 -o kinship_matrix
#mv output/kinship_matrix.* kinship_matrix/
#  /mnt/data2/desmarais/tool_plink/plink --bfile pink_binary_pruned_09/plink_binary.pruned09 --freq


### Make binary plink files for Gemma ###

# mkdir pink_binary
#/mnt/data2/desmarais/tool_plink/plink --file Pruned/plink.data --make-bed --allow-extra-chr --out pink_binary/plink_binary


### Run Gemma
###
# Use the latest installed gemma version - use BIMBAM format
#conda activate gemma-env 

### Only takes one column at a time
#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
#do
#cut -f $i data_BLUE/gemma_pheno.txt > data_BLUE/gemma_pheno_T${i}.txt
#nohup gemma -g data_BLUE/genotypes_for_gemma.txt -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p data_BLUE/gemma_pheno_T${i}.txt -o gwas_results_BLUE_T${i} > gemma_output_BLUE_txt/output_gemma_t${i}.txt &
#done

