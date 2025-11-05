rm(list=ls())
gc()

## load the fam file

fam <- read.table("data/plink_binary.fam")

## load the phenotypes
BLUE1 <- read.table("Breeding_values_NGM.txt",sep='\t',h=TRUE)
BLUE2 <- read.table("Breeding_values_NaCl.txt",sep='\t',h=TRUE)

BV1$traitID=paste0(BV1$traitID,'_NGM')
BV2$traitID=paste0(BV2$traitID,'_NaCl')

all_BV <- rbind(BV1,BV2)

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

all_BV$traitID <- factor(as.character(all_BV$traitID),levels = 
                           paste0(rep(c("SF","SB","FS",'FB',"BS","BF",'area'),2),rep(c("_NGM","_NaCl"),each=7)))
#levels(all_BV$traitID)

### Now we will create the phenotype file for GEMMA
for(i in fam$V1){
  #i="A6140L102"
  
  ## Is there multiple isotypes ?
  temp_iso <- subset(replaced_iso,i==Line | i==Representant_line)
  all_lines_ID <- unique(c(i,temp_iso$Line,temp_iso$Representant_line))
  
  ## Now extract all phenotypic data for these lines
  temp_pheno <- subset(all_BV,pop_label%in%all_lines_ID)
  if(nrow(temp_pheno)==0){
    gemma_pheno=rbind(gemma_pheno,rep(NA,14))
  }else{
    gemma_pheno <- rbind(gemma_pheno,tapply(temp_pheno$BV,temp_pheno$traitID,mean,na.rm=TRUE))
  }
}
head(fam)
dim(gemma_pheno)
dim(fam)

fam$V6[rowSums(is.na(gemma_pheno))<14]=1
#Set to female
fam$V5=2
sort(subset(fam,V6<0)$V1)
write.table(fam,file="data/plink_binary_corrected.fam",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(gemma_pheno,file="data/gemma_pheno.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)


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

#mkdir pink_binary
#/mnt/data2/desmarais/tool_plink/plink --file Pruned/plink.data --make-bed --allow-extra-chr --out pink_binary/plink_binary


### Run Gemma
###
# Use the latest installed gemma version
#conda activate gemma-env 

### Only takes one column at a time
#cut -f 1 phenotypic_data/gemma_pheno.txt > phenotypic_data/gemma_pheno_T1.txt
#nohup gemma -bfile pink_binary/plink_binary -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p phenotypic_data/gemma_pheno_T1.txt -o gwas_results_T1 &
  
#  for i in 2 3 4 5 6 7 8 9 10 11 12 13 14
#do
#cut -f $i phenotypic_data/gemma_pheno.txt > phenotypic_data/gemma_pheno_T${i}.txt
#nohup gemma -bfile pink_binary/plink_binary -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p phenotypic_data/gemma_pheno_T${i}.txt -o gwas_results_T${i} > gemma_output_txt/output_gemma_t${i}.txt &
#  done



######## Could we do permutations ? How long would it take ?
#nohup gemma -bfile pink_binary/plink_binary -k kinship_matrix/kinship_matrix.cXX.txt -lmm 1 -p phenotypic_data/gemma_pheno_T1.txt -o gwas_results_T1 &
  