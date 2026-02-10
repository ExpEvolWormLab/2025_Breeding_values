rm(list=ls())
gc()
library(vcfR)
library(data.table)
library(AGHmatrix)
library(brms)


set.seed(10987345)

vect_P_traits= c("SF","SB","FS","FB","BS","BF","area")
# Load phenotypic data

BLUE1 <- read.table("BLUES_NGM.txt",sep='\t',h=TRUE)
BLUE2 <- read.table("BLUES_NaCl.txt",sep='\t',h=TRUE)

## Create a data frame

phen_NGM = data.frame(line=BLUE1$pop_label,SF=NA,SB=NA,FS=NA,FB=NA,BS=NA,BF=NA,area=NA)
phen_NaCl = data.frame(line=BLUE2$pop_label,SF=NA,SB=NA,FS=NA,FB=NA,BS=NA,BF=NA,area=NA)
for(i in 1:nrow(phen_NGM)){
  phen_NGM[i,2:8]=as.numeric(subset(BLUE1,pop_label==phen_NGM$line[i])$BLUE)
}
for(i in 1:nrow(phen_NaCl)){
phen_NaCl[i,2:8]=as.numeric(subset(BLUE2,pop_label==phen_NaCl$line[i])$BLUE)
}

phen_NGM=unique(phen_NGM); phen_NaCl=unique(phen_NaCl)
dim(phen_NGM)
dim(phen_NaCl)

######## Now estimate the G matrix with brms

# load the genotype data and create the Van Raden matrix
genotype_pruned <- read.vcfR('plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz')
genotype_pruned_matrix <- apply(genotype_pruned@gt[,2:1002],1,function(x){
  as.numeric(substring(x,1,1))
})
str(genotype_pruned_matrix)
genotype_pruned_matrix=2*genotype_pruned_matrix
strainID <- tstrsplit(colnames(genotype_pruned@gt)[2:1002],'_')[[1]]

## Reorder strain names
genotype_pruned_matrix=genotype_pruned_matrix[order(strainID),]
strainID_ordered <- strainID[order(strainID)]

GRM <- Gmatrix((genotype_pruned_matrix), method="VanRaden") # 2112 / 2112
dim(GRM)

## Re-order

GRM_NGM = GRM[strainID_ordered %in% phen_NGM$line,strainID_ordered %in% phen_NGM$line]
GRM_NaCl = GRM[strainID_ordered %in% phen_NaCl$line,strainID_ordered %in% phen_NaCl$line]

phen_NGM_sub = subset(phen_NGM,line%in%strainID)
phen_NaCl_sub = subset(phen_NaCl,line%in%strainID)

strainID_ordered_NGM = strainID_ordered[strainID_ordered %in% phen_NGM$line]
strainID_ordered_NaCl = strainID_ordered[strainID_ordered %in% phen_NaCl$line]

table(phen_NGM_sub$line==strainID_ordered_NGM) # Ok, same order 964
table(phen_NaCl_sub$line==strainID_ordered_NaCl) # Ok, same order 939

rownames(GRM_NGM)=strainID_ordered_NGM
colnames(GRM_NGM)=strainID_ordered_NGM

rownames(GRM_NaCl)=strainID_ordered_NaCl
colnames(GRM_NaCl)=strainID_ordered_NaCl

all(rownames(GRM_NGM) == phen_NGM_sub$line)
all(rownames(GRM_NaCl) == phen_NaCl_sub$line)

## With sommer

library(sommer)

fit.NaCl <- mmer(
  cbind(SF,SB,FS,FB,BS,BF,area) ~ 1,
  random = ~ vsr(line, Gu = GRM_NaCl, Gtc = unsm(7)),
  rcov   = ~ vsr(units, Gtc = unsm(7)),
  data   = phen_NaCl_sub,tolParInv=1e-2
)

Gmat_NaCl <- fit.NaCl$sigma$`line`
Rmat_NaCl <- fit.NaCl$sigma$`units`

save(list=ls(),file='RData/Gmat_sommer_NaCl.RData')
