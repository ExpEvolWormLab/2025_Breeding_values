rm(list=ls());gc()
library(data.table)
library(poolSeq)
library(foreach)
library(matrixStats)
library(Rcpp)
library(stringi)
library(data.table)
library(ggplot2)
library(gplots)
library(matrixStats)
library(limma)
library(qqman)
library(vcfR)

### WARNING HERE ###
setwd("~/Documents/GitHub/2025_Pichon/")

# Code to install the package from source
#install.packages("path_to_file/poolSeq-0.3.5.tar.gz",repos=NULL)

## I will load the sync files with all the sequenced populations and 
## filter it to keep only the CA populations

## Load the population names

a=(read.table('data/all_bam_list_Jan25.txt')$V1)
a=(tstrsplit(a,'_'))[[1]]
a=(tstrsplit(a,'.bam'))[[1]]


## I will now create the list of columns I want to keep in R and run a terminal command from R
## we want all the CA populations and the first three columns
## CA columns  index export ##
vect_index_kept=c(1:length(a))[c(a=="A00" | a=="A6140" | a=="A6G140" | substring(a,1,2)=="CA")]

A0_index <- paste(c(1:3,3+vect_index_kept),collapse = ",")
## Run a terminal command from R
system(paste0("cut -f ",A0_index," data_FM/All_PoolSeq_Jan_2025_polymorphisms_Cemee.sync > data/CApops_PoolSeq_Jan_2025_polymorphisms.sync"))

## Create a table with the corresponding populations
kept.pop <- vect_index_kept
populations=data.frame(population=a[kept.pop],ancestor = "A00")
populations$generations = as.numeric(tstrsplit(populations$population,"G")[[2]])
populations$generations[populations$population=="A00"]=0
populations$generations[populations$population=="A6G140"]=0
populations$replicate = substring(populations$population,1,3)
populations$replicate_num = as.numeric(as.factor(populations$replicate))
rm(kept.pop,a,A0_index)
write.table(populations,file="~/Documents/GitHub/2025_Breeding_values/Poolseq_CA_2025/populations.txt",sep='\t',row.names=FALSE,quote=FALSE)
