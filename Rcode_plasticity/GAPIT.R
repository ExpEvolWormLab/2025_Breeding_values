library("statgenGWAS")

devtools::install_github("jiabowang/GAPIT3")

biocLite("multtest")
install.packages("gplots")
install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d

library("GAPIT")
#nohup gemma -g  -k data_BLUE/kinship.txt -lmm 4 -p data_BLUE/split/gemma_pheno_T${i}.txt -o gwas_results_BLUE_T${i} > gemma_output_BLUE_txt/output_gemma_t${i}.txt &

myY  <- read.table("data_BLUE/gemma_pheno_withLines.txt")
myY=cbind(data.frame(taxa=paste0("L",1:nrow(myY))),myY)
#myG <- read.table("data_BLUE/genotypes_for_gemma.txt", head = FALSE)
myG <- read.table("data_BLUE/genotypes_for_gemma_R08.txt", head = FALSE)
myG.num = t(myG[,4:ncol(myG)])
myG.num=cbind(data.frame(taxa=paste0("L",1:nrow(myY))),myG.num)

GM=data.frame(Name=myG$V1)
GM$Chromosome=as.numeric(factor(tstrsplit(GM$Name,"_")[[1]],levels=c("I","II","III","IV","V","X","MtDNA")))
GM$Position=as.numeric(tstrsplit(GM$Name,"_")[[2]])

myKI <- read.table("data_BLUE/kinship.txt", head = FALSE)
myKI=cbind(data.frame(taxa=paste0("L",1:nrow(myKI))),myKI)

#TEST

#myY  <- read.table("mdp_traits.txt", head = TRUE)
#myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)
#myKI <- read.table("KSN.txt", head = FALSE)

## Issue with the plot, I remove it
#myGAPIT=GAPIT(Y=myY[,1:2],GD=myG.num[,2:ncol(myG.num)],GM=GM,KI=myKI[,2:ncol(myKI)])
myGAPIT=GAPIT(Y=myY[,1:2],GD=myG.num,GM=GM,KI=myKI,  Geno.View.output=FALSE,model="gBLUP")


