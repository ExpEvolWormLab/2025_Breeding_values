rm(list=ls())
gc()

### We need to compute the residuals from the model with only the Kinship matrix

#library(rrBLUP)
#y_hat=NULL
#residuals_pheno=NULL
#for(i in 1:14){
#  print(i)
#  phenotypes <- read.table("data_BLUE/gemma_pheno.txt")
#  Kmat <- as.matrix(read.table("data_BLUE/kinship.txt"))
#  fit_null <- mixed.solve(y = phenotypes[,i], K = Kmat)  
#  u_hat    <- fit_null$u
#  y_hat_temp    <- u_hat + mean(phenotypes[,i], na.rm = TRUE) 
#  y_hat <- cbind(y_hat,y_hat_temp)
#  residuals_pheno <- cbind(residuals_pheno,phenotypes[,i] - y_hat_temp)
#}
#write.table(residuals_pheno,file='data_BLUE/residuals_pheno_for_permutations.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)
#write.table(y_hat,file='data_BLUE/mean_pheno_for_permutations.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)

#system("conda init")
#system("conda activate gemma-env")

## For the first trait:

#j is the number of the trait
for(j in 1:14){
  
  mean_pheno <- read.table("data_BLUE/mean_pheno_for_permutations.txt")
  resid_pheno <- read.table("data_BLUE/residuals_pheno_for_permutations.txt")
  
  run_randomized_GWAS <- function(i){
    
    #generate a rdm number for this sequence
    idx_temp <- abs(round(10000*rnorm(1)))
    # Randomize only the residual phenotype
    pheno <- mean_pheno[,j] + sample(resid_pheno[,j])
    
    # save it
    write.table(pheno,file=paste("temp_perm/pheno",idx_temp,".txt",sep='_'),quote=FALSE,row.names=FALSE,col.names=FALSE)
    gemma_cmd <- paste0("conda run -n gemma-env gemma -g data_BLUE/genotypes_for_gemma.txt -k data_BLUE/kinship.txt -lmm 1 -p ",paste("temp_perm/pheno",idx_temp,".txt",sep='_')," -o gwas_results_",idx_temp)
    system(gemma_cmd)
    p_values <- read.table(file=paste0("output/gwas_results_",idx_temp,".assoc.txt"),sep='\t',h=TRUE)[,c("chr","rs","p_wald")]
    system(paste0("rm ",paste0("output/gwas_results_",idx_temp,".assoc.txt")))
    system(paste0("rm ",paste0("output/gwas_results_",idx_temp,".log.txt")))
    system(paste0("rm ",paste0("temp_perm/pheno_",idx_temp,"_.txt")))
    return(p_values)
  }
  
  library(parallel)
  
  n_cores <- 20
  cl <- makeCluster(n_cores)
  clusterExport(cl, varlist = c("mean_pheno","resid_pheno","j"), envir = environment())
  
  temp_list=list()
  for(i in 1:100) temp_list[[i]]=i
  # Parallel computations
  List_output <-parLapply(cl, temp_list, run_randomized_GWAS)
  stopCluster(cl)
  # Stop the cluster
  
  save(list=ls(),file=paste0('Permutation_Gemma_BLUE_with_residuals/Permutation_gemma_trait_BLUE_',j,'.RData'))
  
}
