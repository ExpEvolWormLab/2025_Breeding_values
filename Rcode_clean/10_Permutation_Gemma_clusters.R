rm(list=ls())
gc()

system("conda activate gemma-env")

## For the first trait:
for(threshold_clust in c(0.2,0.4,0.6,0.8,0.9,1)){
  print(threshold_clust)
  #j is the number of the trait
  
  for(j in 1:14){
  print(j)    
    phenotypes <- read.table("data_BLUE/gemma_pheno.txt")
    
    run_randomized_GWAS <- function(i){
      #generate a rdm number for this sequence
      idx_temp <- abs(round(10000*rnorm(1)))
      # Randomize the phenotype
      pheno <- sample(phenotypes[,j])
      # save it
      write.table(pheno,file=paste("temp_perm/pheno",idx_temp,"_clusters.txt",sep='_'),quote=FALSE,row.names=FALSE,col.names=FALSE)
      gemma_cmd <- paste0("gemma -g data_BLUE_clusters/genotypes_clusters_R",round(10*threshold_clust),".txt -k data_BLUE/kinship.txt -lmm 1 -p ",paste("temp_perm/pheno",idx_temp,"_clusters.txt",sep='_')," -o gwas_results_clusters_BLUE",idx_temp)
      system(gemma_cmd)
      p_values <- read.table(file=paste0("output/gwas_results_clusters_BLUE",idx_temp,".assoc.txt"),sep='\t',h=TRUE)[,c("chr","rs","p_wald")]
      system(paste0("rm ",paste0("output/gwas_results_clusters_BLUE",idx_temp,".assoc.txt")))
      system(paste0("rm ",paste0("output/gwas_results_clusters_BLUE",idx_temp,".log.txt")))
      system(paste0("rm ",paste0("temp_perm/pheno_",idx_temp,"__clusters.txt")))
      return(p_values)
    }
  
    library(parallel)
    
    n_cores <- 10
    cl <- makeCluster(n_cores)
    clusterExport(cl, varlist = c("phenotypes","j"), envir = environment())
    
    temp_list=list()
    for(i in 1:100) temp_list[[i]]=i
    # Parallel computations
    List_output <-parLapply(cl, temp_list, run_randomized_GWAS)
    stopCluster(cl)
    # Stop the cluster
    
    save(list=ls(),file=paste0('Permutation_gemma_trait_clusters_BLUE_R',round(10*threshold_clust),'_',j,'.RData'))
    
  }
}