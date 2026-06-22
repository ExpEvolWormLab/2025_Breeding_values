rm(list=ls());gc()

#threshold_clust=0.2
#load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
#rm(chrID,traitID)
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

# load the poolseq data now
poolseq_env=new.env()
load(file='SaveGLM.RData',envir=poolseq_env)
populations=poolseq_env$populations
all.af.caendr=poolseq_env$all.af
all.cov.caendr=poolseq_env$all.cov
all.counts.caendr=poolseq_env$all.counts
all.counts.alt.caendr=poolseq_env$all.counts.alt
snps_positions.sync.caendr=poolseq_env$snps_positions.sync

## Polarization
polarize = data.frame(poolseq_env$GA.sync@alleles[,1:5])
names(polarize)[1:2]=c("CHROM","POS")
snps_positions.sync.caendr=merge(polarize,snps_positions.sync.caendr)
rm(polarize)
rm(poolseq_env);gc()

## Now we need to ensure the poolseq polarization match the cluster polarization
temp.env=new.env()
load('RData/Genotype_info_Cemee.RData',envir=temp.env)
snp_positions=as.data.frame(temp.env$snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))
names(snp_positions)[3]="rs"
rm(temp.env);gc()

names(snp_positions)[4:5]=c("ref.GWAS","alt.GWAS")
names(snps_positions.sync.caendr)[3:5]=c("ref.SYNC","major.SYNC","minor.SYNC")

merged_pos = merge(snps_positions.sync.caendr,snp_positions)

pol_consistent = ((merged_pos$major.SYNC==merged_pos$ref.GWAS) & (merged_pos$minor.SYNC==merged_pos$alt.GWAS)) # Most are consistent
pol_reversed = ((merged_pos$minor.SYNC==merged_pos$ref.GWAS) & (merged_pos$major.SYNC==merged_pos$alt.GWAS)) # 65k are reversed

incompatible = !(pol_consistent | pol_reversed)
table(incompatible)

merged_pos$pol_status=NA
merged_pos$pol_status[pol_consistent]="consistent"
merged_pos$pol_status[pol_reversed]="reversed"
merged_pos$pol_status[incompatible]="incompatible" # only ~200 positions
rm(incompatible,pol_consistent,pol_reversed)

### remove incompatible and reverse
merged_pos=subset(merged_pos,pol_status!="incompatible")
all.af.caendr=all.af.caendr[rownames(all.af.caendr)%in%merged_pos$ID,]
all.af.caendr[rownames(all.af.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]=1-all.af.caendr[rownames(all.af.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]

all.cov.caendr=all.cov.caendr[rownames(all.cov.caendr)%in%merged_pos$ID,]
all.counts.caendr=all.counts.caendr[rownames(all.counts.caendr)%in%merged_pos$ID,]
all.counts.alt.caendr=all.counts.alt.caendr[rownames(all.counts.alt.caendr)%in%merged_pos$ID,]

all.counts.alt.caendr[rownames(all.counts.alt.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]=
  all.counts.caendr[rownames(all.counts.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]
  
all.counts.caendr[rownames(all.counts.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]=
  all.cov.caendr[rownames(all.cov.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]-
  all.counts.alt.caendr[rownames(all.counts.alt.caendr)%in%subset(merged_pos,pol_status=="reversed")$ID,]

### We can compute a global s and a s per replicate
i=1


library(lme4)
library(parallel)

# Precompute indices once
idx_list <- lapply(3:6, function(j)
  populations$replicate %in% c(1, 2, j))

pop_list <- lapply(idx_list, function(idx)
  populations[idx, ])

library(lme4)
library(parallel)

# Precompute indices once
idx_list <- lapply(3:6, function(j)
  populations$replicate %in% c(1, 2, j))

pop_list <- lapply(idx_list, function(idx)
  populations[idx, ])

compute_s_from_counts <- function(i) {
  
  # Mixed model
  mod1 <- glmer(cbind(all.counts.caendr[i, ],all.counts.alt.caendr[i, ]) ~
      populations$generations +(1 | populations$replicate_fact),family = binomial())
  
  mod2 <- glmer(cbind(all.counts.caendr[i, ],all.counts.alt.caendr[i, ]) ~
                  1 +(1 | populations$replicate_fact),family = binomial())
  
  lrt <- 2 * (logLik(mod1) - logLik(mod2))
  
  out <- c(P = pchisq(lrt, df = 1, lower.tail = FALSE),s = fixef(mod1)[2])
  
  # GA1-GA4 separately
  for (k in seq_along(idx_list)) {
    idx <- idx_list[[k]]
    data <- pop_list[[k]]
    mod1 <- glm(cbind(all.counts.caendr[i, idx],all.counts.alt.caendr[i, idx]) ~ 
                  generations,family = binomial(),data = data)
    
    mod2 <- glm(cbind(all.counts.caendr[i, idx],all.counts.alt.caendr[i, idx]) ~ 1,
      family = binomial(),data = data)
    
    lrt <- 2 * (logLik(mod1) - logLik(mod2))
    out <- c(out,pchisq(lrt, df = 1, lower.tail = FALSE),coef(mod1)[2])
  }
  
  names(out) <- c("P","s","P_GA1","s_GA1","P_GA2","s_GA2","P_GA3","s_GA3","P_GA4","s_GA4")
  return(out)
}

# Use all but one core
ncores <- 20

res <- mclapply(
seq_len(nrow(all.counts.caendr)),compute_s_from_counts,mc.cores = ncores)

results <- do.call(rbind, res)
results <- as.data.frame(results)
save(list=ls(),file='Compute_S_parallel.RData')
