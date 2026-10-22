rm(list=ls());gc()
library(data.table)
library(lme4)
library(parallel)

threshold_clust=0.2
load(paste0("RData_GridLMM_NaCl/Beta_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)
vect_P_traits=c("beta")

# load the poolseq data now
poolseq_env=new.env()
load(file='SaveGLM.RData',envir=poolseq_env)
populations=poolseq_env$populations
all.af.caendr=poolseq_env$all.af
all.cov.caendr=poolseq_env$all.cov
all.counts.caendr=poolseq_env$all.counts
all.counts.alt.caendr=poolseq_env$all.counts.alt
snps_positions.sync.caendr=poolseq_env$snps_positions.sync

polarize_for_mapping=read.table(paste0("data_clustPolarisation_NaCl/Beta_genotypes_clusters_pol_R",round(10*threshold_clust),".txt"),h=TRUE)
head(polarize_for_mapping)

## Polarization
polarize = data.frame(poolseq_env$GA.sync@alleles[,1:5])
names(polarize)[1:2]=c("CHROM","POS")
polarize$pos_labels=paste(polarize$CHROM,polarize$POS,sep='_')
snps_positions.sync.caendr=merge(polarize,snps_positions.sync.caendr)
rm(polarize)
rm(poolseq_env);gc()


names(polarize_for_mapping)[2:3]=c("ref.GWAS","alt.GWAS")
names(snps_positions.sync.caendr)[3:5]=c("ref.SYNC","major.SYNC","minor.SYNC")

merged_pos = merge(snps_positions.sync.caendr,polarize_for_mapping)

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

### We need to add the group information
### We can compute a global s and a s per replicate

clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
merged_pos=merge(merged_pos,clusters_positions[,c("CHROM","POS","group")])

compute_selection <- function(i){
## Exclude GA3
idx_noGA3 = populations$replicate_fact!=5
populations_noGA3=subset(populations,replicate!=5)

tmp.pos = subset(merged_pos,group==i)
if(nrow(tmp.pos)==1){
  count.gp=cbind(all.counts.caendr[row.names(all.counts.caendr)%in%tmp.pos$ID,idx_noGA3],
all.counts.alt.caendr[row.names(all.counts.alt.caendr)%in%tmp.pos$ID,idx_noGA3])
}else{
  count.gp=cbind(colSums(all.counts.caendr[row.names(all.counts.caendr)%in%tmp.pos$ID,idx_noGA3]),
  colSums(all.counts.alt.caendr[row.names(all.counts.alt.caendr)%in%tmp.pos$ID,idx_noGA3]))
}


# Mixed model
mod1 <- glmer(count.gp ~
                generations +(1 | replicate_fact),data=populations_noGA3,family = binomial())

mod2 <- glmer(count.gp ~
                1 +(1 | replicate_fact),data=populations_noGA3,family = binomial())
lrt <- 2 * (logLik(mod1) - logLik(mod2))
out <- c(P = pchisq(lrt, df = 1, lower.tail = FALSE),s = fixef(mod1)[2])


  # GA1/GA2/GA4 separately
for (k in c(3,4,6)){
    idx_kept = populations_noGA3$replicate_fact%in%c(1,2,k)
    mod1 <- glm(count.gp[idx_kept,] ~ populations_noGA3$generations[idx_kept],family = binomial())
    mod2 <- glm(count.gp[idx_kept,] ~ 1,family = binomial())
    
    lrt <- 2 * (logLik(mod1) - logLik(mod2))
    out <- c(out,pchisq(lrt, df = 1, lower.tail = FALSE),coef(mod1)[2])
  }
  out=c("all_gen",i,out)
  names(out) <- c("type","group","P","s","P_GA1","s_GA1","P_GA2","s_GA2","P_GA4","s_GA4")

### Now the same thing but only G10->G50
  ## Exclude GA3 and early generations
  idx_noGA3_no0 = populations$replicate_fact!=5 & populations$generations>=10
  populations_noGA3_no0=subset(populations,replicate!=5 & populations$generations>=10)
  if(nrow(tmp.pos)==1){
    count.gp=cbind(all.counts.caendr[row.names(all.counts.caendr)%in%tmp.pos$ID,idx_noGA3_no0],
                   all.counts.alt.caendr[row.names(all.counts.alt.caendr)%in%tmp.pos$ID,idx_noGA3_no0])
  }else{
    count.gp=cbind(colSums(all.counts.caendr[row.names(all.counts.caendr)%in%tmp.pos$ID,idx_noGA3_no0]),
                   colSums(all.counts.alt.caendr[row.names(all.counts.alt.caendr)%in%tmp.pos$ID,idx_noGA3_no0]))
  }
  
  # Mixed model
  mod1 <- glmer(count.gp ~
                  generations +(1 | replicate_fact),data=populations_noGA3_no0,family = binomial())
  
  mod2 <- glmer(count.gp ~
                  1 +(1 | replicate_fact),data=populations_noGA3_no0,family = binomial())
  lrt <- 2 * (logLik(mod1) - logLik(mod2))
  out2 <- c(P = pchisq(lrt, df = 1, lower.tail = FALSE),s = fixef(mod1)[2])

  # GA1/GA2/GA4 separately
  for (k in c(3,4,6)){
    idx_kept = populations_noGA3_no0$replicate_fact%in%c(1,2,k)
    mod1 <- glm(count.gp[idx_kept,] ~ populations_noGA3_no0$generations[idx_kept],family = binomial())
    mod2 <- glm(count.gp[idx_kept,] ~ 1,family = binomial())
    
    lrt <- 2 * (logLik(mod1) - logLik(mod2))
    out2 <- c(out2,pchisq(lrt, df = 1, lower.tail = FALSE),coef(mod1)[2])
    
  }
  out2=c("no0",i,out2)
  names(out2) <- c("type","group","P","s","P_GA1","s_GA1","P_GA2","s_GA2","P_GA4","s_GA4")
  
  return(rbind(out,out2))
}
  
  # Use all but one core
ncores <- 2
res <- mclapply(sort(unique(merged_pos$group)),compute_selection,mc.cores = ncores)
results <- do.call(rbind, res)
results <- as.data.frame(results)

save(list=ls(),file='RData_GridLMM_NaCl/Beta_compute_S_parallel_group_R02.RData')
write.table(results,file='RData_GridLMM_NaCl/Beta_compute_S_parallel_group_R02.txt',sep='\t',quote=FALSE,row.names=FALSE)

table(results$P<0.05)
table(subset(results,type=="no0")$P<0.05)
table(subset(results,type=="all_gen")$P<0.05)

subset(merged_pos, group==159)
r8.env=new.env()
load('RData_GridLMM_NaCl/Beta_compute_S_parallel_group_R08.RData',envir=r8.env)
sort(unique(subset(r8.env$merged_pos,pos_labels%in%subset(merged_pos, group==159)$pos_labels)$group))
sort(unique(subset(r8.env$merged_pos,pos_labels%in%subset(merged_pos, group==10)$pos_labels)$group))

a=subset(r8.env$results,group%in%(23:37))
for(i in 3:ncol(a)){
  a[,i]=as.numeric(a[,i])
}
# Merged replicates
plot(a$s,-log10(a$P),pch=16,type='n')
text(a$s,-log10(a$P),labels=a$group,col=c("black","red")[as.numeric(as.factor(a$type))])

plot(a$s,-log10(a$P),pch=16,type='n')
text(a$s_GA1,-log10(a$P_GA1),labels=a$group)

plot(a$s,-log10(a$P),pch=16,type='n')
text(a$s_GA2,-log10(a$P_GA2),labels=a$group)

plot(a$s,-log10(a$P),pch=16,type='n')
text(a$s_GA4,-log10(a$P_GA4),labels=a$group,col=c("black","red")[as.numeric(as.factor(a$type))])


results=merge(results,clusters_positions[,c("group","")])




