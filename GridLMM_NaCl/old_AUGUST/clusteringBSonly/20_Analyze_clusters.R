rm(list=ls());gc()
library(data.table)

threshold_clust=0.2
vect_P_traits=c("area.F")

load(paste0("RData_GridLMM_NaCl/Area_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)

## Load the s estimates
s_est.env=new.env()
load('RData_GridLMM_NaCl/Compute_S_parallel.RData',envir=s_est.env)
ls(s_est.env)
merged_pos=s_est.env$merged_pos
s_estimates=s_est.env$results
rm(s_est.env);gc()

## Load the clusters positions
clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Area_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)

## Fix the rs label
s_estimates$rs_dot=s_estimates$rs
s_estimates$rs=paste(tstrsplit(s_estimates$rs_dot,"[.]")[[1]],tstrsplit(s_estimates$rs_dot,"[.]")[[2]],sep='_')

## Add group to the merged_pos
merged_pos_clust=merge(merged_pos,unique(clusters_positions[,c("rs","group")]))
s_estimates_clust=merge(s_estimates,unique(clusters_positions[,c("rs","group")]))
dim(merged_pos_clust)
dim(s_estimates_clust)
vect_s_est=c("P","s","P_GA1","s_GA1","P_GA2","s_GA2","P_GA3","s_GA3","P_GA4","s_GA4")

### Can we compute the s of the clusters ? -> Yes but we need to re-polarize
clust_pol=read.table("data_clustPolarisation_NaCl/Area_genotypes_clusters_pol_R2.txt",h=TRUE)
names(clust_pol)[1]="rs"
names(clust_pol)[2:3]=c("ref.GWAS","alt.GWAS")

merged_pos=merged_pos[,!names(merged_pos)%in%c("ref.GWAS","alt.GWAS")]
merged_pos=merge(merged_pos,clust_pol)

pol_consistent = ((merged_pos$major.SYNC==merged_pos$ref.GWAS) & (merged_pos$minor.SYNC==merged_pos$alt.GWAS))
pol_reversed = ((merged_pos$minor.SYNC==merged_pos$ref.GWAS) & (merged_pos$major.SYNC==merged_pos$alt.GWAS))
incompatible = !(pol_consistent | pol_reversed)
table(incompatible) # all FALSE
table(pol_reversed)

## Reverse some of them
merged_pos$reverse=pol_reversed
s_estimates_clust=merge(s_estimates_clust,merged_pos[,c("rs","reverse")])
for(i in c("s",paste0("s_GA",1:4))){
  s_estimates_clust[,i]=as.numeric(s_estimates_clust$s)
}
## Change reverse
s_estimates_clust$s[s_estimates_clust$reverse]= -s_estimates_clust$s[s_estimates_clust$reverse]
s_estimates_clust$s_GA1[s_estimates_clust$reverse]= -s_estimates_clust$s_GA1[s_estimates_clust$reverse]
s_estimates_clust$s_GA2[s_estimates_clust$reverse]= -s_estimates_clust$s_GA2[s_estimates_clust$reverse]
s_estimates_clust$s_GA3[s_estimates_clust$reverse]= -s_estimates_clust$s_GA3[s_estimates_clust$reverse]
s_estimates_clust$s_GA4[s_estimates_clust$reverse]= -s_estimates_clust$s_GA4[s_estimates_clust$reverse]

### Compute the s of the clusters
clust.s.est=matrix(NA,nrow=length(unique(merged_pos_clust$group)),ncol=length(vect_s_est))
k=0
for(i in vect_s_est){
  k=k+1
  s_estimates_clust[,i]=as.numeric(s_estimates_clust[,i])
  clust.s.est[,k]=as.numeric(tapply(s_estimates_clust[,i],s_estimates_clust$group,mean,na.rm=TRUE))
}
clust.s.est=as.data.frame(clust.s.est)
names(clust.s.est)=vect_s_est
clust.s.est$group=names(tapply(s_estimates_clust[,i],s_estimates_clust$group,mean,na.rm=TRUE))

### Compute the s estimates from the allele frequency changes ?
plot(clust.s.est$s,-log10(clust.s.est$P),pch=16)
plot(s_estimates_clust$s,-log10(s_estimates_clust$P),pch=16)

write.table(clust.s.est,file='Area_clust_s_est_02.txt',sep="\t",row.names=FALSE,quote=FALSE)

### There are some issues - group 92
gp92=subset(s_estimates_clust,group==92)
dim(gp92) # 130 SNPS
plot(gp92$s,-log10(gp92$P),pch=16)

### Can we plot their trajectories ?

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
clust_pol=read.table("data_clustPolarisation_NaCl/Area_genotypes_clusters_pol_R2.txt",h=TRUE)
head(clust_pol)
names(clust_pol)[1]="rs"
names(clust_pol)[2:3]=c("ref.GWAS","alt.GWAS")

snps_positions.sync.caendr$rs=paste(snps_positions.sync.caendr$CHROM,snps_positions.sync.caendr$POS,sep='_')
names(snps_positions.sync.caendr)[3:5]=c("ref.SYNC","major.SYNC","minor.SYNC")

merged_pos = merge(snps_positions.sync.caendr,clust_pol)
dim(merged_pos)

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
sort(table(s_estimates_clust$group))

### OK, can we plot the trajectories of a specific gp  ?
gpXX=subset(s_estimates_clust,group==145)
ID_gpXX = paste(tstrsplit(gpXX$rs,"_")[[1]],tstrsplit(gpXX$rs,"_")[[2]],sep='.')
heatmap(cor(t(all.af.caendr[rownames(all.af.caendr)%in%ID_gpXX,])))
freq_gpXX=all.af.caendr[rownames(all.af.caendr)%in%ID_gpXX,]
dim(freq_gpXX) # 50

plot(freq_gpXX[1,]~populations$generations,pch=16,type='n',ylim=c(0,1))
#for(i in c(3:6)){
i=3
idx= populations$replicate_fact%in%c(1,i)
for(p in 1:nrow(freq_gpXX)){
lines(freq_gpXX[p,idx]~populations$generations[idx])
}
#}





