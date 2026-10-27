rm(list=ls());gc()
threshold_clust=0.2

## Load the selection coefficients for each SNP
select.env=new.env()
load('RData_GridLMM_NaCl/Compute_S_parallel.RData',envir=select.env)
select.results=select.env$results
select.results$P_corrected=p.adjust(select.results$P,method='BH')
table(select.results$P_corrected<0.1) # 68k SNPs
table(select.results$P_corrected<0.01) # 28k SNPs
table(select.results$P_corrected<0.00001) # 3364 SNPs
rm(select.env);gc()

### We need to re-polarize the selection according to the GWAS polarization
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
## First AF and coverage information
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

## Reverse the selection information as well - done before normally ! ->>> REMOVED
#select.results$rs=paste(tstrsplit(select.results$rs,"[.]")[[1]],
#                        tstrsplit(select.results$rs,"[.]")[[2]],sep='_')
#dim(select.results)
#select.results=merge(select.results,merged_pos[,c("rs","pol_status")])
#for(i in c("s",paste0('s_GA',c(1:4)))){
#  select.results[select.results$pol_status=="reversed",i]= -as.numeric(select.results[select.results$pol_status=="reversed",i])
#}

# Load the GWAS cluster positions
clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)
#names(clusters_positions)[names(clusters_positions)=="rs"]="ID"
#clusters_positions$rs=paste(clusters_positions$chr,clusters_positions$ps,sep='.')
#dim(clusters_positions)
dim(merge(clusters_positions,select.results[,c("rs","P_corrected")]))
clusters_positions=merge(clusters_positions,select.results[,c("rs","P_corrected","s")],all.x=TRUE)
head(clusters_positions)
clusters_positions$s=as.numeric(clusters_positions$s)
clusters_positions$P_corrected=as.numeric(clusters_positions$P_corrected)

## Remove GWAS-cluster positions not in pool-seq data
clusters_positions=subset(clusters_positions,!is.na(clusters_positions$s))
dim(clusters_positions)

plot(clusters_positions$s,-log10(clusters_positions$P_corrected))
summary(abs(subset(clusters_positions,P_corrected<0.00001)$s)) # 0.01

## Compute the median/mean/min p-values for the GWAS groups
#clusters_positions$s_abs = abs(clusters_positions$s)
## Min / Max computed
group_selection = data.frame(group=names(tapply(clusters_positions$P_corrected,clusters_positions$group,min,na.rm=TRUE)),
                             minP=as.numeric(tapply(clusters_positions$P_corrected,clusters_positions$group,min,na.rm=TRUE)),
                             medianP=as.numeric(tapply(clusters_positions$P_corrected,clusters_positions$group,median,na.rm=TRUE)),
                             meanP=as.numeric(tapply(clusters_positions$P_corrected,clusters_positions$group,mean,na.rm=TRUE)),
                             mins=as.numeric(tapply(clusters_positions$s,clusters_positions$group,min,na.rm=TRUE)),
                             medians=as.numeric(tapply(clusters_positions$s,clusters_positions$group,median,na.rm=TRUE)),
                             means=as.numeric(tapply(clusters_positions$s,clusters_positions$group,mean,na.rm=TRUE)),
                             maxs=as.numeric(tapply(clusters_positions$s,clusters_positions$group,max,na.rm=TRUE)))


### Alternative way
#group_selection = data.frame(group=names(tapply(clusters_positions$P_corrected,clusters_positions$group,min,na.rm=TRUE)))
#group_selection$maxs = NA; group_selection$minP = NA; group_selection$s_signs=NA
#for(i in 1:nrow(group_selection)){
#  temp_group = subset(clusters_positions, group==group_selection[i,]$group & P_corrected < 0.001)
#  if(nrow(temp_group)>0){
#  group_selection[i,]$s_signs = length(unique(sign(temp_group$s)))
#  group_selection[i,]$maxs = subset(temp_group,abs(temp_group$s)==max(abs(temp_group$s)))$s
#  group_selection[i,]$minP = subset(temp_group,abs(temp_group$s)==max(abs(temp_group$s)))$P_corrected
#  }else{
#  }

#}
### Half of the groups have s in opposite directions... not surprising as most of them are not significantly selected
#table(group_selection$s_signs)



# Load the GWAS cluster results
load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)

vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")
ls()

## Add group to the GWAS results, then P values (min/median/mean)
dim(all_GWAS_clusters)
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters[,c("group","clusterID")])
all_GWAS_clusters=merge(all_GWAS_clusters,group_selection)
dim(all_GWAS_clusters) # 553

plot(all_GWAS_clusters$maxs,-log10(all_GWAS_clusters$minP))
summary(abs(subset(clusters_positions,P_corrected<0.00001)$s)) # 0.01

clusters_positions$ID=paste(clusters_positions$chr,clusters_positions$ps,sep='.')

all.af.caendr.clust = all.af.caendr[rownames(all.af.caendr)%in%clusters_positions$ID,]
all.af.caendr.clust=all.af.caendr.clust[order(rownames(all.af.caendr.clust)),]

table(rownames(all.af.caendr.clust)==clusters_positions$ID)

### Can we compute the frequencies of the clusters / 
clust.af=matrix(NA,nrow=length(unique(clusters_positions$group)),ncol=ncol(all.af.caendr))
for(i in 1:13){
  clust.af[,i]=as.numeric(tapply(all.af.caendr.clust[,i],clusters_positions$group,mean,na.rm=TRUE))
}
clust.af=as.data.frame(clust.af)
colnames(clust.af)=populations$population
clust.af$group=names(tapply(all.af.caendr.clust[,i],clusters_positions$group,mean,na.rm=TRUE))

dim(clust.af)
merged_gp=merge(clust.af,clusters_positions)
dim(merged_gp)

merged_gp$GA4_diff = (merged_gp$GA450-merged_gp$A00)
merged_gp$GA2_diff = (merged_gp$GA250-merged_gp$A00)
merged_gp$GA1_diff = (merged_gp$GA150-merged_gp$A00)

merged_gp$af.diff.mean= (merged_gp$GA1_diff+merged_gp$GA2_diff+merged_gp$GA4_diff)/3

## NB: All of that if for SNPs within clusters

vect_P_traits_bothEnv=paste0(vect_P_traits,"_NaCl")

### Re-set the threshold to be more stringent ?
for(i in vect_P_traits_bothEnv){
  all_GWAS_clusters[,paste0("is_sign_",i)] = all_GWAS_clusters[,paste0("p_value_REML.cor_",i)]<0.05
}

#names(all_GWAS_clusters)[names(all_GWAS_clusters)%in%names(merged_gp)]
#head(all_GWAS_clusters[,c("group","rs","chr","ps")])
#head(merged_gp[,c("group","rs","chr","ps")])
#merged_gp_with_clust=merge(merged_gp,all_GWAS_clusters)


##all_GWAS_clusters$s=all_GWAS_clusters$mins
all_GWAS_clusters$s=all_GWAS_clusters$medians
merged_gp_with_clust=all_GWAS_clusters


pdf("Selection_R02.pdf",h=4)
par(mfrow=c(1,4))

select  <- coef(lm(s~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")

# select  <- coef(lm(s_GA1~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
# plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# 
# select  <- coef(lm(s_GA2~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
# plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# 
# select  <- coef(lm(s_GA4~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
# plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
# plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
# abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")


dev.off()


### We will now load the list of significantly selected loci to highlight them on the plot
selected_GWAS_clusters = read.table(file='Selected_GWAS_clusters.txt')$x

test_group = subset(group_selection,minP<0.00001)#$group
table(test_group$s_signs) # Still more than half with opposite signs... ?

plot(test_group$maxs,-log10(test_group$minP))
table(test_group$group%in%selected_GWAS_clusters)

table(group_selection$minP<0.00001) # 27
table(group_selection$minP<0.0001) # 42
table(group_selection$minP<0.001) # 66
table(group_selection$minP<0.01) # 123

selected_GWAS_clusters=as.numeric(as.character(subset(group_selection,minP<0.001)$group))

merged_gp_with_clust$selected_group = FALSE
merged_gp_with_clust$selected_group[merged_gp_with_clust$group%in%selected_GWAS_clusters] = TRUE
table(merged_gp_with_clust$selected_group) # 40 groups


merged_gp_with_clust$gp_color=as.numeric(as.factor(merged_gp_with_clust$selected_group))
#merged_gp_with_clust$s=merged_gp_with_clust$medians
merged_gp_with_clust$s=merged_gp_with_clust$maxs
summary(abs(merged_gp_with_clust$s))

pdf("Selection_R02_noBeta.pdf",h=4)

select_overall  <- coef(lm(s~beta_SF_NaCl+beta_SB_NaCl+beta_area.F_NaCl,data=merged_gp_with_clust))

par(mfrow=c(1,3))

select  <- coef(lm(s~beta_SF_NaCl,data=merged_gp_with_clust))
select2  <- coef(lm(s~beta_SF_NaCl,data=subset(merged_gp_with_clust,is_sign_SF_NaCl)))

plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3),
     bty="n",las=1,xlab=c("Allelic effect on trait - SF"),ylab="Selection coefficient",type="n")
abline(a=0,b=select2[2],col="darkgreen",lwd=2);abline(h=0,lty=2,col="black");abline(v=0,lty=2,col="black")
#abline(a=0,b=select[2],col="firebrick3",lwd=2);abline(a=0,b=select_overall[2],col="magenta")
temp.sign=subset(merged_gp_with_clust,!is_sign_SF_NaCl)
points(temp.sign$beta_SF_NaCl,temp.sign$s,pch=16,cex=1,col=adjustcolor(c("gray","pink")[temp.sign$gp_color],alpha=.6))
temp.sign=subset(merged_gp_with_clust,is_sign_SF_NaCl)
points(temp.sign$beta_SF_NaCl,temp.sign$s,pch=21,cex=1,bg=c("black","firebrick3")[temp.sign$gp_color])

select  <- coef(lm(s~beta_SB_NaCl,data=merged_gp_with_clust))
select2  <- coef(lm(s~beta_SB_NaCl,data=subset(merged_gp_with_clust,is_sign_SB_NaCl)))
plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3),
     bty="n",las=1,xlab=c("Allelic effect on trait - SB"),ylab="Selection coefficient",type="n")
#abline(a=0,b=select[2],col="firebrick3",lwd=2);abline(a=0,b=select_overall[3],col="magenta")
abline(a=0,b=select2[2],col="darkgreen",lwd=2);abline(h=0,lty=2,col="black");abline(v=0,lty=2,col="black")

temp.sign=subset(merged_gp_with_clust,!is_sign_SB_NaCl)
points(temp.sign$beta_SB_NaCl,temp.sign$s,pch=16,cex=1,col=adjustcolor(c("gray","pink")[temp.sign$gp_color],alpha=.6))
temp.sign=subset(merged_gp_with_clust,is_sign_SB_NaCl)
points(temp.sign$beta_SB_NaCl,temp.sign$s,pch=21,cex=1,bg=c("black","firebrick3")[temp.sign$gp_color])


select  <- coef(lm(s~beta_area.F_NaCl,data=merged_gp_with_clust))
select2  <- coef(lm(s~beta_area.F_NaCl,data=subset(merged_gp_with_clust,is_sign_area.F_NaCl)))
plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3),
     bty="n",las=1,xlab=c("Allelic effect on trait - Body Area"),ylab="Selection coefficient",type="n")
#abline(a=0,b=select[2],col="firebrick3",lwd=2);abline(a=0,b=select_overall[4],col="magenta")
abline(a=0,b=select2[2],col="darkgreen",lwd=2);abline(h=0,lty=2,col="black");abline(v=0,lty=2,col="black")

temp.sign=subset(merged_gp_with_clust,!is_sign_area.F_NaCl)
points(temp.sign$beta_area.F_NaCl,temp.sign$s,pch=16,cex=1,col=adjustcolor(c("gray","pink")[temp.sign$gp_color],alpha=.6))
temp.sign=subset(merged_gp_with_clust,is_sign_area.F_NaCl)
points(temp.sign$beta_area.F_NaCl,temp.sign$s,pch=21,cex=1,bg=c("black","firebrick3")[temp.sign$gp_color])

dev.off()

summary(lm(s~beta_SF_NaCl,data=subset(merged_gp_with_clust,is_sign_SF_NaCl)))
summary(lm(s~beta_SB_NaCl,data=subset(merged_gp_with_clust,is_sign_SB_NaCl)))
summary(lm(s~beta_area.F_NaCl,data=subset(merged_gp_with_clust,is_sign_area.F_NaCl)))

table(merged_gp_with_clust$is_sign_beta_NaCl)


pdf("Selection_R02_BETA.pdf",h=4,w=4)
par(mfrow=c(1,1))
select  <- coef(lm(s~beta_beta_NaCl,data=merged_gp_with_clust))
select2  <- coef(lm(s~beta_beta_NaCl,data=subset(merged_gp_with_clust,is_sign_beta_NaCl)))
temp.sign=subset(merged_gp_with_clust,is_sign_beta_NaCl)
plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s,pch=16,cex=.8,xlim=c(-.3,.3),
     bty="n",las=1,xlab=c("Allelic effect on trait - selection trait"),ylab="Selection coefficient",type="n")
#abline(a=0,b=select[2],col="firebrick3",lwd=2)
abline(a=0,b=select2[2],col="darkgreen",lwd=2);abline(h=0,lty=2,col="black");abline(v=0,lty=2,col="black")
temp.sign=subset(merged_gp_with_clust,!is_sign_beta_NaCl)
points(temp.sign$beta_beta_NaCl,temp.sign$s,pch=16,cex=1,col=adjustcolor(c("gray","pink")[temp.sign$gp_color],alpha=.6))
temp.sign=subset(merged_gp_with_clust,is_sign_beta_NaCl)
points(temp.sign$beta_beta_NaCl,temp.sign$s,pch=21,cex=1,bg=c("black","firebrick3")[temp.sign$gp_color])

summary(lm(s~beta_beta_NaCl,data=subset(merged_gp_with_clust,is_sign_beta_NaCl)))

dev.off()

## Compute the mean change in frequencies per groups
temp.af.diff.per.group=tapply(merged_gp$af.diff.mean,merged_gp$group,mean)
af.diff.per.group=data.frame(group=names(temp.af.diff.per.group),af.diff.mean=as.numeric(temp.af.diff.per.group))
head(af.diff.per.group)

dim(merged_gp_with_clust)
merged_gp_with_clust=merge(merged_gp_with_clust,af.diff.per.group)
dim(merged_gp_with_clust)

#clust.af$diff_GA50 = clust.af$GA150+clust.af$GA250+clust.af$GA410-clust.af$GA110-clust.af$GA210-clust.af$GA410
#dim(merged_gp_with_clust)
#merged_gp_with_clust=merge(merged_gp_with_clust,clust.af[,c("group","diff_GA50")])
#merged_gp_with_clust$contribution_change_beta = merged_gp_with_clust$diff_GA50*merged_gp_with_clust$beta_beta_NaCl
#merged_gp_with_clust$contribution_change_SF = merged_gp_with_clust$diff_GA50*merged_gp_with_clust$beta_SF_NaCl
#merged_gp_with_clust$contribution_change_SB = merged_gp_with_clust$diff_GA50*merged_gp_with_clust$beta_SB_NaCl
#merged_gp_with_clust$contribution_change_area.F = merged_gp_with_clust$diff_GA50*merged_gp_with_clust$beta_area.F_NaCl


merged_gp_with_clust$contribution_change_beta = merged_gp_with_clust$af.diff.mean*merged_gp_with_clust$beta_beta_NaCl
merged_gp_with_clust$contribution_change_SF = merged_gp_with_clust$af.diff.mean*merged_gp_with_clust$beta_SF_NaCl
merged_gp_with_clust$contribution_change_SB = merged_gp_with_clust$af.diff.mean*merged_gp_with_clust$beta_SB_NaCl
merged_gp_with_clust$contribution_change_area.F = merged_gp_with_clust$af.diff.mean*merged_gp_with_clust$beta_area.F_NaCl
temp.sign=subset(merged_gp_with_clust,is_sign_beta_NaCl)

pdf("Expected_change_per_replicates.pdf")
par(mfrow=c(1,4))
plot(abs(temp.sign$beta_SF_NaCl),temp.sign$contribution_change_SF,pch=16,ylim=c(-.06,0.06),bty="n",main="SF",xlab="Absolute GWAS Allelic effect", ylab="Contribution to trait change",
     col=c("black","firebrick3")[temp.sign$gp_color])
abline(h=0,lty=2)
plot(abs(temp.sign$beta_SB_NaCl),temp.sign$contribution_change_SB,pch=16,ylim=c(-.06,0.06),bty="n",main="SB",xlab="Absolute GWAS Allelic effect", ylab="Contribution to trait change",
     col=c("black","firebrick3")[temp.sign$gp_color])
abline(h=0,lty=2)
plot(abs(temp.sign$beta_area.F_NaCl),temp.sign$contribution_change_area.F,pch=16,ylim=c(-.06,0.06),bty="n",main="Body Area",xlab="Absolute GWAS Allelic effect", ylab="Contribution to trait change",
     col=c("black","firebrick3")[temp.sign$gp_color])
abline(h=0,lty=2)
plot(abs(temp.sign$beta_beta_NaCl),temp.sign$contribution_change_beta,pch=16,ylim=c(-.06,0.06),bty="n",main="Composite selected trait",xlab="Absolute GWAS Allelic effect", ylab="Contribution to trait change",
     col=c("black","firebrick3")[temp.sign$gp_color])
abline(h=0,lty=2)
dev.off()

pdf("Expected_change_per_replicate_BETA_only.pdf",h=4,w=4)
plot(abs(temp.sign$beta_beta_NaCl),temp.sign$contribution_change_beta,pch=16,ylim=c(-.06,0.06),bty="n",main="Composite selected trait",xlab="Absolute GWAS Allelic effect", ylab="Contribution to trait change",
     col=c("black","firebrick3")[temp.sign$gp_color])
abline(h=0,lty=2)
dev.off()


betas_TRUE = read.table("txt_BLUEs/betas_NaCl.txt",sep='\t')$V2

sum(as.numeric(subset(merged_gp_with_clust,group=="354")[,paste0('beta_',vect_P_traits[1:7],"_NaCl")])*betas_TRUE)
subset(merged_gp_with_clust,group=="354")$beta_beta_NaCl
subset(merged_gp_with_clust,group=="354")$diff_GA50
subset(merged_gp_with_clust,group=="354")$s
subset(merged_gp_with_clust,group=="354")$rs
subset(clust.af,group=="354")





