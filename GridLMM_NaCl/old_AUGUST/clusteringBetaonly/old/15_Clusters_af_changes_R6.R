rm(list=ls());gc()
threshold_clust=0.6
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

## Polarization
polarize = data.frame(poolseq_env$GA.sync@alleles[,1:5])
names(polarize)[1:2]=c("CHROM","POS")
snps_positions.sync.caendr=merge(polarize,snps_positions.sync.caendr)
rm(polarize)
rm(poolseq_env);gc()

## Now we need to ensure the poolseq polarization match the cluster polarization
clust_pol=read.table("data_clustPolarisation_NaCl/Beta_genotypes_clusters_pol_R6.txt",h=TRUE)
head(clust_pol)
names(clust_pol)[1]="rs"
names(clust_pol)[2:3]=c("ref.GWAS","alt.GWAS")

snps_positions.sync.caendr$rs=paste(snps_positions.sync.caendr$CHROM,snps_positions.sync.caendr$POS,sep='_')
names(snps_positions.sync.caendr)[3:5]=c("ref.SYNC","major.SYNC","minor.SYNC")

merged_pos = merge(snps_positions.sync.caendr,clust_pol)
dim(merged_pos)

pol_consistent = ((merged_pos$major.SYNC==merged_pos$ref.GWAS) & (merged_pos$minor.SYNC==merged_pos$alt.GWAS)) 
pol_reversed = ((merged_pos$minor.SYNC==merged_pos$ref.GWAS) & (merged_pos$major.SYNC==merged_pos$alt.GWAS))

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

clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Beta_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)

## Add group to the merged_pos
merged_pos_clust=merge(merged_pos,unique(clusters_positions[,c("rs","group")]))
dim(merged_pos_clust) # 7.7k
length(unique(merged_pos_clust$group)) # 197
length(unique(clusters_positions$group)) # 601

all.af.caendr.clust = all.af.caendr[rownames(all.af.caendr)%in%merged_pos_clust$ID,]
all.af.caendr.clust=all.af.caendr.clust[order(rownames(all.af.caendr.clust)),]
table(rownames(all.af.caendr.clust)==merged_pos_clust$ID)

### Can we compute the frequencies of the clusters / 
clust.af=matrix(NA,nrow=length(unique(merged_pos_clust$group)),ncol=ncol(all.af.caendr))
for(i in 1:13){
  clust.af[,i]=as.numeric(tapply(all.af.caendr.clust[,i],merged_pos_clust$group,mean,na.rm=TRUE))
}
clust.af=as.data.frame(clust.af)
colnames(clust.af)=populations$population
clust.af$group=names(tapply(all.af.caendr.clust[,i],merged_pos_clust$group,mean,na.rm=TRUE))

clust_s_est=read.table("Beta_clust_s_est_06.txt",h=TRUE,sep='\t')
dim(clust_s_est)
dim(clust.af)

merged_gp <- (merge(clust.af,clust_s_est))
merged_gp$GA4_diff = (merged_gp$GA450-merged_gp$A00)
merged_gp$GA2_diff = (merged_gp$GA250-merged_gp$A00)
merged_gp$GA1_diff = (merged_gp$GA150-merged_gp$A00)

merged_gp$af.diff.mean= (merged_gp$GA1_diff+merged_gp$GA2_diff+merged_gp$GA4_diff)/3

plot(merged_gp$s,merged_gp$af.diff.mean)
subset(merged_gp,af.diff.mean>.3)

GA4_diff = (clust.af$GA450-clust.af$A00)
GA2_diff = (clust.af$GA250-clust.af$A00)
GA1_diff = (clust.af$GA150-clust.af$A00)

##############################################################################
hist(abs(GA1_diff))
table(GA1_diff>0.2 & GA2_diff>0.2 & GA4_diff>0.2) # 9
table(-GA1_diff>0.2 & -GA2_diff>0.2 & -GA4_diff>0.2) # 10
## Some clusters are selected in all the replicates in the same direction ! ##
#############################################################################
group_selected = (GA1_diff>0.2 & GA2_diff>0.2 & GA4_diff>0.2) | (-GA1_diff>0.2 & -GA2_diff>0.2 & -GA4_diff>0.2)
group_selected <- clust.af[group_selected,]$group

all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters[,c("group","rs")])
vect_P_traits_bothEnv=paste0(vect_P_traits,"_NaCl")
### Re-set the threshold to be more stringent ?
for(i in vect_P_traits_bothEnv){
  all_GWAS_clusters[,paste0("is_sign_",i)] = all_GWAS_clusters[,paste0("p_value_REML.cor_",i)]<0.05
}

merged_gp_with_clust=merge(merged_gp,all_GWAS_clusters)

plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s,pch=16,cex=.8)
abline(h=0,lty=2,col="red",lwd=2)
abline(lm(s~beta_beta_NaCl,data=merged_gp_with_clust),lwd=3,col="red")
summary(lm(s~beta_beta_NaCl,data=merged_gp_with_clust))

temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_beta_NaCl","is_sign_beta_NaCl")])
temp$GA4_diff = (temp$GA450-temp$A00)
temp$GA2_diff = (temp$GA250-temp$A00)
temp$GA1_diff = (temp$GA150-temp$A00)

sum(temp$beta_area.F_NaCl*temp$GA4_diff)
sum(temp$beta_area.F_NaCl*temp$GA2_diff)
sum(temp$beta_area.F_NaCl*temp$GA1_diff)

#plot(temp$GA4_diff~temp$beta_area.F_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
#abline(h=0,lty=2,col="orange",lwd=2)
#points(temp$GA4_diff[!temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[!temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
#points(temp$GA4_diff[temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

#plot(temp$GA2_diff~temp$beta_area.F_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
#abline(h=0,lty=2,col="orange",lwd=2)
#points(temp$GA2_diff[!temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[!temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
#points(temp$GA2_diff[temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

#plot(temp$GA1_diff~temp$beta_area.F_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
#abline(h=0,lty=2,col="orange",lwd=2)
#points(temp$GA1_diff[!temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[!temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
#points(temp$GA1_diff[temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

# aa=rep(0,10000);aa2=rep(0,10000)
# for(i in 1:10000){
#   aa[i]=sum(temp$beta_area.F_NaCl*sample(temp$GA4_diff))
#   aa2[i]=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*sample(temp$GA4_diff)[temp$is_sign_area.F_NaCl])
# }
# 
# hist(aa,xlim=c(-2.5,2.5))
# par(new=T);hist(aa2,xlim=c(-2.5,2.5),col="blue")
# abline(v=sum(temp$beta_area.F_NaCl*temp$GA4_diff),col="red",lwd=2)
# abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA4_diff[temp$is_sign_area.F_NaCl]),col="blue",lwd=2)
# 
# hist(aa2,xlim=c(-.6,.6))
# abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA4_diff[temp$is_sign_area.F_NaCl]),lwd=2,col="red")
# 
# bb=rep(0,10000)
# bb2=rep(0,10000)
# for(i in 1:10000){
#   bb[i]=sum(temp$beta_area.F_NaCl*sample(temp$GA1_diff))
#   bb2[i]=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*sample(temp$GA1_diff)[temp$is_sign_area.F_NaCl])
# }
# hist(bb,xlim=c(-2.5,2.5))
# par(new=T);hist(bb2,xlim=c(-2.5,2.5),col="blue")
# abline(v=sum(temp$beta_area.F_NaCl*temp$GA1_diff),col="red",lwd=2)
# abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA1_diff[temp$is_sign_area.F_NaCl]),col="blue",lwd=2)
# 
# 
# cc=rep(0,10000)
# cc2=rep(0,10000)
# for(i in 1:10000){
#   cc[i]=sum(temp$beta_area.F_NaCl*sample(temp$GA2_diff))
#   cc2[i]=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*sample(temp$GA2_diff)[temp$is_sign_area.F_NaCl])
# }
# hist(cc,xlim=c(-2.5,2.5))
# par(new=T);hist(cc2,xlim=c(-2.5,2.5),col="blue")
# abline(v=sum(temp$beta_area.F_NaCl*temp$GA1_diff),col="red",lwd=2)
# abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA2_diff[temp$is_sign_area.F_NaCl]),col="blue",lwd=2)
# 
# 
# pdf("Hist_change_GAs.pdf",h=5)
# 
# hist(bb,main="",xlim=c(-2,2))
# abline(v=sum(temp$beta_area.F_NaCl*temp$GA1_diff),col="red",lwd=2)
# 
# hist(cc,main="",xlim=c(-2,2))
# abline(v=sum(temp$beta_area.F_NaCl*temp$GA2_diff),col="red",lwd=2)
# 
# hist(aa,main="",xlim=c(-2,2))
# abline(v=sum(temp$beta_area.F_NaCl*temp$GA4_diff),col="red",lwd=2)
# 
# dev.off()
# 
# table(sum(temp$beta_area.F_NaCl*temp$GA1_diff)>bb)
# table(sum(temp$beta_area.F_NaCl*temp$GA2_diff)>cc)
# table(sum(temp$beta_area.F_NaCl*temp$GA4_diff)>aa)
# 
# table(sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA1_diff[temp$is_sign_area.F_NaCl])>bb2)
# table(sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA2_diff[temp$is_sign_area.F_NaCl])>cc2)
# table(sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA4_diff[temp$is_sign_area.F_NaCl])>aa2)

#### Now with the trait beta ####

temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_beta_NaCl","is_sign_beta_NaCl")])
temp$GA4_diff = (temp$GA450-temp$A00)
temp$GA2_diff = (temp$GA250-temp$A00)
temp$GA1_diff = (temp$GA150-temp$A00)

sum(temp$beta_beta_NaCl*temp$GA4_diff)
sum(temp$beta_beta_NaCl*temp$GA2_diff)
sum(temp$beta_beta_NaCl*temp$GA1_diff)


plot(temp$GA1_diff~temp$beta_beta_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$GA1_diff[!temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[!temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$GA1_diff[temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

plot(temp$GA2_diff~temp$beta_beta_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$GA2_diff[!temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[!temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$GA2_diff[temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

plot(temp$GA4_diff~temp$beta_beta_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$GA4_diff[!temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[!temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$GA4_diff[temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

temp$A00_v1=temp$A00
temp$A00_v1[temp$A00_v1>0.5]=1-temp$A00_v1[temp$A00_v1>0.5]
plot(temp$A00_v1,abs(temp$GA4_diff),col="grey",pch=16)
points(temp$A00_v1[temp$is_sign_beta_NaCl],abs(temp$GA4_diff)[temp$is_sign_beta_NaCl],col="red",pch=16)

aa=rep(0,10000);aa2=rep(0,10000)
for(i in 1:10000){
  aa[i]=sum(temp$beta_beta_NaCl*sample(temp$GA4_diff))
  aa2[i]=sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*sample(temp$GA4_diff)[temp$is_sign_beta_NaCl])
}
hist(aa,xlim=c(-2.5,2.5))
par(new=T);hist(aa2,xlim=c(-2.5,2.5),col="blue")
abline(v=sum(temp$beta_beta_NaCl*temp$GA4_diff),col="red",lwd=2)
abline(v=sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*temp$GA4_diff[temp$is_sign_beta_NaCl]),col="blue",lwd=2)

bb=rep(0,10000)
bb2=rep(0,10000)
for(i in 1:10000){
  bb[i]=sum(temp$beta_beta_NaCl*sample(temp$GA1_diff))
  bb2[i]=sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*sample(temp$GA1_diff)[temp$is_sign_beta_NaCl])
}
hist(bb,xlim=c(-2.5,2.5))
par(new=T);hist(bb2,xlim=c(-2.5,2.5),col="blue")
abline(v=sum(temp$beta_beta_NaCl*temp$GA1_diff),col="red",lwd=2)
abline(v=sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*temp$GA1_diff[temp$is_sign_beta_NaCl]),col="blue",lwd=2)

cc=rep(0,10000)
cc2=rep(0,10000)
for(i in 1:10000){
  cc[i]=sum(temp$beta_beta_NaCl*sample(temp$GA2_diff))
  cc2[i]=sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*sample(temp$GA2_diff)[temp$is_sign_beta_NaCl])
}
hist(cc,xlim=c(-2.5,2.5))
par(new=T);hist(cc2,xlim=c(-2.5,2.5),col="blue")
abline(v=sum(temp$beta_beta_NaCl*temp$GA1_diff),col="red",lwd=2)
abline(v=sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*temp$GA2_diff[temp$is_sign_beta_NaCl]),col="blue",lwd=2)

table(sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*temp$GA4_diff[temp$is_sign_beta_NaCl])>aa2)
table(sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*temp$GA1_diff[temp$is_sign_beta_NaCl])>bb2)
table(sum(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]*temp$GA2_diff[temp$is_sign_beta_NaCl])>cc2)

## Distribution of alleles 

### Contributions of large vs small effect alleles ??
vect_P_traits_bothEnv=paste0(vect_P_traits,"_NaCl")
#temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_beta_NaCl","is_sign_beta_NaCl","beta_beta_NGM","is_sign_beta_NGM")])
temp=merge(clust.af,all_GWAS_clusters[,c("group",paste0("beta_",vect_P_traits_bothEnv),paste0("is_sign_",vect_P_traits_bothEnv))])
temp$pq_A00 = temp$A00*(1-temp$A00)


plot(temp$pq_A00[temp$is_sign_beta_NaCl],I(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]^2),pch=16)

diff_GA1_10=temp$GA110-temp$A00
diff_GA1_36=temp$GA1G36-temp$GA110
diff_GA1_50=temp$GA150-temp$GA1G36

diff_GA2_10=temp$GA210-temp$A00
diff_GA2_36=temp$GA2G36-temp$GA210
diff_GA2_50=temp$GA250-temp$GA2G36

diff_GA4_10=temp$GA410-temp$A00
diff_GA4_36=temp$GA4G36-temp$GA110
diff_GA4_50=temp$GA450-temp$GA1G36

contributions_GA110 = (diff_GA1_10 * temp$beta_beta_NaCl)
contributions_GA136 = (diff_GA1_36 * temp$beta_beta_NaCl)
contributions_GA150 = (diff_GA1_50 * temp$beta_beta_NaCl)

contributions_GA210 = (diff_GA2_10 * temp$beta_beta_NaCl)
contributions_GA236 = (diff_GA2_36 * temp$beta_beta_NaCl)
contributions_GA250 = (diff_GA2_50 * temp$beta_beta_NaCl)

contributions_GA410 = (diff_GA4_10 * temp$beta_beta_NaCl)
contributions_GA436 = (diff_GA4_36 * temp$beta_beta_NaCl)
contributions_GA450 = (diff_GA4_50 * temp$beta_beta_NaCl)

plot(abs(diff_GA1_10),temp$beta_beta_NaCl)
plot((diff_GA1_10[temp$is_sign_beta_NaCl]),temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
plot((diff_GA1_36[temp$is_sign_beta_NaCl]),temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
plot((diff_GA1_50[temp$is_sign_beta_NaCl]),temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])

avg_change=(diff_GA1_50+diff_GA1_10+diff_GA1_36)+(diff_GA2_50+diff_GA2_10+diff_GA2_36)+(diff_GA4_50+diff_GA4_10+diff_GA4_36)
avg_change=avg_change/3
#pdf("Beta_change_NaCl_R6.pdf",h=4,w=5)
plot((temp$beta_beta_NaCl),(diff_GA1_50+diff_GA2_50+diff_GA4_50)/3,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),((diff_GA1_50+diff_GA2_50+diff_GA4_50)/3)[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
#dev.off()

#pdf("Beta_change_NaCl_per_replicate_R6.pdf",h=4,w=5)
plot((temp$beta_beta_NaCl),diff_GA1_50,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),diff_GA1_50[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")

plot((temp$beta_beta_NaCl),diff_GA2_50,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),diff_GA2_50[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")

plot((temp$beta_beta_NaCl),diff_GA4_50,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),diff_GA4_50[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")

#dev.off()

plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),contributions_GA110[temp$is_sign_beta_NaCl],pch=16,cex=.8)
abline(h=0,lty=2,col="red")

plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),contributions_GA136[temp$is_sign_beta_NaCl],pch=16,cex=.8)
abline(h=0,lty=2,col="red")

plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),contributions_GA150[temp$is_sign_beta_NaCl],pch=16,cex=.8)
abline(h=0,lty=2,col="red")

## Color with the starting frequencies ?
temp$pq=temp$A00*(1-temp$A00)
summary(temp$pq)
colfunc <- colorRampPalette(c("blue", "red"))
#pdf("Neda_seminar_Contribution_allGA_R6.pdf",w=5)
#par(mfrow=c(3,1))
plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),(contributions_GA110+contributions_GA136+contributions_GA150)[temp$is_sign_beta_NaCl],pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,bty="n",xlab="Absolute value of allelic effect on the trait",ylab="Contribution of allele effect on trait variation",
     main="GA1",ylim=c(-0.035,0.11))
abline(h=0,lty=2,col="red")
abline(h=c(-0.01,0.01),lty=2,col="black")

plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),(contributions_GA210+contributions_GA236+contributions_GA250)[temp$is_sign_beta_NaCl],pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,bty="n",xlab="Absolute value of allelic effect on the trait",ylab="Contribution of allele effect on trait variation",
     main="GA2",ylim=c(-0.035,0.11))
abline(h=0,lty=2,col="red")
abline(h=c(-0.01,0.01),lty=2,col="black")

plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),(contributions_GA410+contributions_GA436+contributions_GA450)[temp$is_sign_beta_NaCl],pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,bty="n",xlab="Absolute value of allelic effect on the trait",ylab="Contribution of allele effect on trait variation",
     main="GA4",ylim=c(-0.035,0.11))
abline(h=0,lty=2,col="red")
abline(h=c(-0.01,0.01),lty=2,col="black")

#dev.off()

# Averaged accross 3 replicates ?
avg_contribution=(contributions_GA110+contributions_GA136+contributions_GA150)[temp$is_sign_beta_NaCl]+
  (contributions_GA210+contributions_GA236+contributions_GA250)[temp$is_sign_beta_NaCl]+
(contributions_GA410+contributions_GA436+contributions_GA450)[temp$is_sign_beta_NaCl]
avg_contribution=avg_contribution/3

#pdf("Average_contribution_R6.pdf",h=5)
plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),avg_contribution,pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,bty="n",xlab="Absolute value of allelic effect on the trait",ylab="Contribution of allele effect on trait variation")
abline(h=0,lty=2,col="red")
abline(h=c(-0.01,0.01),lty=2,col="black")
#dev.off()

### Are the ones that are selected more/less pleiotropic ?
### In which environment ?
temp$group=as.numeric(as.character(temp$group))
selected_loci_idx = avg_contribution>0.015 & abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])>0.15
### ALT
#selected_loci_idx = avg_contribution>0.01 & abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])>0.1
selected_loci = (temp$group[temp$is_sign_beta_NaCl])[selected_loci_idx]
selected_QTL=subset(temp,group%in%selected_loci)
dim(selected_QTL) # Only 2

#plot(diff_GA1_50[temp$is_sign_beta_NaCl][selected_loci_idx],diff_GA4_50[temp$is_sign_beta_NaCl][selected_loci_idx],asp=1)
#abline(v=0);abline(h=0);abline(a=0,b=1)

top10_diff_GA1_50=selected_QTL$GA150-selected_QTL$GA1G36
top10_diff_GA2_50=selected_QTL$GA250-selected_QTL$GA2G36
top10_diff_GA4_50=selected_QTL$GA450-selected_QTL$GA4G36

dim(subset(clusters_positions,group%in%selected_QTL$group))

subset(clusters_positions,group%in%selected_QTL$group[1])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[2])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[3])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[4])[,1:3] # CHR V 6Mb


subset(clusters_positions,group%in%selected_QTL$group[5])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[6])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[7])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[8])[,1:3] # CHR V 6Mb
subset(clusters_positions,group%in%selected_QTL$group[9])[,1:3] # CHR V 6Mb












plot(temp$pq[temp$is_sign_beta_NaCl],temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]^2,pch=16)
tot_cont_High=(contributions_GA110+contributions_GA136+contributions_GA150)[temp$is_sign_beta_NaCl]>0.15

plot(temp$pq[temp$is_sign_beta_NGM],temp$beta_beta_NGM[temp$is_sign_beta_NGM]^2,pch=16)









###### OK so are these largest effect alleles reflecting the G matrix ?
### THIS COULD EXPLAIN WHY WE OBSERVE DIFFERENCES BETWEEN SELECTED/NON SELECTED ENVIRONMENTS ###

temp=merge(clust.af,all_GWAS_clusters)
temp_lowBeta = subset(temp,abs(beta_beta_NaCl)<0.15 & is_sign_beta_NaCl)
temp_highBeta = subset(temp,abs(beta_beta_NaCl)>0.15 & is_sign_beta_NaCl)

plot(cor(temp_highBeta[,paste0("beta_",vect_P_traits_bothEnv[c(1:7,16)])])[1:7,8],cor(temp_lowBeta[,paste0("beta_",vect_P_traits_bothEnv[c(1:7,16)])])[1:7,8],pch=16,asp=1)
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)


hist(cor(temp_highBeta[,paste0("beta_",vect_P_traits_bothEnv[c(1:7,16)])])[1:7,8]/cor(temp_lowBeta[,paste0("beta_",vect_P_traits_bothEnv[c(1:7,16)])])[1:7,8])

## Try to build the G matrix over time ?

#temp=merge(clust.af,all_GWAS_clusters)
#genotypes <- read.table(file=paste0("data_BLUE_clusters_CHR//genotypes_matrix_clusters_R",round(10*threshold_clust),".txt"))
#genotypes=genotypes[rownames(genotypes)%in%temp$clusterID,]
#genotypes=genotypes[order(rownames(genotypes)),]
#table(rownames(genotypes)==temp$clusterID)
  
#pred_traits_lowB=NULL
#genotypes_lowB=genotypes[abs(temp$beta_beta_NaCl)<0.15,]
#for(i in 1:16){
#  pred_traits_lowB = cbind(pred_traits_lowB,unlist(apply(genotypes_lowB,2,function(x){
#  sum(x[temp_lowBeta[,paste0("is_sign_",vect_P_traits_bothEnv[i])]]*temp_lowBeta[,paste0("beta_",vect_P_traits_bothEnv[i])][temp_lowBeta[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])})))
#}

#pred_traits_highB=NULL
#genotypes_highB=genotypes[abs(temp$beta_beta_NaCl)>0.15,]
#for(i in 1:16){
#  pred_traits_highB = cbind(pred_traits_highB,unlist(apply(genotypes_highB,2,function(x){
#    sum(x[temp_highBeta[,paste0("is_sign_",vect_P_traits_bothEnv[i])]]*temp_highBeta[,paste0("beta_",vect_P_traits_bothEnv[i])][temp_highBeta[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])})))
#}



## Predicted mean change over time

a=sum(diff_GA1_10[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
b=sum(diff_GA1_36[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
c=sum(diff_GA1_50[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])

a2=sum(diff_GA2_10[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
b2=sum(diff_GA2_36[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
c2=sum(diff_GA2_50[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])

a4=sum(diff_GA4_10[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
b4=sum(diff_GA4_36[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])
c4=sum(diff_GA4_50[temp$is_sign_beta_NaCl]*temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])

plot(c(0,10,36,50),c(0,a/10,b/26,c/14),type="b",las=1,bty="n",ylab="Predicted change per generation",xlab="generations",ylim=c(-0.01,0.1),pch=16)
points(c(0,10,36,50),c(0,a2/10,b2/26,c2/14),col="red",type="b",pch=16)
points(c(0,10,36,50),c(0,a4/10,b4/26,c4/14),col="blue",type="b",pch=16)

plot(c(0,10,36,50),c(0,a,a+b,a+b+c),type="b",las=1,bty="n",ylab="Cumulative predicted change",xlab="generations",ylim=c(0,2.6),pch=16)
points(c(0,10,36,50),c(0,a2,a2+b2,a2+b2+c2),col="red",type="b",pch=16)
points(c(0,10,36,50),c(0,a4,a4+b4,a4+b4+c4),col="blue",type="b",pch=16)


## Area
temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_area.F_NaCl","is_sign_area.F_NaCl","beta_area.F_NGM","is_sign_area.F_NGM","beta_SF_NaCl","is_sign_SF_NaCl","beta_SF_NGM","is_sign_SF_NGM")])
temp$pq_A00 = temp$A00*(1-temp$A00)

pdf("Area_effect_size_distribution.pdf",h=5)
plot(temp$pq_A00[temp$is_sign_area.F_NaCl],I(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]^2),pch=16,
     bty="n",las=1,ylab="",xlab="")
dev.off()

pdf("Area_effect_size_distribution_NGM.pdf",h=5)
plot(temp$pq_A00[temp$is_sign_area.F_NGM],I(temp$beta_area.F_NGM[temp$is_sign_area.F_NGM]^2),pch=16,
     bty="n",las=1,ylab="",xlab="")
dev.off()

pdf("Area_effect_size_distribution_BOTH.pdf",h=4,w=6)
plot(temp$pq_A00[temp$is_sign_area.F_NGM],I(temp$beta_area.F_NGM[temp$is_sign_area.F_NGM]^2),pch=16,
     bty="n",las=1,ylab="",col=adjustcolor("cornflowerblue",alpha=.6),xlab="")
points(temp$pq_A00[temp$is_sign_area.F_NaCl],I(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]^2),pch=16,
       bty="n",las=1,ylab="",xlab="",col=adjustcolor("orange",alpha=.6),ylim=c(0,0.02))
legend(0,0.02,c("Low Salt","High Salt"),col=c("cornflowerblue","orange"),pch=16)
dev.off()

contributions_GA110 = (diff_GA1_10 * temp$beta_area.F_NaCl)
contributions_GA136 = (diff_GA1_36 * temp$beta_area.F_NaCl)
contributions_GA150 = (diff_GA1_50 * temp$beta_area.F_NaCl)

contributions_GA210 = (diff_GA2_10 * temp$beta_area.F_NaCl)
contributions_GA236 = (diff_GA2_36 * temp$beta_area.F_NaCl)
contributions_GA250 = (diff_GA2_50 * temp$beta_area.F_NaCl)

contributions_GA410 = (diff_GA4_10 * temp$beta_area.F_NaCl)
contributions_GA436 = (diff_GA4_36 * temp$beta_area.F_NaCl)
contributions_GA450 = (diff_GA4_50 * temp$beta_area.F_NaCl)

avg_contribution=(contributions_GA110+contributions_GA136+contributions_GA150)[temp$is_sign_area.F_NaCl]+
  (contributions_GA210+contributions_GA236+contributions_GA250)[temp$is_sign_area.F_NaCl]+
  (contributions_GA410+contributions_GA436+contributions_GA450)[temp$is_sign_area.F_NaCl]
avg_contribution=avg_contribution/3

pdf("Average_contribution_Area.pdf",h=5)
plot(abs(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]),avg_contribution,pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,
     bty="n",xlab="Absolute value of allelic effect on the trait",
     ylab="Contribution of allele effect on trait variation")
abline(h=0,lty=2,col="red")
abline(h=c(-0.007,0.007),lty=2,col="black")
dev.off()

avg_change=(diff_GA1_50+diff_GA1_10+diff_GA1_36)+(diff_GA2_50+diff_GA2_10+diff_GA2_36)+(diff_GA4_50+diff_GA4_10+diff_GA4_36)
avg_change=avg_change/3
pdf("change_NaCl.pdf",h=4,w=5)
plot((temp$beta_area.F_NaCl),avg_change,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]),avg_change[temp$is_sign_area.F_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
dev.off()

pdf("change_NGM.pdf",h=4,w=5)
plot((temp$beta_area.F_NGM),avg_change,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_area.F_NGM[temp$is_sign_area.F_NGM]),avg_change[temp$is_sign_area.F_NGM],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
dev.off()

plot(temp$beta_area.F_NaCl^2,contributions_GA110,pch=16,cex=.8)
plot(abs(diff_GA1_10),temp$beta_area.F_NaCl)

plot((diff_GA1_10[temp$is_sign_area.F_NaCl]),temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])


## Predicted mean change over time

a=sum(diff_GA1_10[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])
b=sum(diff_GA1_36[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])
c=sum(diff_GA1_50[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])

a2=sum(diff_GA2_10[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])
b2=sum(diff_GA2_36[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])
c2=sum(diff_GA2_50[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])

a4=sum(diff_GA4_10[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])
b4=sum(diff_GA4_36[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])
c4=sum(diff_GA4_50[temp$is_sign_area.F_NaCl]*temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl])

plot(c(0,10,36,50),c(0,a/10,b/26,c/14),type="b",las=1,bty="n",ylab="Predicted change per generation",xlab="generations",ylim=c(-0.01,0.1),pch=16)
points(c(0,10,36,50),c(0,a2/10,b2/26,c2/14),col="red",type="b",pch=16)
points(c(0,10,36,50),c(0,a4/10,b4/26,c4/14),col="blue",type="b",pch=16)

plot(c(0,10,36,50),c(0,a,a+b,a+b+c),type="b",las=1,bty="n",ylab="Cumulative predicted change",xlab="generations",ylim=c(0,2.6),pch=16)
points(c(0,10,36,50),c(0,a2,a2+b2,a2+b2+c2),col="red",type="b",pch=16)
points(c(0,10,36,50),c(0,a4,a4+b4,a4+b4+c4),col="blue",type="b",pch=16)


## Try to build the G matrix over time ?
temp=merge(clust.af,all_GWAS_clusters[,c("group",paste0("beta_",vect_P_traits_bothEnv),paste0("is_sign_",vect_P_traits_bothEnv))])

pred_traits_A00=rep(0,16)
pred_traits_GA110=rep(0,16)
pred_traits_GA136=rep(0,16)
pred_traits_GA150=rep(0,16)

for(i in 1:16){
  pred_traits_A00[i]=sum((temp$A00*temp[,paste0("beta_",vect_P_traits_bothEnv[i])])[temp[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])
  pred_traits_GA110[i]=sum((temp$GA110*temp[,paste0("beta_",vect_P_traits_bothEnv[i])])[temp[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])  
  pred_traits_GA136[i]=sum((temp$GA1G36*temp[,paste0("beta_",vect_P_traits_bothEnv[i])])[temp[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])
  pred_traits_GA150[i]=sum((temp$GA150*temp[,paste0("beta_",vect_P_traits_bothEnv[i])])[temp[,paste0("is_sign_",vect_P_traits_bothEnv[i])]])
}
pred_GA1=cbind(pred_traits_A00,pred_traits_GA110,pred_traits_GA136,pred_traits_GA150)
plot(c(0,10,36,50),pred_GA1[15,])
