rm(list=ls());gc()
library(data.table)
threshold_clust=0.2
load(paste0("RData_GridLMM_NaCl/Area_GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
rm(chrID,traitID)
vect_P_traits=c("area.F")

s_est.env=new.env()
load('RData_GridLMM_NaCl/Compute_S_parallel.RData',envir=s_est.env)
ls(s_est.env)
merged_pos=s_est.env$merged_pos
s_estimates=s_est.env$results
rm(s_est.env);gc()
clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/Area_all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)

## Fix the rs label
s_estimates$rs_dot=s_estimates$rs
s_estimates$rs=paste(tstrsplit(s_estimates$rs_dot,"[.]")[[1]],tstrsplit(s_estimates$rs_dot,"[.]")[[2]],sep='_')
## Add group to the merged_pos
merged_pos_clust=merge(merged_pos,unique(clusters_positions[,c("rs","group")]))
s_estimates_clust=merge(s_estimates,unique(clusters_positions[,c("rs","group")]))

s_estimates$s=as.numeric(s_estimates$s)
s_estimates$P=as.numeric(s_estimates$P)
s2=subset(s_estimates,rs%in%all_GWAS$rs)
s2$chr=tstrsplit(s2$rs_dot,"[.]")[[1]]
s2$chr=factor(s2$chr,levels=c("I","II","III","IV","V","X"))
s2$ps=tstrsplit(s2$rs_dot,"[.]")[[2]]
### Try a plot
plot(subset(s2,abs(s2$s)<1)$s,pch=16,cex=.4,col=adjustcolor("gray",alpha=.4))
plot(-log10(subset(s2,abs(s2$s)<1)$P),pch=16,cex=.4,col=adjustcolor("gray",alpha=.4))

s3=subset(s2,abs(s2$s)<1)
data_manhattan = data.frame(CHR=as.numeric(s3$chr),BP=as.numeric(s3$ps),SNP=c(1:nrow(s3)),
                            P=-log10(p.adjust(s3$P,method="BH")))
library(qqman)
manhattan(subset(data_manhattan,P>1),chr="CHR",bp="BP",p="P",logp=FALSE)


#s_estimates_clust = s_estimates[merged_pos$rs%in%merged_pos_clust$rs,]
dim(merged_pos_clust) # 31k
dim(s_estimates_clust) # 7590
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

merged_pos$reverse=pol_reversed

s_estimates_clust=merge(s_estimates_clust,merged_pos[,c("rs","reverse")])
for(i in c("s",paste0("s_GA",1:4))){
  s_estimates_clust[,i]=as.numeric(s_estimates_clust[,i])
}
## Change reverse
s_estimates_clust$s[s_estimates_clust$reverse]= -s_estimates_clust$s[s_estimates_clust$reverse]
s_estimates_clust$s_GA1[s_estimates_clust$reverse]= -s_estimates_clust$s_GA1[s_estimates_clust$reverse]
s_estimates_clust$s_GA2[s_estimates_clust$reverse]= -s_estimates_clust$s_GA2[s_estimates_clust$reverse]
s_estimates_clust$s_GA3[s_estimates_clust$reverse]= -s_estimates_clust$s_GA3[s_estimates_clust$reverse]
s_estimates_clust$s_GA4[s_estimates_clust$reverse]= -s_estimates_clust$s_GA4[s_estimates_clust$reverse]

### Can we compute the s of the clusters ?
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
head(clust.s.est)

### Compute the s estimates from the allele frequency changes ?
plot(clust.s.est$s,-log10(clust.s.est$P),pch=16)
write.table(clust.s.est,file='Area_clust_s_est_02.txt',sep="\t",row.names=FALSE,quote=FALSE)

##############################################################################
## Some clusters are selected in all the replicates in the same direction ! ##
#############################################################################

all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters[,c("group","rs")])
vect_P_traits_bothEnv=paste0(vect_P_traits,"_NaCl")
### Re-set the threshold to be more stringent ?
for(i in vect_P_traits_bothEnv){
  all_GWAS_clusters[,paste0("is_sign_",i)] = all_GWAS_clusters[,paste0("p_value_REML.cor_",i)]<0.05
}

temp=merge(clust.s.est,all_GWAS_clusters[,c("group","beta_area.F_NaCl","is_sign_area.F_NaCl")])
head(temp$s)

plot(temp$s_GA1~temp$beta_area.F_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n',ylim=c(-.2,.2))
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$s_GA1[!temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[!temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$s_GA1[temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

plot(temp$s_GA2~temp$beta_area.F_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n',ylim=c(-.2,.2))
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$s_GA2[!temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[!temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$s_GA2[temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

plot(temp$s_GA4~temp$beta_area.F_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n',ylim=c(-.2,.2))
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$s_GA4[!temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[!temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$s_GA4[temp$is_sign_area.F_NaCl]~temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

#### Now with the trait beta ####

temp=merge(clust.s.est,all_GWAS_clusters[,c("group","beta_beta_NaCl","is_sign_beta_NaCl")])

plot(temp$s_GA4~temp$beta_beta_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$s_GA4[!temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[!temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$s_GA4[temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

plot(temp$s_GA1~temp$beta_beta_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$s_GA1[!temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[!temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$s_GA1[temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

plot(temp$s_GA2~temp$beta_beta_NaCl,pch=16,cex=.8,col=adjustcolor("grey",alpha=.8),las=1,bty="n",xlab="Effect on the trait",ylab="Allele frequency change",type='n')
abline(h=0,lty=2,col="orange",lwd=2)
points(temp$s_GA2[!temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[!temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("grey",alpha=.8))
points(temp$s_GA2[temp$is_sign_beta_NaCl]~temp$beta_beta_NaCl[temp$is_sign_beta_NaCl],pch=16,cex=.8,col=adjustcolor("red",alpha=.8))

### Contributions of large vs small effect alleles ??
vect_P_traits_bothEnv=paste0(vect_P_traits,"_NaCl")
#temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_beta_NaCl","is_sign_beta_NaCl","beta_beta_NGM","is_sign_beta_NGM")])
temp=merge(clust.s.est,all_GWAS_clusters[,c("group",paste0("beta_",vect_P_traits_bothEnv),paste0("is_sign_",vect_P_traits_bothEnv))])

#pdf("Beta_change_NaCl_R2.pdf",h=4,w=5)
plot((temp$beta_beta_NaCl),abs(temp$s),pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="S estimate")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),temp$s[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
#dev.off()

plot((temp$beta_area.F_NaCl),temp$s,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]),temp$s[temp$is_sign_area.F_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")


#pdf("Area_change_NaCl_R2.pdf",h=4,w=5)
plot((temp$beta_area.F_NaCl),(diff_GA1_50+diff_GA2_50+diff_GA4_50)/3,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]),((diff_GA1_50+diff_GA2_50+diff_GA4_50)/3)[temp$is_sign_area.F_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
#dev.off()

#pdf("Beta_change_NaCl_per_replicate_R2.pdf",h=4,w=5)
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
#pdf("Neda_seminar_Contribution_allGA_R2.pdf",w=5)
par(mfrow=c(3,1))
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

#pdf("Average_contribution_R2.pdf",h=5)
plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),avg_contribution,pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,bty="n",xlab="Absolute value of allelic effect on the trait",ylab="Contribution of allele effect on trait variation")
abline(h=0,lty=2,col="red")
abline(h=c(-0.01,0.01),lty=2,col="black")
#dev.off()

### Are the ones that are selected more/less pleiotropic ?
### In which environment ?
temp$group=as.numeric(as.character(temp$group))
selected_loci_idx = avg_contribution>0.015 & abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])>0.15
selected_loci = (temp$group[temp$is_sign_beta_NaCl])[selected_loci_idx]
selected_QTL=subset(temp,group%in%selected_loci)
dim(selected_QTL) # 7

diff_GA1_50[selected_loci_idx]

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


#pdf("Top_Beta_Area_change_NaCl_R2.pdf",h=4,w=5)
plot((temp$beta_area.F_NaCl),(diff_GA1_50+diff_GA2_50+diff_GA4_50)/3,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((selected_QTL$beta_area.F_NaCl),((top10_diff_GA1_50+top10_diff_GA2_50+top10_diff_GA4_50)/3),col="orange",pch=16,cex=1)
points((selected_QTL$beta_area.F_NaCl[selected_QTL$is_sign_area.F_NaCl]),((top10_diff_GA1_50+top10_diff_GA2_50+top10_diff_GA4_50)/3)[selected_QTL$is_sign_area.F_NaCl],col="red",pch=16)
abline(h=0,lty=2,col="red")
#dev.off()



pleio_NaCl = rowSums(selected_QTL[,paste0('is_sign_',vect_P_traits_bothEnv[1:7])])
pleio_NaCl_all = rowSums(temp[,paste0('is_sign_',vect_P_traits_bothEnv[1:7])])
table(pleio_NaCl)

#pdf("Pleiotropy_selection_NaCl_R2.pdf",h=5)
plot(c(0:7),as.numeric(table(pleio_NaCl_all))/nrow(temp),type="b",ylim=c(0,1),xlim=c(0,7),col="orange",ylab="Number of QTLs",xlab="Number of traits affected")
#dev.off()


#pdf("Pleiotropy_selection_R2.pdf",h=5)
plot(c(0:7),as.numeric(table(pleio_NaCl_all))/nrow(temp),type="b",ylim=c(0,1),xlim=c(0,7),col="orange",ylab="Number of QTLs",xlab="Number of traits affected")
points(c(1:5),as.numeric(table(pleio_NaCl))/nrow(selected_QTL),type="b",ylim=c(0,1),xlim=c(0,7),col="orange",pch=16)
#dev.off()









