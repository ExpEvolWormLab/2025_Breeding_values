rm(list=ls());gc()

threshold_clust=0.2
load(paste0("RData_GridLMM_NaCl//GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))

rm(chrID,traitID)
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

clusters_positions <- read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep='\t',h=TRUE)

## Add group to the merged_pos
merged_pos_clust=merge(merged_pos,unique(clusters_positions[,c("rs","group")]))
dim(merged_pos_clust) # 31k
length(unique(merged_pos_clust$group)) # 553
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
clust_s_est=read.table("clust_s_est_02.txt",h=TRUE,sep='\t')

merged_gp <- (merge(clust.af,clust_s_est))
merged_gp$GA4_diff = (merged_gp$GA450-merged_gp$A00)
merged_gp$GA2_diff = (merged_gp$GA250-merged_gp$A00)
merged_gp$GA1_diff = (merged_gp$GA150-merged_gp$A00)

merged_gp$af.diff.mean= (merged_gp$GA1_diff+merged_gp$GA2_diff+merged_gp$GA4_diff)/3

all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters[,c("group","rs")])
vect_P_traits_bothEnv=paste0(vect_P_traits,"_NaCl")
### Re-set the threshold to be more stringent ?
for(i in vect_P_traits_bothEnv){
  all_GWAS_clusters[,paste0("is_sign_",i)] = all_GWAS_clusters[,paste0("p_value_REML.cor_",i)]<0.05
}

merged_gp_with_clust=merge(merged_gp,all_GWAS_clusters)


summary(lm(s~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_area.F_NaCl,data=subset(merged_gp_with_clust,is_sign_area.F_NaCl)))
summary(lm(s~beta_SF_NaCl,data=subset(merged_gp_with_clust,is_sign_SF_NaCl)))
summary(lm(s~beta_SB_NaCl,data=subset(merged_gp_with_clust,is_sign_SB_NaCl)))
summary(lm(s~beta_FS_NaCl,data=subset(merged_gp_with_clust,is_sign_FS_NaCl)))
summary(lm(s~beta_FB_NaCl,data=subset(merged_gp_with_clust,is_sign_FB_NaCl)))
summary(lm(s~beta_BF_NaCl,data=subset(merged_gp_with_clust,is_sign_BF_NaCl)))
summary(lm(s~beta_BS_NaCl,data=subset(merged_gp_with_clust,is_sign_BS_NaCl)))
summary(lm(s~beta_beta_NaCl,data=subset(merged_gp_with_clust,is_sign_beta_NaCl)))
summary(lm(s~beta_beta_NaCl,data=subset(merged_gp_with_clust,is_sign_beta_NaCl)))
summary(lm(s~beta_area.F_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_SF_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_SB_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_beta_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_FS_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_FB_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_BF_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_SB_NaCl,data=merged_gp_with_clust))
summary(lm(s~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl+beta_FS_NaCl+beta_FB_NaCl+
             beta_BS_NaCl+beta_BF_NaCl,data=merged_gp_with_clust))
summary(lm(s~(beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl)^2,data=merged_gp_with_clust))
drop1(lm(s~(beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl),data=merged_gp_with_clust),test='Chisq')


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

select  <- coef(lm(s_GA1~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s_GA1,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")

select  <- coef(lm(s_GA2~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s_GA2,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")

select  <- coef(lm(s_GA4~beta_area.F_NaCl+beta_SF_NaCl+beta_SB_NaCl,data=merged_gp_with_clust))
plot(merged_gp_with_clust$beta_area.F_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[2]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SF_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[3]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_SB_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")
plot(merged_gp_with_clust$beta_beta_NaCl,merged_gp_with_clust$s_GA4,pch=16,cex=.8,xlim=c(-.3,.3))
abline(a=0,b=select[4]);abline(h=0,lty=2,col="firebrick3");abline(v=0,lty=2,col="firebrick3")


dev.off()


### We will now load the list of significantly selected loci to highlight them on the plot
selected_GWAS_clusters = read.table(file='Selected_GWAS_clusters.txt')$x

merged_gp_with_clust$selected_group = FALSE
merged_gp_with_clust$selected_group[merged_gp_with_clust$group%in%selected_GWAS_clusters] = TRUE
table(merged_gp_with_clust$selected_group) # 40 groups, OK

merged_gp_with_clust$gp_color=as.numeric(as.factor(merged_gp_with_clust$selected_group))

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
select2  <- coef(lm(s~beta_SF_NaCl,data=subset(merged_gp_with_clust,is_sign_area.F_NaCl)))
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
abline(a=0,b=select[2],col="firebrick3",lwd=2);abline(a=0,b=select2[2],col="darkgreen",lwd=2);abline(h=0,lty=2,col="black");abline(v=0,lty=2,col="black")
temp.sign=subset(merged_gp_with_clust,!is_sign_beta_NaCl)
points(temp.sign$beta_beta_NaCl,temp.sign$s,pch=16,cex=1,col=adjustcolor(c("gray","pink")[temp.sign$gp_color],alpha=.6))
temp.sign=subset(merged_gp_with_clust,is_sign_beta_NaCl)
points(temp.sign$beta_beta_NaCl,temp.sign$s,pch=21,cex=1,bg=c("black","firebrick3")[temp.sign$gp_color])

summary(lm(s~beta_beta_NaCl,data=subset(merged_gp_with_clust,is_sign_beta_NaCl)))

dev.off()

names(merged_gp_with_clust)


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

betas_TRUE = read.table("txt_BLUEs/betas_NaCl.txt",sep='\t')$V2

sum(as.numeric(subset(merged_gp_with_clust,group=="354")[,paste0('beta_',vect_P_traits[1:7],"_NaCl")])*betas_TRUE)
subset(merged_gp_with_clust,group=="354")$beta_beta_NaCl
subset(merged_gp_with_clust,group=="354")$diff_GA50
subset(merged_gp_with_clust,group=="354")$s
subset(merged_gp_with_clust,group=="354")$rs
subset(clust.af,group=="354")

temp.env=new.env()
load('RData/Genotype_info_Cemee.RData',envir=temp.env)

table(temp.env$genotype_info[temp.env$snp_positions[,3]=="IV_2017218",])


## OK there is likely a problem of polarization here ? How do we do that ??


subset(all_GWAS_clusters,group%in%group_selected)$is_sign_area.F_NaCl

colSums(subset(all_GWAS_clusters,group%in%group_selected)[,paste0("is_sign_",paste0(vect_P_traits,"_NaCl"))])
colSums(all_GWAS_clusters[,paste0("is_sign_",paste0(vect_P_traits,"_NaCl"))])

temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_area.F_NaCl","is_sign_area.F_NaCl")])
temp$GA4_diff = (temp$GA450-temp$A00)
temp$GA2_diff = (temp$GA250-temp$A00)
temp$GA1_diff = (temp$GA150-temp$A00)

sum(temp$beta_area.F_NaCl*temp$GA4_diff)
sum(temp$beta_area.F_NaCl*temp$GA2_diff)
sum(temp$beta_area.F_NaCl*temp$GA1_diff)

aa=rep(0,10000);aa2=rep(0,10000)
for(i in 1:10000){
  aa[i]=sum(temp$beta_area.F_NaCl*sample(temp$GA4_diff))
  aa2[i]=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*sample(temp$GA4_diff)[temp$is_sign_area.F_NaCl])
}

hist(aa,xlim=c(-2.5,2.5))
par(new=T);hist(aa2,xlim=c(-2.5,2.5),col="blue")
abline(v=sum(temp$beta_area.F_NaCl*temp$GA4_diff),col="red",lwd=2)
abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA4_diff[temp$is_sign_area.F_NaCl]),col="blue",lwd=2)

hist(aa2,xlim=c(-.6,.6))
abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA4_diff[temp$is_sign_area.F_NaCl]),lwd=2,col="red")

bb=rep(0,10000)
bb2=rep(0,10000)
for(i in 1:10000){
  bb[i]=sum(temp$beta_area.F_NaCl*sample(temp$GA1_diff))
  bb2[i]=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*sample(temp$GA1_diff)[temp$is_sign_area.F_NaCl])
}
hist(bb,xlim=c(-2.5,2.5))
par(new=T);hist(bb2,xlim=c(-2.5,2.5),col="blue")
abline(v=sum(temp$beta_area.F_NaCl*temp$GA1_diff),col="red",lwd=2)
abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA1_diff[temp$is_sign_area.F_NaCl]),col="blue",lwd=2)


cc=rep(0,10000)
cc2=rep(0,10000)
for(i in 1:10000){
  cc[i]=sum(temp$beta_area.F_NaCl*sample(temp$GA2_diff))
  cc2[i]=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*sample(temp$GA2_diff)[temp$is_sign_area.F_NaCl])
}
hist(cc,xlim=c(-2.5,2.5))
par(new=T);hist(cc2,xlim=c(-2.5,2.5),col="blue")
abline(v=sum(temp$beta_area.F_NaCl*temp$GA1_diff),col="red",lwd=2)
abline(v=sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA2_diff[temp$is_sign_area.F_NaCl]),col="blue",lwd=2)


pdf("Hist_change_GAs.pdf",h=5)

hist(bb,main="",xlim=c(-2,2))
abline(v=sum(temp$beta_area.F_NaCl*temp$GA1_diff),col="red",lwd=2)

hist(cc,main="",xlim=c(-2,2))
abline(v=sum(temp$beta_area.F_NaCl*temp$GA2_diff),col="red",lwd=2)

hist(aa,main="",xlim=c(-2,2))
abline(v=sum(temp$beta_area.F_NaCl*temp$GA4_diff),col="red",lwd=2)

dev.off()

table(sum(temp$beta_area.F_NaCl*temp$GA1_diff)>bb)
table(sum(temp$beta_area.F_NaCl*temp$GA2_diff)>cc)
table(sum(temp$beta_area.F_NaCl*temp$GA4_diff)>aa)

table(sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA1_diff[temp$is_sign_area.F_NaCl])>bb2)
table(sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA2_diff[temp$is_sign_area.F_NaCl])>cc2)
table(sum(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]*temp$GA4_diff[temp$is_sign_area.F_NaCl])>aa2)

#### Now with the trait beta ####

temp=merge(clust.af,all_GWAS_clusters[,c("group","beta_beta_NaCl","is_sign_beta_NaCl")])
temp$GA4_diff = (temp$GA450-temp$A00)
temp$GA2_diff = (temp$GA250-temp$A00)
temp$GA1_diff = (temp$GA150-temp$A00)

sum(temp$beta_beta_NaCl*temp$GA4_diff)
sum(temp$beta_beta_NaCl*temp$GA2_diff)
sum(temp$beta_beta_NaCl*temp$GA1_diff)

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

plot(temp$pq_A00[temp$is_sign_area.F_NaCl],I(temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]^2),pch=16)

plot(temp$pq_A00[temp$is_sign_SF_NaCl],I(temp$beta_SF_NaCl[temp$is_sign_SF_NaCl]^2),pch=16)
plot(temp$pq_A00[temp$is_sign_SB_NaCl],I(temp$beta_SB_NaCl[temp$is_sign_SB_NaCl]^2),pch=16)
plot(temp$pq_A00[temp$is_sign_FS_NaCl],I(temp$beta_SB_NaCl[temp$is_sign_FS_NaCl]^2),pch=16)
plot(temp$pq_A00[temp$is_sign_FB_NaCl],I(temp$beta_SB_NaCl[temp$is_sign_FB_NaCl]^2),pch=16)
plot(temp$pq_A00[temp$is_sign_BS_NaCl],I(temp$beta_SB_NaCl[temp$is_sign_BS_NaCl]^2),pch=16)
plot(temp$pq_A00[temp$is_sign_BF_NaCl],I(temp$beta_SB_NaCl[temp$is_sign_BF_NaCl]^2),pch=16)

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
pdf("Beta_change_NaCl_R2.pdf",h=4,w=5)
plot((temp$beta_beta_NaCl),(diff_GA1_50+diff_GA2_50+diff_GA4_50)/3,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),((diff_GA1_50+diff_GA2_50+diff_GA4_50)/3)[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
dev.off()

pdf("Area_change_NaCl_R2.pdf",h=4,w=5)
plot((temp$beta_area.F_NaCl),(diff_GA1_50+diff_GA2_50+diff_GA4_50)/3,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_area.F_NaCl[temp$is_sign_area.F_NaCl]),((diff_GA1_50+diff_GA2_50+diff_GA4_50)/3)[temp$is_sign_area.F_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")
dev.off()

pdf("Beta_change_NaCl_per_replicate_R2.pdf",h=4,w=5)
plot((temp$beta_beta_NaCl),diff_GA1_50,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),diff_GA1_50[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")

plot((temp$beta_beta_NaCl),diff_GA2_50,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),diff_GA2_50[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")

plot((temp$beta_beta_NaCl),diff_GA4_50,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
points((temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),diff_GA4_50[temp$is_sign_beta_NaCl],col="red",pch=16,cex=.8)
abline(h=0,lty=2,col="red")

dev.off()


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
pdf("Neda_seminar_Contribution_allGA_R2.pdf",w=5)
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

dev.off()

# Averaged accross 3 replicates ?
avg_contribution=(contributions_GA110+contributions_GA136+contributions_GA150)[temp$is_sign_beta_NaCl]+
  (contributions_GA210+contributions_GA236+contributions_GA250)[temp$is_sign_beta_NaCl]+
(contributions_GA410+contributions_GA436+contributions_GA450)[temp$is_sign_beta_NaCl]
avg_contribution=avg_contribution/3

pdf("Average_contribution_R2.pdf",h=5)
plot(abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl]),avg_contribution,pch=16,cex=.8,
     col=colfunc(26)[1+floor(100*temp$pq)],las=1,bty="n",xlab="Absolute value of allelic effect on the trait",ylab="Contribution of allele effect on trait variation")
abline(h=0,lty=2,col="red")
abline(h=c(-0.01,0.01),lty=2,col="black")
dev.off()

### Are the ones that are selected more/less pleiotropic ?
### In which environment ?
temp$group=as.numeric(as.character(temp$group))
selected_loci_idx = avg_contribution>0.015 & abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])>0.15
### ALT
#selected_loci_idx = avg_contribution>0.01 & abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])>0.1
selected_loci = (temp$group[temp$is_sign_beta_NaCl])[selected_loci_idx]
selected_QTL=subset(temp,group%in%selected_loci)
dim(selected_QTL) # 7

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


pdf("Top_Beta_Area_change_NaCl_R2.pdf",h=4,w=5)
plot((temp$beta_area.F_NaCl),(diff_GA1_50+diff_GA2_50+diff_GA4_50)/3,pch=16,cex=.8,bty="n",las=1,xlab="cluster effect on the trait",ylab="Allele frequency change")
#points((selected_QTL$beta_area.F_NaCl)[1],((top10_diff_GA1_50+top10_diff_GA2_50+top10_diff_GA4_50)[1]/3),col="orange",pch=16,cex=1)
#points((selected_QTL$beta_area.F_NaCl)[9],((top10_diff_GA1_50+top10_diff_GA2_50+top10_diff_GA4_50)[9]/3),col="orange",pch=16,cex=1)
points((selected_QTL$beta_area.F_NaCl),((top10_diff_GA1_50+top10_diff_GA2_50+top10_diff_GA4_50)/3),col="orange",pch=16,cex=1)
points((selected_QTL$beta_area.F_NaCl[selected_QTL$is_sign_area.F_NaCl]),((top10_diff_GA1_50+top10_diff_GA2_50+top10_diff_GA4_50)/3)[selected_QTL$is_sign_area.F_NaCl],col="red",pch=16)
abline(h=0,lty=2,col="red")
dev.off()

pleio_NaCl = rowSums(selected_QTL[,paste0('is_sign_',vect_P_traits_bothEnv[1:7])])
pleio_NaCl_all = rowSums(temp[,paste0('is_sign_',vect_P_traits_bothEnv[1:7])])
table(pleio_NaCl)

pdf("Pleiotropy_selection_NaCl_R2.pdf",h=5)
plot(c(0:7),as.numeric(table(pleio_NaCl_all))/nrow(temp),type="b",ylim=c(0,1),xlim=c(0,7),col="orange",ylab="Number of QTLs",xlab="Number of traits affected")
dev.off()


pdf("Pleiotropy_selection_R2.pdf",h=5)
plot(c(0:7),as.numeric(table(pleio_NaCl_all))/nrow(temp),type="b",ylim=c(0,1),xlim=c(0,7),col="orange",ylab="Number of QTLs",xlab="Number of traits affected")
points(c(1:5),as.numeric(table(pleio_NaCl))/nrow(selected_QTL),type="b",ylim=c(0,1),xlim=c(0,7),col="orange",pch=16)
dev.off()











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
