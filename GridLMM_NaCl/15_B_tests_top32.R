

selected_loci_idx = avg_contribution>0.01 & abs(temp$beta_beta_NaCl[temp$is_sign_beta_NaCl])>0.1
selected_loci = (temp$group[temp$is_sign_beta_NaCl])[selected_loci_idx]
selected_QTL=subset(temp,group%in%selected_loci)
dim(selected_QTL) # 21 QTLs now 

gp_gmat = tstrsplit(rownames(genotyped_matrix),"_")[[1]]

heatmap(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,])))

## 1 is alone
sort(round(100*abs(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,]))[,1])))
## 2 and 3 are clustered
sort(round(100*abs(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,]))[,2])))
sort(round(100*abs(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,]))[,3])))
## All 7 also correlated ? / 3 neg / 3 pos
sort(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,]))[,4])
sort(round(100*abs(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,]))[,5])))
sort(round(100*abs(cor(t(genotyped_matrix[as.numeric(gp_gmat)%in%selected_QTL$group,]))[,6])))


