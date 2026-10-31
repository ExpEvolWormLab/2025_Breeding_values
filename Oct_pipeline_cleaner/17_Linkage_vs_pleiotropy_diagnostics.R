## 17_Linkage_vs_pleiotropy_diagnostics.R
## Within-cluster trait-specific lead SNP comparison and conditional tests.
## Important: shared GWAS association in LD DOES NOT prove causal pleiotropy.
## This screen detects resolvable *distinct* linked signals; it cannot establish
## that one molecular variant affects multiple traits.
## Conditional p-values are exploratory (lead SNP selected using same data).

traits <- c('SF','SB','FS','FB','BS','BF','area.F')
label <- 'R5'
lead_ld_r2_unresolved <- .95
min_RIL <- 40L
out_dir <- 'pleiotropy_results/R5/linkage_diagnostics'
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

p <- new.env(parent=globalenv())
load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData', envir=p)
G <- p$genotype_matrix
if(!is.matrix(G))G<-as.matrix(G)
storage.mode(G)<-'numeric'
ph <- p$phen_NaCl
ph$pop_label<-as.character(ph$pop_label)
ids <- intersect(rownames(G),ph$pop_label)
G<-G[ids,,drop=FALSE]
Y<-as.matrix(ph[match(ids,ph$pop_label),traits,drop=FALSE])
storage.mode(Y)<-'numeric'
rm(p,ph)

e <- new.env(parent=globalenv())
load('GridLMM_results_NaCl/All_Gwas.RData',envir=e)
S <- e$all_GWAS
needed <- c('rs',paste0('pvalue_T',seq_along(traits)))
if(!all(needed %in% names(S)))stop('Unexpected SNP GWAS table columns')
if(anyDuplicated(S$rs))stop('Duplicated SNP GWAS identifiers')
rownames(S) <- as.character(S$rs)

c <- new.env(parent=globalenv())
load(paste0('RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_',label,
            '_with_residuals.RData'), envir=c)
Q <- c$all_GWAS_clusters
sigs<-paste0('is_sign_',traits,'_NaCl')
if(!all(c('group',sigs) %in% names(Q)))stop('Unexpected cluster GWAS columns')
Q$group <- as.character(Q$group)
if(anyDuplicated(Q$group))stop('Duplicated LD cluster groups')

M<-read.table(paste0('txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_',label,'.txt'),
              header=TRUE,sep='\t',stringsAsFactors=FALSE,check.names=FALSE)
if(!all(c('group','rs') %in% names(M)))stop('Cluster member table needs group and rs')
M$group<-as.character(M$group)
members<-split(as.character(M$rs),M$group)

partial_p <- function(y,g1,g2) {
  ok<-is.finite(y)&is.finite(g1)&is.finite(g2)
  if(sum(ok)<min_RIL)return(c(NA_real_,NA_real_))
  x<-cbind(1,g1[ok],g2[ok]); yy<-y[ok]
  if(qr(x)$rank<3L)return(c(NA_real_,NA_real_))
  fit<-lm.fit(x,yy)
  df<-sum(ok)-3L
  if(df<=0)return(c(NA_real_,NA_real_))
  se<-sqrt(sum(fit$residuals^2)/df*diag(chol2inv(qr.R(qr(x)))))
  pp<-2*pt(-abs(fit$coefficients[2:3]/se[2:3]),df)
  as.numeric(pp)
}

result<-NULL
for(i in seq_len(nrow(Q))){
  affected<-traits[which(as.logical(unlist(Q[i,sigs],use.names=FALSE)))]
  if(length(affected)<2)next
  ss<-intersect(members[[Q$group[i]]], intersect(rownames(S),colnames(G)))
  if(length(ss)<2)next
  leads<-sapply(affected,function(tr){
    idx<-match(tr,traits)
    pv<-S[ss,paste0('pvalue_T',idx)]
    if(all(!is.finite(pv)))return(NA_character_)
    ss[which.min(replace(pv,!is.finite(pv),Inf))]
  })
  for(pair in combn(affected,2,simplify=FALSE)) {
    a<-pair[1];b<-pair[2]; s1<-leads[a];s2<-leads[b]
    if(is.na(s1)||is.na(s2))next
    x<-G[,s1];z<-G[,s2]
    ok<-is.finite(x)&is.finite(z)
    r2<-if(sum(ok)>10&&sd(x[ok])>0&&sd(z[ok])>0)cor(x[ok],z[ok])^2 else NA_real_
    ps1<-partial_p(Y[,a],x,z)
    ps2<-partial_p(Y[,b],x,z)
    chr1<-strsplit(s1,'_',fixed=TRUE)[[1]]
    chr2<-strsplit(s2,'_',fixed=TRUE)[[1]]
    distance<-if(length(chr1)>=2&&length(chr2)>=2&&chr1[1]==chr2[1])
      abs(as.numeric(chr1[2])-as.numeric(chr2[2])) else NA_real_
    status<-if(s1==s2 || (is.finite(r2) && r2>=lead_ld_r2_unresolved))
      'unresolved_high_LD_or_shared_lead' else if(
      all(is.finite(c(ps1[1],ps2[2])))&& ps1[1]<.01 && ps2[2]<.01)
      'candidate_distinct_linked_signals' else 'unresolved_no_clear_separation'
    result<-rbind(result,data.frame(group=Q$group[i],chr=as.character(Q$chr[i]),
      trait_A=a,trait_B=b,lead_A=s1,lead_B=s2,
      lead_distance_bp=distance,lead_genotype_r2=r2,
      p_A_given_B=ps1[1],p_B_given_A=ps2[2],
      n_cluster_SNPs=length(ss),classification=status,
      stringsAsFactors=FALSE))
  }
  if(i%%150==0)cat('Processed',i,'clusters\n')
}
if(!is.null(result)) {
  write.csv(result,file.path(out_dir,'17_Within_cluster_lead_comparisons.csv'),row.names=FALSE)
  tab<-as.data.frame(table(result$classification),stringsAsFactors=FALSE)
  names(tab)<-c('classification','n_trait_pairs')
  write.csv(tab,file.path(out_dir,'17_Classification_summary.csv'),row.names=FALSE)
  pdf(file.path(out_dir,'17_Lead_SNP_LD_by_class.pdf'),width=7,height=5)
  boxplot(lead_genotype_r2~classification,data=result,las=2,
          ylab='r2 between trait-specific lead SNPs',
          main='Linked signals versus unresolved pleiotropy',cex.axis=.7)
  dev.off()
  print(tab)
} else warning('No clusters with two traits and two available mapped SNPs')
cat('\nCAUTION: These are exploratory locus-resolution diagnostics.\n',
 'Distinct lead SNPs in LD are not proof of two causal variants.\n',
 'Shared peaks are not proof of true molecular pleiotropy.\n')
