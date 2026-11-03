## 57 — Chromosome V local RIL haplotypes and trait associations
## Run from the project root, after scripts 54–56.
## Only homozygous RIL marker calls are treated as resolved 4-marker haplotypes.
## This is NOT a phased-founder ancestry inference or a causal fine-map.
options(stringsAsFactors=FALSE)
set.seed(20261102L)
out <- 'selection_results/chrV_deep_dive/57_haplotypes'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_out <- function(x,n)write.csv(x,file.path(out,n),row.names=FALSE,na='')
first <- function(paths){p<-paths[file.exists(paths)];if(length(p))p[[1]]else NA_character_}
required <- function(paths,desc){z<-first(paths);if(is.na(z))stop('Missing ',desc,': ',paste(paths,collapse=' ; '),call.=FALSE);message(desc,': ',z);z}
read_tab <- function(p)read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
trait_names <- c('SF','SB','FS','FB','BS','BF','area.F')
focal_groups <- c('848','851','853','854')
min_class_n <- 10L
n_background_pcs <- 5L

A <- read_tab(required(c('selection_results/chrV_deep_dive/54_V3_12_cluster_prediction_audit.csv'),
                       'script 54 prediction audit'))
S <- read_tab(required(c('selection_results/chrV_deep_dive/56_Exact_SNP_ALT_and_phase_aligned_trajectories.csv'),
                       'script 56 allele counts'))
L <- read_tab(required(c('selection_results/chrV_deep_dive/55_Selected_RIL_markers.csv'),
                       'script 55 marker list'))
prep <- required(c('RData/Preparatory_GLMM_NaCl_only_withBeta.RData'),
                 'RIL phenotypes and genotypes')
if(!all(c('group','neutral_test_tag') %in% names(A)))stop('54 audit missing group/neutral_test_tag')
if(!all(c('SNP_ID','group','sample_id','ALT_frequency') %in% names(S)))
  stop('56 table missing required columns')
if(!all(c('rs','group','phase')%in%names(L)))stop('55 marker list malformed')
A$group <- as.character(A$group); L$group<-as.character(L$group); S$group<-as.character(S$group)
if(anyDuplicated(A$group))stop('Duplicate cluster IDs in 54 audit')
ix <- match(focal_groups,A$group)
if(anyNA(ix))stop('One of groups 848/851/853/854 is absent from 54 audit')
tags <- data.frame(group=focal_groups,SNP_ID=as.character(A$neutral_test_tag[ix]))
if(anyNA(tags$SNP_ID)||anyDuplicated(tags$SNP_ID))stop('Missing or repeated focal tags')
tags$phase <- as.numeric(L$phase[match(tags$SNP_ID,L$rs)])
if(any(!tags$phase%in%c(-1,1)))stop('No unambiguous phase for focal tag')
base <- S[S$sample_id=='A00',c('SNP_ID','ALT_frequency'),drop=FALSE]
if(anyDuplicated(base$SNP_ID))stop('More than one common A00 baseline per site')
tags$A00_ALT_freq <- as.numeric(base$ALT_frequency[match(tags$SNP_ID,base$SNP_ID)])
if(any(!is.finite(tags$A00_ALT_freq))||any(tags$A00_ALT_freq==0.5))
  stop('Cannot determine ancestral minor state for all tags')
## Biological allele orientation: 1 = allele rare in A00 (GWAS ALT or REF).
tags$rare_allele <- ifelse(tags$A00_ALT_freq<0.5,'ALT','REF')
tags$rarity_vs_phase <- ifelse((tags$rare_allele=='ALT' & tags$phase==1) |
                                  (tags$rare_allele=='REF' & tags$phase==-1),
                              'same','opposite')
write_out(tags,'57_Four_anchor_marker_and_ancestral_orientation_QC.csv')

P <- new.env(parent=globalenv());load(prep,envir=P)
if(!all(c('genotype_matrix','phen_NaCl')%in%ls(P)))stop('Preparation RData missing genotype_matrix/phen_NaCl')
G <- P$genotype_matrix
if(is.data.frame(G))G <- as.matrix(G)
if(!is.matrix(G)||is.null(rownames(G))||is.null(colnames(G)))
  stop('genotype_matrix must be RIL x GWAS-ALT SNP with row and column IDs')
ph <- P$phen_NaCl
if(!all(c('pop_label',trait_names)%in%names(ph)))stop('phen_NaCl lacks expected traits')
ph$pop_label <- as.character(ph$pop_label)
if(anyDuplicated(rownames(G))||anyDuplicated(ph$pop_label))stop('Duplicate RIL phenotype/genotype IDs')
ids <- intersect(rownames(G),ph$pop_label)
if(length(ids)<100L)stop('Fewer than 100 RILs match phenotype and genotype records')
if(!all(tags$SNP_ID%in%colnames(G)))stop('Missing anchor markers: ',paste(setdiff(tags$SNP_ID,colnames(G)),collapse=', '))
## Extract only relevant genotype columns; do not convert the genome-wide matrix.
focus <- unique(c(tags$SNP_ID,'V_6300718'))
focus <- focus[focus%in%colnames(G)]
Hdosage <- G[ids,focus,drop=FALSE];storage.mode(Hdosage)<-'double'
if(any(is.finite(Hdosage)&!(Hdosage%in%c(0,1,2))))stop('Non-diploid dosages in focal SNPs')
Y <- as.matrix(ph[match(ids,ph$pop_label),trait_names,drop=FALSE]);storage.mode(Y)<-'double'
rownames(Y)<-ids

## Kinship PCs for descriptive adjustment. They do NOT replace a full genomic
## mixed model; note that a GRM may also absorb some local haplotype signal.
pcs <- NULL
if('GRM_sub'%in%ls(P)){
  K <- P$GRM_sub
  if(is.matrix(K)&&!is.null(rownames(K))&&all(ids%in%rownames(K))&&
     all(ids%in%colnames(K))){
    K <- K[ids,ids,drop=FALSE]
    if(max(abs(K-t(K)),na.rm=TRUE)<1e-6&&all(is.finite(K))){
      eg <- eigen((K+t(K))/2,symmetric=TRUE)
      kk <- min(n_background_pcs,sum(eg$values>1e-8))
      if(kk>0){
        pcs <- eg$vectors[,seq_len(kk),drop=FALSE]
        pcs <- sweep(pcs,2,pmax(apply(pcs,2,sd),1e-10),'/')
        colnames(pcs) <- paste0('kinship_PC',seq_len(kk))
        rownames(pcs) <- ids
      }
    }
  }
}
if(is.null(pcs))warning('No usable GRM_sub: only unadjusted trait models will be reported')
## Optional broad local RIL backgrounds from ONE representative per R5 group
## in V:6--8 Mb. This avoids overweighting 844/858 simply because they
## contain more association SNPs. These are UNSUPERVISED genotype patterns,
## not phased founding haplotypes or estimated local ancestry probabilities.
L <- L[!duplicated(L$rs),,drop=FALSE]
reps <- L[!duplicated(L$group),c('group','rs'),drop=FALSE]
for(k in seq_len(nrow(tags))){
  j<-match(tags$group[k],reps$group)
  if(!is.na(j))reps$rs[j]<-tags$SNP_ID[k]
}
reps <- reps[reps$rs%in%colnames(G),,drop=FALSE]
background <- rep(NA_character_,length(ids))
if(nrow(reps)>=8L){
  ZZ <- G[ids,as.character(reps$rs),drop=FALSE];storage.mode(ZZ)<-'double'
  if(any(is.finite(ZZ)&!(ZZ%in%c(0,1,2))))stop('Unexpected dosage in broad-background panel')
  ZZ[is.finite(ZZ)&ZZ==1]<-NA_real_ # heterozygotes are unphased
  ZZ[is.finite(ZZ)&ZZ==2]<-1
  missing_frac<-rowMeans(!is.finite(ZZ))
  valid <- missing_frac<=0.2
  if(sum(valid)>=80L){
    local <- ZZ[valid,,drop=FALSE]
    for(j in seq_len(ncol(local))){
      p <- mean(local[,j],na.rm=TRUE)
      if(!is.finite(p))stop('No called genotypes at background representative')
      local[!is.finite(local[,j]),j]<-as.integer(p>0.5)
    }
    hc <- hclust(dist(local,method='manhattan')/ncol(local),method='average')
    cl <- cutree(hc,h=0.15)  # fixed pre-analysis Hamming distance threshold
    background[valid] <- sprintf('BG%03d',cl)
    write_out(reps,'57_Broad_local_background_representative_SNPs.csv')
    write_out(data.frame(pop_label=ids,local_background_proxy=background,
                         missing_fraction_representatives=missing_frac),
              '57_Broad_local_background_RIL_assignments.csv')
  }
}
rm(G,P);invisible(gc())

M <- Hdosage[,tags$SNP_ID,drop=FALSE]
## Homozygous dosage 2 = ALT present; 0 = REF present. Exclude heterozygotes.
B <- matrix(NA_integer_,nrow=nrow(M),ncol=ncol(M),dimnames=dimnames(M))
for(j in seq_len(ncol(M))){
  valid <- is.finite(M[,j]) & M[,j]%in%c(0,2)
  B[valid,j] <- if(tags$rare_allele[j]=='ALT')as.integer(M[valid,j]==2)
                 else as.integer(M[valid,j]==0)
}
bit_sig <- apply(B,1,function(z)if(anyNA(z))NA_character_ else paste0(z,collapse=''))
pattern <- data.frame(pop_label=ids,rare_haplotype_4tag=bit_sig,
                       complete_4tag=!is.na(bit_sig),local_background_proxy=background,
                       stringsAsFactors=FALSE)
for(j in seq_len(ncol(B)))pattern[[paste0('rare_',tags$group[j])]] <- B[,j]
write_out(pattern,'57_RIL_four_marker_haplotype_assignments.csv')
pattern_counts <- sort(table(bit_sig,useNA='no'),decreasing=TRUE)
counts <- data.frame(rare_haplotype_4tag=names(pattern_counts),n_RIL=as.integer(pattern_counts),
                     frequency_among_resolved=as.numeric(pattern_counts)/sum(pattern_counts))
write_out(counts,'57_RIL_four_marker_haplotype_frequencies.csv')
if(any(!is.na(background))){
  cross <- as.data.frame(table(local_background_proxy=background,
                               rare_haplotype_4tag=bit_sig,useNA='no'))
  names(cross)[3] <- 'n_RIL'
  write_out(cross,'57_Background_vs_four_tag_pattern_counts.csv')
}
if(sum(!is.na(bit_sig))<80L)stop('Too few fully resolved RIL four-marker haplotypes')

## Regional LD from ALL homozygous RILs, not independent selected sites.
ld <- do.call(rbind,lapply(seq_len(ncol(B)-1L),function(i){
  do.call(rbind,lapply((i+1L):ncol(B),function(j){
    ok<-is.finite(B[,i])&is.finite(B[,j]);n<-sum(ok)
    rr<-if(n>25L&&sd(B[ok,i])>0&&sd(B[ok,j])>0)cor(B[ok,i],B[ok,j])else NA_real_
    data.frame(group_A=tags$group[i],group_B=tags$group[j],n_complete=n,
               correlation_ancestral_rare=rr,r2=rr^2)
  }))
}))
write_out(ld,'57_Anchor_marker_pairwise_LD.csv')

## Test each of seven *existing* adjusted RIL BLUEs, not individual recordings.
## Haplotype classes occurring <10 times are grouped as OTHER; do not interpret
## pooled OTHER as a biological founder haplotype.
count_by_id <- setNames(counts$n_RIL,counts$rare_haplotype_4tag)
classes <- ifelse(is.na(bit_sig),NA_character_,
                  ifelse(count_by_id[bit_sig]>=min_class_n,bit_sig,'OTHER'))
if(length(unique(classes[!is.na(classes)]))<2L)
  warning('Only one common multilocus class; trait association tests may be unavailable')
all_trait_data <- data.frame(pop_label=ids,raw_haplotype=bit_sig,haplotype_class=classes,Y,
                             stringsAsFactors=FALSE)
if(!is.null(pcs))all_trait_data <- cbind(all_trait_data,pcs)
write_out(all_trait_data,'57_RIL_haplotype_and_trait_BLUEs.csv')

compare <- function(data,trait,use_pcs=FALSE){
  d<-data[is.finite(data[[trait]])&!is.na(data$haplotype_class),,drop=FALSE]
  d$haplotype_class<-droplevels(factor(d$haplotype_class))
  if(nrow(d)<80||nlevels(d$haplotype_class)<2)
    return(data.frame(n_RIL=nrow(d),n_classes=nlevels(d$haplotype_class),
      F=NA_real_,p=NA_real_,delta_adj_R2=NA_real_))
  cc<-if(use_pcs&&any(grepl('^kinship_PC',names(d))))
    grep('^kinship_PC',names(d),value=TRUE) else character()
  rhs0<-if(length(cc))paste(cc,collapse=' + ') else '1'
  fit0<-lm(as.formula(paste('`',trait,'` ~ ',rhs0,sep='')),data=d)
  fit1<-lm(as.formula(paste('`',trait,'` ~ ',rhs0,' + haplotype_class',sep='')),data=d)
  av<-anova(fit0,fit1)
  data.frame(n_RIL=nrow(d),n_classes=nlevels(d$haplotype_class),
             F=as.numeric(av$F[2]),p=as.numeric(av$`Pr(>F)`[2]),
             delta_adj_R2=summary(fit1)$adj.r.squared-summary(fit0)$adj.r.squared)
}
res <- do.call(rbind,lapply(trait_names,function(tr){
  a<-compare(all_trait_data,tr,FALSE)
  b<-compare(all_trait_data,tr,TRUE)
  data.frame(trait=tr,n_RIL=a$n_RIL,n_haplotype_classes=a$n_classes,
    p_unadjusted=a$p,delta_adj_R2_unadjusted=a$delta_adj_R2,
    p_with_kinship_PCs=if(is.null(pcs))NA_real_ else b$p,
    delta_adj_R2_with_PCs=if(is.null(pcs))NA_real_ else b$delta_adj_R2)
}))
res$BH_q_unadjusted<-p.adjust(res$p_unadjusted,method='BH')
res$BH_q_with_kinship_PCs<-p.adjust(res$p_with_kinship_PCs,method='BH')
write_out(res,'57_Haplotype_association_seven_trait_BLUEs_EXPLORATORY.csv')

## Exploratory *conditional* test of the two nominated linked SNPs from 17.
## Not new confirmatory GWAS: variants were chosen after seeing phenotype data.
if('V_6300718'%in%colnames(Hdosage)){
  a<-as.numeric(Hdosage[,'V_6322080']);b<-as.numeric(Hdosage[,'V_6300718'])
  test_rows<-list()
  for(tr in c('SF','FB','BS')){
    df<-data.frame(y=Y[,tr],A=a,B=b)
    if(!is.null(pcs))df<-cbind(df,pcs)
    df<-df[complete.cases(df),,drop=FALSE]
    pc_names<-grep('^kinship_PC',names(df),value=TRUE)
    base_rhs<-if(length(pc_names))paste(pc_names,collapse=' + ')else '1'
    fullfit<-lm(as.formula(paste('y ~ ',base_rhs,' + A + B')),data=df)
    fA<-lm(as.formula(paste('y ~ ',base_rhs,' + B')),data=df)
    fB<-lm(as.formula(paste('y ~ ',base_rhs,' + A')),data=df)
    pA<-anova(fA,fullfit)$`Pr(>F)`[2]; pB<-anova(fB,fullfit)$`Pr(>F)`[2]
    test_rows[[length(test_rows)+1L]]<-data.frame(trait=tr,n_RIL=nrow(df),
      SNP_A='V_6322080',SNP_B='V_6300718',p_A_given_B=pA,p_B_given_A=pB,
      adjusted_for_kinship_PCs=length(pc_names)>0)
  }
  write_out(do.call(rbind,test_rows),'57_Within848_conditional_leads_EXPLORATORY.csv')
}

pdf(file.path(out,'57_Four_tag_haplotype_frequencies.pdf'),width=9,height=5.5)
par(mar=c(6,5,4,1))
barplot(counts$n_RIL,names.arg=counts$rare_haplotype_4tag,las=2,
        col='grey55',ylab='Number of RILs',xlab='Initially-rare-allele state: 848 / 851 / 853 / 854',
        main='Chromosome V four-tag RIL multilocus patterns')
mtext('RIL multilocus genotype ≠ reconstructed founder haplotype',side=3,line=.2,cex=.7)
dev.off()
writeLines(c('RIL four-tag patterns are unphased only at homozygous sites.',
  'Heterozygous or missing calls are excluded from complete haplotype assignments.',
  'Rare-allele orientation uses exactly one shared A00 Pool-seq baseline.',
  'Common patterns in inbred lines can proxy haplotypes; recombination and phase remain possible.',
  'Trait association tests are exploratory on adjusted BLUEs, not causal or mixed-model fine-mapping.',
  'Ancestry is NOT assigned; external founder genotype comparisons would be required.',
  paste('resolved_RIL',sum(!is.na(bit_sig))),paste('common_haplotype_classes',length(unique(classes[!is.na(classes)])))),
  file.path(out,'57_Interpretation_and_QC.txt'))
cat('Script 57 completed. Output directory: ',out,'\n')
