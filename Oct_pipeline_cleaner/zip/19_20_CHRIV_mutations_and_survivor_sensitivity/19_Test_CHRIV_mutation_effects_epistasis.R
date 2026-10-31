## 19_Test_CHRIV_mutation_effects_epistasis.R
## Targeted RIL genotype/phenotype association for two chromosome IV variants.
## Source from the project root after script 00. This script does not change GWAS outputs.
## Trait BLUEs are line-level adjusted phenotypes, so this is a corroborating
## analysis, not an identical refit of the original replicate-level GridLMM model.
## Kinship excludes ALL chromosome IV markers (LOCO).

positions <- c(8892663L, 10499135L)
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
preparation_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
origin_file <- '19_RIL_source_populations.csv' # optional, columns pop_label,source_population
out_dir <- 'mutation_effects_chrIV'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
alpha <- 0.05
power_effect_sizes_SD <- c(.10, .25, .50, 1.0)

suppressPackageStartupMessages(library(AGHmatrix))
p <- new.env(parent = globalenv())
load(preparation_file, envir = p)
need <- c('genotype_matrix','genotype_pruned_matrix','phen_NaCl')
if (!all(need %in% ls(p))) stop('Preparation file lacks: ',
                                 paste(setdiff(need, ls(p)),collapse=', '))
G <- p$genotype_matrix
P <- p$genotype_pruned_matrix
Y <- p$phen_NaCl
Y$pop_label <- as.character(Y$pop_label)
ids <- intersect(rownames(G), intersect(rownames(P), Y$pop_label))
if (anyDuplicated(ids) || anyDuplicated(Y$pop_label)) stop('Duplicated RIL labels')
P <- P[ids,,drop=FALSE]
Y <- Y[match(ids,Y$pop_label),,drop=FALSE]

## Require a unique exact variant identifier. 00 uses IV_POSITION.
find_marker <- function(pos, marker_names) {
  expected <- paste0('IV_',pos)
  if (expected %in% marker_names) return(expected)
  matches <- grep(paste0('^IV[_:.]', pos, '($|[_:.])'),marker_names,value=TRUE)
  if (length(matches)!=1L)
    stop('Expected exactly one genotype column for IV:', pos,
         '; found ', length(matches), '. Confirm original VCF coordinates.')
  matches
}
marker_cols <- vapply(positions,find_marker,character(1),marker_names=colnames(G))
## Extract only two genotype columns: the full G matrix may be several GB.
g_raw <- G[ids,marker_cols,drop=FALSE]
rm(G)
p$genotype_matrix <- NULL
colnames(g_raw) <- c('mut_8892663','mut_10499135')
if (any(is.finite(g_raw) & (g_raw < -1e-6 | g_raw > 2+1e-6)))
  stop('Mutation genotype dosages outside [0,2]; check ALT-allele encoding')
if (any(is.finite(g_raw) & abs(g_raw-round(g_raw))>1e-5))
  warning('Fractional/imputed genotype dosages detected. The association uses ',
          'the continuous dosages; genotype-class counts are rounded solely ',
          'for descriptive QC.')
## ALT-homozygote contrast is exactly 1 unit, 0/2 homozygote dosage.
g <- g_raw / 2

## Optional genotype-cohort annotation. Its absence must never be substituted
## with cohort labels guessed from RIL names or mutation genotype.
origins_available <- file.exists(origin_file)
if (origins_available) {
  org <- read.csv(origin_file, stringsAsFactors=FALSE)
  if (!all(c('pop_label','source_population') %in% names(org)))
    stop('Source-population CSV must have pop_label,source_population')
  if (anyDuplicated(org$pop_label)) stop('Duplicate population-source annotation')
  origin <- org$source_population[match(ids,org$pop_label)]
  origin <- ifelse(is.na(origin) | !nzchar(origin), 'unknown', origin)
  cat('RIL source annotation: ',paste(table(origin),collapse=', '),'\n')
} else {
  origin <- rep('unknown',length(ids))
  warning('No RIL source-population annotation supplied. Mutation effects may ',
          'be confounded by ancestral/evolved RIL origin or genetic background.')
}

## Descriptive genotype distribution, including heterozygotes if present.
qc <- data.frame(pop_label=ids,source_population=origin,
                 dosage_8892663=as.numeric(g_raw[,1]),
                 dosage_10499135=as.numeric(g_raw[,2]),
                 stringsAsFactors=FALSE)
write.csv(qc,file.path(out_dir,'19_RIL_mutation_genotypes.csv'),row.names=FALSE)
hap <- as.data.frame(table(factor(round(g_raw[,1]),levels=0:2),
                           factor(round(g_raw[,2]),levels=0:2)),stringsAsFactors=FALSE)
names(hap) <- c('dosage_8892663','dosage_10499135','n_RILs')
write.csv(hap,file.path(out_dir,'19_Two_locus_genotype_counts.csv'),row.names=FALSE)
valid_g <- complete.cases(g)
r2 <- if(sum(valid_g)>4 && all(apply(g[valid_g,,drop=FALSE],2,sd)>0))
  cor(g[valid_g,1],g[valid_g,2])^2 else NA_real_
cat('RILs with both SNP calls:',sum(valid_g),'\n',
    'Two-locus genotype r2:',r2,'\n')
print(hap,row.names=FALSE)

## Match script 01's LOCO relatedness; NEVER condition on chromosome IV
## markers in the K matrix when assessing IV mutation effects.
keep_pruned <- !grepl('^IV_',colnames(P))
if (sum(keep_pruned)<100L) stop('Too few non-IV SNPs for LOCO kinship')
P_other <- P[,keep_pruned,drop=FALSE]
cat('Building chromosome-IV LOCO kinship from ',ncol(P_other),' pruned SNPs\n')
K <- AGHmatrix::Gmatrix(P_other,method='VanRaden')
K <- as.matrix(K[ids,ids,drop=FALSE])
K <- (K+t(K))/2
rm(P,P_other)
p$genotype_pruned_matrix <- NULL
if(any(!is.finite(K)))stop('LOCO kinship has nonfinite entries')

## An EMMAX-style variance-component score approximation: estimate a single
## variance ratio under the null for each phenotype, then hold it fixed for
## marginal, two-locus additive and interaction GLS tests. All p-values are
## approximate, and observed BLUE estimation error is not modeled separately.
whiten_reml <- function(y,K,X0) {
  ev <- eigen(K,symmetric=TRUE)
  lam <- pmax(0,ev$values)
  U <- ev$vectors
  yr <- as.vector(crossprod(U,y)); Xr <- crossprod(U,X0)
  n <- length(y); q <- ncol(X0)
  profile <- function(log_delta) {
    w <- 1/(lam+exp(log_delta))
    Xt <- crossprod(Xr,Xr*w)
    if (rcond(Xt)<1e-12) return(Inf)
    b <- solve(Xt,crossprod(Xr,w*yr))
    e <- yr-drop(Xr%*%b)
    rss <- sum(w*e^2)
    if (!is.finite(rss) || rss<=0) return(Inf)
    ld <- determinant(Xt,logarithm=TRUE)
    if (ld$sign<=0) return(Inf)
    (n-q)*log(rss/(n-q))+sum(log(lam+exp(log_delta)))+as.numeric(ld$modulus)
  }
  opt <- optimize(profile,interval=c(-9,9))
  delta <- exp(opt$minimum)
  w <- 1/sqrt(lam+delta)
  list(U=U,whiten=function(M) sweep(crossprod(U,M),1,w,'*'),
       ystar=yr*w,delta=delta,profile_value=opt$objective)
}

fit_gls <- function(ystar, Xstar) {
  n <- length(ystar);p <- ncol(Xstar)
  qrX <- qr(Xstar,tol=1e-9)
  if (qrX$rank<p || n<=p+3) return(NULL)
  b <- as.vector(qr.coef(qrX,ystar))
  e <- ystar-drop(Xstar%*%b)
  rss <- sum(e^2); sig <- rss/(n-p)
  xx <- crossprod(Xstar)
  vi <- try(solve(xx),silent=TRUE)
  if(inherits(vi,'try-error')) return(NULL)
  se <- sqrt(pmax(0,diag(vi)*sig))
  list(beta=b,SE=se,vi=vi,rss=rss,sigma2=sig,df=n-p,p=p)
}

nested_p <- function(simple,complex) {
  if(is.null(simple)||is.null(complex))return(NA_real_)
  df1<-complex$p-simple$p
  if(df1<=0)return(NA_real_)
  f<-pmax(0,(simple$rss-complex$rss)/df1)/(complex$rss/complex$df)
  pf(f,df1,complex$df,lower.tail=FALSE)
}

results <- list(); powers <- list(); main_powers <- list()
for (trait in traits) {
  yy <- as.numeric(Y[[trait]])
  good <- is.finite(yy) & valid_g
  n <- sum(good)
  if (n<30L) { warning('Too few RILs for ',trait);next }
  y <- yy[good];g1 <- as.numeric(g[good,1]);g2 <- as.numeric(g[good,2])
  o <- factor(origin[good])
  ## Unknown is a single level -> intercept-only null.
  X0 <- if(length(levels(o))>1L) model.matrix(~o) else matrix(1,n,1,
                                                               dimnames=list(NULL,'(Intercept)'))
  if (qr(X0)$rank<ncol(X0))stop('Source-population design is rank deficient')
  whiten <- whiten_reml(y,K[good,good,drop=FALSE],X0)
  W <- function(X) whiten$whiten(as.matrix(X))
  ys <- whiten$ystar
  f0 <- fit_gls(ys,W(X0))
  X1 <- cbind(X0,mutation_8892663=g1)
  X2 <- cbind(X0,mutation_10499135=g2)
  XA <- cbind(X0,mutation_8892663=g1,mutation_10499135=g2)
  XI <- cbind(XA,interaction=g1*g2)
  f1 <- fit_gls(ys,W(X1));f2 <- fit_gls(ys,W(X2))
  fa <- fit_gls(ys,W(XA));fi <- fit_gls(ys,W(XI))
  read_coef <- function(f,idx) if (!is.null(f)) f$beta[idx] else NA_real_
  read_se <- function(f,idx) if (!is.null(f)) f$SE[idx] else NA_real_
  ## In the joint model, the last 2 columns always correspond to the mutations.
  ## Separate per-mutation effects remain available when collinear variants
  ## prevent estimating both together.
  ma <- ncol(X0)+1L;mb <- ma+1L;mi <- mb+1L
  beta1 <- read_coef(fa,ma); beta2 <- read_coef(fa,mb)
  se1 <- read_se(fa,ma);se2 <- read_se(fa,mb)
  betaI <- read_coef(fi,mi);seI <- read_se(fi,mi)
  p_coef <- function(b,se,df) if(is.finite(b)&&is.finite(se)&&se>0&&is.finite(df))
    2*pt(abs(b/se),df,lower.tail=FALSE) else NA_real_
  results[[trait]] <- data.frame(trait=trait,n_RILs=n,
    n_alt_8892663=sum(g1>0),n_alt_10499135=sum(g2>0),
    n_double_alt=sum(g1>0 & g2>0),
    n_single_8892663=sum(g1>0 & g2==0),
    n_single_10499135=sum(g1==0 & g2>0),
    two_locus_r2=if(sd(g1)>0 && sd(g2)>0)cor(g1,g2)^2 else NA_real_,
    n_origin_groups=length(unique(o)),
    origin_adjusted=length(unique(o))>1L,
    residual_delta=whiten$delta,
    marginal_beta_8892663=read_coef(f1,ma),
    marginal_beta_10499135=read_coef(f2,ma),
    marginal_p_8892663=nested_p(f0,f1),
    marginal_p_10499135=nested_p(f0,f2),
    joint_additive_estimable=!is.null(fa),
    joint_p_both=nested_p(f0,fa),
    conditional_beta_8892663=beta1,
    conditional_SE_8892663=se1,
    conditional_p_8892663=if(!is.null(fa))p_coef(beta1,se1,fa$df)else NA_real_,
    conditional_beta_10499135=beta2,
    conditional_SE_10499135=se2,
    conditional_p_10499135=if(!is.null(fa))p_coef(beta2,se2,fa$df)else NA_real_,
    epistasis_estimable=!is.null(fi),
    interaction_design_condition_number=tryCatch(kappa(W(XI)),error=function(e) NA_real_),
    epistasis_beta=betaI,epistasis_SE=seI,
    epistasis_p=nested_p(fa,fi),
    stringsAsFactors=FALSE)
  ## Conditional 80%-power threshold, accounting for observed two-locus LD,
  ## kinship, phenotype spread and OPTIONAL source-population covariates.
  ## This approximates power for a fixed dosage-scale interaction contrast.
  if(!is.null(fi)) {
    effect_se_SD <- seI/sd(y)
    critical <- qt(1-alpha/2,fi$df)
    thresh <- (critical + qt(.80,fi$df))*effect_se_SD
    power_at <- vapply(power_effect_sizes_SD,function(d) {
      ncp <- d/effect_se_SD
      pt(-critical,fi$df,ncp=ncp) + pt(critical,fi$df,ncp=ncp,lower.tail=FALSE)
    },numeric(1))
  } else {
    effect_se_SD<-thresh<-NA_real_
    power_at<-rep(NA_real_,length(power_effect_sizes_SD))
  }
  powers[[trait]] <- data.frame(trait=trait,
    epistasis_SE_in_trait_SD=effect_se_SD,
    minimum_abs_epistasis_for_80pct_power_SD=thresh,
    effect_sizes_SD=paste(power_effect_sizes_SD,collapse=';'),
    powers=paste(signif(power_at,3),collapse=';'),
    stringsAsFactors=FALSE)
  ## Main-effect power is presented BOTH marginally and conditional on the
  ## second mutation. Conditional power may be unidentifiable in strong LD.
  for (target in c('8892663','10499135')) {
    single <- if (target=='8892663') f1 else f2
    for (model_type in c('marginal','joint_conditional')) {
      ff <- if (model_type=='marginal') single else fa
      idx <- if (model_type=='marginal') ma else
        if (target=='8892663') ma else mb
      zz <- if(!is.null(ff)) ff$SE[idx]/sd(y) else NA_real_
      df <- if(!is.null(ff)) ff$df else NA_real_
      crit <- if(is.finite(df)) qt(1-alpha/2,df) else NA_real_
      min80 <- if(is.finite(zz) && zz>0)
        (crit+qt(.80,df))*zz else NA_real_
      pp <- if(is.finite(zz) && zz>0) vapply(power_effect_sizes_SD,
        function(delta) {
          nc <- delta/zz
          pt(-crit,df,ncp=nc)+pt(crit,df,ncp=nc,lower.tail=FALSE)
        },numeric(1)) else rep(NA_real_,length(power_effect_sizes_SD))
      main_powers[[length(main_powers)+1L]] <- data.frame(trait=trait,
        variant=paste0('IV:',target),model=model_type,
        SE_per_trait_SD=zz,minimum_abs_effect_for_80pct_power_SD=min80,
        effect_sizes_SD=paste(power_effect_sizes_SD,collapse=';'),
        powers=paste(signif(pp,3),collapse=';'))
    }
  }
  cat('Trait',trait,'n',n,' add estimable',!is.null(fa),
      ' interaction estimable',!is.null(fi),'\n')
}
if(!length(results))stop('No phenotype could be analyzed')
out <- do.call(rbind,results)
power <- do.call(rbind,powers)
## Adjust tests across traits; do not present nominal p-values as GWAS hits.
for (nm in c('marginal_p_8892663','marginal_p_10499135',
             'conditional_p_8892663','conditional_p_10499135','epistasis_p')) {
  out[[paste0(nm,'_BH')]] <- p.adjust(out[[nm]],method='BH')
}
write.csv(out,file.path(out_dir,'19_Mutation_effects_and_epistasis.csv'),row.names=FALSE)
write.csv(power,file.path(out_dir,'19_Epistasis_power.csv'),row.names=FALSE)
write.csv(do.call(rbind,main_powers),
          file.path(out_dir,'19_Main_effect_power.csv'),row.names=FALSE)

## Plot SB two-locus genotype classes directly; these are adjusted BLUEs and
## NOT fitted interaction estimates or causal gene perturbation effects.
sb <- as.numeric(Y$SB)
ix <- is.finite(sb) & valid_g
d <- data.frame(SB=sb[ix],
  haplotype=paste0('(',g_raw[ix,1],',',g_raw[ix,2],')'))
pdf(file.path(out_dir,'19_SB_by_CHRIV_genotypes.pdf'),width=7,height=4.6)
boxplot(SB~haplotype,data=d,xlab='ALT dosages (IV:8892663, IV:10499135)',
        ylab='Adjusted RIL SB BLUE',main='Chromosome IV candidate genotypes',
        outline=FALSE,col='grey90')
stripchart(SB~haplotype,data=d,vertical=TRUE,method='jitter',
           pch=16,cex=.35,col=adjustcolor('steelblue4',alpha.f=.3),add=TRUE)
dev.off()
cat('\nSaved outputs in: ',out_dir,'\n',sep='')
cat('Interpretation: kinship-adjusted, BLUE-level associations; RIL origin ',
    'adjustment requires documented source-population IDs.\n')
cat('If the interaction design is rank deficient, epistasis is unidentifiable ',
    'with these two segregating loci, regardless of nominal sample size.\n')
