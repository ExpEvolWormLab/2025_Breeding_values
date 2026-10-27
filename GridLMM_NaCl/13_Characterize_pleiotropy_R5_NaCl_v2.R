rm(list=ls());gc()
library(data.table)
library(lme4)

threshold_clust = 0.5
load(paste0("RData_GridLMM_NaCl/GWAS_cluster_results_BLUE_concatenated_R",round(10*threshold_clust),"_with_residuals.RData"))
load(file="GridLMM_results_NaCl/All_Gwas.RData")

vect_P_traits=paste0(c("SF",'SB','FS','FB','BS','BF',"area.F","beta"),"_NaCl")

all_GWAS_clusters$sign_overall <- rowSums(all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits)])
dim(all_GWAS_clusters)
table(all_GWAS_clusters$is_sign_area.F_NaCl)
table(all_GWAS_clusters$sign_overall)

## We will merge the group ID from the clusters using the mean POS
all_cluster_merged = read.table(paste0("txt_GridLMM_NaCl/all_sign_GWAS_with_clusters_R",round(10*threshold_clust),".txt"),sep="\t",header=TRUE)
  #read.table(file=paste0('txt_BLUEs//all_sign_GWAS_with_clusters.txt'),sep='\t',h=TRUE)
all_cluster_merged$CHROM=factor(all_cluster_merged$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
all_cluster_merged$POS=as.numeric(all_cluster_merged$POS)
all_cluster_merged = all_cluster_merged[order(all_cluster_merged$CHROM,all_cluster_merged$POS),]

df_clusters=read.table(paste0("data_BLUE_GridLMM_NaCl/data_clusters_BLUE_R",round(10*threshold_clust),".txt"),h=TRUE)
all_GWAS_clusters=merge(all_GWAS_clusters,df_clusters)
dim(all_GWAS_clusters) # 1058 clusters for R2
### Add the tip/arm/center tag
rkmap <- read.table("Caendr_filtering/rockman_kruglya_map_liftOver.txt",h=TRUE)
all_GWAS$CHROM=tstrsplit(all_GWAS$rs,"_")[[1]]
table(all_GWAS$CHROM)
all_GWAS$chr_type=NA
for(i in 1:nrow(rkmap)){
  if(sum(all_GWAS$CHROM==rkmap$CHROM[i] & all_GWAS$ps>=rkmap$Start[i] & all_GWAS$ps<=rkmap$End[i])>0){
  all_GWAS$chr_type[all_GWAS$CHROM==rkmap$CHROM[i] & all_GWAS$ps>=rkmap$Start[i] & all_GWAS$ps<=rkmap$End[i]]=rkmap$type[i]
  }
}
table(all_GWAS$chr_type)
table(is.na(all_GWAS$chr_type))
summary(subset(all_GWAS,is.na(chr_type) & CHROM=='V')[,c("CHROM","ps")]$ps) # Some tips on CHR V
all_GWAS$chr_type[all_GWAS$CHROM=="V" & all_GWAS$ps>=20341612]="right_tip"
table(subset(all_GWAS,is.na(chr_type))$CHROM) # Only MtDNA left

## Now also add the group ID to the GWAS
dim(all_GWAS)
dim(all_cluster_merged)
all_GWAS_sign <- (merge(all_GWAS,all_cluster_merged[,c("rs","group","sign_overall","POS")]))
dim(all_GWAS_sign)
names(all_cluster_merged)[!names(all_cluster_merged)%in%names(all_GWAS_sign)]

### How much is in Caendr  - Founders
#library(vcfR)
#founders_Caendr <- read.vcfR("Caendr_filtering/WI.20250625.hard-filter.isotype.CemeeFounders_poly.vcf.gz")
#caendr_pos <- as.data.frame(founders_Caendr@fix[,1:2])

#sign_for_Caendr <- all_GWAS_sign
#names(sign_for_Caendr)[2:5]=c("CHROM","POS","ALT","REF")
#dim(merge(sign_for_Caendr,data.frame(founders_Caendr@fix)))
#dim(merge(sign_for_Caendr,data.frame(founders_Caendr@fix[,1:2])))

gp_nbSNPs=tapply(all_GWAS_sign$group,all_GWAS_sign$group,length)
gp_min=tapply(all_GWAS_sign$ps,all_GWAS_sign$group,min)
gp_max=tapply(all_GWAS_sign$ps,all_GWAS_sign$group,max)
gp_type=tapply(all_GWAS_sign$chr_type,all_GWAS_sign$group,unique)
table(as.numeric(lapply(gp_type,length)))

#gp_type[as.numeric(lapply(gp_type,length))==2] # We decide center>arm>tip
for(i in c(1:length(gp_type))[as.numeric(lapply(gp_type,length))==2]){
  if(gp_type[[i]][1]=="left_tip"){
    gp_type[[i]]='left_arm'
  }else if(gp_type[[i]][1]=="left_arm"){
    gp_type[[i]]='center'
  }else{
    gp_type[[i]]=gp_type[[i]][[2]]
  }
}
table(as.numeric(lapply(gp_type,length))) # OK
gp_type=unlist(gp_type)

gp_length=data.frame(group=names(gp_min),clust_length=as.numeric(gp_max-gp_min),nb_SNPs=as.numeric(gp_nbSNPs),chr_type=gp_type)
head(gp_length)

pdf("Cluster_description_R5.pdf",h=5,w=5)
plot(gp_length$nb_SNPs,gp_length$clust_length,pch=16,cex=.8,col=c("red","cornflowerblue","lightgreen","lightgreen","cornflowerblue")[as.numeric(as.factor(gp_length$chr_type))],
     log="xy",xlab="Number of SNPs in QTL",ylab="QTL length (bp)")
legend(50,100,c("center","arm","tip"),pch=16,col=c("red","cornflowerblue","lightgreen"))
dev.off()
pdf("Cluster_nbSNPs_R5.pdf",h=3,w=4)
hist(log10(gp_length$nb_SNPs),xlab="Nb SNPs (log10)",main='')
dev.off()
pdf("Cluster_length_R5.pdf",h=3,w=4)
hist(log10(gp_length$clust_length),xlab="Cluster length (log10 bp)",main='')
dev.off()

## ClusterID & group are the same thing
all_GWAS_clusters=merge(all_GWAS_clusters,gp_length)
head(all_GWAS_clusters)

chr_type2=all_GWAS_clusters$chr_type
all_GWAS_clusters$chr_type2=c("center","arm","tip","tip","arm")[as.numeric(as.factor(all_GWAS_clusters$chr_type))]

hist(log(all_GWAS_clusters$nb_SNPs))
model1 <- glm(sign_overall~log(nb_SNPs)+clust_length+chr_type2,
              family="quasipoisson",data=subset(all_GWAS_clusters,sign_overall>0))
summary(model1) # Yes, the more SNPs, the more pleiotropic

boxplot(all_GWAS_clusters$clust_length+1~all_GWAS_clusters$sign_overall,log="y")

### Add the pleiotropy
#dim(merge(all_GWAS_sign,unique(all_GWAS_clusters[,c("group","sign_overall","sign_NaCl","sign_NGM")])))

#pheno_NGM=read.table("Breeding_values_NGM.txt",sep='\t',h=TRUE)
#pheno_NaCl=read.table("Breeding_values_NaCl.txt",sep='\t',h=TRUE)
#gemma_pheno=read.table('data_BLUE//gemma_pheno.txt',sep='\t',h=FALSE)

# Clusters
sign_all_df_clust <- all_GWAS_clusters[,names(all_GWAS_clusters)%in%paste0("is_sign_",vect_P_traits)]
table(rowSums(sign_all_df_clust)>=1) # ~737 Clusters / 10 are no longer significant
all_GWAS_clusters$sign_overall <- rowSums(sign_all_df_clust)>=1

for(i in 1:8){
  all_GWAS[,paste0("is_sign_",vect_P_traits[i])] = p.adjust(all_GWAS[,paste0("pvalue_T",i)],method="BH")<0.2
  print(table(all_GWAS[,paste0("is_sign_",vect_P_traits[i])]))
}

# SNPs
sign_all_df <- all_GWAS[,names(all_GWAS)%in%paste0("is_sign_",vect_P_traits)]
table(rowSums(sign_all_df)>=1) # 33k SNPs

all_GWAS$sign_overall <- rowSums(sign_all_df)>=1

rowS_NaCl_clust <-table(rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits)]))
rowS_NaCl <-table(rowSums(all_GWAS[,paste0("is_sign_",vect_P_traits)]))

vect_col=c(adjustcolor("grey",alpha=0.3),adjustcolor("firebrick3",alpha=1))
i=7
plot(-log10(all_GWAS_clusters[,paste0("p_value_REML.cor_",vect_P_traits[i])])~all_GWAS_clusters[,paste0("beta_",vect_P_traits[i])],
     pch=16,xlab='Beta coefficient',ylab="-log(p)",
     col=vect_col[1+as.numeric(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits[i])])],
     main=paste0('Trait: ',vect_P_traits[i]))

#plot((all_GWAS_clusters$beta_SF_NaCl),(all_GWAS_clusters$beta_area.F_NaCl),pch=16,xlim=c(-.1,.1),ylim=c(-.1,.1))


### Heatmap plots between beta/sign ?
heatmap(cor(na.omit(all_GWAS_clusters[,paste0("beta_",vect_P_traits)])),
        labRow=vect_P_traits,
        labCol=vect_P_traits)

## With the non-significant
plot(rowS_NaCl/sum(rowS_NaCl),type='b',col='cornflowerblue',ylim=c(.00000001,1),log="y",xlim=c(0,8))
points(rowS_NaCl_clust/sum(rowS_NaCl_clust),type='b',col='darkorange',lty=2)


pdf("Lille_Pleio_R5_clusters_NaCl.pdf",w=4,h=4)
plot(rowS_NaCl,type='b',col=adjustcolor('darkorange',alpha=.4),ylim=c(1,1000000),log="y",xlim=c(0,8),las=1,bty="n",xlab="",xaxt="n",ylab="")
axis(1,0:8)
points(rowS_NaCl_clust,type='b',col='darkorange',lty=2)
legend(4,1000000,lty=c(1,2),lwd=3,col=c("darkorange","darkorange"),c('SNPs',"LD blocks"))
dev.off()

library(gplots)
library(limma) 
vennDiagram(vennCounts(
  cbind(as.numeric(all_GWAS_clusters$is_sign_SF_NaCl),as.numeric(all_GWAS_clusters$is_sign_area.F_NaCl))),
  names=c("SF NaCl","Area NaCl"),main='All')

load('RData/Gmatrices_P.RData')
heatmap(G_mat_all)

####### Reproduce the Wagner/Frachon pipeline 

# RIL pleiotropic scaling: Wagner et al. (2008) covariance correction
# and Frachon et al. (2017) Manhattan / Euclidean exponents.

all_GWAS_clusters_sign_NaCl <- na.omit(all_GWAS_clusters[rowSums(all_GWAS_clusters[,paste0("is_sign_",vect_P_traits)])>0,])
nb_sign_clust_NaCl <- nrow(all_GWAS_clusters_sign_NaCl)
all_GWAS_clusters_sign=all_GWAS_clusters_sign_NaCl

vect_P_traits=c("SF_NaCl","SB_NaCl","FS_NaCl","FB_NaCl","BS_NaCl","BF_NaCl","area.F_NaCl","beta_NaCl")
traits <- as.character(vect_P_traits)
qtl <- all_GWAS_clusters_sign
beta_cols <- paste0("beta_", traits)
sig_cols  <- paste0("is_sign_", traits)

glmm.env=new.env()
load('RData/Preparatory_GLMM.RData',envir=glmm.env)
data_NaCl=glmm.env$all_BLUEs_genotyped_NaCl
names(data_NaCl)[1:7]=paste0(names(data_NaCl)[1:7],"_NaCl")
names(data_NaCl)[15]=paste0(names(data_NaCl)[15],"_NaCl")
rm(glmm.env)

is_arm <- all_GWAS_clusters_sign$chr_type2=="arm"
is_center <- all_GWAS_clusters_sign$chr_type2=="center"


stopifnot(length(traits) >= 2L,
          !anyDuplicated(traits),
          all(c(traits, "pop_label", "date_str", "session_id",
                "logD", "rel_humidity", "temperature") %in% names(data_NaCl)),
          all(c(beta_cols, sig_cols) %in% names(qtl)),
          length(is_arm) == nrow(qtl),
          length(is_center) == nrow(qtl))

# 1. Get adjusted RIL means AND their covariance of estimation error.
# Treat RIL as FIXED while estimating adjusted means. This is analogous
# to the user's original model and makes replication-error correction
# possible. This estimates among-line covariance, not a full multivariate
# genetic random-effect covariance model.

#For each trait, the error contribution is calculated from the covariance matrix of the adjusted RIL means.
#This accounts for correlations in estimation errors between RIL means for the same trait.
#It does not fully account for correlated measurement errors between different traits.
#A multivariate mixed model would be needed for that more complete correction.

fit_one_trait <- function(trait) {
  
  required <- c(
    trait, "pop_label", "logD",
    "rel_humidity", "temperature",
    "date_str", "session_id"
  )
  
  dat <- data_NaCl[, required, drop = FALSE]
  dat <- droplevels(
    dat[complete.cases(dat), , drop = FALSE]
  )
  
  dat$pop_label <- factor(dat$pop_label)
  dat$date_str <- factor(dat$date_str)
  dat$session_id <- factor(dat$session_id)
  dat$Y <- dat[[trait]]
  
  # 1. Fit the same mixed model
  mod <- lme4::lmer(
    Y ~ logD + rel_humidity + temperature +
      pop_label +
      (1 | date_str) + (1 | session_id),
    data = dat,
    REML = TRUE
  )
  
  # 2. Create prediction data for all RILs
  newdata <- data.frame(
    pop_label = factor(
      levels(dat$pop_label),
      levels = levels(dat$pop_label)
    ),
    logD = 0,
    rel_humidity = 0,
    temperature = 0
  )
  
  # 3. Construct fixed-effects design matrix
  fixed_formula <- lme4::nobars(formula(mod))
  
  X <- model.matrix(
    delete.response(terms(fixed_formula)),
    data = newdata,
    contrasts.arg = attr(
      lme4::getME(mod, "X"), "contrasts"
    )
  )
  
  beta <- lme4::fixef(mod)
  Vbeta <- as.matrix(vcov(mod))
  
  # Match matrix columns to fixed effects
  X <- X[, names(beta), drop = FALSE]
  Vbeta <- Vbeta[names(beta), names(beta),
                 drop = FALSE]
  
  # 4. Adjusted means (equivalent to fixed-only
  # predict() at the specified covariates)
  mu <- drop(X %*% beta)
  
  names(mu) <- as.character(newdata$pop_label)
  rownames(X) <- names(mu)
  
  # Retain X and Vbeta to calculate uncertainty
  # later, after determining common RILs
  list(
    model = mod,
    means = mu,
    X = X,
    Vbeta = Vbeta,
    n_replicates = table(dat$pop_label)
  )
}


fits <- setNames(lapply(traits, fit_one_trait), traits)

# 2. Keep only RILs with estimable means for EVERY trait.
lines <- Reduce(intersect, lapply(fits, function(z) names(z$means)))
lines <- sort(lines)
L <- length(lines)
p <- length(traits)
if (L <= p)
  stop("Need more common RILs than traits to estimate a stable covariance: ",
       L, " RILs, ", p, " traits")

M <- vapply(fits, function(z) unname(z$means[lines]), numeric(L))
rownames(M) <- lines
colnames(M) <- traits
G_observed <- stats::cov(M)

# 3. Wagner-style diagonal correction using the sampling-error variance
# of adjusted RIL means. For each trait, V is the full covariance across
# estimated means. Centering matrix H = I - J/L yields
#   E[ sample variance of errors ] = tr(H V) / (L - 1).
# This is better than blindly averaging squared SEs when mean estimates
# share uncertainty from common environmental regression coefficients.

err_mean <- vapply(fits, function(z) {
  
  # Only use RILs common to all traits
  X <- z$X[lines, , drop = FALSE]
  
  # Center the prediction design matrix
  Xc <- sweep(X, 2, colMeans(X), "-")
  
  # Equivalent to trace(H V) / (L - 1),
  # where V = X Vbeta X'
  sum((Xc %*% z$Vbeta) * Xc) /
    (nrow(Xc) - 1)
  
}, numeric(1))



# Function to apply a fraction of the correction
correct_G <- function(alpha) {
  G <- G_observed
  diag(G) <- diag(G_observed) - alpha * err_mean
  G
}

# Small numerical tolerance
tol <- 1e-10 * max(diag(G_observed))

# Check whether G is positive definite
valid_G <- function(alpha) {
  G <- correct_G(alpha)
  
  if (any(!is.finite(G)) || any(diag(G) <= tol))
    return(FALSE)
  
  min(eigen(G, symmetric = TRUE,
            only.values = TRUE)$values) > tol
}

if (!valid_G(0))
  stop("Even the uncorrected G is not positive definite")

if (valid_G(1)) {
  alpha <- 1
} else {
  # Find the largest admissible correction
  lower <- 0
  upper <- 1
  
  for (i in 1:60) {
    mid <- (lower + upper)/2
    
    if (valid_G(mid)) {
      lower <- mid
    } else {
      upper <- mid
    }
  }
  
  # Slight margin away from singularity
  alpha <- lower * (1 - 1e-6)
}

G_corrected <- correct_G(alpha)

cat("Fraction of correction retained:", alpha, "\n")
cat("Eigenvalues of corrected G:\n")
print(eigen(G_corrected, symmetric = TRUE)$values)


R_corrected <- stats::cov2cor(G_corrected)
R_observed <- stats::cov2cor(G_observed)

cat("Common RILs:", L, "\n")
print(data.frame(trait = traits,
                 observed_var = diag(G_observed),
                 error_var = unname(err_mean),
                 corrected_var = diag(G_corrected)), row.names = FALSE)
print(round(R_corrected, 3))

# 4. Effective number of affected traits (Wagner et al.).
# Use POPULATION eigenvalue variance: mean((lambda - 1)^2), not var(lambda).
neff_from_corr <- function(affected, R) {
  n <- length(affected)
  if (n == 0L) return(NA_real_)
  if (n == 1L) return(1)
  eig <- eigen(R[affected, affected, drop = FALSE],
               symmetric = TRUE, only.values = TRUE)$values
  ne <- n - mean((eig - 1)^2)
  if (!is.finite(ne) || ne < 1 - 1e-6 || ne > n + 1e-6)
    stop("Neff outside [1,N]: check covariance/correlation matrix")
  max(1, min(n, ne))
}

# Unit checks: fully independent -> N; perfectly correlated -> 1
I4 <- diag(4)
dimnames(I4) <- list(letters[1:4], letters[1:4])
J4 <- matrix(1, 4, 4, dimnames = dimnames(I4))
stopifnot(isTRUE(all.equal(neff_from_corr(letters[1:4], I4), 4)),
          isTRUE(all.equal(neff_from_corr(letters[1:4], J4), 1)))

# 5. Effect distances. Do NOT standardize already SD-standardized GWAS betas.
B <- as.matrix(qtl[, beta_cols, drop = FALSE])
if (!is.numeric(B)) stop("QTL beta columns must be numeric")

S <- vapply(sig_cols, function(nm) {
  s <- qtl[[nm]]
  if (is.logical(s)) return(replace(s, is.na(s), FALSE))
  if (is.numeric(s) && all(na.omit(s) %in% c(0, 1)))
    return(replace(as.logical(s), is.na(s), FALSE))
  stop("Significance columns must contain TRUE/FALSE or 1/0: ", nm)
}, logical(nrow(qtl)))
colnames(S) <- traits

B[!S] <- 0
if (any(!is.finite(B))) stop("Non-finite GWAS beta for significant QTL trait")
T_Manhattan <- rowSums(abs(B))
T_Euclidean <- sqrt(rowSums(B^2))
Ntraits <- rowSums(S)

Neff <- vapply(seq_len(nrow(qtl)), function(i) {
  neff_from_corr(traits[S[i, ]], R_corrected)
}, numeric(1))
Neff_uncorrected <- vapply(seq_len(nrow(qtl)), function(i) {
  neff_from_corr(traits[S[i, ]], R_observed)
}, numeric(1))

stopifnot(all(T_Manhattan >= T_Euclidean - 1e-10),
          all(Neff[is.finite(Neff)] <= Ntraits[is.finite(Neff)] + 1e-6))

temp.NaCl=rep(TRUE,length(arm))

# *** Required: all selectors must refer to rows of qtl in SAME ORDER. ***
arm <- is_arm %in% TRUE
center <- is_center %in% TRUE
if (any(arm & center)) stop("Overlapping arm and center QTL assignments")
keep <- (temp.NaCl %in% TRUE) & (arm | center) &
  is.finite(Neff) & (Neff > 0) &
  is.finite(T_Manhattan) & (T_Manhattan > 0) &
  is.finite(T_Euclidean) & (T_Euclidean > 0)

d <- data.frame(
  qtl_row = which(keep),
  group = factor(ifelse(arm[keep], "Arm", "Center"),
                 levels = c("Arm", "Center")),
  Ntraits = Ntraits[keep],
  Neff = Neff[keep],
  Neff_uncorrected = Neff_uncorrected[keep],
  Manhattan = T_Manhattan[keep],
  Euclidean = T_Euclidean[keep]
)

print(table(d$group))
print(tapply(d$Neff, d$group, function(x) length(unique(x))))
if (nrow(d) < 10 || any(table(d$group) < 4) ||
    any(tapply(d$Neff, d$group, function(x) length(unique(x))) < 2))
  stop("Too few QTLs or distinct Neff values for arm/center regressions")

# 6. Original SD-scale GWAS effects + log-linear power-law fitting.
# A_g*Neff^c_g, with multiplicative error.
# Manhattan no-synergy c0 = 1; Euclidean no-synergy c0 = 0.5.
fit_scaling <- function(df, metric, c0) {
  dat <- df
  dat$lx <- log(dat$Neff)
  dat$z <- log(dat[[metric]]) - c0 * dat$lx
  m0 <- lm(z ~ group, data = dat)        # c_A = c_C = c0
  mc <- lm(z ~ group + lx, data = dat)  # common exponent
  mf <- lm(z ~ group * lx, data = dat) # group-specific exponents
  if (anyNA(coef(mf)) || df.residual(mf) < 1)
    stop("Scaling model is rank-deficient for ", metric)
  
  b <- coef(mf)
  V <- vcov(mf)
  tt <- function(v) {
    excess <- sum(v * b)
    se <- sqrt(as.numeric(t(v) %*% V %*% v))
    tstat <- excess / se
    c(exponent = c0 + excess, SE = se,
      lower95_one_sided = c0 + excess - qt(.95, df.residual(mf)) * se,
      p_one_sided = pt(tstat, df.residual(mf), lower.tail = FALSE))
  }
  La <- setNames(rep(0, length(b)), names(b))
  La["lx"] <- 1
  Lc <- La
  Lc["groupCenter:lx"] <- 1
  groups <- rbind(Arm = tt(La), Center = tt(Lc))
  
  # Difference: center c - arm c = interaction coefficient
  int_p <- anova(mc, mf)[["Pr(>F)"]][2]
  int_delta <- unname(b["groupCenter:lx"])
  int_CI <- confint(mf)["groupCenter:lx", ]
  
  # Common-exponent test (conditional on shared-exponent assumption)
  sc <- coef(summary(mc))["lx", ]
  common <- c(exponent = c0 + sc["Estimate"],
              p_one_sided = pt(sc["t value"], df.residual(mc),
                               lower.tail = FALSE))
  list(data = dat, null = m0, common_model = mc, full = mf,
       groups = groups, interaction_p = int_p,
       delta_center_minus_arm = int_delta, delta_CI = int_CI,
       common = common, c0 = c0, metric = metric)
}

resM <- fit_scaling(d, "Manhattan", 1)
resE <- fit_scaling(d, "Euclidean", 0.5)

slope_tests <- rbind(
  data.frame(metric = "Manhattan", group = rownames(resM$groups),
             resM$groups, row.names = NULL),
  data.frame(metric = "Euclidean", group = rownames(resE$groups),
             resE$groups, row.names = NULL)
)
slope_tests$p_Holm <- p.adjust(slope_tests$p_one_sided, method = "holm")

region_tests <- data.frame(
  metric = c("Manhattan", "Euclidean"),
  center_minus_arm = c(resM$delta_center_minus_arm,
                       resE$delta_center_minus_arm),
  CI_low = c(resM$delta_CI[1], resE$delta_CI[1]),
  CI_high = c(resM$delta_CI[2], resE$delta_CI[2]),
  p = c(resM$interaction_p, resE$interaction_p)
)
region_tests$p_Holm <- p.adjust(region_tests$p, method = "holm")

cat("\nGroup-specific one-sided tests (c > theoretical null):\n")
print(slope_tests, row.names = FALSE, digits = 4)
cat("\nRegion differences (center - arm):\n")
print(region_tests, row.names = FALSE, digits = 4)
cat("\nCommon exponent (use only if biologically/statistically justified):\n")
print(rbind(Manhattan = resM$common, Euclidean = resE$common))

# 7. Figure: fitted powers, pointwise 95% CIs, constrained null curves.
plot_scaling <- function(res) {
  df <- res$data
  cols <- c(Arm = "forestgreen", Center = "deepskyblue3")
  plot(df$Neff, df[[res$metric]], type = "n", las = 1, bty = "n",
       xlab = "Effective number of traits (Neff)",
       ylab = paste(res$metric, "total QTL effect"),
       main = paste0(res$metric, " (null exponent = ", res$c0, ")"))
  
  for (g in levels(df$group)) {
    sub <- df[df$group == g, ]
    grid <- seq(min(sub$Neff), max(sub$Neff), length.out = 200)
    nd <- data.frame(lx = log(grid),
                     group = factor(rep(g, length(grid)),
                                    levels = levels(df$group)))
    pp <- predict(res$full, newdata = nd, interval = "confidence")
    # These are CI bands for the conditional median, NOT observations.
    mid <- exp(pp[, "fit"] + res$c0 * log(grid))
    low <- exp(pp[, "lwr"] + res$c0 * log(grid))
    high <- exp(pp[, "upr"] + res$c0 * log(grid))
    polygon(c(grid, rev(grid)), c(low, rev(high)), border = NA,
            col = adjustcolor(cols[g], alpha.f = 0.15))
    lines(grid, mid, col = cols[g], lwd = 3)
    null <- exp(predict(res$null, newdata = nd) + res$c0 * log(grid))
    lines(grid, null, col = cols[g], lwd = 2, lty = 2)
  }
  points(df$Neff, df[[res$metric]], pch = 21, col = "grey30",
         bg = adjustcolor(cols[as.character(df$group)], alpha.f = .65))
  legend("topleft", c("Arms fitted", "Centers fitted",
                      "Arms null", "Centers null"),
         col = rep(cols, 2), lty = c(1, 1, 2, 2), lwd = c(3, 3, 2, 2),
         bty = "n", cex = 0.8)
}

pdf("NaCl_pleiotropy_Wagner_Frachon.pdf", width = 11, height = 5)
oldpar <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
plot_scaling(resM)
plot_scaling(resE)
par(oldpar)
dev.off()

write.csv(d, "NaCl_pleiotropy_QTL_data.csv", row.names = FALSE)
write.csv(slope_tests, "NaCl_pleiotropy_exponent_tests.csv", row.names = FALSE)
write.csv(region_tests, "NaCl_pleiotropy_region_tests.csv", row.names = FALSE)
write.csv(G_corrected, "NaCl_G_corrected.csv")

# 8. Residual diagnostics (read interactively, or save to another PDF).
# par(mfrow = c(2, 2)); plot(resM$full); par(mfrow = c(2, 2)); plot(resE$full)
# Alternative sensitivity check: repeat the Neff calculation with
# R_observed / Neff_uncorrected and refit the models on that predictor.



#### WARNING : 5. Two limitations to retain in your final interpretation
#First, Wagner et al.'s approach provides a correction for trait correlations and measurement error, but the effective pleiotropy estimate and the GWAS effect estimates still contain statistical uncertainty. The p-values in this script are conditional on those estimates.
#Second, because your QTLs were selected using significant GWAS associations, ascertainment bias remains possible: larger-effect loci may be detected for more traits. A simulation reproducing your GWAS detection process would be the strongest additional check before claiming a biological mechanism of synergistic pleiotropy.
#The script uses the log-linear representation of the power law, so its main tests assume approximately normal residuals on the log-effect scale. This is not identical to Wagner et al.'s original-scale regression comparisons.




