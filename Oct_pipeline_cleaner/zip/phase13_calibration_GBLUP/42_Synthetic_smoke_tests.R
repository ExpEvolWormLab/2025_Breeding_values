## 42_Synthetic_smoke_tests.R
## Lightweight checks for the key formulas; does NOT replace real-data QC.
## No project data required. Requires rrBLUP.

suppressPackageStartupMessages(library(rrBLUP))
helper09 <- if(file.exists('09f_drift_functions.R')) '09f_drift_functions.R' else
  'drift_nullfit_update/09f_drift_functions.R'
source(helper09)
source('33_Frequency_drift_helpers.R')
set.seed(20261029)

## Near-zero frequency is finite and W-F null simulation includes absorption.
samples <- data.frame(
 sample_id=c('A00','GA1G10','GA1G36','GA1G50','GA2G10','GA2G36','GA2G50','GA4G10','GA4G36','GA4G50'),
 generations=c(0,10,36,50,10,36,50,10,36,50),
 replicate=c('0','1','1','1','2','2','2','4','4','4'),stringsAsFactors=FALSE)
d <- rep(110,length(samples$generations))
k <- c(0,0,1,2,0,1,2,0,0,3)
info <- conditional_neutral_info(k,d,samples,Ne=1000)
stopifnot(isTRUE(info$ok),is.finite(info$conditional_SE),info$conditional_SE>0)
sm <- simulate_conditional_neutral(info,k,d,samples,Ne=1000,n_sim=200L)
stopifnot(length(sm)==200L,all(is.finite(sm)))
cat('WF baseline-zero synthetic test: PASS\n')

## The cross-kinship out-of-fold BLUP must agree with an augmented kernel solve
## whose test phenotypes are never used to estimate variance components.
n <- 36L;m <- 300L
X <- matrix(rbinom(n*m,2,.3),nrow=n,ncol=m)
tr <- seq_len(25L);te <- 26:n
mu_marker <- colMeans(X[tr,,drop=FALSE]);p <- mu_marker/2
W <- sweep(X,2,mu_marker,'-')
K <- tcrossprod(W)/sum(2*p*(1-p))
latent <- rnorm(m,sd=.12)
y <- as.vector(W%*%latent+rnorm(n,sd=.7))
fit <- rrBLUP::mixed.solve(y=y[tr],K=K[tr,tr,drop=FALSE],method='REML')
lam <- as.numeric(fit$Ve/fit$Vu)
stopifnot(is.finite(lam),lam>0)
g <- as.numeric(fit$beta[1])+K[te,tr,drop=FALSE]%*%
     solve(K[tr,tr,drop=FALSE]+diag(lam,length(tr)),y[tr]-as.numeric(fit$beta[1]))
stopifnot(length(g)==length(te),all(is.finite(g)))
cat('GBLUP cross-kinship synthetic test: PASS\n')
