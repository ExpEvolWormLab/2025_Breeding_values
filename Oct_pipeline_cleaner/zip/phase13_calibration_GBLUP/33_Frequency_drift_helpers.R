## Shared helpers for frequency-stratified, baseline-conditional neutral testing.
## This is a posterior predictive neutral diagnostic, not a selection coefficient.
## A00 is the only pre-registered baseline. A6G140 is excluded by design.

if(!exists('fit_drift_site_null', mode='function'))
  stop('First source the original 09f_drift_functions.R')

frequency_class <- function(p) {
  maf <- pmin(p, 1-p)
  cut(maf, breaks=c(-Inf, .01, .05, .10, .25, Inf), right=FALSE,
      labels=c('MAF<0.01','0.01-0.05','0.05-0.10','0.10-0.25','0.25-0.50'))
}

## No unknown individual-pool size parameter is supported by this conjugate
## baseline posterior. The direct read model is the pre-specified baseline.
conditional_neutral_info <- function(k, d, pop, Ne=1000L) {
  if(length(k)!=nrow(pop)||length(d)!=nrow(pop))stop('Pool length mismatch')
  fit <- fit_drift_site_null(k,d,pop,Ne=Ne,pool_individuals=NA_real_)
  blank <- list(ok=FALSE,p0=NA_real_,p0_lo=NA_real_,p0_hi=NA_real_,
                a0=NA_real_,b0=NA_real_,mean_slope=NA_real_,
                slope=NA_real_,conditional_SE=NA_real_,conditional_z=NA_real_,
                conditional_p_normal=NA_real_,n=fit$n_pools,n_G0=fit$n_G0,
                n_evolved_reps=fit$n_evolved_reps,
                baseline_ALT_reads=NA_real_,baseline_depth=NA_real_,
                weights=fit$weights)
  if(!is.finite(fit$slope)||!all(is.finite(fit$weights)))return(blank)
  ok <- is.finite(k)&is.finite(d)&d>0&k>=0&k<=d
  g0 <- which(ok & pop$generations==0)
  if(length(g0)!=1L)return(blank)
  a <- k[g0]+.5; b <- d[g0]-k[g0]+.5
  mean_p <- a/(a+b)
  var_p <- a*b/((a+b)^2*(a+b+1))
  expected_hetero <- mean_p*(1-mean_p)-var_p
  idx <- which(ok & pop$generations>0)
  if(length(idx)<4)return(blank)
  w <- fit$weights
  t <- as.numeric(pop$generations[idx]); r <- as.character(pop$replicate[idx]); ww <- w[idx]
  f <- 1-(1-1/(2*Ne))^outer(t,t,pmin)
  C <- var_p + expected_hetero*f*outer(r,r,'==')
  ## Conditional expected read variance after drift; read sampling is direct.
  diag(C) <- diag(C) + expected_hetero*(1-diag(f))/d[idx]
  predicted <- sum(w[g0]*k[g0]/d[g0]) + sum(ww)*mean_p
  variance <- drop(crossprod(ww, C %*% ww))
  if(!is.finite(variance)||variance<=0)return(blank)
  se <- sqrt(variance)
  z <- (fit$slope-predicted)/se
  list(ok=TRUE,p0=mean_p,p0_lo=qbeta(.025,a,b),p0_hi=qbeta(.975,a,b),
       a0=a,b0=b,mean_slope=predicted,slope=fit$slope,
       conditional_SE=se,conditional_z=z,
       conditional_p_normal=2*pnorm(abs(z),lower.tail=FALSE),
       n=fit$n_pools,n_G0=fit$n_G0,n_evolved_reps=fit$n_evolved_reps,
       baseline_ALT_reads=k[g0],baseline_depth=d[g0],weights=w)
}

## Condition on the observed A00 read count. Integrate unknown true p0 over
## Beta(k0+1/2, d0-k0+1/2), then simulate diploid WF drift and Pool-seq reads.
## Independent populations share p0 and drift independently afterward.
## Absorption (loss/fixation), all sequence generations, and read counts are
## preserved. No simulations use observed evolved frequencies as p0.
simulate_conditional_neutral <- function(info, k, d, pop, Ne=1000L,
                                         n_sim=4000L) {
  if(!isTRUE(info$ok))return(numeric(0))
  n_sim <- as.integer(n_sim)
  if(n_sim<50L)stop('n_sim must be >=50')
  evolved <- unique(as.character(pop$replicate[pop$generations>0]))
  tt <- as.integer(pop$generations)
  B <- n_sim
  p0 <- rbeta(B, info$a0, info$b0)
  P <- matrix(rep(p0,length(evolved)),nrow=B,ncol=length(evolved))
  w <- info$weights
  y0 <- sum(w[tt==0]*k[tt==0]/d[tt==0])
  ## The original single G0 is fixed, so only evolved read outcomes contribute.
  slope_sim <- rep(y0,B)
  valid <- is.finite(d)&d>0&is.finite(k)&k>=0&k<=d
  evolve_at <- split(which(valid&tt>0),tt[valid&tt>0])
  for(gen in seq_len(max(tt))){
    P[] <- rbinom(length(P),size=2L*as.integer(Ne),prob=as.vector(P))/(2*Ne)
    ind <- evolve_at[[as.character(gen)]]
    if(length(ind))for(j in ind){
      rr <- match(as.character(pop$replicate[j]),evolved)
      reads <- rbinom(B,size=as.integer(round(d[j])),prob=P[,rr])
      slope_sim <- slope_sim + w[j]*reads/d[j]
    }
  }
  slope_sim-info$mean_slope
}

mc_neutral_test <- function(info, k, d, pop, Ne=1000L, n_sim=4000L){
  if(!isTRUE(info$ok))return(list(p=NA_real_, exceed=NA_integer_, B=0L,
    neutral_sd=NA_real_,normal_fp01=NA_real_,obs=NA_real_))
  draws <- simulate_conditional_neutral(info,k,d,pop,Ne=Ne,n_sim=n_sim)
  observed <- info$slope-info$mean_slope
  exceed <- sum(abs(draws)>=abs(observed))
  p <- (exceed+1)/(length(draws)+1)
  list(p=p,exceed=exceed,B=length(draws),neutral_sd=sd(draws),
       normal_fp01=mean(2*pnorm(abs(draws/info$conditional_SE),lower.tail=FALSE)<.01),
       obs=observed,sim=draws)
}
