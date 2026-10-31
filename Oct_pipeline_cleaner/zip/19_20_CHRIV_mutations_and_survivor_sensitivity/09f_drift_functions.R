## Shared utilities for scripts 09c and 09d.
## Alleles are always oriented to the GWAS ALT allele.

read_drift_inputs <- function(pool_file, selection_file, exclude_G0 = character()) {
  p <- new.env(parent = baseenv())
  load(pool_file, envir = p)
  sel <- read.csv(selection_file, stringsAsFactors = FALSE)
  pop <- p$populations
  pop$sample_id <- as.character(pop$sample_id)
  pop$generations <- as.integer(pop$generations)
  pop$replicate <- as.character(pop$replicate)

  if (nrow(pop) != 13L || sum(pop$generations == 0) != 2L)
    stop("Expected the verified 13-pool design with TWO G0 samples")
  if (anyDuplicated(pop$sample_id)) stop("Duplicate sample IDs")

  # Exclusion is justified by independent experimental history, not p-values.
  # Original sync columns and GWAS allele orientation are never modified.
  if (length(exclude_G0)) {
    if (!all(exclude_G0 %in% pop$sample_id))
      stop("Excluded sample not found in manifest: ",
           paste(setdiff(exclude_G0, pop$sample_id), collapse = ", "))
    if (any(pop$generations[match(exclude_G0, pop$sample_id)] != 0L))
      stop("Only generation-zero samples may be excluded here")
  }
  remaining <- pop$sample_id[!pop$sample_id %in% exclude_G0]
  pop <- pop[match(remaining, pop$sample_id), , drop = FALSE]
  if (sum(pop$generations == 0L) != 1L && length(exclude_G0))
    stop("Expected exactly one retained generation-zero sample")

  minor <- p$all.counts.caendr
  coverage <- p$all.cov.caendr
  if (!identical(dim(minor), dim(coverage)) ||
      !identical(rownames(minor), rownames(coverage)))
    stop("Pool-seq count and coverage matrices are not aligned")
  if (!all(pop$sample_id %in% colnames(minor)))
    stop("Pool-seq matrix lacks retained samples")
  minor <- minor[, pop$sample_id, drop = FALSE]
  coverage <- coverage[, pop$sample_id, drop = FALSE]

  if (anyDuplicated(sel$rs) || anyDuplicated(sel$ID))
    stop("Selection results contain duplicate SNP IDs")
  ix <- match(sel$rs, rownames(minor))
  if (anyNA(ix)) stop("Some script-09 SNPs are missing from script 08")
  minor <- minor[ix, , drop = FALSE]
  coverage <- coverage[ix, , drop = FALSE]

  if (!all(sel$orientation %in% c("minor_is_alt", "minor_is_ref")))
    stop("Unresolved GWAS-ALT allele orientation")
  k_alt <- minor
  rev <- which(sel$orientation == "minor_is_ref")
  if (length(rev))
    k_alt[rev, ] <- coverage[rev, , drop = FALSE] - minor[rev, , drop = FALSE]

  if (any(k_alt[is.finite(k_alt)] < 0) ||
      any(k_alt[is.finite(k_alt) & is.finite(coverage)] >
          coverage[is.finite(k_alt) & is.finite(coverage)]))
    stop("Some ALT counts are outside [0, depth]")

  list(sel = sel, pop = pop, k_alt = k_alt, depth = coverage, excluded_G0 = exclude_G0)
}

## Returns a site-specific GLS slope test of a common linear allele-frequency
## trend, with neutral Wright-Fisher covariance and read-sampling uncertainty.
## A linear frequency slope is a trajectory statistic, NOT a direct estimate of s.
## Null-score approximation: p(1-p) uses the pooled GLS intercept under H0.
fit_drift_site_null <- function(k, d, pop, Ne = 1000,
                           pool_individuals = NA_real_) {
  t <- pop$generations
  r <- pop$replicate
  good <- is.finite(k) & is.finite(d) & d > 0 & k >= 0 & k <= d
  keep <- which(good)

  empty <- list(p0 = NA_real_, p0_null = NA_real_, slope = NA_real_, SE = NA_real_,
                z = NA_real_, p = NA_real_, n_pools = sum(good),
                n_G0 = sum(good & t == 0), n_evolved_reps = NA_integer_,
                baseline_diff = NA_real_, baseline_Z = NA_real_,
                weights = rep(NA_real_, length(k)),
                null_weights = rep(NA_real_, length(k)),
                slope_scale = NA_real_)
  if (sum(good & t == 0) == 0L || sum(good & t > 0) < 4L ||
      length(unique(t[good])) < 3L ||
      length(unique(r[good & t > 0])) < 2L)
    return(empty)
  if (!is.finite(Ne) || Ne <= 1) stop("Ne must be >1 diploid individuals")

  initial <- good & t == 0
  ## Jeffreys smoothing prevents a zero estimated covariance at boundaries.
  p0 <- (sum(k[initial]) + .5) / (sum(d[initial]) + 1)
  y <- k[keep] / d[keep]
  tk <- t[keep]
  rk <- r[keep]
  depths <- d[keep]

  ## Neutral drift covariance: Cov(p_s,p_t) ~ p0(1-p0) F_min(s,t),
  ## for times from the SAME evolving replicate; G0 measurement errors
  ## are independent conditional on their shared unknown ancestral p0.
  fmin <- 1 - (1 - 1 / (2 * Ne)) ^ outer(tk, tk, pmin)
  same_rep <- outer(rk, rk, "==")
  C <- fmin * same_rep
  diag(C) <- diag(C) + 1 / depths

  ## Two-stage sampling: optionally include sampling N individuals into a pool.
  ## This is an approximation (independent samples, without finite-pop. corr.).
  if (is.finite(pool_individuals)) {
    if (pool_individuals < 2) stop("pool_individuals must be >= 2")
    diag(C) <- diag(C) + (1 - 1 / depths) / (2 * pool_individuals)
  }

  X <- cbind(1, tk)
  R <- try(chol(C), silent = TRUE)
  if (inherits(R, "try-error")) return(empty)
  rhs <- cbind(X, y)
  sol <- backsolve(R, forwardsolve(t(R), rhs))
  S <- crossprod(X, sol[, 1:2, drop = FALSE])
  v <- crossprod(X, sol[, 3])
  detS <- S[1, 1] * S[2, 2] - S[1, 2]^2
  if (!is.finite(detS) || detS <= 0) return(empty)

  slope <- as.numeric((S[1, 1] * v[2] - S[1, 2] * v[1]) / detS)
  ## Restricted (neutral-null) ancestral-frequency estimate.
  ## The null hypothesis is a shared constant true frequency across all pools.
  ## GLS covariance C has been factored as p0*(1-p0)*C, so the GLS mean
  ## estimate under H0 does not depend on that unknown scale factor.
  ## This uses ALL pools to stabilize the plug-in null SE, unlike the
  ## G0-only plug-in which produced an anticonservative extreme tail.
  p0_null_raw <- sum(sol[, 1] * y) / sum(sol[, 1])
  lower <- 0.5 / (sum(depths) + 1)
  p0_null <- max(lower, min(1 - lower, p0_null_raw))
  se <- sqrt(p0_null * (1 - p0_null) * S[1, 1] / detS)
  if (!is.finite(se) || se <= 0) return(empty)

  ## Slope linear weights: b = w' y. Used in exact simulations without
  ## refitting every replicate. Zero weights for unobserved samples.
  Sinv <- matrix(c(S[2, 2], -S[1, 2], -S[1, 2], S[1, 1]),
                 2, 2) / detS
  w <- drop(sol[, 1:2, drop = FALSE] %*% Sinv[, 2])
  weights <- rep(0, length(k))
  weights[keep] <- w
  u <- drop(sol[, 1]) / sum(sol[, 1])
  null_weights <- rep(0, length(k))
  null_weights[keep] <- u
  slope_scale <- sqrt(S[1, 1] / detS)

  g0 <- which(good & t == 0)
  bdiff <- if (length(g0) == 2L)
    k[g0[2]] / d[g0[2]] - k[g0[1]] / d[g0[1]] else NA_real_

  bse <- NA_real_
  if (length(g0) == 2L) {
    bvar <- p0 * (1 - p0) * sum(1 / d[g0])
    if (is.finite(pool_individuals))
      bvar <- bvar + p0 * (1 - p0) *
        sum((1 - 1 / d[g0]) / (2 * pool_individuals))
    if (is.finite(bvar) && bvar > 0) bse <- bdiff / sqrt(bvar)
  }

  z <- slope / se
  list(p0 = p0, p0_null = p0_null, slope = slope, SE = se, z = z,
       p = 2 * pnorm(abs(z), lower.tail = FALSE),
       n_pools = length(keep), n_G0 = length(g0),
       n_evolved_reps = length(unique(r[good & t > 0])),
       baseline_diff = bdiff, baseline_Z = bse, weights = weights,
       null_weights = null_weights, slope_scale = slope_scale)
}

## Exact discrete-generation *neutral* diploid Wright-Fisher simulation.
## Any retained G0 pool is independently read-sampled from the starting p0;
## the four evolved populations drift independently conditional on p0.
## Coverage is held at the observed SNP-by-pool depth.
## If pool_individuals is set, sample diploid chromosomes into each library
## before drawing sequencing reads (otherwise assumes effectively large pools).
simulate_wf_site <- function(p0, depths, pop, Ne, n_sim = 500L,
                             pool_individuals = NA_real_) {
  t <- pop$generations
  r <- pop$replicate
  evolve <- unique(r[t > 0])
  B <- as.integer(n_sim)
  Y <- matrix(0, nrow = B, ncol = length(t))
  good <- is.finite(depths) & depths > 0

  sample_reads <- function(true_p, depth) {
    p <- true_p
    if (is.finite(pool_individuals)) {
      p <- rbinom(length(p), size = 2L * as.integer(pool_individuals),
                  prob = p) / (2L * as.integer(pool_individuals))
    }
    rbinom(length(p), size = as.integer(round(depth)), prob = p) / depth
  }

  for (j in which(t == 0 & good))
    Y[, j] <- sample_reads(rep(p0, B), depths[j])

  P <- matrix(p0, nrow = B, ncol = length(evolve))
  for (gen in seq_len(max(t))) {
    P[] <- rbinom(length(P), size = 2L * as.integer(Ne),
                   prob = as.vector(P)) / (2 * Ne)
    jj <- which(t == gen & good)
    if (length(jj)) for (j in jj) {
      rr <- match(r[j], evolve)
      Y[, j] <- sample_reads(P[, rr], depths[j])
    }
  }
  Y
}
