MCMCglmm multivariate among-RIL covariance: fix and recovery
============================================================

WHAT FAILED
Your original script 07a fitted three MCMCglmm chains and then called
parse_block(), which incorrectly expected 28 covariance terms for 7 traits.
MCMCglmm supplied a full 7 x 7 matrix (49 named columns). The error is
post-processing; it is not evidence that the fitted chains were invalid.

ROUTE A: ORIGINAL R SESSION STILL OPEN AFTER THE ERROR
1. Check: length(chains_VCV)    # expected 3
2. Set the working directory to the same project root as your earlier run.
3. Copy 07a_RESCUE_completed_chains_FROM_CURRENT_R_SESSION.R there.
4. In THAT SAME R SESSION run:
   source("07a_RESCUE_completed_chains_FROM_CURRENT_R_SESSION.R")
5. It should save:
   RData_G_NaCl/07a_final_unstructured_chain01.rds
   RData_G_NaCl/07a_final_unstructured_chain02.rds
   RData_G_NaCl/07a_final_unstructured_chain03.rds
6. Copy 07a_Estimate_multivariate_RIL_G_FIXED.R into the project root.
7. Change `postprocess_only <- FALSE` to TRUE (near the start) for safety.
8. Run:
   source("07a_Estimate_multivariate_RIL_G_FIXED.R")
   The script should load all three chains and generate output without refitting.

ROUTE B: ORIGINAL SESSION CLOSED BUT .RData WORKSPACE BACKUP EXISTS
If you saved an R workspace after the error, load it first and check
that chains_VCV and the other original model variables are present.
Then follow Route A. Do not source the original script again.

ROUTE C: NO IN-MEMORY CHAINS AND NO CHECKPOINTS/WORKSPACE BACKUP
The original script did not save completed chains to disk before failing.
If no recoverable workspace exists, MCMC must be refit. Run the fixed 07a
script with `postprocess_only <- FALSE` (its default). It checkpoints each
completed chain atomically. On restart, it reuses compatible checkpoints
and fits only missing chains.

VALIDATION
Inspect RData_G_NaCl/07a_MCMC_diagnostics.csv:
 - Rhat ideally <1.05 for all relevant parameters, and ESS >=200.
 - Inspect trace / autocorrelation plots as well; thresholds alone do not
   establish convergence, identifiability or prior robustness.
 - RIL covariance here is among-line covariance; NOT automatically an
   additive genetic covariance matrix G_A (no pedigree/GRM used).
 - Inspect 07a_G_correlation_credible_intervals.csv and
   07a_G_variance_credible_intervals.csv.
 - The output 07a_multivariate_G_final.RData is consumed by 07b.

After satisfactory diagnostics, run the supplied stricter
07b_Refit_pleiotropy_with_multivariate_G_FIXED.R.

These scripts have been statically checked but NOT run against your data:
this environment has neither the phenotype files nor an R executable.
