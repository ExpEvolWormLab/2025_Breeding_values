# Follow-up: correcting the one-G0 neutral tail

Your neutral WF calibration gave an approximately nominal 0.05 tail but a strongly inflated 0.001 tail (for Ne=1000, 0.00615 instead of 0.001). **Therefore previous 09c BH discoveries should not be claimed as selected SNPs.** The issue is present for all Ne settings. With only a single G0, its observed frequency is a noisy plug-in estimate for the variance `p0*(1-p0)`; particularly at the boundaries this can badly distort extreme Z values.

This is a **candidate fix requiring simulation validation**, not a validated replacement. It fits the same GLS slope and drift/read covariance but evaluates the null variance at a GLS estimate of the common allele frequency using **all pools under the neutral null**, instead of the G0-only frequency. This is a standard restricted-null (score-test-like) variance estimate. It makes no claim that later samples actually experienced no selection: estimating the variance *under the null* is how this hypothesis test is defined.

The simulation script re-estimates both the slope and the null-frequency estimator in each simulated dataset. Since only genotype frequencies change and the covariance geometry is fixed by coverage and generations, this refit is algebraically vectorized without repeating the same 12x12 Cholesky decomposition thousands of times.

## Files and order

Place `09f_drift_functions.R`, `09f_Drift_GLS_nullfit.R`, and `09g_WF_calibrate_nullfit.R` together in the project root. Run:

```r
source("09f_Drift_GLS_nullfit.R")
source("09g_WF_calibrate_nullfit.R")
```

Check **both** files retain the correct exclusion (currently `exclude_G0 <- "A6G140"`; change if A00 was the bottlenecked sample). The old 09c/09d results remain unchanged. Output directory: `selection_results/drift_oneG0_nullfit/`.

Read:

- `09g_WF_null_calibration_summary.csv` — Ne-wide false-positive rates at actual p thresholds of 0.05, 0.01, 0.001, 0.0001.
- `09g_WF_calibration_by_frequency.csv` — **most important**, because rare-frequency loci can remain anti-conservative even if the global mean/SD looks fine. Do not aggregate away frequency-specific behavior.
- `09f_SNP_drift_Ne1000.csv` — estimated effects and approximate p-values; its `p_drift_BH_common` gives a *prespecified common-G0-frequency sensitivity subset* (0.10–0.90), but these calls are **provisional until their null calibration passes**. `p_drift_BH_all` refers to all sites and should not be treated as reliable while rare-allele tails are inflated.

The pilot samples up to 12 loci per ancestral-frequency bin and simulates 500 neutral trajectories per locus and Ne, totaling around 100 loci × 500 × 3 Ne. Increase `n_sim` and `random_per_freq_bin` only if needed; at p<0.0001 pilot counts are still noisy. This does **not** validate the very small probabilities required for genome-wide BH; a more rigorous tail model or calibrated empirical genome-wide null would still be needed for definitive selection calls.

## Assumptions and outstanding limitations

- Diploid WF Ne=500/1000/2000, four independent evolving lines sharing a common ancestor; conditional on one retained G0 Pool-seq observation.
- Sequencing-read binomial sampling; unspecified number of individual animals contributing to each library (`pool_individuals=NA_real_`). This is **not** a beta-binomial or full two-stage sampling likelihood.
- G0 frequency is noisy; simulation uses the observed G0 estimate as the assumed *true neutral frequency* for each sampled locus. This conditions on estimates and does not model uncertain ancestry beyond its resampling of G0 counts in every simulation.
- The SNP subset was filtered in old script 08 with 13 samples; for final analysis, rebuild filtering excluding bottlenecked G0 if necessary.
- Linkage makes nearby SNP tests non-independent; count genomic regions rather than calling each SNP an independent selection event.
- Analytic GLS assumes linear *frequency* trends, not a precise biological selection coefficient s.
- No R interpreter or private biological matrices are available in the creation environment: code is **statically inspected, not executed on your Pool-seq dataset**.
