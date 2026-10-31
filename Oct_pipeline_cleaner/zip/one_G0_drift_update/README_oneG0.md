# Drift-aware Pool-seq analysis with one retained G0 sample

## Scientific decision

Your independent experimental history suggests one G0 sample underwent a bottleneck. If it does not represent the common ancestral population of the four evolved replicates, excluding that sample from the drift analysis is justified *independently of SNP significance*. The previous G0 comparison was overdispersed: SD(Z)=1.58, 21.4% beyond ±1.96, with no major sequencing-depth trend. The earlier rho=0.0166 is **not** a valid estimate of technical dispersion because it may reflect the real bottleneck; do not insert it as an error correction.

**Verify the excluded sample ID.** This package defaults to `exclude_G0 <- "A6G140"` because the former PCA code explicitly omitted that sample. If the bottlenecked sample was `A00`, change `exclude_G0` to `"A00"` in **both** scripts before running. The retained G0 must represent the ancestor of the four evolved replicates.

## Files

- `09c_drift_functions.R`: modified shared input loader that drops the specified G0 **after** reading and aligning the prepared 13-sample matrices; the neutral GLS code itself is unchanged and already handles one baseline.
- `09c_Drift_GLS_oneG0.R`: genome-wide drift-aware GLS at primary Ne=1000 (optional 500 and 2000). Results go to **`selection_results/drift_oneG0/`** rather than overwriting your existing two-G0 analysis. The legacy logit/binomial slopes are retained under names ending `_twoG0`, because they have *not* been recomputed for the one-G0 design.
- `09d_WF_neutral_calibration_oneG0.R`: frequency- and depth-stratified neutral WF simulation, including the rare-frequency bins. Re-estimates the ancestral frequency and refits the full GLS statistic on **every simulated dataset**, unlike the earlier fixed-weight shortcut. This improves the calibration check near frequency boundaries but increases the number of small GLS fits. Pilot default: 200 simulations per locus, 2 random loci per stratum, 10 top candidates and Ne=500/1000/2000.

## Run from project root

Place all three `.R` files in the same directory, then:

```r
source("09c_Drift_GLS_oneG0.R")
source("09d_WF_neutral_calibration_oneG0.R")
```

The code reads the original `08_PoolSeq_prepared.RData` and `09_Selection_SNP_results.csv` because those contain the verified SNP/allele matches. It does **not** use the old two-G0 selection estimates for its new drift test. The prepared SNP universe was filtered using all 13 samples in script 08; this is acceptable as an initial controlled comparison, but for a fully one-G0 filtering policy, rerun script 08 with filtering based only on retained samples, while maintaining the original sync indices and allele orientation metadata.

## Expected verification

```r
res <- read.csv("selection_results/drift_oneG0/09c_SNP_drift_Ne1000.csv")
print(table(res$n_G0, useNA="ifany"))  # finite sites should use 1
print(unique(res$retained_G0))          # A00 unless you changed exclusion
print(sum(is.finite(res$p_drift)))
print(c(BH_001=sum(res$p_drift_BH < 0.01, na.rm=TRUE),
        BH_005=sum(res$p_drift_BH < 0.05, na.rm=TRUE)))
print(read.csv("selection_results/drift_oneG0/09d_WF_null_calibration_summary.csv"))
print(read.csv("selection_results/drift_oneG0/09d_WF_calibration_by_stratum.csv"))
```

`baseline_Z` and `baseline_freq_difference` are correctly **NA** with only one G0 sample. There is no direct G0-versus-G0 dispersion estimate after exclusion. The neutral simulations are conditional on the retained estimated ancestral frequency, with that frequency re-estimated for the simulated test; their results assess model calibration, not the biological identity of the retained ancestral sample.

**Do not yet use the number of BH-significant SNPs as definitive selection calls.** Check the frequency-specific simulated null (especially <0.02 and >0.98), sensitivity to Ne, and availability of the actual number of individuals contributing DNA to each pool. The simulation pilot has too little extreme-tail resolution for formal genome-wide Monte Carlo FDR. Linkage also makes neighboring SNP detections non-independent; interpret genomic regions or clusters rather than treating all significant SNPs as separate selection events.

**No need to rerun the GWAS, phase-aligned QTL clustering, cluster mapping, or MCMCglmm for this change.** Scripts 10–11 will need to reference the new `drift_oneG0` output after calibration is acceptable.

## Limitations

- Original 08 coverage filtering still used 13 pools; see above.
- This remains an approximate Gaussian GLS using a Wright–Fisher covariance; it is not a joint likelihood for s, Ne, technical overdispersion, and correlated founder sampling.
- A single G0 improves biological consistency if the other G0 is unrepresentative, but reduces precision for initial frequency. Its read/pool-sampling uncertainty remains included.
- If the bottlenecked sample actually founded the evolved populations, simply discarding it would be biologically incorrect. Confirm experimental provenance.
- Files were statically checked but not executed against your R data, and R is not available in this environment.
