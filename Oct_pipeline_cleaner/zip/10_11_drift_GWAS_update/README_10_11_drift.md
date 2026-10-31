# Link 09f drift-aware Pool-seq trajectories to R5 GWAS clusters

Run after `09f_Drift_GLS_nullfit.R` (and optionally `09g_WF_calibrate_nullfit.R`) from the project root:

```r
source("10_Selection_GWAS_overlap_09f.R")
source("11_GWAS_selection_cluster_associations_09f.R")
```

Both scripts read `selection_results/drift_oneG0_nullfit/09f_SNP_drift_Ne1000.csv` and write to `selection_results/drift_oneG0_nullfit/R5/`. They do not overwrite scripts 09, 09c, 09f, or earlier 10/11 output. They presume `threshold_clust = 0.5` and `Ne = 1000` by default.

**Interpretation:** The observed 157 common-frequency genome-wide 09f BH q<0.01 SNPs are *provisional*. Neutral simulation calibration is encouraging for ancestral frequency 0.10–0.90 at p = 0.01 down to ~0.0001, but cannot validate much smaller genome-wide tail p-values. The output deliberately calls these *provisional* rather than definitive selected loci. The selected-GWAS-cluster overlap is a descriptive count of preidentified GWAS-cluster groups containing one or more provisional SNPs; SNPs in LD are not independent.

**Matching and orientation:** GWAS SNP `rs` is `CHROM_POS`, and 09f `ID` is `CHROM_POS`, while 09f `rs` is `CHROM.POS`. Scripts explicitly match to `ID` and never use row order. 09f's `slope_freq` is a *linear ALT-frequency change per generation*; NOT a logit coefficient or an exact selection coefficient. 03's cluster dosage is ALT dosage for `phase=+1`, and `2 - ALT dosage` for `phase=-1`. Thus 11 multiplies ALT-frequency slope by the existing ±1 phase. Each cluster is summarized by the median of available phase-aligned SNP slopes, primarily those with ancestral frequencies in [0.1,0.9]. All-frequency summaries are retained solely as a sensitivity comparison.

**Files:**

- `10_Overlap_summary.csv`: total genome-wide discovery counts; how many GWAS SNPs and GWAS clusters have coverage, a common-frequency test, or a provisional discovery.
- `10_GWAS_cluster_overlap.csv`, `10_GWAS_SNPs_with_09f.csv`, `10_Selected_common_SNPs_q01_PROVISIONAL.csv`: detailed descriptive join.
- `11_Clusters_selection_and_GWAS_09f.csv`, `11_SNP_phase_aligned_09f.csv`: phase-corrected continuous selection trajectories with all coverage and frequency-subset flags.
- `11_GWAS_selection_associations_COMMON_exploratory.csv`: eight standardized GWAS effect vs. median continuous common-frequency slope associations, using cluster as the analysis unit.
- `11_GWAS_selection_associations_ALL_sensitivity.csv`: comparison using all frequency classes, including poorly calibrated boundary alleles.
- `11_GWAS_vs_drift_frequency_slopes_COMMON.pdf`: eight trait panels; red points indicate clusters with at least one provisional common-frequency SNP at BH q<0.01.

**Association inference:** Ordinary lm p-values are descriptive (GWAS cluster selection, LD across clusters, unequal cluster information). A supplementary chromosome-wise circular-shift association test rotates frequency-slope vectors within chromosomes, retaining each chromosome's ordering/autocorrelation. This is a *sensitivity/exploratory null*, not a calibrated formal proof: the wraparound, marker spacing, genome-wide allele sharing, and GWAS ascertainment may violate exchangeability. Holm adjustment across eight phenotypes is provided, but does not cure those limitations. Interpret direction, effect size, coverage, consistency across chromosomes and LD thresholds; do not claim causal overlap or independently selected loci.

**Design checks:** Confirm that A6G140 is the bottlenecked G0 sample. The script logs a warning if 09f retained/excluded G0 metadata differ from A00/A6G140. The initial script 08 SNP coverage filtering was performed with both G0 samples and has not been rebuilt; for definitive reporting, rerun coverage filtering on the 12-sample design, as noted in previous analysis.

The scripts were statically checked but cannot be executed on your private Pool-seq/GWAS results here, and no R installation is available in this environment.
