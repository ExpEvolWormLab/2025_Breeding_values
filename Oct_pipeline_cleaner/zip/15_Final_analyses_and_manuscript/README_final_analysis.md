# Final analysis: script 15

## What to run

Copy `15_Final_prediction_chrIV_analysis.R` to the CeMEE project root and run:

```r
source("15_Final_prediction_chrIV_analysis.R")
```

**Dependencies:** requires the corrected Pool-seq alignment and already completed scripts 09f, 11 and 13. Base R only. It reads:

- `selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Cluster_predictors_and_response.csv`
- `selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/13_Leave_chromosome_out_fitted.csv`

It doesn't run another GWAS, selection fit or simulation. The new area.F-only predictor is calibrated on five chromosomes and evaluated on the held-out chromosome. Existing LOCO predictions are not refitted.

**Outputs** in `selection_results/drift_oneG0_nullfit/R5/multivariate_prediction/final_chrIV/`:

- `15_Model_comparison.csv`: descriptive in-sample correlations and held-out predictive R² for all clusters, non-IV clusters, and IV alone (subset evaluation; not re-trained on excluded regions).
- `15_Heldout_block_bootstrap.csv`: 2-Mb genomic-block bootstrap for paired MSE differences, particularly full minus three-trait, full minus SB and full minus area.F.
- `15_By_chromosome_paired_errors.csv`: one row per chromosome and model comparison.
- `15_CHRIV_cluster_prediction_errors.csv`: errors per GWAS cluster, with genomic positions.
- `15_Chromosome_IV_prediction_diagnostics.pdf`: four-panel summary including IV spatial residuals and overall prediction.
- `15_Analysis_QC.txt`: concise run summary.

## Two optional chromosome-IV mutations

After the first run, fill the templates `15_chrIV_mutations.csv` and `15_chrIV_mutation_counts.csv` created in the project root. **Alternatively edit the included templates before the first run.**

`15_chrIV_mutations.csv` must contain `mutation_id,chromosome,position_bp`, with positions in the same reference genome assembly as the GWAS SNP positions; optional `sb_association_note` and `origin_evidence` columns are informative only. No positions were supplied by the user, so both position cells are intentionally empty. Re-run script 15 after entering them. The script will then output `15_CHRIV_mutation_window_sensitivity.csv` and show the positions in the figure. The analysis compares windows ±0.25, ±0.5, ±1 and ±2 Mb to the effect of excluding random chromosome-IV-centered windows, an exploratory contextual comparison not a valid causal test or a full LD-matched null.

`15_chrIV_mutation_counts.csv` contains one row per mutation and sequenced pool. Supply `mutant_reads`, `depth` for the exact same mutant allele at every time point (preferred), or `allele_frequency` directly. Keep `sample_id`, `generations`, and `replicate` consistent with the checked sample manifest; the bottlenecked generation-zero sample can be omitted if scientifically justified. Once populated, `15_CHRIV_mutation_trajectories.pdf` is produced. For a mutation not present in the founder-SNP reference, obtain raw information separately; the script-08 filtered Pool-seq matrices **cannot** serve as evidence for it.

## Interpretation and stopping rule

The negative chromosome-IV MSE contrast in the current analysis is small (−3.53e−9 per-cluster squared-frequency-slope units). **Do not conclude the mutations explain that difference solely because both occur on IV.** Report a regional effect only if the independent mutation coordinates localize the excess prediction error and the comparison is robust to window width; otherwise report an interesting chromosome-IV exception that remains unexplained.

The full seven-trait score had held-out R² 0.051; the reduced three-trait score 0.035 and SB alone 0.054. A positive predictive effect is supported, but not superiority over SB alone. The current selection gradient intervals exclude zero only for area.F; other gradients can contribute as uncertain weighted effects but are not individually demonstrated selection targets.

## Check before manuscript submission

- Original draft refers to **three** evolved populations. The manifest labels four evolved trajectories (one missing generation 50). Reconcile with experiment log.
- Verify whether `A6G140` is the independently established bottlenecked G0 sample, as configured in script 09f. Do not infer this only from a PCA plot.
- Distinguish new mutations arising during the salt-selection experiment from low-frequency standing alleles or variants that arose during prior domestication.
- Manuscript Results/Methods draft is in `Methods_and_Results_draft.md` and `.docx`, with unverified claims clearly marked.

The R script has been statically reviewed but **not executed on the user's private datasets** in this environment.
