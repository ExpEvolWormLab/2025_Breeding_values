# Pleiotropy manuscript: finalization analyses 59–62

These scripts extend the **frozen**, frequency-inclusive 802-cluster R5 analysis. They do not revise the GWAS, R5 cluster boundaries, sample manifest, phase alignment, ancestry labels, or main selection-gradient estimates. Run them **from the project root** with the existing `Oct_pipeline_cleaner` directory present for script 62.

## Run order and intended decisions

| Script | Main question | Input requirements | Typical output |
| --- | --- | --- | --- |
| **59** `59_Nested_chromosome_trait_comparison.R` | Does an *a priori* seven-trait beta perform better than a single trait selected without using the held-out chromosome? | Script-11 R5 effects; script-13 weights; script-45 802-cluster outputs | LOCO predictions, nested fold trait choices, block-bootstrap MSE intervals, sensitivity with V_3 excluded |
| **60** `60_Paired_RIL_pheno_prediction.R` | Does GBLUP exceed a fixed 75-marker GWAS predictor for the *same held-out RILs*? | Individual-level `51_Complete_three_method_RIL_predictions.csv` (not merely metrics); optionally script-40 predictions | Paired predictive R², MSE, confidence intervals; publication plot |
| **61** `61_Joint_gradient_uncertainty_802.R` | Are genome-wide beta predictions and the V_3 directional predictions stable to *joint* gradient uncertainty? | Script-00 prep RData including `model_beta`, selection-gradient CSV, scripts 11, 13, 45 | Draw-specific 802-cluster LOCO R², gradient intervals, per-cluster directional stability |
| **62** `62_Fixed_p0_WF_validation_resumable.R` | Is the published conditional WF neutral test calibrated when the true ancestral frequency is FIXED but A00 reads are sampled? | Original 43b inputs: Pool-seq matrices, raw sync, R5 tags, 43 config/helpers, 09f and 33 functions | Restartable simulation checkpoints, empirical type-I curves, exact binomial CIs and multiplicity-adjusted inflation diagnostics |

## Calls

```r
source('Oct_pipeline_cleaner/59_Nested_chromosome_trait_comparison.R')
source('Oct_pipeline_cleaner/60_Paired_RIL_pheno_prediction.R')
source('Oct_pipeline_cleaner/61_Joint_gradient_uncertainty_802.R')
# 62: FIRST run pilot, inspect outputs, then set CONFIG$mode='final'
# Script 62 automatically uses separate checkpoint directories for pilot/final.
source('Oct_pipeline_cleaner/62_Fixed_p0_WF_validation_resumable.R')
```

Scripts can be run independently (different input requirements). No MCMCglmm refit is needed. Script 62 is the most computationally demanding. To extend it to Ne sensitivity, change `CONFIG$Ne_grid`; use a **new output directory** when changing settings.

## Evidence and interpretation

- **59** selects each single trait *inside* the outer training chromosomes via inner chromosome CV; calibration still uses evolutionary data from the other chromosomes. Retrospective selection of SB across the entire dataset is deliberately NOT the fair comparison. Bootstrap intervals are conditional on fitted predictions and do not refit nested model selection in each bootstrap draw.
- **60** compares only held-out predictors. The historical in-sample QTL reconstruction is explicitly excluded. Bootstrap source-population groups when available; otherwise RIL bootstrap intervals are potentially optimistic due to relatedness.
- **61** uses the original fitness-model coefficient covariance. These intervals are conditional on frozen GWAS effects and evolved observations: they do not propagate GWAS discovery, RIL sampling, evolving-population variance, or LD dependencies. It does **not** use the MCMCglmm G posterior as selection-gradient uncertainty.
- **62** mirrors the **original** conditional WF test, including finite A00 read sampling, diploid Wright–Fisher simulation, and original Pool-seq depths. It explicitly counts non-estimable null datasets as non-rejections and flags possible anti-conservatism. The pilot/final outer simulation probes nominal 0.01 and 0.05 tails ONLY; even a non-rejection of an inflation test **does not validate the approximately 1e-5 tail needed for 802-test BY selection significance**. Thus the 12 candidate selection clusters remain conditional until proper extreme-tail/null-calibration evidence is available. This file does not replace `43b_Calibration_interpretation_gate.txt` and cannot authorize 44 significance claims.

## Previously established baseline

- 1,089 R5 LD association groups, 802 with testable phase-aligned evolutionary response.
- Full seven-trait LOCO R² ≈ 0.0796; SB ≈ 0.0716; covariance PC1 ≈ 0.0474.
- Major chromosome V block V_3 is influential and LD-correlated; a four-tag haplotype comparison and GA omission analysis were already done in scripts 54–58. Do not treat its neighboring tags as four independent selection events.
- MCMCglmm among-RIL covariance finished; the pleiotropy-scaling excess over a naive N_eff null is not yet demonstrated to be caused by molecular synergistic pleiotropy.

## Manuscript stop/go recommendation

After scripts 59–61, freeze genome-wide prediction conclusions and revise the figures. If script 62 reveals inflation or cannot establish biologically relevant extreme tails, remove definitive BY-selected-locus claims from the main results. The paper's strongest defensible contribution remains **partial, independently measured trait-effect-based prediction of allele-frequency evolution**, with chromosome I/IV/V as regional illustrations of linkage, founder ancestry and domestication-origin alleles.

## Missing data (local project only)

The uploaded summary archive does not contain `51_Complete_three_method_RIL_predictions.csv`, `40_Out_of_fold_GBLUP_predictions.csv`, `RData/Preparatory_GLMM_NaCl_only_withBeta.RData`, or the raw Pool-seq inputs. Scripts 60–62 therefore require those **already-existing local project files**. We intentionally fail with a clear missing-file message rather than silently change the analysis.

## Validation status

Script 59 was independently re-implemented and checked against script 45 (max deviation < 1e-9) on all 802 clusters. Scripts 60–62 are designed from the supplied original R code but **could not be run** in this environment because R / required local raw inputs are unavailable. Check their generated QC and diagnostic output before accepting results.
