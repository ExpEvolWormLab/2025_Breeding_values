# Pleiotropy analysis summary

**4 November 2026 — for scientific discussion**

The primary result is reproducible: ancestral fitness-weighted GWAS effects predict a modest component of subsequent genomic frequency change across 802 R5 clusters. Multivariate prediction is not demonstrably superior to the single-trait comparators. Exact focal trajectories are now complete, while confirmatory selection-outlier and causal-variant claims remain unsupported.

## Modifications

New scripts **69–74** freeze the R5 registry, SNP membership and signed phase, shared-A00 sample manifest, GWAS effects and gradient weights; reproduce chromosome-held-out predictions; audit neutral calibration; assess regional genotype identifiability; extract exact Pool-seq counts; and generate four main figures, eight supplementary figures and a manuscript consistency audit. Historical scripts and numerical outputs were preserved.

The regional analysis retains the original phenotype inclusion, source-population covariates, founder proxies and chromosome-excluded kinship. Added conditional maximum-likelihood comparisons hold the baseline covariance fixed; they do not compare REML likelihoods across different fixed effects. New multiple-testing families are labelled separately from historical corrections.

## Main numerical results

| Analysis | Result |
|---|---:|
| R5 catalogue / covered / eligible clusters | 1,089 / 804 / 802 |
| Raw fitness-weighted score–response correlation | r = 0.304510 |
| Historical 2-Mb block 95% interval for r | 0.199220–0.431139 |
| Seven-trait β LOCO predictive R² | 0.079606 |
| Fixed SB / nested single-trait / PC1 R² | 0.071579 / 0.053430 / 0.047407 |
| β R² excluding V_3, with refitting | 0.053385 |
| Matched RIL β prediction: GBLUP / 75-marker R² | 0.103350 / 0.052386 |

These reproduce the historical benchmarks. Paired block-bootstrap intervals for β's advantage over SB and the nested comparator include zero. The response is a **phase-aligned GLS allele-frequency slope**, not a selection coefficient; LOCO calibration uses evolutionary responses from the other chromosomes.

All eight focal sites were recovered across 13 samples: **104 observed records**, including both chromosome IV variants, with no missing or zero-depth records. All REF alleles match WS283. The primary subset contains 80 records, with A00 represented once per site.

For chromosome V, the β-model cest-24/scd-2 genotype classes 00/01/10/11 contain **863/13/10/54 RILs**. Focal conditional coefficients are statistically estimable, but the small discordant classes do not establish causal separation from linked variation. The scd-2 FB/β associations remain nonsignificant after the new focal-family correction (q ≈ 0.109; historical broader-family q ≈ 0.151).

## Tests performed

The analysis run passed both existing script-42 smoke checks and **18 lightweight checks**, including key/row alignment, duplicate detection, phase and coordinate consistency, matrix dimensions, shared-A00 handling, singular-model diagnostics, frozen hashes and exact-site validation. Held-out-response perturbation tests found no leakage into the corresponding β, SB, nested or PC1 predictions.

An independent complete run reproduced **68 numerical tables byte-for-byte**. Historical regional sample counts and nested F-test p-values reproduced exactly. Completed MCMC checkpoint metadata and input hashes were verified: three chains and 6,000 retained draws. No full GWAS, MCMC, phenotype-CV fit or raw-count GLS catalogue was rerun; saved phenotype-bootstrap and gradient-uncertainty summaries were reused.

## Methodological issues and manuscript corrections

- **Neutral tails remain INDETERMINATE.** The 2,500-inner-draw calibration pilot cannot resolve the approximately 8.6 × 10⁻⁶ rank-1 BY threshold. A conservative nuisance-uniform bound yielded no BY discoveries at Ne = 1,000, but does not validate the original test or demonstrate neutrality. Retain the 12 historical hits only as conditional/exploratory candidates.
- **Correct ascertainment and uncertainty wording.** Of 1,089 clusters, 875 have component-trait discovery support and 214 are supported only through the eighth mapped phenotype, the β composite. Historical block resampling was global, not chromosome-stratified. Distinguish posterior means of correlations from correlations of the posterior mean covariance matrix.
- **Preserve interpretive limits.** Covariance is among-RIL, not necessarily additive G; pleiotropy scaling lacks a calibrated genetic-architecture null. Linked regional associations do not identify independent causal alleles. Historical sa935-to-WS283 coordinate identity and movie/RIL sample-flow and experimental provenance still require reconciliation.

Analysis artifacts and detailed evidence are in `/Users/mallard/Documents/ChatGPT/Pleiotropy/manuscript_finalization/`, particularly `FINAL_REPORT.md`, `VALIDATION_SUMMARY.md`, and `reports/74_Manuscript_consistency_audit.md`. This summary adds documentation only; it does not change analyses or establish publication readiness.
