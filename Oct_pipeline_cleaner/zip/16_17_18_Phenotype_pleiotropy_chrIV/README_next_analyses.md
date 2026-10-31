# Next analyses: phenotypic prediction, pleiotropy versus linkage, and chromosome IV

Run from the same project root used for scripts 00–15. These scripts do not alter existing analyses.

## 16 — Cross-validated SNP prediction of RIL BLUEs

`source("16_Crossvalidated_phenotype_prediction.R")`

- **Primary source:** the complete `genotype_matrix` saved by script 00; it deliberately does NOT select features from GWAS results or LD clusters defined using all phenotypes.
- Five held-out folds of RIL identifiers; within each training fold and for each trait separately, rank SNPs by marginal association; greedily prune markers with training-set genotype r² >= 0.8; refit a ridge model for top K markers, with shrinkage selected by inner validation; predict never-seen RILs.
- Outputs: `phenotype_prediction/CV_RILs/16_Prediction_vs_number_of_SNPs.csv`, plot, folds, SNP ranks, out-of-fold predictions, trait-correlation comparison.
- Computes held-out R² = 1 - SSE/SST. These can be negative. This is more defensible than in-sample phenotype reconstruction from GWAS-selected loci.
- **Caveat:** RIL BLUEs were estimated using the entire phenotype data; folds are random rather than completely independent ancestry groups. This tests prediction between held-out RIL identifiers within the same panel, not transfer to a distinct evolution population. Feature selection is within fold, but upstream BLUE estimation is not refitted for each fold.
- Predicted-trait correlations reflect genetic prediction from a common marker panel, and are not direct estimates of the among-line genetic correlation matrix G.
- Script can be computationally demanding when the full SNP matrix is large: `max_leads`, `K_values`, and `n_outer` can be reduced for a pilot.

## 17 — What the LD clusters can and cannot tell us about pleiotropy

`source("17_Linkage_vs_pleiotropy_diagnostics.R")`

- Within each cluster significant for at least two traits, selects the SNP with the smallest original mixed-model GWAS p-value for each trait, then measures LD between lead variants and performs descriptive two-SNP conditional tests on RIL BLUEs.
- Reports `candidate_distinct_linked_signals` only when leads are distinguishable and each has a conditional association for its corresponding phenotype. Many clusters will be unresolved because high LD makes two possible causal configurations nearly observationally equivalent.
- **Not a formal causal pleiotropy proof.** The conditional p-values are post-selection and must not be reported as confirmatory. A full multi-trait fine-mapping method with locus-wide LD and joint trait likelihoods, more recombinants or perturbation experiments would be required to distinguish shared causal variants from tight linkage more convincingly.

## 18 — Mutation candidates on chromosome IV

`source("18_CHRIV_founder_and_Poolseq_variants.R")`

- Set `vcf_file` to the *indexed* bgzip VCF that contains the 16 founders. Provide `18_founder_sample_ids.txt` containing the 16 exact VCF sample names, one per line. Set `sync_file` if different.
- Automatically queries positions **IV:8892663** and **IV:10499135** from the VCF and reads exact rows from the raw Pool-seq SYNC file; it does not depend on the founder-variant-filtered script-08 matrix.
- Exports VCF founder genotype calls, sample-specific ALT frequencies and coverage, trajectory plots, and a comparison of script-15 prediction errors near candidate coordinates.
- Zero observed ALT reads at G0 is a detection limit, not proof of strictly zero population frequency. Founder 0/0 calls establish absence from *sampled* haplotypes only (provided calls are not missing). Confirm assembly and allele orientation.
- The already computed IV-error proximity check is descriptive and based on very few clusters. It may CONTRADICT the explanation that these mutations drive worse full-index prediction on IV.

## Suggested stopping point

1. Run script 16 and determine whether held-out phenotype R² increases then plateaus as K grows, and whether the matrix of predicted-trait correlations recovers the overall pattern of BLUE correlations.
2. Run script 17; report the fraction of multi-trait clusters that are unresolved and focus any fine-mapping or experimental validation on a few informative loci. Do not equate multiple associations within one LD cluster with proven causal pleiotropy.
3. Run script 18 to document mutation coverage, founder genotypes and trajectories, then review whether chromosome-IV prediction errors truly co-localize.
4. Put results of 16 and 17 into the manuscript only after seeing the outputs, with phenotypic predictions properly labeled as held-out RIL prediction and variant sharing as unresolved without stronger evidence.
