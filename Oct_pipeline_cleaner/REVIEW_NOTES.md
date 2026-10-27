# First-pass refactor: scripts 00 and 01

The scripts preserve the NaCl-only preparation, the seven phenotypes, the composite selection phenotype, the GridLMM mixed-model formula, leave-one-chromosome-out kinship, and all downstream output names. They are **not yet validated on your private datasets**.

## 00_Prepare_data_reviewed.R

- Removed unused exploratory imports, the `population`/`population2` variables, print-only consistency checks, manually built `bimbam` table, and the positional per-observation isotype replacement loop.
- Applied existing population, block, and NaCl filters, environmental covariate scaling, and representative-isotype assignments in the same order as the legacy script.
- Joined genotyped RILs using explicit sample identifiers, not ordering assumptions.
- Reads *all* sample columns in the pruned VCF rather than hard-coded sample columns `2:1002`.
- Parses alternate-allele dosage from the complete GT field; supports phased and unphased diploids, homozygous haploid-coded calls; rejects unsupported alleles/missing values instead of silently making genotypes up. In particular, original `2 * first GT character` mishandled heterozygotes.
- Reorders kinship rows **and columns** with explicit matching to the full genotype matrix, never simply overwriting dimnames.
- Loads the existing selection-gradient coefficients from `GridLMM_results_mixed/beta_NaCl.txt`, checks trait order/length, and constructs the same weighted composite from seven standardized phenotypes. With a legacy unnamed `x` column, exact trait order remains an assumption. Prefer a file with `trait` and `beta` columns.
- Deliberately *removed* `Selection_gradient.pdf` generation: the old 7-rows-per-RIL BLUE reshaping was position-dependent; the locally fitted coefficients were not used to construct the composite phenotype. Verify which method actually produced `beta_NaCl.txt` before reintroducing this plot.
- Saves `RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData`, exactly the file expected by 01 and 09; adds basic manifest and `sessionInfo()`.

## 01_Mapping_final_reviewed.R

- Retains the original eight tested traits (seven ecological traits plus the derived `beta`), model formula, REML setting, and GWAS output filenames.
- Switches nested loops: constructs LOCO kinship **once per chromosome** (six times rather than 6 × 8 = 48 times), then maps all traits.
- Matches RILs by names across full and pruned matrices, verifies kinship rows and columns and SNP chromosome naming, and excludes markers with any missing genotype as before.
- Preserves the repeated standardization of the phenotype at mapping time. This includes standardizing `beta`, just as in the original script.
- Preserves results at `GridLMM_results_NaCl/gwas_results_<trait>_CHR_<chr>_NaCl.txt`. Other scripts do not need output-path modifications.

## Run and validation

1. Run from the project root where `data/`, `RData/`, `plink_pruning/`, and `GridLMM_results_mixed/` are present.
2. Save the old generated `.RData.RData` and GWAS output directory if you want to compare results before overwriting.
3. Run `source('00_Prepare_data_reviewed.R')` and inspect reported numbers of observations, RILs, SNPs and errors. Compare sample genotypes, VCF sample labels, and `GRM_sub` to your previous objects using explicit sample names.
4. Run `source('01_Mapping_final_reviewed.R')` (preferably start with `chromosomes_to_run <- "I"` and `traits_to_run <- "SF"`) and compare marker-by-marker p-values and effect estimates with the old run.
5. **Stop and investigate discrepancies**, especially if the VCF contains heterozygotes, sample-order differences, missing genotypes, or the old `beta` coefficients were mismatched by trait order.
6. Install package dependencies in the same R environment as the original project (`vcfR`, `AGHmatrix`, `GridLMM`). Use the new sessionInfo outputs to record version differences. No random sampling is introduced here, so setting a random seed is not part of the scripts.

## Outstanding decisions

- The selection-gradient file is an upstream input, but its provenance is currently external to script 00. Confirm how it is estimated and that its units/order match the seven standardized traits. Script `00b_Compute_betas.R` looks exploratory and should be reviewed separately.
- Original code permits any missingness in the full genotype matrix and drops incomplete GWAS SNPs at mapping; the revised script retains that behavior. Missing values in the pruned GRM input now trigger an error and require an explicit policy rather than implicit treatment by `Gmatrix`.
- The RIL phenotype replication structure and potential in-sample use of selection-gradient weights are scientific questions, not automatically changed by this programming refactor.
