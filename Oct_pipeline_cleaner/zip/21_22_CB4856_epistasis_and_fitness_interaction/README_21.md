# Chromosome IV mutations, CB4856 ancestry, epistasis and the beta index

## Scope

Scripts **21a** (optional, tentative CB4856 local founder-ancestry proxy), **21** (targeted association analysis), and **22** (nonlinear phenotypic-selection test). Script 21 fits the **seven original NaCl traits and `beta`**, with a chromosome-IV-excluded kinship matrix. It does not overwrite scripts 19 or 20. No RIL or Pool-seq data are included.

The **19/20 results were not present in the latest conversation upload**; they are required to interpret the prior mutation and survivor-only conclusions. This package lets you extend 19 regardless. Send the resulting `19_*` and `20_*` output CSV files as well as the new `21_*` output CSVs.

## Background ancestry matters

The two candidate variants are **IV:8,892,663 (C>T)** and **IV:10,499,135 (A>G)**. Absence from the 16 founder samples is compatible with mutations arising during 140 generations of lab domestication before the salt-evolution experiment; RIL segregation must be established empirically.

CB4856 **local founder ancestry is NOT the same variable as the new allele**. Because the mutations may occur only on a CB4856-derived background, `mutation × CB4856` can be **strictly unidentifiable**: `M × H = M` if every mutation carrier has `H = 1`. More RILs do not solve perfect collinearity; recombinant or independently occurring backgrounds are required. The two local ancestry variables are different because recombination could occur between the sites.

### Preferred input: verified local founder ancestry

Create `21_CB4856_local_ancestry.csv` in the project root with:

```csv
pop_label,CB4856_8892663,CB4856_10499135
RIL_001,1,1
RIL_002,0,0
RIL_003,1,0
```

Use your documented local-ancestry assignment at each mutation position (0/1 for predominantly inbred RILs, or probability 0..1 if genuinely inferred). Missing/uncertain calls should be left blank. Do not derive the CB4856 origin from presence of either mutation; that would make the analysis circular.

If founder ancestry calls are **not** available, optional script `21a_Infer_CB4856_ancestry_PROXY.R` can produce a *conservative similarity-based proxy* using the original VCF, all 16 founder sample IDs and the RIL genotype matrix. Set `vcf_file`, `cb_founder_name` and the founder list path before running. **This is not a phased founder-ancestry HMM**. Check its confidence margins, sample sizes and ambiguous calls. If the proxy has limited discrimination, obtain/estimate proper founder ancestry rather than forcing 0/1 assignments.

## Run

1. Update, if available, `19_RIL_source_populations.csv` with real, **not guessed**, population origins (`pop_label,source_population`).
2. If founder ancestry is known, save `21_CB4856_local_ancestry.csv`. If not, run `21a_Infer_CB4856_ancestry_PROXY.R` after configuring the paths (or let script 21 run mutation-only).
3. Run `source('21_Mutations_CB4856_background_all8.R')`.

## Statistical analyses

Outcome traits: `SF`, `SB`, `FS`, `FB`, `BS`, `BF`, `area.F`, `beta`. Script 21 prefers `phen_NaCl`, the **one-row-per-RIL adjusted BLUEs** saved by the recomputed script 00; if `beta` is absent from the BLUE table, it reconstructs it from the saved index weights and phenotype standardization. This avoids the previous script-19 incompatibility with alternate phenotype-object names. If line-level BLUEs do not exist, it explicitly warns and aggregates replicate observations as a lower-quality fallback.

Genotypes: SNP ALT dosages divided by two (inbred homozygote contrast); models include optional categorical RIL source-population fixed effects and chromosome-IV-excluded kinship (LOCO).

Within each trait, script 21 fits the following nested **BLUE-level approximate GLS** models with fixed null-fitted variance ratio:

- `mutation_only`: M1, M2 separately and jointly, then M1 × M2.
- `with_CB4856`: local ancestry H1/H2; mutation effects conditional on H1/H2; M1 × M2; and each mutation × ancestry interaction if identifiable.

The file `21_All_eight_trait_effects.csv` contains effect, standard error, approximate p-value, BH within eight traits by term, and minimum ~80%-power detectable absolute effect in **trait SD units**. `21_Model_identifiability_and_power_QC.csv` includes the number of carriers, recombinants, mutation–ancestry correlations, model rank/estimability, and nested model tests. Rank-deficient effects are `NA`, not arbitrary estimates. Counts and origins appear in separate QC tables.

The file `21_Beta_trait_contributions.csv`, when possible, decomposes the **same-model** beta association/interaction into the seven original weighted trait contributions using the exact `beta` weights and shared GLS whitening. This checks whether the SB contribution and area.F contribution oppose or reinforce one another. It **does not prove causal selection on SB interacts with selection on area.F**, which requires independent fitness or manipulation information.

## Optional direct test: does area.F modify selection on SB? (script 22)

This is a **different hypothesis** from genetic epistasis. Script `22_Test_SBxArea_fitness_interaction.R` loads the original `model_beta` saved by the recomputed script 00, and tests:

`relative_fitness ~ all seven standardized BLUE traits + SB^2 + area.F^2 + SB:area.F`

It reports the nested model comparison, conventional and HC3 robust inference, RIL bootstrap intervals for the SB:area.F interaction, and the estimated conditional association between SB and fitness at body-area values -1, 0 and +1 SD. This is an exploratory **correlational-selection** analysis, not evidence of causality, epistasis between SNPs, or direct interaction between separate selection regimes. It requires the recomputed preparation file to contain `model_beta`.

```r
source('22_Test_SBxArea_fitness_interaction.R')
```

Main outputs in `selection_results/phenotypic_selection_interactions/`:
`22_Selection_regression_coefs.csv`, `22_RIL_bootstrap_interaction.csv`, `22_SB_gradient_by_area.csv` and `22_Nested_fitness_models.csv`.

### Interpretation safeguards

- The `beta` index is a **derived linear combination of traits**; evidence for an effect on `beta` is not independent of component-trait tests.
- ``M1 × M2`` is statistical epistasis at the phenotype BLUE level, not necessarily molecular epistasis.
- Mutation effects may be confounded with local ancestry, linked variants or RIL origin; controlling for local haplotype helps only when the required independent genotype variation exists.
- BLUE-level p-values are approximate, and the original replicate-level adjustment and line-specific BLUE uncertainty are not propagated. Strong effects should be corroborated in a relatedness-adjusted replicate-level model if feasible.
- Power calculations are conditional on observed residual variation, allele frequencies and available genotype contrasts; they are not guarantees of detection.

## Upload for interpretation

- `mutation_effects_chrIV/CB4856_all8/21_All_eight_trait_effects.csv`
- `mutation_effects_chrIV/CB4856_all8/21_Model_identifiability_and_power_QC.csv`
- `mutation_effects_chrIV/CB4856_all8/21_Mutation_by_CB4856_cells.csv` (if ancestry provided)
- `mutation_effects_chrIV/CB4856_all8/21_Beta_trait_contributions.csv` (if produced)
- `mutation_effects_chrIV/19_Mutation_effects_and_epistasis.csv`
- `selection_results/survivors_only_Ne1000/20_Matched_cluster_associations.csv`
- `selection_results/survivors_only_Ne1000/20_Matched_LOCO_prediction.csv`
- `selection_results/phenotypic_selection_interactions/22_Selection_regression_coefs.csv`
- `selection_results/phenotypic_selection_interactions/22_RIL_bootstrap_interaction.csv`

## Validation status

Code has been statically inspected against the original saved-object conventions, but **R and the user's private genotype matrices are unavailable here**, so results and runtime behavior cannot yet be certified.
