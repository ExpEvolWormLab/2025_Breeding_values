# Local QTL competition: extended-region analysis

**Source:** user-supplied `deep_dive_extended(1).zip`, including scripts 26–32 outputs. Exact focal variant trajectories are additionally drawn from the preceding user-supplied `deep_dive.zip` raw-sync results for V, I and IV, which were not included in the new archive's ordinary Pool-seq SNP matrix. Analyses below are descriptive and exploratory unless otherwise noted.

## Executive interpretation

The highest-scoring local-pair effects represent contrasts between genomic marker backgrounds. The evidence **does not currently establish genetic epistasis or selection-induced competition between causally independent loci**. R5 beta is a constructed linear combination of seven traits; its inclusion as an eighth coordinate in the original competition scatterplot partly duplicates the other coordinates. Marker effects can differ in direction because of allele coding, LD phase and conditional-model specifications. Evolutionary slope orientation must be matched to the same marker before contrasting it with the conditional beta coefficient.

Across the extended regional epistasis screen, 432 proposed interaction contrasts were recorded, 256 were judged estimable and none survived BH correction (minimum q ≈ 0.272). These tests are screening results, not negative evidence for all moderate interactions.

## Ranked examples

### 1. Chromosome V beta region: best alternative-haplotype example

- The three historical SNPs are V:6,322,080 (ALT absent from all 16 founders), V:6,637,695 (CB4507-private ALT) and V:6,826,421 (CB4507-private ALT).
- In the relatedness-adjusted **joint marker** beta model (n=940), V:6,637,695 has effect +0.1614 (SE 0.0627, p=0.0102, BH across 8 traits for that marker q=0.0409). R5-844 tag has effect -0.0402 (SE 0.0145, p=0.00565, within-marker eight-trait q=0.0452). These do not survive a separate BH adjustment over *all* regional terms (q=0.1513 and 0.1125, respectively).
- Their seven/eight-trait effect vectors are nearly opposite: original scatterplot correlation -0.894, with 8/8 opposite sign coordinates. **This is a heuristic statistic, not a significance test**.
- The two genotype indicators have 590 double-reference, 283 R5-844-only, 67 V:6,637,695-only and **zero double-positive** beta-measured RILs. Thus the two main associations can be contrasted with the double-negative genotype, but their interaction cannot be estimated. The pattern is compatible with alternative local haplotypes rather than two recombining independent causal alleles.
- The three focal V SNPs rise from approximately 7–9% in A00 to approximately 47–50% in GA4 at generation 50, with weaker changes in GA1 and GA2; these are exact read-count trajectories, **not** formal q-valued selection tests for those sites.
- R5 cluster 848 has 0/1 SNPs with provisional drift q<0.01 across all/retained GA, and cluster 854 has 0/2. R5 851 and 853 have no common-frequency tests. CB4507 founder ancestry proxy failed the pre-set coverage criteria. Consequently independent CB4507-background effects are **not established**.
- **Suggested manuscript role:** alternative associated haplotypes with contrasting multi-trait effects, potentially both responding to selection, but no demonstrated interaction.

### 2. Chromosome I beta region: most interpretable founder-allele composition

- I:7,167,738 ALT is carried by PX174 and RC301. I:7,491,656 and I:7,596,359 ALT are RC301-private. The R5 53 shared-background proxy was admitted in n≈558 RILs; RC301-specific R5 54 ancestry proxy did not pass the threshold.
- Original R5 GWAS beta q: 0.0113 for group 53, 0.000357 for group 54. Common-frequency SNPs with provisional selection q<0.01: 7/7 in R5 53 and 2/2 in R5 54 (full/retained GA counts).
- Joint marker additions improve SB as a whole after founder-proxy adjustment (model p=4.9e-5), but **do not improve beta convincingly** (model p=0.151). Individual beta terms have wide confidence intervals. The shared site I:7,167,738 is associated with multiple behaviors, with signs whose weighted fitness contributions partly cancel.
- For I:7,167,738 versus I:7,491,656, the beta-complete genotype cells are 298, 3, 160, 83: one cell has only three RILs. The interaction was flagged unidentifiable. For I:7,491,656 vs R5 54, 448,10,7,79 cells exist, but effects remain highly linked.
- **Suggested manuscript role:** linked founder-private vs shared variation with independent evolutionary observations; conditional effects on beta unresolved.

### 3. Chromosome IV: clearest mutation-associated behavioral signal, not established local competition

- Varying local founder-genotype similarity and two post-founder mutant markers are both present. The CB4856 proxy H730 is associated with SB before adding focal variants (estimate -0.1014, p≈0.000727), but not among mutation noncarriers (estimate -0.0352, p≈0.280). Thus a standalone CB4856 background effect was **not** confirmed.
- IV:10,499,135 has conditional SF +0.3418 and SB +0.3371 with q≈1.1e-5 for each; its beta association is +0.0598 (p≈0.132), not significant. Its beta contribution via SB (+0.087) is offset substantially by SF (-0.064). IV:8,892,663 has no compelling joint beta effect.
- R5 725 and 730 each show strong region-level drift-aware signals (15/28 and 77/90 selected SNPs under all/retained GA).
- Interactions are not significant after correction. IV:10,499,135 vs R5 725 is nearly nested: cells 465,267,4,92, so a reliable interaction is not estimable. IV:8,892,663 vs IV:10,499,135 has all four adequately represented in n=828 beta-complete RILs, but beta interaction p≈0.556.
- **Suggested manuscript role:** multiple behavioral effects on a local founder-derived haplotype, net fitness cancellation and limits of separating linked mutation from background.

### 4. Chromosome III beta clusters 504/506: compelling evolution, weak conditional independence

- Clusters R5 504 and R5 506 have opposite beta-GWAS effect directions and opposite phase-aligned frequency slopes. Their common-frequency selected SNP counts are 19/10 and 32/38 (all/retained).
- In the founder-adjusted regional model, **the R5 506 tag was dropped as collinear with founder proxy H506**; R5 504's beta effect is -0.0418, p=0.104, and H506's beta association is p=0.136.
- A plot comparing H506 with R5 504 is therefore **not** evidence that the two original clusters are independently favored.
- **Suggested manuscript role:** phase-opposed genomic response in overlapping LD regions, requiring phased haplotype-level modeling.

### 5. Chromosome III historical body-area region

- Joint model has area.F p≈2.70e-4 and beta p≈0.00538, but no single conditional tag passes within-term eight-trait q<0.05.
- The original competition plot's R5 590 vs R5 598 pair has effect-profile correlation -0.796 and opposite signs in seven/eight coordinates, but neither individual marker has strong conditional beta evidence.
- R5 598 has one selected SNP q<0.01 in full GA and zero retained, so its selection signal lacks the robustness of chromosome I/IV.

### 6. Chromosome X historical body-area region

- Joint phenotype associations are strong for area.F (p≈3.39e-5) and beta (p≈9.30e-4). The conditional R5 1074 tag has beta effect +0.0528 (p≈0.000268, within-term eight-trait q≈0.00214) and area.F +0.000921 (q≈0.0325).
- Cluster R5 1087 has two SNPs below drift q<0.01 in all GA, but zero retained and a negative cluster slope. Its **conditional** area.F effect is not detectable (p≈0.788). Consequently this is better a **discordant GWAS/evolution region** than demonstrated competition among causal components.

### 7. Other nominated loci: valuable controls

- **X R5 1010**: conditional area.F effect +0.001678 (q≈0.00116), beta effect +0.06710 (q≈0.00116), but phase-aligned cluster frequency slope is negative in both full (-0.000811) and retained (-0.000875) analyses; one selected common SNP in both. This is a useful *predictive failure/discordance* case.
- **II R5 249**: significant SF, SB, FS and BS effects (within-marker q<0.05) but virtually zero beta (+/-0.005), with 25/17 selected common SNPs in full/retained GA. This is a promising *multi-trait cancellation* example even though it does not represent competing QTLs.
- **V R5 886**: 78/54 selected common SNPs in full/retained GA, strong original beta GWAS, but its conditional tag and MY16 proxy are not significant in the reduced founder-classified subset of ≈605 beta-measured RILs. Good robust evolution case, not a convincing founder-separated effect.
- **II historical body-area regions**: body-area joint associations persist, but no common-frequency SNP reaches provisional selection q<0.01 in either analysis; no strong local-competition case.

## What the figures show

`Regional_QTL_multitrait_selection_atlas.pdf` contains 11 pages. Each region has (A) conditional seven-trait effects normalized by the matched trait SD (dots mark within-term q<0.05 over eight outcomes), (B) beta effects with nominal 95% CI, (C) phase-aligned R5 median allele-frequency slopes in all vs retained GA (SNP counts q<0.01 shown alongside cluster names), and (D) model sample size and ancestry/interaction QC.

`Key_case_studies_genotype_support.pdf` highlights chromosomes V, I, IV and III 504/506 and prints observed genotype cells. `Exact_site_evolutionary_trajectories.pdf` shows the ancestral and salt-evolved allele-frequency read proportions for V, I, IV variants, using raw-sync extraction from the *previous* user-supplied deep-dive.

## Critical interpretation safeguards

1. The original competition scatterplot is **a heuristic on estimates**, not an inferential test for competing loci. It includes beta as an eighth coordinate even though beta is derived from the other seven traits. Correlations of effect profiles are particularly unstable when coefficients are weak/noisy.
2. A marker's orientation may be ALT/reference, but a cluster's evolutionary slope is **phase-aligned**. Direct sign comparisons between arbitrary tag genotype effects and phase-aligned cluster slopes require explicit reorientation.
3. Founder proximity assignments are genotype-similarity **proxies**, not confirmed ancestry. Some key proxies (especially CB4507 in chromosome V) were too poorly classified for inclusion.
4. Global and term-wise multiple-testing corrections differ. The chromosome V beta marker p-values are suggestive within their eight-trait tests, but not q<0.05 over *all* regional terms.
5. A high count of selected SNPs within a broad LD cluster does not represent independent targets of selection. Genome-wide extreme tail calibration for the drift-aware test is provisional.
6. No current data demonstrate direct antagonistic selection between two causal alleles, nor validated epistasis. Stronger tests need verified haplotype assignments, their joint phenotype/genotype counts and matched phase-oriented Pool-seq trajectories.
