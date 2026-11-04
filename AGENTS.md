# Repository Guidelines

# Research coding instructions

This repository contains R analyses for evolutionary
genetics research.

## Scientific integrity

- Never change a statistical model without explaining
  the scientific consequences.
- Distinguish coding errors from methodological problems.
- Preserve validated analyses and numerical results.
- Never invent results or claim unexecuted code was tested.


## R development

- Preserve existing data structures when possible.
- Make minimal, targeted changes.
- Add clear error messages and input validation.
- Document new dependencies.
- Use fixed random seeds for reproducible tests.
- Prefer synthetic or reduced datasets for initial tests.

## Verification

- Check syntax after modifying a script.
- Test relevant functions independently when possible.
- Compare new and previous outputs.
- Report which tests passed and which could not run.
- Do not run expensive MCMC or GWAS jobs without approval.

## Git

- Do not commit or push without explicit approval.
- Do not overwrite unrelated user changes.
- Show the diff and summarize all modified files.


## Project Structure & Module Organization

This R research repository studies GWAS, pleiotropy, and Pool-seq evolution. `Oct_pipeline_cleaner/` contains the current numbered pipeline and helper scripts. Preparation and mapping (00–05) feed covariance analysis (07a–07b), Pool-seq analysis (08 onward), and unified R5 analysis (43–48). Manuscript extension directories contain scripts 59–68.

`data/`, `plink_pruning/`, and `Poolseq_GA_2026/` supply inputs. Generated objects, tables, and figures live in `RData*/`, `selection_results/`, `pleiotropy_results/`, and `phenotype_prediction/`. Many inputs and outputs are ignored by Git. Treat `Rcodes_BLUEs/`, `Rcode_clean*/`, and older `GridLMM*` implementations as historical unless explicitly required.

## Development & Validation Commands

Run commands from the repository root in a clean R session:

```sh
git status --short
Rscript --vanilla Oct_pipeline_cleaner/42_Synthetic_smoke_tests.R
Rscript --vanilla Oct_pipeline_cleaner/00_Prepare_data_reviewed.R
```

These inspect changes, run lightweight synthetic checks, and regenerate preparation outputs, respectively. Preparation requires local datasets and writes files. There is no centralized build system or dependency lockfile. Dependencies include GridLMM, AGHmatrix, lme4, vcfR, igraph, poolSeq, MCMCglmm, coda, rrBLUP, and PLINK.

Script `48_Run_unified_main_analysis.R` orchestrates stages 43–46 after upstream inputs exist. It can launch expensive simulations; do not run it as a routine check. Verify trait and LD-threshold settings first: current mapping defaults do not regenerate the complete R5 workflow.

## Coding Style & Naming Conventions

Use two-space indentation and explicit package namespaces where helpful. Preserve numbered descriptive filenames, existing output contracts, and `R5` threshold labels. Match samples, traits, and alleles by identifiers rather than row positions. Check duplicate keys, missing values, dimensions, and REF/ALT orientation before modeling. No repository-wide formatter or linter is configured.

## Testing Guidelines

Validation uses standalone R scripts and `stopifnot()` assertions; no formal coverage target is defined. Extend synthetic checks for changed formulas and alignment logic. Record seeds, parameters, input provenance, and `sessionInfo()`. Synthetic success does not establish biological calibration or extreme-tail validity.

## Commit & Pull Request Guidelines

Recent commits use short messages such as `update MCMCglmm`; no strict convention is evident. Prefer specific imperative summaries. PRs should describe affected stages, scientific assumptions, output changes, validation performed, and checks omitted. Link relevant issues and include figures when visual outputs change. Review `.gitignore` before adding datasets or generated artifacts.
