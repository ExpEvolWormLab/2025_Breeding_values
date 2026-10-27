
## 01_Mapping_final.R - reviewed version
## Run from project root, after 00_Prepare_data_reviewed.R.

library(AGHmatrix)
library(GridLMM)

input_file <- "RData/Preparatory_GLMM_NaCl_only_withBeta.RData"
stopifnot(file.exists(input_file))
prep <- new.env(parent = globalenv())
load(input_file, envir = prep)

required_objects <- c("data_for_GWAS_NaCl_genotyped_NaCl", "genotype_matrix",
                      "genotype_pruned_matrix")
if (!all(required_objects %in% ls(prep)))
  stop("Missing input objects: ",
       paste(setdiff(required_objects, ls(prep)), collapse = ", "))

data_for_GWAS_NaCl_genotyped_NaCl <- prep$data_for_GWAS_NaCl_genotyped_NaCl
genotype_matrix <- prep$genotype_matrix
genotype_pruned_matrix <- prep$genotype_pruned_matrix

vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
chromosomes <- c("I", "II", "III", "IV", "V", "X")
allowed_chromosomes <- c(chromosomes, "MtDNA")
chromosomes_to_run <- chromosomes
traits_to_run <- vect_P_traits
mc_cores <- 10
overwrite_results <- TRUE

## Use names as keys
ril_ids <- rownames(genotype_matrix)
if (is.null(ril_ids) || anyDuplicated(ril_ids) ||
    !all(unique(data_for_GWAS_NaCl_genotyped_NaCl$pop_label) %in% ril_ids))
  stop("Phenotype/genotype RIL IDs are inconsistent")
if (is.null(rownames(genotype_pruned_matrix)) ||
    !setequal(ril_ids, rownames(genotype_pruned_matrix)))
  stop("Full and pruned genotype matrices have different RIL sets")

genotype_pruned_matrix <- genotype_pruned_matrix[ril_ids, , drop = FALSE]

CHR_genotype_matrix <- sub("_.*$", "", colnames(genotype_matrix))
CHR_genotype_matrix_pruned <- sub("_.*$", "", colnames(genotype_pruned_matrix))

output_dir <- "results_mapping"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


for (chrID in chromosomes_to_run) {

  cat("\nMapping chromosome", chrID, "\n")
  genotype_matrix_CHR <- genotype_matrix[,
    CHR_genotype_matrix == chrID, drop = FALSE]
  genotype_pruned_otherCHR <- genotype_pruned_matrix[,
    CHR_genotype_matrix_pruned != chrID, drop = FALSE]

  if (!ncol(genotype_matrix_CHR) || !ncol(genotype_pruned_otherCHR))
    stop("No target or non-target SNPs for chromosome ", chrID)

  ## Only complete SNPs are tested, as in the original mapping script
  is.NA <- colSums(is.na(genotype_matrix_CHR))
  X_CHR <- genotype_matrix_CHR[, is.NA == 0, drop = FALSE]
  if (!ncol(X_CHR)) stop("No complete GWAS markers on chromosome ", chrID)

  GRM_CHR <- AGHmatrix::Gmatrix(genotype_pruned_otherCHR,
                               method = "VanRaden")
  if (is.null(rownames(GRM_CHR)) || is.null(colnames(GRM_CHR)) ||
      !all(ril_ids %in% rownames(GRM_CHR)) ||
      !all(ril_ids %in% colnames(GRM_CHR)))
    stop("LOCO kinship has incorrect/missing RIL identifiers for ", chrID)
  GRM_CHR <- GRM_CHR[ril_ids, ril_ids, drop = FALSE]
  stopifnot(identical(rownames(X_CHR), rownames(GRM_CHR)))

  for (traitID in traits_to_run) {

    cat("  Trait:", traitID, "\n")
    y <- data_for_GWAS_NaCl_genotyped_NaCl[[traitID]]
    if (sum(is.finite(y)) < 3 || !is.finite(sd(y, na.rm = TRUE)) ||
        sd(y, na.rm = TRUE) == 0)
      stop("Invalid phenotype distribution for ", traitID)

    data_for_GWAS_NaCl_genotyped_NaCl$y <- as.numeric(scale(y))

    output_file <- file.path(output_dir,
      paste0("gwas_results_", traitID, "_CHR_", chrID, "_NaCl.txt"))
    if (!overwrite_results && file.exists(output_file))
      stop("Existing GWAS output, set overwrite_results=TRUE: ", output_file)

    gwas.NaCl.CHR <- GridLMM::GridLMM_GWAS(
      formula = y ~ 1 + logD + rel_humidity + temperature +
        (1 | date_str) + (1 | pop_label) + (1 | session_id),
      test_formula = ~1, reduced_formula = ~0,
      data = data_for_GWAS_NaCl_genotyped_NaCl,
      X = X_CHR,
      X_ID = "pop_label",
      relmat = list(pop_label = GRM_CHR),
      mc.cores = mc_cores,
      method = "REML"
    )

    NaCl_results <- gwas.NaCl.CHR$results
    if (is.null(NaCl_results) || !nrow(NaCl_results))
      stop("GWAS returned no results for ", traitID, " on ", chrID)
    write.table(NaCl_results, output_file,
                quote = FALSE, row.names = FALSE)
  }
}

capture.output(sessionInfo(), file = file.path(output_dir, "01_sessionInfo.txt"))
cat("GWAS completed; files in:", output_dir, "\n")
