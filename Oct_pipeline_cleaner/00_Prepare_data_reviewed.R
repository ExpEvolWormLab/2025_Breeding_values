library(vcfR)
library(AGHmatrix)

## ------------------------------------------------------------
## 1. Phenotypes and environmental covariates
## ------------------------------------------------------------
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")
vect_P_traits <- traits

phen_file <- "data/all_phenotypes_from_movies_fixed.txt"

phen_data <- read.delim(phen_file, stringsAsFactors = FALSE,
                        check.names = FALSE)

phen_data$pop_label <- as.character(phen_data$pop_label)
phen_data$block_id <- as.character(phen_data$block_id)
phen_data$date_str <- as.character(phen_data$date_str)

## Define the session info
phen_data$session_id <- substr(phen_data$date_str, 1, 4)
phen_data$session_id[!is.na(phen_data$session_id) &
                     phen_data$session_id == "2014"] <- "2013"
phen_data$session_id[!is.na(phen_data$block_id) &
                     substr(phen_data$block_id, 1, 1) == "F"] <- "F"
phen_data$session_id[!is.na(phen_data$block_id) &
                     substr(phen_data$block_id, 1, 2) == "SB"] <- "SB"

## Exclude populations and non-CeMEE lines from phenotypic data
phen_data <- subset(phen_data,
                    !is.na(pop_label) & nzchar(pop_label) &
                    !is.na(block_id) &
                    !(block_id %in% paste0("B", 300:350)) &
                    !(pop_label %in% c("A00", "CX12311L0", "LSJ2L0")))

## ## Convert and scale
env_covariates <- c("temperature", "rel_humidity", "logD")
for (i in env_covariates) {
  phen_data[[i]] <- as.numeric(as.character(phen_data[[i]]))
}
phen_data <- phen_data[complete.cases(phen_data[, env_covariates]), , drop = FALSE]

for (i in env_covariates) {
  phen_data[[i]] <- as.numeric(scale(phen_data[[i]]))
}

kept_column <- c(traits, env_covariates, "session_id", "date_str", "pop_label")
data_for_GWAS_NaCl <- phen_data[phen_data$env_label == "NaCl",
                                kept_column, drop = FALSE]
if (!nrow(data_for_GWAS_NaCl)) stop("No NaCl phenotypes after filtering")
if (anyNA(data_for_GWAS_NaCl[, c("session_id", "date_str", "pop_label")]))
  stop("Missing session/date/RIL identifiers in NaCl phenotype data")

data_for_GWAS_NaCl$env_label <- "NaCl"
data_for_GWAS_NaCl$date_str <- factor(data_for_GWAS_NaCl$date_str)

cat("NaCl observations before genotype matching:", nrow(data_for_GWAS_NaCl), "\n")

## ------------------------------------------------------------
## 2. Map duplicate/isotype line names to representative lines
## ------------------------------------------------------------
removed_file <- "data/removed_isotypes_Cemee_2025.txt"
corr_file <- "data/correspondance_isotypes.csv"
stopifnot(file.exists(removed_file), file.exists(corr_file))

removed_isotypes <- read.table(removed_file, stringsAsFactors = FALSE)[[1]]
corr_isotypes <- read.csv(corr_file, stringsAsFactors = FALSE)

if (!all(c("Line", "Representant_line") %in% names(corr_isotypes)))
  stop("Isotype correspondence table needs Line and Representant_line")

clean_iso <- function(x) gsub("_sorted", "", gsub("CeMee", "", x, fixed = TRUE),
                              fixed = TRUE)
removed_isotypes <- unique(clean_iso(removed_isotypes))
corr_isotypes$Line <- clean_iso(corr_isotypes$Line)
corr_isotypes$Representant_line <- clean_iso(corr_isotypes$Representant_line)
corr_isotypes <- unique(corr_isotypes[, c("Line", "Representant_line")])
replaced_iso <- subset(corr_isotypes,
                       Line %in% removed_isotypes & Line != Representant_line)

## Ambiguous mappings must not be settled by row order
if (anyDuplicated(replaced_iso$Line))
  stop("Some removed isotypes map to multiple representatives")

old_labels <- data_for_GWAS_NaCl$pop_label
ix <- match(old_labels, replaced_iso$Line)
replace <- !is.na(ix)
data_for_GWAS_NaCl$pop_label[replace] <- replaced_iso$Representant_line[ix[replace]]
cat("Isotype labels replaced:", sum(replace), "observations\n")

## ------------------------------------------------------------
## 3. Full genotypes: keep the same SNP and line identifiers
## ------------------------------------------------------------
geno_file <- "RData/Genotype_info_Cemee.RData"
stopifnot(file.exists(geno_file))
geno_env <- new.env(parent = globalenv())
load(geno_file, envir = geno_env)
if (!all(c("genotype_info", "snp_positions") %in% ls(geno_env)))
  stop("Genotype RData must contain genotype_info and snp_positions")

genotype_info <- geno_env$genotype_info
snp_positions <- as.data.frame(geno_env$snp_positions)
stopifnot(all(c("CHROM", "POS", "ID") %in% names(snp_positions)))
if (nrow(genotype_info) != nrow(snp_positions))
  stop("Genotype SNP rows and SNP metadata do not match")

snp_positions$CHROM <- factor(as.character(snp_positions$CHROM),
                              levels = c("I", "II", "III", "IV", "V", "X", "MtDNA"))
snp_positions$POS <- as.numeric(as.character(snp_positions$POS))
if (anyNA(snp_positions$CHROM) || anyNA(snp_positions$POS) ||
    anyNA(snp_positions$ID) || anyDuplicated(snp_positions$ID))
  stop("Missing/duplicate SNP identifiers or invalid chromosome/position")

ord <- order(snp_positions$CHROM, snp_positions$POS)
genotype_info <- genotype_info[ord, , drop = FALSE]
snp_positions <- snp_positions[ord, , drop = FALSE]

if (is.null(colnames(genotype_info)) || anyDuplicated(colnames(genotype_info)))
  stop("Full genotype matrix needs unique RIL column names")

genotyped_lines <- colnames(genotype_info)
data_for_GWAS_NaCl_genotyped <- subset(data_for_GWAS_NaCl, pop_label %in% genotyped_lines)
if (!nrow(data_for_GWAS_NaCl_genotyped)) stop("No phenotype/genotype line overlap")

## Match explicitly using identifiers, never positional assumptions
phenotyped_lines <- unique(data_for_GWAS_NaCl_genotyped$pop_label)
genotype_info_kept <- genotype_info[, match(phenotyped_lines, genotyped_lines),
                                     drop = FALSE]
genotype_matrix <- t(as.matrix(genotype_info_kept))
storage.mode(genotype_matrix) <- "numeric"
rownames(genotype_matrix) <- phenotyped_lines
colnames(genotype_matrix) <- as.character(snp_positions$ID)
if (any(!is.finite(genotype_matrix[!is.na(genotype_matrix)])))
  stop("Invalid numeric values in full genotype matrix")

cat("Genotyped NaCl observations:", nrow(data_for_GWAS_NaCl_genotyped), "\n")
cat("Genotyped RILs:", nrow(genotype_matrix), "\n")

## ------------------------------------------------------------
## 4. Pruned genotypes and kinship: explicit RIL alignment
## ------------------------------------------------------------
pruned_file <- "plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz"
stopifnot(file.exists(pruned_file))
genotype_pruned <- vcfR::read.vcfR(pruned_file, verbose = FALSE)

## vcfR::extract.gt returns SNP rows x sample columns
GT <- vcfR::extract.gt(genotype_pruned, element = "GT", as.numeric = FALSE)
if (is.null(dim(GT)) || ncol(GT) == 0)
  stop("No genotype calls found in pruned VCF")

## Original VCF samples have a suffix after an underscore
strainID <- sub("_.*$", "", colnames(GT))
if (anyDuplicated(strainID))
  stop("VCF sample IDs collide after removing suffixes; check sample naming")
colnames(GT) <- strainID

## Convert GT to alternative-allele count (0/1/2).
## This handles phased/unphased heterozygotes correctly
gt_string <- sub(":.*$", "", as.vector(GT))
dosage <- rep(NA_real_, length(gt_string))
dosage[gt_string %in% c("0|0", "0/0", "0")] <- 0
dosage[gt_string %in% c("0|1", "1|0", "0/1", "1/0")] <- 1
dosage[gt_string %in% c("1|1", "1/1", "1")] <- 2
unknown <- !(gt_string %in% c("0|0", "0/0", "0", "0|1", "1|0",
                              "0/1", "1/0", "1|1", "1/1", "1", ".",
                              "./.", ".|.")) & !is.na(gt_string)
if (any(unknown))
  stop("Unsupported/multiallelic GT values in pruned VCF: ",
       paste(head(unique(gt_string[unknown])), collapse = ", "))

GT_dosage <- matrix(dosage, nrow = nrow(GT), ncol = ncol(GT),
                    dimnames = dimnames(GT))

if (!all(phenotyped_lines %in% colnames(GT_dosage)))
  stop("Some phenotyped RILs are missing from the pruned VCF: ",
       paste(head(setdiff(phenotyped_lines, colnames(GT_dosage))), collapse = ", "))

genotype_pruned_matrix <- t(GT_dosage[, phenotyped_lines, drop = FALSE])

GRM_sub <- AGHmatrix::Gmatrix(genotype_pruned_matrix, method = "VanRaden")
GRM_sub <- GRM_sub[phenotyped_lines, phenotyped_lines, drop = FALSE]
stopifnot(identical(rownames(GRM_sub), rownames(genotype_matrix)),
          identical(colnames(GRM_sub), rownames(genotype_matrix)))

data_for_GWAS_NaCl_genotyped_NaCl <- data_for_GWAS_NaCl_genotyped
is.NA <- colSums(is.na(genotype_matrix))
cat("Unusable full-genotype SNPs (any NA):", sum(is.NA > 0), "\n")

## ------------------------------------------------------------
## 5. Composite selection phenotype (beta): existing weights
## ------------------------------------------------------------

beta_file <- "GridLMM_results_mixed/beta_NaCl.txt"
beta_tab <- read.table(beta_file, header = TRUE, stringsAsFactors = FALSE)

if (all(c("trait", "beta") %in% names(beta_tab))) {
  if (anyDuplicated(beta_tab$trait) ||
      !setequal(as.character(beta_tab$trait), traits))
    stop("Selection-gradient table must name the seven traits exactly once")
  beta_NaCl <- setNames(as.numeric(beta_tab$beta), as.character(beta_tab$trait))[traits]
} else if ("x" %in% names(beta_tab)) {
  ## Backward compatible with the original unnamed 'x' column
  if (nrow(beta_tab) != length(traits))
    stop("Expected exactly seven selection gradients in beta_NaCl.txt")
  beta_NaCl <- setNames(as.numeric(beta_tab$x), traits)
  message("Selection gradients loaded from 'x' in trait order: ",
          paste(traits, collapse = ", "))
} else {
  stop("Selection gradient file needs trait/beta columns or legacy x column")
}

for (i in traits) {
  if (sum(is.finite(data_for_GWAS_NaCl_genotyped_NaCl[[i]])) < 2)
    stop("Insufficient phenotypes for: ", i)
  data_for_GWAS_NaCl_genotyped_NaCl[[i]] <- as.numeric(scale(data_for_GWAS_NaCl_genotyped_NaCl[[i]]))
}

data_for_GWAS_NaCl_genotyped_NaCl$beta <- as.numeric(
  as.matrix(data_for_GWAS_NaCl_genotyped_NaCl[, traits, drop = FALSE]) %*% beta_NaCl
)

## ------------------------------------------------------------
## 6. Save only objects consumed by downstream scripts
## ------------------------------------------------------------

dir.create("RData", showWarnings = FALSE, recursive = TRUE)
output_file <- "RData/Preparatory_GLMM_NaCl_only_withBeta.RData"
save(data_for_GWAS_NaCl_genotyped_NaCl, genotype_matrix, genotype_pruned_matrix,
     is.NA, GRM_sub, file = output_file)
writeLines(c("Preparation script: 00_Prepare_data_reviewed.R",
             paste("NaCl observations:", nrow(data_for_GWAS_NaCl_genotyped_NaCl)),
             paste("RILs:", nrow(genotype_matrix)),
             paste("Full SNPs:", ncol(genotype_matrix)),
             paste("Pruned SNPs:", ncol(genotype_pruned_matrix)),
             paste("Selection beta traits:", paste(names(beta_NaCl), collapse = ", "))),
           con = "RData/00_preparation_manifest.txt")
capture.output(sessionInfo(), file = "RData/00_sessionInfo.txt")
cat("Saved:", output_file, "\n")

