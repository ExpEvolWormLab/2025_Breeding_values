## 00_Prepare_data.R - reviewed version
## Run from the project root. Retains the output objects/filename used by 01 and 09.
## Scientific choices deliberately retained: NaCl only; original RIL exclusions;
## environmental covariates scaled over the filtered full phenotype table;
## seven traits, RIL BLUEs fitted from replicate phenotypes, and selection gradients
## and log absolute fecundity (data/fertility.rda).

library(vcfR)
library(AGHmatrix)
library(lme4)

## ------------------------------------------------------------
## 1. Phenotypes and environmental covariates
## ------------------------------------------------------------
traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F")
vect_P_traits <- traits

phen_file <- "data/all_phenotypes_from_movies_fixed.txt"
stopifnot(file.exists(phen_file))
phen_data <- read.delim(phen_file, stringsAsFactors = FALSE,
                        check.names = FALSE)

required <- c(traits, "env_label", "pop_label", "date_str", "block_id",
              "temperature", "rel_humidity", "logD")
missing_cols <- setdiff(required, names(phen_data))
if (length(missing_cols))
  stop("Missing phenotype columns: ", paste(missing_cols, collapse = ", "))

phen_data$pop_label <- as.character(phen_data$pop_label)
phen_data$block_id <- as.character(phen_data$block_id)
phen_data$date_str <- as.character(phen_data$date_str)

## Retain the original session definitions
phen_data$session_id <- substr(phen_data$date_str, 1, 4)
phen_data$session_id[!is.na(phen_data$session_id) &
                     phen_data$session_id == "2014"] <- "2013"
phen_data$session_id[!is.na(phen_data$block_id) &
                     substr(phen_data$block_id, 1, 1) == "F"] <- "F"
phen_data$session_id[!is.na(phen_data$block_id) &
                     substr(phen_data$block_id, 1, 2) == "SB"] <- "SB"

## Retain the original exclusion criteria
phen_data <- subset(phen_data,
                    !is.na(pop_label) & nzchar(pop_label) &
                    !is.na(block_id) &
                    !(block_id %in% paste0("B", 300:350)) &
                    !(pop_label %in% c("A00", "CX12311L0", "LSJ2L0")))

## Convert explicitly before scaling (as.numeric(factor) is unsafe)
env_covariates <- c("temperature", "rel_humidity", "logD")
for (i in env_covariates) {
  phen_data[[i]] <- suppressWarnings(as.numeric(as.character(phen_data[[i]])))
}
phen_data <- phen_data[complete.cases(phen_data[, env_covariates]), , drop = FALSE]

for (i in env_covariates) {
  if (!is.finite(sd(phen_data[[i]])) || sd(phen_data[[i]]) == 0)
    stop("Cannot standardize environmental covariate: ", i)
  phen_data[[i]] <- as.numeric(scale(phen_data[[i]]))
}

kept_column <- c(traits, env_covariates, "session_id", "date_str", "pop_label")
data_for_BLUE_NaCl <- phen_data[phen_data$env_label == "NaCl",
                                kept_column, drop = FALSE]
if (!nrow(data_for_BLUE_NaCl)) stop("No NaCl phenotypes after filtering")
if (anyNA(data_for_BLUE_NaCl[, c("session_id", "date_str", "pop_label")]))
  stop("Missing session/date/RIL identifiers in NaCl phenotype data")
data_for_BLUE_NaCl$env_label <- "NaCl"
all_BLUEs <- data_for_BLUE_NaCl  # Historical object name: these are replicate observations
all_BLUEs$date_str <- factor(all_BLUEs$date_str)  # Retain original model coding

cat("NaCl observations before genotype matching:", nrow(all_BLUEs), "\n")

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

old_labels <- all_BLUEs$pop_label
ix <- match(old_labels, replaced_iso$Line)
replace <- !is.na(ix)
all_BLUEs$pop_label[replace] <- replaced_iso$Representant_line[ix[replace]]
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
all_BLUEs_genotyped <- subset(all_BLUEs, pop_label %in% genotyped_lines)
if (!nrow(all_BLUEs_genotyped)) stop("No phenotype/genotype line overlap")

## Match explicitly using identifiers, never positional assumptions
phenotyped_lines <- unique(all_BLUEs_genotyped$pop_label)
genotype_info_kept <- genotype_info[, match(phenotyped_lines, genotyped_lines),
                                     drop = FALSE]
genotype_matrix <- t(as.matrix(genotype_info_kept))
storage.mode(genotype_matrix) <- "numeric"
rownames(genotype_matrix) <- phenotyped_lines
colnames(genotype_matrix) <- as.character(snp_positions$ID)
if (any(!is.finite(genotype_matrix[!is.na(genotype_matrix)])))
  stop("Invalid numeric values in full genotype matrix")

cat("Genotyped NaCl observations:", nrow(all_BLUEs_genotyped), "\n")
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
## This handles phased/unphased heterozygotes correctly, unlike 2*first allele.
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
pruned_ids <- as.character(genotype_pruned@fix[, "ID"])
if (anyNA(pruned_ids) || any(pruned_ids %in% c("", ".")) ||
    anyDuplicated(pruned_ids))
  stop("Pruned VCF requires unique nonmissing SNP IDs")
colnames(genotype_pruned_matrix) <- pruned_ids
if (anyNA(genotype_pruned_matrix))
  stop("Missing pruned genotypes: decide on an explicit imputation/filtering rule")

GRM_sub <- AGHmatrix::Gmatrix(genotype_pruned_matrix, method = "VanRaden")
if (is.null(rownames(GRM_sub)) || is.null(colnames(GRM_sub)))
  stop("Gmatrix returned unnamed lines; cannot verify matrix alignment")
if (!all(phenotyped_lines %in% rownames(GRM_sub)) ||
    !all(phenotyped_lines %in% colnames(GRM_sub)))
  stop("Kinship RIL identifiers do not match phenotype/genotype identifiers")
GRM_sub <- GRM_sub[phenotyped_lines, phenotyped_lines, drop = FALSE]
stopifnot(identical(rownames(GRM_sub), rownames(genotype_matrix)),
          identical(colnames(GRM_sub), rownames(genotype_matrix)))

## Preserve this output object for scripts that use NaCl observations
all_BLUEs_genotyped_NaCl <- all_BLUEs_genotyped
is.NA <- colSums(is.na(genotype_matrix))
cat("Unusable full-genotype SNPs (any NA):", sum(is.NA > 0), "\n")

## ------------------------------------------------------------
## 5. Estimate adjusted RIL BLUEs directly from replicate phenotypes
## ------------------------------------------------------------
## Each NaCl phenotype is modeled on its ORIGINAL trait scale.
## pop_label is FIXED (BLUEs), while date and session are random effects.
## Covariates were standardized across the filtered phenotype table in section 1.
## Predict at the reference value 0 (not necessarily the NaCl-only mean)
## and exclude date/session random effects (re.form = NA).
##
## Use all eligible NaCl observations, including non-genotyped RILs, to
## estimate the experimental nuisance effects. The fitness regression below
## uses only RILs that also have a fecundity estimate.

fit_BLUE <- function(trait) {

  cols <- c(trait, "pop_label", "logD", "rel_humidity", "temperature",
            "date_str", "session_id")
  dat <- all_BLUEs[, cols, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (!nrow(dat)) stop("No observations for trait: ", trait)

  dat$y <- as.numeric(as.character(dat[[trait]]))
  if (any(!is.finite(dat$y)))
    stop("Non-finite phenotype values for trait: ", trait)

  dat$pop_label <- factor(dat$pop_label)
  dat$date_str <- factor(dat$date_str)
  dat$session_id <- factor(dat$session_id)
  if (nlevels(dat$pop_label) < 2L ||
      nlevels(dat$date_str) < 2L ||
      nlevels(dat$session_id) < 2L)
    stop("Insufficient levels for the mixed model of trait: ", trait)

  mod <- lme4::lmer(
    y ~ logD + rel_humidity + temperature + pop_label +
      (1 | date_str) + (1 | session_id),
    data = dat, REML = TRUE
  )

  ## Rank deficiency can make some line means non-identifiable.
  if (length(attr(lme4::getME(mod, "X"), "col.dropped")))
    stop("Rank-deficient fixed effects for trait: ", trait)

  line_names <- levels(dat$pop_label)
  newdata <- data.frame(
    pop_label = factor(line_names, levels = line_names),
    logD = 0, rel_humidity = 0, temperature = 0
  )

  blue <- as.numeric(predict(mod, newdata = newdata, re.form = NA))
  if (length(blue) != length(line_names) || any(!is.finite(blue)))
    stop("Invalid adjusted RIL predictions for trait: ", trait)

  cat(trait, ": ", nrow(dat), " measurements, ", length(line_names),
      " RILs", if (lme4::isSingular(mod)) " (singular random effects)" else "",
      "\n", sep = "")

  data.frame(pop_label = line_names, trait = trait, BLUE = blue,
             stringsAsFactors = FALSE)
}

## Compute each trait once; predictions are line means adjusted to the
## same reference environmental conditions.
BLUE_by_trait <- lapply(traits, fit_BLUE)
BLUE_NaCl <- do.call(rbind, BLUE_by_trait)
rownames(BLUE_NaCl) <- NULL
if (anyDuplicated(BLUE_NaCl[, c("pop_label", "trait")]))
  stop("Duplicate RIL/trait combinations in computed BLUEs")

## Write an explicit long-format table for auditability and compatibility.
## Unlike the historical file, trait identity must not rely on row order.
dir.create("RData_BLUEs", recursive = TRUE, showWarnings = FALSE)
write.table(BLUE_NaCl, "RData_BLUEs/BLUES_NaCl.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)

## One row per RIL, one column per phenotype (NA if that RIL lacks a trait).
phen_NaCl <- reshape(BLUE_NaCl, idvar = "pop_label", timevar = "trait",
                     direction = "wide")
blue_names <- paste0("BLUE.", traits)
if (!all(blue_names %in% names(phen_NaCl)))
  stop("One or more traits missing from computed BLUEs")
phen_NaCl <- phen_NaCl[, c("pop_label", blue_names), drop = FALSE]
names(phen_NaCl)[-1] <- traits
if (anyDuplicated(phen_NaCl$pop_label))
  stop("Duplicate RIL after constructing BLUE phenotype table")
write.csv(phen_NaCl, "RData_BLUEs/BLUES_NaCl_wide.csv", row.names = FALSE)

## ------------------------------------------------------------
## 5a. Estimate selection gradients from RIL BLUEs and fecundity
## ------------------------------------------------------------
## Fecundity is log absolute fecundity. Fit relative fitness against the
## seven BLUE-standardized traits, using one observation per RIL.
fitness_file <- "data/fertility.rda"
stopifnot(file.exists(fitness_file))

## Fitness identifiers need the same isotype normalization as phenotype data.
normalize_RIL <- function(x) {
  x <- clean_iso(as.character(x))
  ix <- match(x, replaced_iso$Line)
  x[!is.na(ix)] <- replaced_iso$Representant_line[ix[!is.na(ix)]]
  x
}

fitness_env <- new.env(parent = baseenv())
load(fitness_file, envir = fitness_env)
if (!exists("ecoefs", envir = fitness_env, inherits = FALSE))
  stop("fertility.rda must contain 'ecoefs'")
ecoef <- fitness_env$ecoefs
if (!all(c("line", "env", "fertility") %in% names(ecoef)))
  stop("ecoefs must contain line, env, and fertility")
fertility_NaCl <- ecoef[!is.na(ecoef$env) & ecoef$env == "NaCl",
                       c("line", "fertility"), drop = FALSE]
names(fertility_NaCl)[1] <- "pop_label"
fertility_NaCl$pop_label <- normalize_RIL(fertility_NaCl$pop_label)
fertility_NaCl$fertility <- suppressWarnings(
  as.numeric(as.character(fertility_NaCl$fertility)))
if (anyDuplicated(fertility_NaCl$pop_label))
  stop("More than one NaCl fecundity estimate per RIL; resolve before merging")

## One row per RIL; don't duplicate fitness across phenotype replicates.
phen_NaCl_with_fitness <- merge(phen_NaCl, fertility_NaCl,
                               by = "pop_label", sort = FALSE)
phen_NaCl_with_fitness <- phen_NaCl_with_fitness[
  complete.cases(phen_NaCl_with_fitness[, c(traits, "fertility")]) &
    apply(phen_NaCl_with_fitness[, c(traits, "fertility"), drop = FALSE],
          1, function(row) all(is.finite(as.numeric(row)))),
  , drop = FALSE]
if (nrow(phen_NaCl_with_fitness) <= length(traits) + 2)
  stop("Insufficient complete RILs for seven-trait selection gradients")
cat("RILs used to estimate selection gradients:",
    nrow(phen_NaCl_with_fitness), "\n")

## Stable back-transformation of log absolute fecundity to relative fitness.
lf <- phen_NaCl_with_fitness$fertility
absolute_scaled <- exp(lf - max(lf))
phen_NaCl_with_fitness$rel_fit <- absolute_scaled / mean(absolute_scaled)
stopifnot(isTRUE(all.equal(mean(phen_NaCl_with_fitness$rel_fit), 1,
                           tolerance = 1e-10)))

## Standardize each RIL BLUE across the COMPLETE fitness-matched RILs.
## Preserve the raw BLUE mean/SD for the GWAS-index scale conversion.
BLUE_mean <- vapply(phen_NaCl_with_fitness[, traits, drop = FALSE],
                    mean, numeric(1))
BLUE_sd <- vapply(phen_NaCl_with_fitness[, traits, drop = FALSE],
                  sd, numeric(1))
if (any(!is.finite(BLUE_sd) | BLUE_sd <= 0))
  stop("Cannot standardize BLUEs: zero or invalid among-RIL SD")
for (i in traits)
  phen_NaCl_with_fitness[[i]] <-
    (phen_NaCl_with_fitness[[i]] - BLUE_mean[i]) / BLUE_sd[i]

## Lande-Arnold multivariate directional selection gradients.
model_beta <- lm(reformulate(traits, response = "rel_fit"),
                 data = phen_NaCl_with_fitness)
if (model_beta$rank != length(traits) + 1L)
  stop("Selection-gradient model is rank deficient: correlated traits")
beta_NaCl <- coef(model_beta)[traits]
stopifnot(length(beta_NaCl) == length(traits),
          all(is.finite(beta_NaCl)))

## ------------------------------------------------------------
## 5b. Construct the selection-index phenotype for GWAS
## ------------------------------------------------------------
## The seven GWAS phenotypes were historically standardized across
## REPLICATE observations, not across adjusted RIL BLUEs. Convert the
## regression gradients to that scale; otherwise a one-unit increase
## in z_replicate is mistaken for one SD of among-RIL BLUE variation.
replicate_mean <- vapply(all_BLUEs_genotyped_NaCl[, traits, drop = FALSE],
                         function(x) mean(x, na.rm = TRUE), numeric(1))
replicate_sd <- vapply(all_BLUEs_genotyped_NaCl[, traits, drop = FALSE],
                       function(x) sd(x, na.rm = TRUE), numeric(1))
if (any(!is.finite(replicate_sd) | replicate_sd <= 0))
  stop("Cannot standardize GWAS replicate phenotypes")

## These are the weights to use on replicate-standardized phenotypes.
## Equivalent (up to an additive constant) to the linear predictor
## formed from the original trait values and unstandardized gradients.
beta_for_mapping <- beta_NaCl * replicate_sd[traits] / BLUE_sd[traits]

for (i in traits)
  all_BLUEs_genotyped_NaCl[[i]] <-
    (all_BLUEs_genotyped_NaCl[[i]] - replicate_mean[i]) / replicate_sd[i]
all_BLUEs_genotyped_NaCl$beta <- drop(
  as.matrix(all_BLUEs_genotyped_NaCl[, traits, drop = FALSE]) %*%
    beta_for_mapping[traits]
)

## ------------------------------------------------------------
## 5c. Save gradients, uncertainty, and a selection-gradient figure
## ------------------------------------------------------------
beta_summary <- coef(summary(model_beta))[traits, , drop = FALSE]
beta_CI <- confint(model_beta)[traits, , drop = FALSE]
selection_gradients <- data.frame(
  trait = traits,
  beta = as.numeric(beta_NaCl[traits]),
  SE = as.numeric(beta_summary[, "Std. Error"]),
  lower95 = as.numeric(beta_CI[, 1]),
  upper95 = as.numeric(beta_CI[, 2]),
  p_value = as.numeric(beta_summary[, "Pr(>|t|)"]),
  BLUE_sd = as.numeric(BLUE_sd[traits]),
  replicate_sd = as.numeric(replicate_sd[traits]),
  beta_for_mapping = as.numeric(beta_for_mapping[traits]),
  stringsAsFactors = FALSE
)

dir.create("RData", recursive = TRUE, showWarnings = FALSE)
write.csv(selection_gradients, "RData/selection_gradients_NaCl.csv",
          row.names = FALSE)

pdf("Selection_gradient.pdf", width = 5.5, height = 4.5)
ci_range <- range(c(selection_gradients$lower95,
                    selection_gradients$upper95, 0), finite = TRUE)
plot(seq_along(traits), selection_gradients$beta, type = "n",
     xaxt = "n", xlab = "Trait", ylab = "Directional selection gradient",
     ylim = ci_range, bty = "n", las = 1)
axis(1, at = seq_along(traits), labels = traits)
abline(h = 0, col = "firebrick3", lty = 2)
arrows(seq_along(traits), selection_gradients$lower95,
       seq_along(traits), selection_gradients$upper95,
       code = 3, angle = 90, length = 0.05)
points(seq_along(traits), selection_gradients$beta,
       pch = 21, bg = "grey75")
dev.off()

cat("Selection gradients and 95% CIs written to RData/selection_gradients_NaCl.csv\n")
cat("Selection-gradient plot written to Selection_gradient.pdf\n")

## ------------------------------------------------------------
## 6. Save only objects consumed by downstream scripts
## ------------------------------------------------------------
dir.create("RData", showWarnings = FALSE, recursive = TRUE)
output_file <- "RData/Preparatory_GLMM_NaCl_only_withBeta.RData.RData"
save(all_BLUEs_genotyped_NaCl, genotype_matrix, genotype_pruned_matrix,
     is.NA, GRM_sub, beta_NaCl, beta_for_mapping, selection_gradients,
     model_beta, BLUE_mean, BLUE_sd, replicate_mean, replicate_sd,
     BLUE_NaCl, phen_NaCl,
     file = output_file)
writeLines(c("Preparation script: 00_Prepare_data_recomputed_BLUEs.R",
             paste("NaCl observations:", nrow(all_BLUEs_genotyped_NaCl)),
             paste("RILs:", nrow(genotype_matrix)),
             paste("Full SNPs:", ncol(genotype_matrix)),
             paste("Pruned SNPs:", ncol(genotype_pruned_matrix)),
             paste("Selection gradient RILs:", nrow(phen_NaCl_with_fitness)),
             paste("BLUE generation: seven univariate lmer models, pop_label fixed"),
             "BLUE predictions: logD=rel_humidity=temperature=0, re.form=NA",
             "BLUE input: all eligible NaCl replicate phenotypes, original trait scales",
             paste("Fitness scale: exp(log-fecundity) / mean(exp(log-fecundity))"),
             paste("Composite index: BLUE-standardized gradients converted to replicate scale"),
             paste("Selection beta traits:", paste(names(beta_NaCl), collapse = ", "))),
           con = "RData/00_preparation_manifest.txt")
capture.output(sessionInfo(), file = "RData/00_sessionInfo.txt")
cat("Saved:", output_file, "\n")

## BLUEs and selection gradients are now estimated here and used directly for GWAS.
