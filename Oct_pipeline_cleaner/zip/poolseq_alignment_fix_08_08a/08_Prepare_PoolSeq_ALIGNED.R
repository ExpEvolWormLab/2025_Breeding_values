## 08_Prepare_PoolSeq.R
## Refactors the ACTIVE data-preparation section of old 14_Load_PoolSeq.R.
## Creates a reusable, explicitly ordered Pool-seq dataset for selection analysis.
## Does NOT estimate selection coefficients; that is the next script.

library(poolSeq)
library(matrixStats)

sync_file <- "Poolseq_GA_2026/GA_populations_polymorphisms_v2.sync"
bam_file <- "Poolseq_GA_2026/all_bam_list.txt"
founders_file <- "data/caendr_SNPs.txt"
samples_file <- "Poolseq_GA_2026/08_PoolSeq_sample_manifest.csv"
out_dir <- "RData_GridLMM_NaCl"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Original coverage filter thresholds (relative to sample mean coverage).
max_relative_coverage <- 4
min_relative_median <- .15
max_relative_median <- 5

## ------------------------------------------------------------
## 1. Sample metadata and EXACT sync-column order
## ------------------------------------------------------------
## The original script used 13 of 15 BAM-list entries, in fixed order.
## This is preserved but made explicit in a reusable sample manifest.
## IMPORTANT: Check 08_PoolSeq_sample_manifest.csv against sync columns once.
if (file.exists(samples_file)) {
  populations <- read.csv(samples_file, stringsAsFactors = FALSE)
} else {
  paths <- read.table(bam_file, stringsAsFactors = FALSE)[[1]]
  if (length(paths) < 15L)
    stop("Expected at least 15 BAM-list entries for the original 13-pool design")
  a <- vapply(strsplit(basename(paths), "_", fixed = TRUE), `[`,
              character(1), 1)
  a <- sub("\\.bam$", "", a)
  a[1:3] <- c("A00_early", "A00", "A00_late")
  a <- a[c(2, 4:15)]
  populations <- data.frame(
    sample_id = a,
    sync_column = seq.int(4L, length.out = length(a)),
    generations = c(0, 0, 10, 50, 36, 10, 50, 36, 10, 36, 10, 50, 36),
    replicate = c(1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6),
    stringsAsFactors = FALSE
  )
  write.csv(populations, samples_file, row.names = FALSE)
  message("Created ", samples_file, ". Verify the sync-column/sample mapping.")
}

required <- c("sample_id", "sync_column", "generations", "replicate")
if (!all(required %in% names(populations)) ||
    anyDuplicated(populations$sample_id) ||
    !identical(as.integer(populations$sync_column),
               seq.int(4L, length.out = nrow(populations))))
  stop("Sample manifest must have unique IDs and consecutive sync columns 4:... ")

populations$population <- populations$sample_id
populations$replicate_fact <- factor(populations$replicate)

write.table(data.frame(
    sample_id = populations$sample_id,
    sync_file = basename(sync_file),
    sync_column = populations$sync_column,
    generation = populations$generations,
    replicate = populations$replicate),
  file = file.path(out_dir, "08_pools.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE)

## ------------------------------------------------------------
## 2. Read sync with EXPLICIT minor-allele polarization
## ------------------------------------------------------------
## poolSeq::read.sync() defaults to polarization='minor'. State it explicitly
## so that downstream counts are not mistakenly called reference counts.
GA.sync <- poolSeq::read.sync(sync_file,
  gen = populations$generations, repl = populations$replicate,
  polarization = "minor")
all.af <- as.matrix(poolSeq::af(GA.sync,
                     gen = populations$generations, repl = populations$replicate))
all.cov <- as.matrix(poolSeq::coverage(GA.sync,
                     gen = populations$generations, repl = populations$replicate))

if (ncol(all.af) != nrow(populations) ||
    !identical(dim(all.af), dim(all.cov)) ||
    !identical(rownames(all.af), rownames(all.cov)))
  stop("Pool-seq frequency/coverage matrices do not align to sample manifest")

## CRITICAL: poolSeq::af() and coverage() SORT the extracted columns by
## replicate and generation, even when the sync columns / manifest are in a
## different order. The original 14_Load_PoolSeq.R sorted populations to
## match. Here we instead explicitly REORDER both matrices to the verified
## manifest, before assigning sample names. Never blindly relabel columns.
expected_af <- paste0("F", populations$generations,
                      ".R", populations$replicate, ".freq")
expected_cov <- paste0("F", populations$generations,
                       ".R", populations$replicate, ".cov")
if (anyDuplicated(expected_af) || anyDuplicated(expected_cov))
  stop("Duplicate replicate-generation pairs in the manifest")
if (!setequal(colnames(all.af), expected_af) ||
    !setequal(colnames(all.cov), expected_cov))
  stop("poolSeq column labels differ from manifest replicate-generation pairs; inspect native colnames")
frequency_native_columns <- colnames(all.af)
coverage_native_columns <- colnames(all.cov)
idx_af <- match(expected_af, frequency_native_columns)
idx_cov <- match(expected_cov, coverage_native_columns)
if (anyNA(idx_af) || anyNA(idx_cov))
  stop("Failed to map poolSeq columns to named samples")
write.csv(data.frame(
  sample_id = populations$sample_id,
  sync_column = populations$sync_column,
  generations = populations$generations,
  replicate = populations$replicate,
  original_af_column_position = idx_af,
  original_af_label = frequency_native_columns[idx_af],
  original_cov_column_position = idx_cov,
  original_cov_label = coverage_native_columns[idx_cov]),
  file.path(out_dir, "08_verified_sample_column_mapping.csv"),
  row.names = FALSE)
all.af <- all.af[, idx_af, drop = FALSE]
all.cov <- all.cov[, idx_cov, drop = FALSE]
colnames(all.af) <- colnames(all.cov) <- populations$sample_id
sample_alignment_validated <- TRUE
sample_alignment_method <- "poolSeq_native_replicate_generation_column_match"

## The index used by poolSeq is CHROM.POS.
ids <- rownames(all.af)
split_ids <- strsplit(ids, ".", fixed = TRUE)
if (any(lengths(split_ids) != 2L))
  stop("Unexpected sync marker IDs (expected CHROM.POS)")
snps_positions.sync <- data.frame(
  ID = ids,
  CHROM = vapply(split_ids, `[`, character(1), 1),
  POS = as.numeric(vapply(split_ids, `[`, character(1), 2)),
  stringsAsFactors = FALSE
)

## ------------------------------------------------------------
## 3. Coverage filter, retaining the original scientific thresholds
## ------------------------------------------------------------
mean_cov <- colMeans(all.cov, na.rm = TRUE)
if (any(!is.finite(mean_cov) | mean_cov <= 0))
  stop("Invalid average coverage in at least one pool")
rel_cov <- sweep(all.cov, 2, mean_cov, "/")
snps_positions.sync$minCov <- matrixStats::rowMins(rel_cov, na.rm = TRUE)
snps_positions.sync$maxCov <- matrixStats::rowMaxs(rel_cov, na.rm = TRUE)
snps_positions.sync$medianCov <- matrixStats::rowMedians(rel_cov,
                                                         na.rm = TRUE)
kept.alleles <- snps_positions.sync$maxCov <= max_relative_coverage &
  snps_positions.sync$medianCov < max_relative_median &
  snps_positions.sync$medianCov > min_relative_median
kept.alleles[is.na(kept.alleles)] <- FALSE

## ------------------------------------------------------------
## 4. Restrict to founder SNPs (match by CHROM and POS)
## ------------------------------------------------------------
founders <- read.table(founders_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
if (!all(c("CHROM", "POS") %in% names(founders)))
  stop("Founder SNP table needs CHROM and POS columns")
founder_ids <- unique(paste(founders$CHROM, founders$POS, sep = "."))
keep <- kept.alleles & snps_positions.sync$ID %in% founder_ids

all.af.caendr <- all.af[keep, , drop = FALSE]
all.cov.caendr <- all.cov[keep, , drop = FALSE]
snps_positions.sync.caendr <- snps_positions.sync[keep, , drop = FALSE]

## Convert frequency to reads only ONCE; complementary counts sum to depth.
## 'all.counts.caendr' represents the MINOR allele, not the REF allele.
all.counts.caendr <- round(all.af.caendr * all.cov.caendr)
all.counts.alt.caendr <- all.cov.caendr - all.counts.caendr
rownames(snps_positions.sync.caendr) <- NULL

## ------------------------------------------------------------
## 5. Allele metadata (the allele frequency is for the MINOR allele)
## ------------------------------------------------------------
## poolSeq::alleles() can return loci in a different order. Match by ID.
allele_table <- as.data.frame(poolSeq::alleles(GA.sync),
                              stringsAsFactors = FALSE)
if (!all(c("chr", "pos", "ref", "major", "minor") %in%
         names(allele_table)))
  stop("Unexpected poolSeq::alleles() columns; inspect the installed version")
allele_table$ID <- paste(allele_table$chr, allele_table$pos, sep = ".")
if (anyDuplicated(allele_table$ID))
  stop("Ambiguous allele annotations at repeated CHROM.POS")
k <- match(snps_positions.sync.caendr$ID, allele_table$ID)
if (anyNA(k)) stop("Missing allele annotations for filtered SNPs")
allele_metadata <- data.frame(
  ID = snps_positions.sync.caendr$ID,
  CHROM = snps_positions.sync.caendr$CHROM,
  POS = snps_positions.sync.caendr$POS,
  ref.SYNC = as.character(allele_table$ref[k]),
  major.SYNC = as.character(allele_table$major[k]),
  minor.SYNC = as.character(allele_table$minor[k]),
  stringsAsFactors = FALSE
)

## Match the previous column names, which old script 16 expected.
snps_positions.sync.caendr <- merge(snps_positions.sync.caendr,
  allele_metadata[, c("ID", "ref.SYNC", "major.SYNC", "minor.SYNC")],
  by = "ID", sort = FALSE)
snps_positions.sync.caendr <- snps_positions.sync.caendr[
  match(rownames(all.af.caendr), snps_positions.sync.caendr$ID), , drop = FALSE]
rownames(snps_positions.sync.caendr) <- NULL

## Strictly verify correspondence once: a frequency at row i is for its
## named SNP, with the MINOR allele defined in allele_metadata.
if (!identical(rownames(all.af.caendr), snps_positions.sync.caendr$ID) ||
    !identical(rownames(all.af.caendr), rownames(all.counts.caendr)))
  stop("Pool-seq SNP alignment failed")

## ------------------------------------------------------------
## 6. Save the exact starting point for the selection analysis
## ------------------------------------------------------------
coverage_qc <- data.frame(
  sample_id = populations$sample_id,
  generations = populations$generations,
  replicate = populations$replicate,
  mean_depth = mean_cov,
  median_depth = matrixStats::colMedians(all.cov, na.rm = TRUE)
)
write.csv(coverage_qc, file.path(out_dir, "08_coverage_QC.csv"),
          row.names = FALSE)
write.csv(data.frame(
  step = c("read.sync", "coverage filter", "founder SNP intersection"),
  n_snps = c(nrow(all.af), sum(kept.alleles), sum(keep))),
  file.path(out_dir, "08_SNP_filter_QC.csv"), row.names = FALSE)
write.csv(allele_metadata, file.path(out_dir, "08_PoolSeq_alleles.csv"),
          row.names = FALSE)

## We save the filtered counts and alleles, not the entire several-GB sync.
## The GWAS-effect-allele polarization will be done ONCE in the next script.
save(populations, snps_positions.sync.caendr, allele_metadata,
     sample_alignment_validated, sample_alignment_method,
     all.af.caendr, all.cov.caendr,
     all.counts.caendr, all.counts.alt.caendr,
     file = file.path(out_dir, "08_PoolSeq_prepared.RData"),
     compress = "gzip")
capture.output(sessionInfo(), file = file.path(out_dir, "08_sessionInfo.txt"))
cat("Saved", nrow(all.af.caendr), "filtered SNPs for",
    ncol(all.af.caendr), "pools in", out_dir, "\n")
cat("Next: match minor alleles to GWAS effect alleles BEFORE estimating s.\n")
