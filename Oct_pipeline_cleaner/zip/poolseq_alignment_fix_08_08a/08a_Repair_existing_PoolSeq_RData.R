## Repair the sample-column relabeling bug in the ORIGINAL revised script 08
## WITHOUT rereading the large sync file. Run ONCE from the project root.
## Original poolSeq::af()/coverage() matrices were sorted by replicate then
## generation, but script 08 labeled columns with UNSORTED sample manifest IDs.
## This script restores the correct mapping in each of the saved sample matrices.

pool_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared.RData'
backup_file <- 'RData_GridLMM_NaCl/08_PoolSeq_prepared_BEFORE_ALIGNMENT_FIX.RData'
report_file <- 'RData_GridLMM_NaCl/08_sample_alignment_repair_report.csv'

if (!file.exists(pool_file)) stop('Missing script-08 prepared Pool-seq data')
p <- new.env(parent = globalenv())
load(pool_file, envir = p)
if (isTRUE(p$sample_alignment_validated))
  stop('This Pool-seq RData is already marked sample-alignment validated; no repair was performed')
if (file.exists(backup_file))
  stop('A pre-repair backup already exists; do not reapply this column permutation')

required <- c('populations','all.af.caendr','all.cov.caendr',
              'all.counts.caendr','all.counts.alt.caendr')
if (!all(required %in% ls(p)))
  stop('Incomplete script-08 RData object set')
pop <- p$populations
if (!all(c('sample_id','generations','replicate') %in% names(pop)))
  stop('Missing manifest sample ID, generation, or replicate')
ids <- as.character(pop$sample_id)
times <- as.numeric(pop$generations)
rep_ids <- as.numeric(as.character(pop$replicate))
if (anyNA(times) || anyNA(rep_ids) || anyDuplicated(ids))
  stop('Invalid sample identifiers, generations, or replicate IDs')

## Recover the actual sample identity of each stored column based on the
## native poolSeq output order, NOT the (incorrect) stored column label.
native_order <- order(rep_ids, times)
true_identity_of_saved_column <- ids[native_order]
if (identical(true_identity_of_saved_column, ids))
  stop('The manifest is already in poolSeq output order; cannot establish that the known bug applies')
if (anyDuplicated(paste(rep_ids, times)))
  stop('Replicate-generation pairs not unique')

matrix_names <- c('all.af.caendr','all.cov.caendr',
                  'all.counts.caendr','all.counts.alt.caendr')
for (nm in matrix_names) {
  mat <- p[[nm]]
  if (!is.matrix(mat) || ncol(mat) != length(ids) ||
      !identical(colnames(mat), ids))
    stop('Unexpected column layout in ', nm, '; no changes saved')
  if (!identical(rownames(mat), rownames(p$all.af.caendr)))
    stop('SNP row order differs in ', nm)
}

## Apply the same bijective column permutation to every matrix.
for (nm in matrix_names) {
  mat <- p[[nm]]
  colnames(mat) <- true_identity_of_saved_column
  p[[nm]] <- mat[, ids, drop = FALSE]
}

## Recompute depth QC with corrected columns.
if (requireNamespace('matrixStats', quietly = TRUE)) {
  p$sample_alignment_validated <- TRUE
  p$sample_alignment_method <- 'reconstructed_poolSeq_native_replicate_generation_order'
  cov <- p$all.cov.caendr
  depth_qc <- data.frame(sample_id = ids, generations = times,
                         replicate = rep_ids,
                         mean_depth = colMeans(cov, na.rm = TRUE),
                         median_depth = matrixStats::colMedians(cov, na.rm = TRUE))
} else {
  stop('matrixStats is required for the coverage QC output')
}

## Basic internal integrity check: minor + complementary counts = coverage.
## Restrict the check to finite entries to accommodate missing observations.
minor <- p$all.counts.caendr
other <- p$all.counts.alt.caendr
cov <- p$all.cov.caendr
finite <- is.finite(minor) & is.finite(other) & is.finite(cov)
if (!all((minor[finite] + other[finite]) == cov[finite]))
  stop('Count/complement/coverage inconsistency; nothing has been saved')

report <- data.frame(
  stored_column_position = seq_along(ids),
  previous_INCORRECT_label = ids,
  true_sample_in_stored_column = true_identity_of_saved_column,
  actual_generation = times[native_order],
  actual_replicate = rep_ids[native_order],
  label_was_wrong = ids != true_identity_of_saved_column,
  stringsAsFactors = FALSE)
cat('Incorrect sample labels in old script-08 RData:', sum(report$label_was_wrong), '\n')
print(report[report$label_was_wrong, , drop = FALSE], row.names = FALSE)

## Save a backup before modifying the project input.
if (!file.copy(pool_file, backup_file, overwrite = FALSE))
  stop('Could not make a backup of the original Pool-seq RData')
temp <- tempfile(pattern = '08_PoolSeq_repaired_', tmpdir = dirname(pool_file),
                 fileext = '.RData')
save(list = ls(p, all.names = TRUE), envir = p,
     file = temp, compress = 'gzip')
if (!file.copy(temp, pool_file, overwrite = TRUE))
  stop('Could not write the repaired Pool-seq data; recover from backup: ', backup_file)
unlink(temp)
write.csv(report, report_file, row.names = FALSE)
write.csv(depth_qc, 'RData_GridLMM_NaCl/08_coverage_QC_REPAIRED_FILTERED_SNPS.csv', row.names = FALSE)
cat('\nRepaired and validated', pool_file, '\n')
cat('Backup:', backup_file, '\n')
cat('Column audit:', report_file, '\n')
cat('Original 08_coverage_QC.csv is not corrected by this fast repair (it was computed on unfiltered SNPs).\n')
cat('IMPORTANT: rerun 09f, 09g, 10, 11, 12/12b before interpreting selection results.\n')
