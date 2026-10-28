## 03_Create_matrix_for_cluster_mapping_CHR.R
## Build cluster genotypes from the significant SNPs assigned in script 02.
## Uses the same genotype_matrix as the SNP GWAS, not the old GEMMA text matrix.
## SNP dosages are oriented using signed LD BEFORE averaging within clusters.
## Run from the project root after 02_Clustering_from_LD_phase_aligned.R.

library(igraph)

thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
threshold_label <- function(x) paste0("R", round(10 * x))
cluster_dir <- "txt_GridLMM_NaCl"
output_dir <- "data_BLUE_GridLMM_NaCl"
ld_file <- "data_BLUE_GridLMM_NaCl/ld_window_100.hap.ld"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## 1. Load SNP genotypes used for the GWAS (RILs in rows, SNPs in columns)
prep <- new.env(parent = globalenv())
load("RData/Preparatory_GLMM_NaCl_only_withBeta.RData", envir = prep)
genotype_matrix <- prep$genotype_matrix
rm(prep)

genotype_matrix <- as.matrix(genotype_matrix)
storage.mode(genotype_matrix) <- "numeric"
if (anyDuplicated(rownames(genotype_matrix)) ||
    anyDuplicated(colnames(genotype_matrix)))
  stop("Duplicate RIL or SNP identifiers in genotype_matrix")

## Genotype coding must be 0/1 or 0/1/2. RIL diploid dosage is typically 0/2.
observed_values <- sort(unique(as.vector(genotype_matrix)))
observed_values <- observed_values[!is.na(observed_values)]
if (!length(observed_values) || any(!observed_values %in% c(0, 1, 2)))
  stop("Inspect genotype_matrix: expected 0/1 or 0/1/2 alternative-allele dosages")

dosage_max <- if (2 %in% observed_values) 2 else 1
cat("Genotype dosage range:", paste(observed_values, collapse = ", "),
    "; opposite allele dosage =", dosage_max, "- dosage\n")

## 2. Load signed LD. The D/genotype-correlation sign agreement was
## validated during script 02; keep the same LD file and SNP universe.
ld_snps <- read.table(ld_file, header = TRUE, check.names = TRUE)
if (!all(c("CHR", "POS1", "POS2", "R.2", "D") %in% names(ld_snps)))
  stop("LD file needs CHR, POS1, POS2, R.2 and signed D")
ld_snps$R.2 <- as.numeric(as.character(ld_snps$R.2))
ld_snps$D <- as.numeric(as.character(ld_snps$D))
ld_snps$ID1 <- paste(ld_snps$CHR, ld_snps$POS1, sep = "_")
ld_snps$ID2 <- paste(ld_snps$CHR, ld_snps$POS2, sep = "_")

## Same strongest-LD spanning-tree rule as script 02. Phase is defined
## relative to the first SNP in genomic order within each LD component.
phase_one_component <- function(ids, ed) {
  if (length(ids) == 1L) return(setNames(1L, ids))

  gr <- igraph::graph_from_data_frame(
    data.frame(from = ed$ID1, to = ed$ID2,
               cost = 1 - ed$R.2, edge_sign = sign(ed$D)),
    directed = FALSE, vertices = data.frame(name = ids))
  if (igraph::components(gr)$no != 1L)
    stop("Missing signed LD edges within component containing ", ids[1])

  tree <- igraph::mst(gr, weights = igraph::E(gr)$cost)
  te <- igraph::as_data_frame(tree, what = "edges")
  adj <- setNames(vector("list", length(ids)), ids)
  for (k in seq_len(nrow(te))) {
    a <- as.character(te$from[k]); b <- as.character(te$to[k])
    s <- as.integer(te$edge_sign[k])
    tmp <- adj[[a]]; tmp[b] <- s; adj[[a]] <- tmp
    tmp <- adj[[b]]; tmp[a] <- s; adj[[b]] <- tmp
  }

  phase <- setNames(rep(NA_integer_, length(ids)), ids)
  phase[ids[1]] <- 1L
  queue <- character(length(ids)); head <- 1L; tail <- 1L
  queue[tail] <- ids[1]
  while (head <= tail) {
    current <- queue[head]; head <- head + 1L
    neighbors <- adj[[current]]
    for (neighbor in names(neighbors)) {
      if (!is.na(phase[neighbor])) next
      phase[neighbor] <- phase[current] * unname(neighbors[neighbor])
      tail <- tail + 1L
      queue[tail] <- neighbor
    }
  }
  if (anyNA(phase)) stop("Unassigned phase in component ", ids[1])
  phase
}

## 3. Construct one cluster-genotype matrix per LD threshold
for (threshold_clust in thresholds_clust) {
  label <- threshold_label(threshold_clust)
  cat("\nConstructing cluster genotypes at R2 >", threshold_clust, "\n")

  f <- file.path(cluster_dir,
                 paste0("all_sign_GWAS_with_clusters_", label, ".txt"))
  all_cluster_merged <- read.table(f, header = TRUE, sep = "\t",
                                    stringsAsFactors = FALSE)
  all_cluster_merged <- all_cluster_merged[
    order(match(all_cluster_merged$CHROM, c("I", "II", "III", "IV", "V", "X")),
          as.numeric(all_cluster_merged$POS)), , drop = FALSE]
  if (anyDuplicated(all_cluster_merged$rs) ||
      !all(all_cluster_merged$rs %in% colnames(genotype_matrix)))
    stop("Missing or duplicated GWAS SNP genotype IDs at ", label)
  if (!all(c("group", "cluster", "CHROM", "POS") %in% names(all_cluster_merged)))
    stop("Unexpected columns in cluster table at ", label)

  ids <- as.character(all_cluster_merged$rs)
  phase <- setNames(rep(1L, length(ids)), ids)
  component <- all_cluster_merged$cluster
  names(component) <- ids

  ## Include all edges of each original LD component, not just edges within
  ## the final sign-split group: sign-split members may lack a direct LD edge.
  ed <- ld_snps[is.finite(ld_snps$R.2) & ld_snps$R.2 > threshold_clust &
                  is.finite(ld_snps$D) & ld_snps$D != 0 &
                  ld_snps$ID1 %in% ids & ld_snps$ID2 %in% ids,
                c("ID1", "ID2", "R.2", "D"), drop = FALSE]
  ed$component <- unname(component[ed$ID1])
  comp2 <- unname(component[ed$ID2])
  ed <- ed[!is.na(ed$component) & ed$component == comp2, , drop = FALSE]
  edge_rows <- split(seq_len(nrow(ed)), as.character(ed$component))
  comp_rows <- split(which(!is.na(component)), as.character(component[!is.na(component)]))

  for (k in names(comp_rows)) {
    snps <- ids[comp_rows[[k]]]
    eidx <- edge_rows[[k]]
    if (is.null(eidx)) stop("Missing LD edges for component ", k)
    phase[snps] <- phase_one_component(snps, ed[eidx, , drop = FALSE])[snps]
  }

  ## Independent validation: compare inferred orientation with every LD edge.
  n_conflict <- sum(phase[ed$ID1] * phase[ed$ID2] != sign(ed$D))
  cat("  Signed LD conflicts:", n_conflict, "of", nrow(ed), "edges\n")
  if (threshold_clust >= 0.4 && n_conflict > 0)
    stop("Phase check differs from script 02 at ", label)

  groups <- sort(unique(as.integer(all_cluster_merged$group)))
  group_rows <- split(seq_len(nrow(all_cluster_merged)),
                      as.character(all_cluster_merged$group))
  ril_ids <- rownames(genotype_matrix)
  geno_matrix_groups <- matrix(NA_real_, nrow = length(groups),
                               ncol = length(ril_ids),
                               dimnames = list(NULL, ril_ids))

  df_clusters <- data.frame(group = groups, CHROM = NA_character_,
                            POS = NA_integer_, clusterID = NA_character_,
                            stringsAsFactors = FALSE)
  n_snp <- integer(length(groups))
  for (j in seq_along(groups)) {
    idx <- group_rows[[as.character(groups[j])]]
    snps <- ids[idx]
    chr <- unique(as.character(all_cluster_merged$CHROM[idx]))
    if (length(chr) != 1L) stop("Cluster spans chromosomes: ", groups[j])

    pos <- floor(mean(as.numeric(all_cluster_merged$POS[idx])))
    clusterID <- paste(groups[j], chr, pos, sep = "_")
    df_clusters$CHROM[j] <- chr
    df_clusters$POS[j] <- pos
    df_clusters$clusterID[j] <- clusterID
    n_snp[j] <- length(snps)

    ## Match SNPs by identifier; do NOT assign groups positionally.
    dosage <- t(genotype_matrix[, snps, drop = FALSE])
    reversed <- phase[snps] == -1L
    if (any(reversed))
      dosage[reversed, ] <- dosage_max - dosage[reversed, , drop = FALSE]

    ## Do not code missing genotypes as 3: this is not a biological dosage.
    avg <- colMeans(dosage, na.rm = TRUE)
    avg[is.nan(avg)] <- NA_real_
    geno_matrix_groups[j, ] <- avg
  }
  rownames(geno_matrix_groups) <- df_clusters$clusterID

  ## Keep legacy filenames + metadata columns required by later scripts.
  write.table(df_clusters,
              file.path(output_dir, paste0("data_clusters_BLUE_", label, ".txt")),
              row.names = FALSE, quote = FALSE)
  write.table(geno_matrix_groups,
              file.path(output_dir, paste0("genotypes_matrix_clusters_", label, ".txt")),
              quote = FALSE)
  write.csv(data.frame(rs = ids,
                       group = all_cluster_merged$group,
                       component = all_cluster_merged$cluster,
                       phase = unname(phase[ids])),
            file.path(output_dir, paste0("SNP_phase_clusters_", label, ".csv")),
            row.names = FALSE)

  missing <- rowSums(is.na(geno_matrix_groups))
  variation <- apply(geno_matrix_groups, 1, sd, na.rm = TRUE)
  qc <- data.frame(clusterID = df_clusters$clusterID, group = groups,
                   CHROM = df_clusters$CHROM, n_SNPs = n_snp,
                   missing_RILs = missing,
                   missing_fraction = missing / ncol(geno_matrix_groups),
                   sd_genotype = variation)
  write.csv(qc, file.path(output_dir, paste0("Cluster_genotypes_QC_", label, ".csv")),
            row.names = FALSE)
  cat("  QTL clusters:", nrow(geno_matrix_groups),
      "; with missing RILs:", sum(missing > 0),
      "; with no dosage variation:", sum(!is.na(variation) & variation == 0), "\n")
  rm(ed, edge_rows, comp_rows, geno_matrix_groups)
  invisible(gc())
}

capture.output(sessionInfo(), file = file.path(output_dir, "03_sessionInfo.txt"))
cat("Finished cluster genotypes in", output_dir, "\n")
