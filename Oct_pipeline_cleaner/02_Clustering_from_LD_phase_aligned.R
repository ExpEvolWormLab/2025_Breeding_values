## 02_Clustering_from_LD.R
## SNP-level GWAS results -> LD components -> sign-compatible QTL groups.
## Run after 01_Mapping_final_reviewed.R, from the project root.
## Same GWAS significance and LD thresholds; phase-align SNP effect signs
## using signed LD (D) only during the sign-compatibility split.

library(igraph)

vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")
chromosomes <- c("I", "II", "III", "IV", "V", "X")
env_label <- "NaCl"
thresholds_clust <- c(0.2, 0.4, 0.5, 0.6, 0.8, 0.9, 0.99)
threshold_label <- function(x) paste0("R", round(10 * x))

mapping_dir <- "results_mapping"       # Output directory of script 01
all_gwas_dir <- "GridLMM_results_NaCl"
clusters_dir <- "txt_GridLMM_NaCl"
plot_dir <- "plots_GridLMM_NaCl"
ld_dir <- "data_BLUE_GridLMM_NaCl"
ld_prefix <- file.path(ld_dir, "ld_window_100")
ld_file <- paste0(ld_prefix, ".hap.ld")
positions_file <- file.path(ld_dir, "sign_positions.txt")
ld_ids_file <- file.path(ld_dir, "ld_window_100_SNP_IDs.txt")

## Set TRUE when mapping results or significant SNP set have changed.
## You should regenerate LD after the newly launched GWAS.
recompute_ld <- TRUE
vcf_file <- "/Users/mallard/Nextcloud/CeMEE_vcf/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz"
make_manhattan_plots <- TRUE
split_cores <- 6L  # Use >1 on macOS/Linux; Windows uses sequential lapply
phase_align_effects <- TRUE  # FALSE reproduces original sign splitting
phase_check_pairs <- 500L   # check D sign versus actual GridLMM genotype dosages

for (d in c(all_gwas_dir, clusters_dir, plot_dir, ld_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------
## 1. Gather SNP-level GWAS results from six chromosomes
## ------------------------------------------------------------

all_GWAS <- NULL

for (i in seq_along(vect_P_traits)) {
  traitJ <- vect_P_traits[i]
  cat("Loading GWAS:", traitJ, "\n")

  all_GWAS.temp <- do.call(rbind, lapply(chromosomes, function(chr) {
    f <- file.path(mapping_dir,
                   paste0("gwas_results_", traitJ, "_CHR_", chr, "_", env_label, ".txt"))
    tab <- read.table(f, header = TRUE, check.names = FALSE)
    tab <- tab[, c("X_ID", "beta.5", "p_value_REML"), drop = FALSE]
    names(tab) <- c("rs", paste0("beta_T", i), paste0("pvalue_T", i))
    tab
  }))
  rownames(all_GWAS.temp) <- NULL

  if (anyDuplicated(all_GWAS.temp$rs))
    stop("Duplicated SNP IDs in GWAS results for trait ", traitJ)

  ## Keep the common SNP universe if a trait contains fewer tested SNPs.
  if (is.null(all_GWAS)) {
    all_GWAS <- all_GWAS.temp
  } else {
    n_before <- nrow(all_GWAS)
    all_GWAS <- merge(all_GWAS, all_GWAS.temp, by = "rs")
    if (nrow(all_GWAS) < n_before)
      warning("GWAS SNP intersection shrank when merging trait ", traitJ)
  }
}

## Parse the usual CHR_POS SNP names produced in script 00.
split_ID <- strsplit(as.character(all_GWAS$rs), "_", fixed = TRUE)
all_GWAS$chr <- factor(vapply(split_ID, `[`, character(1), 1), levels = chromosomes)
all_GWAS$ps <- as.numeric(vapply(split_ID, `[`, character(1), 2))
if (anyNA(all_GWAS$chr) || anyNA(all_GWAS$ps))
  stop("Cannot interpret GWAS SNP IDs as CHR_POS")
all_GWAS <- all_GWAS[order(all_GWAS$chr, all_GWAS$ps), , drop = FALSE]
rownames(all_GWAS) <- NULL
all_GWAS$CHROM <- all_GWAS$chr

## SNP-wise Benjamini-Hochberg adjustment across all chromosomes for each trait.
for (i in seq_along(vect_P_traits)) {
  pv <- all_GWAS[[paste0("pvalue_T", i)]]
  qv <- p.adjust(pv, method = "BH")
  all_GWAS[[paste0("is_sign_T", i)]] <- !is.na(qv) & qv < 0.2

  if (make_manhattan_plots) {
    plot_data <- data.frame(
      CHR = as.integer(all_GWAS$chr),
      BP = all_GWAS$ps,
      SNP = as.character(all_GWAS$rs),
      P = -log10(pmax(qv, .Machine$double.xmin))
    )
    pdf(file.path(plot_dir, paste0("Manhattan_GWAS_Trait_", vect_P_traits[i],
                                   "_", env_label, ".pdf")), width = 10)
    qqman::manhattan(plot_data, chr = "CHR", bp = "BP", snp = "SNP", p = "P",
                     logp = FALSE, cex = 0.6,
                     highlight = plot_data$SNP[all_GWAS[[paste0("is_sign_T", i)]]],
                     suggestiveline = -log10(0.2), genomewideline = -log10(0.05))
    dev.off()
  }
  cat("  BH q < 0.2:", sum(all_GWAS[[paste0("is_sign_T", i)]]), "SNPs\n")
}

sig_cols <- paste0("is_sign_T", seq_along(vect_P_traits))
all_GWAS$sign_overall <- rowSums(all_GWAS[, sig_cols, drop = FALSE]) >= 1
all_GWAS_sign <- all_GWAS[all_GWAS$sign_overall, , drop = FALSE]
cat("Significant SNPs in at least one trait:", nrow(all_GWAS_sign), "\n")

## Object used later by 13_Characterize_pleiotropy_R5_NaCl_v2.R
save(all_GWAS, file = file.path(all_gwas_dir, "All_Gwas.RData"))

## ------------------------------------------------------------
## 2. Calculate/load SNP-level LD once (NOT once per R2 cutoff)
## ------------------------------------------------------------

sign_positions <- unique(data.frame(CHR = as.character(all_GWAS_sign$chr),
                                    POS = all_GWAS_sign$ps))
write.table(sign_positions, positions_file, sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)

## The old LD file cannot be assumed current after a complete new GWAS.
if (recompute_ld) {
  if (!file.exists(vcf_file)) stop("Set vcf_file to your filtered VCF: ", vcf_file)
  sign_vcf <- file.path(ld_dir, "Cemee2025.sign.GWAS_BLUE.filt.vcf.gz")

  status <- system2("bcftools", c("view", "-T", shQuote(positions_file),
                                 "-Oz", "-o", shQuote(sign_vcf), shQuote(vcf_file)))
  if (status != 0) stop("bcftools extraction failed")

  ## Verify that selected VCF records cover the GWAS SNP universe.
  extracted <- suppressWarnings(system2(
    "bcftools", c("query", "-f", shQuote("%CHROM\\t%POS\\n"), shQuote(sign_vcf)),
    stdout = TRUE))
  positions <- strsplit(extracted, "\t", fixed = TRUE)
  extracted_ids <- vapply(positions, function(z) paste(z, collapse = "_"), character(1))
  missing_ids <- setdiff(all_GWAS_sign$rs, extracted_ids)
  if (length(missing_ids))
    stop(length(missing_ids), " significant GWAS SNPs missing from extracted VCF")

  status <- system2("vcftools", c("--gzvcf", shQuote(sign_vcf), "--hap-r2",
                                  "--ld-window", "100", "--out", shQuote(ld_prefix)))
  if (status != 0 || !file.exists(ld_file)) stop("vcftools LD calculation failed")
  writeLines(all_GWAS_sign$rs, ld_ids_file)
} else {
  if (!file.exists(ld_file)) stop("LD file absent; set recompute_ld = TRUE")
  if (file.exists(ld_ids_file)) {
    old_ids <- readLines(ld_ids_file)
    if (!setequal(old_ids, all_GWAS_sign$rs))
      stop("LD was calculated on a different significant SNP set; recompute LD")
  } else {
    warning("LD input SNP list is undocumented: set recompute_ld = TRUE before final analysis")
  }
}

ld_snps <- read.table(ld_file, header = TRUE, check.names = TRUE)
if (!all(c("CHR", "POS1", "POS2", "R.2") %in% names(ld_snps)))
  stop("LD file needs CHR, POS1, POS2, and R.2 columns")
ld_snps$R.2 <- suppressWarnings(as.numeric(as.character(ld_snps$R.2)))
ld_snps$ID1 <- paste(ld_snps$CHR, ld_snps$POS1, sep = "_")
ld_snps$ID2 <- paste(ld_snps$CHR, ld_snps$POS2, sep = "_")

## ------------------------------------------------------------
## 2b. Validate that signed D uses the same allele coding as GWAS.
## ------------------------------------------------------------

if (phase_align_effects) {
  if (!'D' %in% names(ld_snps))
    stop('This LD file has no D column; --hap-r2 is required')

  input_prep <- new.env(parent = globalenv())
  load('RData/Preparatory_GLMM_NaCl_only_withBeta.RData', envir = input_prep)
  X_check <- input_prep$genotype_matrix

  check_edges <- ld_snps[
    is.finite(ld_snps$R.2) & ld_snps$R.2 > 0.5 &
      ld_snps$ID1 %in% colnames(X_check) &
      ld_snps$ID2 %in% colnames(X_check), , drop = FALSE]
  check_edges$D <- as.numeric(as.character(check_edges$D))
  check_edges <- check_edges[is.finite(check_edges$D) &
                               check_edges$D != 0, , drop = FALSE]
  if (!nrow(check_edges)) stop('No eligible edges for signed D check')

  set.seed(123)
  check_edges <- check_edges[sample.int(nrow(check_edges),
                 min(phase_check_pairs, nrow(check_edges))), , drop = FALSE]
  check_edges$r_genotype <- vapply(seq_len(nrow(check_edges)), function(k) {
    suppressWarnings(cor(X_check[, check_edges$ID1[k]],
                         X_check[, check_edges$ID2[k]],
                         use = 'pairwise.complete.obs'))
  }, numeric(1))
  check_edges <- check_edges[is.finite(check_edges$r_genotype), , drop = FALSE]
  if (!nrow(check_edges)) stop('No valid genotype correlations for D check')
  match_fraction <- mean(sign(check_edges$D) ==
                           sign(check_edges$r_genotype))
  cat('Signed D / GWAS genotype correlation agreement:',
      round(match_fraction, 4), '\n')
  cat('Median difference between LD R2 and genotype r2:',
      median(abs(check_edges$R.2 - check_edges$r_genotype^2)), '\n')
  if (match_fraction < 0.99)
    stop('Signed D disagrees with GWAS genotype coding: do not phase-align')
  rm(X_check, input_prep)
  invisible(gc())
}

align_component_signs <- function(component_data, component_edges, beta_cols) {
  ids <- as.character(component_data$rs)
  if (length(ids) <= 1L) {
    return(list(data = component_data, phase = setNames(1L, ids),
                conflicts = 0L, n_edges = 0L))
  }
  ed <- component_edges[
    component_edges$ID1 %in% ids & component_edges$ID2 %in% ids,
    c('ID1', 'ID2', 'R.2', 'D'), drop = FALSE]
  ed$D <- as.numeric(as.character(ed$D))
  ed$R.2 <- as.numeric(as.character(ed$R.2))
  ed <- ed[is.finite(ed$R.2) & is.finite(ed$D) & ed$D != 0, , drop = FALSE]
  if (!nrow(ed)) stop('No signed LD edges in component')

  graph <- igraph::graph_from_data_frame(
    data.frame(from = ed$ID1, to = ed$ID2,
               cost = 1 - ed$R.2, edge_sign = sign(ed$D)),
    directed = FALSE, vertices = data.frame(name = ids))
  if (igraph::components(graph)$no != 1L)
    stop('Finite signed LD edges do not connect the entire component')

  # Minimum cost = maximum R2 spanning tree
  tree <- igraph::mst(graph, weights = igraph::E(graph)$cost)
  te <- igraph::as_data_frame(tree, what = 'edges')

  adj <- setNames(vector('list', length(ids)), ids)
  for (k in seq_len(nrow(te))) {
    a <- as.character(te$from[k]); b <- as.character(te$to[k])
    s <- as.integer(te$edge_sign[k])
    tmp <- adj[[a]]; tmp[b] <- s; adj[[a]] <- tmp
    tmp <- adj[[b]]; tmp[a] <- s; adj[[b]] <- tmp
  }

  phase <- setNames(rep(NA_integer_, length(ids)), ids)
  root <- ids[1L]
  phase[root] <- 1L
  queue <- character(length(ids))
  head <- 1L; tail <- 1L; queue[tail] <- root
  while (head <= tail) {
    current <- queue[head]; head <- head + 1L
    neighbors <- adj[[current]]
    for (neighbor in names(neighbors)) {
      if (!is.na(phase[neighbor])) next
      phase[neighbor] <- phase[current] * unname(neighbors[neighbor])
      tail <- tail + 1L; queue[tail] <- neighbor
    }
  }
  if (anyNA(phase)) stop('Phase could not be assigned to every SNP')

  bad <- phase[ed$ID1] * phase[ed$ID2] != sign(ed$D)
  out <- component_data
  # Only a temporary copy of SNP betas for sign splitting is changed!
  for (nm in beta_cols)
    out[[nm]] <- out[[nm]] * unname(phase[out$rs])

  list(data = out, phase = phase,
       conflicts = sum(bad), n_edges = nrow(ed),
       weighted_conflict_fraction = sum(ed$R.2[bad]) / sum(ed$R.2))
}


## ------------------------------------------------------------
## 3. Split one LD component by compatible effect directions
##     Preserve your original greedy algorithm and SNP order.
## ------------------------------------------------------------

split_clusters <- function(sign_all) {
  beta_cols <- paste0("beta_T", seq_along(vect_P_traits))
  sig_cols <- paste0("is_sign_T", seq_along(vect_P_traits))

  directions <- as.matrix(sign_all[, beta_cols, drop = FALSE])
  sig <- as.matrix(sign_all[, sig_cols, drop = FALSE])
  directions[!sig] <- 0
  directions <- sign(directions)

  ## Original rule: larger significant-effect profile goes first.
  n_nonzero <- rowSums(directions != 0)
  ord <- order(n_nonzero, decreasing = TRUE)
  directions <- directions[ord, , drop = FALSE]
  rs <- sign_all$rs[ord]

  ## Identical sign patterns must receive identical subgroup membership.
  patterns <- apply(directions, 1, paste, collapse = ":")
  keep <- !duplicated(patterns)
  unique_directions <- directions[keep, , drop = FALSE]
  unique_patterns <- patterns[keep]
  group <- rep(NA_integer_, nrow(unique_directions))
  compatible <- function(a, b) all(a == b | a == 0 | b == 0)

  current_group <- 0L
  for (i in seq_len(nrow(unique_directions))) {
    if (!is.na(group[i])) next
    current_group <- current_group + 1L
    group[i] <- current_group
    if (i == nrow(unique_directions)) next

    for (j in seq.int(i + 1L, nrow(unique_directions))) {
      if (!is.na(group[j])) next
      old <- which(group == current_group)
      if (all(vapply(old, function(k)
        compatible(unique_directions[k, ], unique_directions[j, ]), logical(1))))
        group[j] <- current_group
    }
  }
  data.frame(rs = rs, subgroup = group[match(patterns, unique_patterns)],
             stringsAsFactors = FALSE)
}

## ------------------------------------------------------------
## 4. LD components and final clusters for each threshold
## ------------------------------------------------------------

for (threshold_clust in thresholds_clust) {
  cat("\nClustering at R2 >", threshold_clust, "\n")
  edges_clust <- ld_snps[is.finite(ld_snps$R.2) &
                           ld_snps$R.2 > threshold_clust, c("ID1", "ID2"), drop = FALSE]

  ## Construct the graph with the full SNP universe: singleton SNPs included.
  edges_clust <- as.data.frame(edges_clust)
  edges_clust <- edges_clust[edges_clust$ID1 %in% all_GWAS_sign$rs &
                               edges_clust$ID2 %in% all_GWAS_sign$rs &
                               edges_clust$ID1 != edges_clust$ID2, , drop = FALSE]
  g <- igraph::graph_from_data_frame(edges_clust, directed = FALSE,
                 vertices = data.frame(name = all_GWAS_sign$rs))
  membership <- igraph::components(g)$membership
  names(membership) <- igraph::V(g)$name
  degree <- igraph::degree(g)
  names(degree) <- igraph::V(g)$name
  cluster <- unname(membership[all_GWAS_sign$rs])
  cluster[unname(degree[all_GWAS_sign$rs]) == 0] <- NA_integer_

  all_GWAS_sign_with_clust <- all_GWAS_sign
  all_GWAS_sign_with_clust$POS <- all_GWAS_sign_with_clust$ps
  all_GWAS_sign_with_clust$cluster <- cluster

  ## Keep old scientific definition: LD component, further split by sign.
  in_component <- which(!is.na(cluster))
  by_component <- split(in_component, cluster[in_component])
  if (phase_align_effects) {
    signed_edges <- ld_snps[is.finite(ld_snps$R.2) &
                               ld_snps$R.2 > threshold_clust,
                             c('ID1', 'ID2', 'R.2', 'D'), drop = FALSE]
    eid <- match(signed_edges$ID1, all_GWAS_sign$rs)
    group_index <- cluster[eid]
    signed_edges <- split(signed_edges, group_index)
  }
  split_one <- function(idx) {
    this_comp <- cluster[idx[1]]
    component_data <- all_GWAS_sign_with_clust[idx, , drop = FALSE]
    if (phase_align_effects) {
      aligned <- align_component_signs(
        component_data, signed_edges[[as.character(this_comp)]],
        paste0('beta_T', seq_along(vect_P_traits)))
      # Original beta_T* columns remain untouched in the exported SNP table.
      component_data <- aligned$data
      phase_stats <- data.frame(component = this_comp,
                                n_SNPs = length(idx),
                                conflicts = aligned$conflicts,
                                n_edges = aligned$n_edges,
                                weighted_conflict_fraction =
                                  aligned$weighted_conflict_fraction)
    } else {
      phase_stats <- NULL
    }
    out <- split_clusters(component_data)
    out$key <- paste0('component_', this_comp, '_', out$subgroup)
    list(assignments = out, phase_stats = phase_stats)
  }
  if (.Platform$OS.type == "unix" && split_cores > 1L) {
    component_splits <- parallel::mclapply(by_component, split_one,
                                           mc.cores = split_cores)
  } else {
    component_splits <- lapply(by_component, split_one)
  }
  split_assignments <- if (length(component_splits))
    do.call(rbind, lapply(component_splits, `[[`, 'assignments'))
  else data.frame(rs = character(0), subgroup = integer(0), key = character(0))

  if (phase_align_effects && length(component_splits)) {
    phase_summary <- do.call(rbind, lapply(component_splits, `[[`, 'phase_stats'))
    write.csv(phase_summary,
              file.path(clusters_dir, paste0('phase_alignment_',
                                            threshold_label(threshold_clust), '.csv')),
              row.names = FALSE)
    cat('  LD-sign conflicts outside spanning trees:',
        sum(phase_summary$conflicts), 'of', sum(phase_summary$n_edges), 'edges\n')
    if (any(phase_summary$weighted_conflict_fraction > 0.01))
      warning('Some LD components have contradictory signed phase links; ',
              'inspect phase_alignment_*.csv')
  }

  group_key <- split_assignments$key[match(all_GWAS_sign_with_clust$rs,
                                            split_assignments$rs)]
  singles <- is.na(cluster)
  group_key[singles] <- paste0("singleton_", all_GWAS_sign_with_clust$rs[singles])

  ## Preserve downstream SNP-level columns and numeric group identifiers.
  all_cluster_merged <- all_GWAS_sign_with_clust
  all_cluster_merged$group_key <- group_key
  all_cluster_merged <- all_cluster_merged[
    order(all_cluster_merged$CHROM, all_cluster_merged$POS), , drop = FALSE]
  all_cluster_merged$group <- match(all_cluster_merged$group_key,
                                    unique(all_cluster_merged$group_key))
  all_cluster_merged$group_key <- NULL
  rownames(all_cluster_merged) <- NULL

  if (anyNA(all_cluster_merged$group) ||
      anyDuplicated(all_cluster_merged$rs) ||
      nrow(all_cluster_merged) != nrow(all_GWAS_sign))
    stop("Clustering did not assign every significant SNP exactly once")

  out_file <- file.path(clusters_dir, paste0("all_sign_GWAS_with_clusters_",
                                            threshold_label(threshold_clust), ".txt"))
  write.table(all_cluster_merged, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Final clusters:", length(unique(all_cluster_merged$group)),
      "; singleton SNPs:", sum(tabulate(all_cluster_merged$group) == 1), "\n")
}

capture.output(sessionInfo(), file = file.path(clusters_dir, "02_sessionInfo.txt"))
cat("Done. GWAS/cluster outputs in", all_gwas_dir, "and", clusters_dir, "\n")
