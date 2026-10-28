## OPTIONAL: Align SNP effect signs within one LD connected component.
## Use only AFTER confirming sign(D) matches signed correlations from the
## exact GWAS genotype_matrix (see 02_LD_sign_diagnostic.R).
## Uses the strongest-LD spanning tree so chains need not all correlate
## directly with one lead SNP. Reports phase inconsistencies on extra edges.

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
