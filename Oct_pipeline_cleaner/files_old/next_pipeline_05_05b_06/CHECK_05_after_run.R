## CHECK_05_after_run.R
## Run after 05_Concatenate_cluster_GWAS.R for the R5 primary analysis.

label <- "R5"
f <- file.path("RData_GridLMM_NaCl",
                paste0("GWAS_cluster_results_BLUE_concatenated_", label,
                       "_with_residuals.RData"))
z <- new.env(parent = globalenv())
load(f, envir = z)
a <- z$all_GWAS_clusters
meta <- read.table(file.path("data_BLUE_GridLMM_NaCl",
                   paste0("data_clusters_BLUE_", label, ".txt")),
                   header = TRUE, stringsAsFactors = FALSE)

stopifnot(nrow(a) == nrow(meta),
          !anyDuplicated(a$rs),
          setequal(as.character(a$rs), as.character(meta$clusterID)),
          identical(as.character(a$chr), as.character(a$CHROM)),
          identical(as.numeric(a$ps), as.numeric(a$POS)))

for (trait in c("SF", "SB", "FS", "FB", "BS", "BF", "area.F", "beta")) {
  tag <- paste0(trait, "_NaCl")
  p <- a[[paste0("p_value_REML_", tag)]]
  q <- a[[paste0("p_value_REML.cor_", tag)]]
  sig <- a[[paste0("is_sign_", tag)]]
  observed <- which(is.finite(p))
  stopifnot(isTRUE(all.equal(q[observed], p.adjust(p[observed], "BH"))),
            identical(sig, !is.na(q) & q < 0.05),
            all(is.na(q[!is.finite(p)])))
  cat(trait, ": tested", length(observed), "; q < 0.05:", sum(sig), "\n")
}
cat("PASS: R5 cluster-GWAS metadata, p-values and ID alignment validated.\n")
