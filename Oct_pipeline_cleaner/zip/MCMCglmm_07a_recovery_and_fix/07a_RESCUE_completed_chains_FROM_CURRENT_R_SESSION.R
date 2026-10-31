## 07a_RESCUE_completed_chains_FROM_CURRENT_R_SESSION.R
## Run this in the SAME live R session in which the original script stopped.
## It NEVER refits a model. It saves completed VCV draws as resumption files.
## Do this BEFORE clearing the workspace or closing R.

required <- c("chains_VCV", "traits", "phen", "trait_mean", "trait_sd",
              "pair_n", "settings", "run_mode", "residual_structure",
              "n_chains", "fixed_form", "random_form", "rcov_form", "prior")
missing <- required[!vapply(required,
  function(nm) exists(nm, envir = .GlobalEnv, inherits = FALSE), logical(1))]
if (length(missing)) stop("Cannot rescue; objects missing in current R session: ",
                          paste(missing, collapse = ", "))
stopifnot(length(chains_VCV) == n_chains,
          all(vapply(chains_VCV, is.matrix, logical(1))))

output_dir <- "RData_G_NaCl"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_paths <- c("data/all_phenotypes_from_movies_fixed.txt",
                 "data/removed_isotypes_Cemee_2025.txt",
                 "data/correspondance_isotypes.csv")
input_hashes <- tools::md5sum(input_paths)
if (anyNA(input_hashes)) stop("Missing input files; run from project root")
checkpoint_meta <- list(
  traits = traits,
  input_hashes = input_hashes,
  trait_mean = trait_mean[traits], trait_sd = trait_sd[traits],
  pair_n = pair_n, n_records = nrow(phen),
  n_RIL = nlevels(phen$pop_label),
  n_dates = nlevels(phen$date_str),
  n_sessions = nlevels(phen$session_id),
  run_mode = run_mode, residual_structure = residual_structure,
  settings = settings, n_chains = n_chains,
  fixed = deparse(fixed_form), random = deparse(random_form),
  rcov = deparse(rcov_form), prior = prior
)
for (chain in seq_len(n_chains)) {
  path <- file.path(output_dir, sprintf("07a_%s_%s_chain%02d.rds",
                                           run_mode, residual_structure, chain))
  if (file.exists(path)) stop("Already exists: ", path,
                               " (do not overwrite blindly)")
  tmp <- tempfile("07a_rescue_", tmpdir = output_dir)
  saveRDS(list(VCV = chains_VCV[[chain]], meta = checkpoint_meta),
          tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Unable to save checkpoint: ", path)
  cat("Recovered chain", chain, "->", path, "\n")
}
cat("All completed chains recovered. You can now run the FIXED 07a script.\n")
