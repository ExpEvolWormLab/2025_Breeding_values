## 30_Build_extended_region_tables.R
## Extend the regional deep-dive to include the historical examples plus
## a small set of newly suggested R5 QTLs for direct competition analyses.
## Run from the project root.

existing_group_file <- '26_Regions_and_R5_groups.csv'
existing_focal_file <- '26_Focal_variants.csv'
novel_file <- 'deep_dive/discovery/29_Novel_QTL_shortlist_spatially_thinned.csv'
## Optional fallbacks if you keep the previous outputs in another folder.
if(!file.exists(existing_group_file))
  existing_group_file <- 'phase10_region_deepdive/26_Regions_and_R5_groups.csv'
if(!file.exists(existing_focal_file))
  existing_focal_file <- 'phase10_region_deepdive/26_Focal_variants.csv'
if(!file.exists(novel_file))
  novel_file <- 'historical_QTL_R5/deep_dive/discovery/29_Novel_QTL_shortlist_spatially_thinned.csv'

out_group_file <- '30_Regions_and_R5_groups_extended.csv'
out_focal_file <- '30_Focal_variants_extended.csv'
out_manifest <- '30_Extended_region_manifest.csv'

seed_candidates <- data.frame(
  region = c('beta_III_pair','beta_III_pair','beta_V_886','beta_area_X_1010','SB_II_249'),
  group  = c('504','506','886','1010','249'),
  stringsAsFactors = FALSE
)
## You can add more rows manually if you want to expand the panel.

stopifnot(file.exists(existing_group_file), file.exists(existing_focal_file), file.exists(novel_file))
old_groups <- read.csv(existing_group_file, stringsAsFactors = FALSE)
old_focal  <- read.csv(existing_focal_file, stringsAsFactors = FALSE)
novel <- read.csv(novel_file, stringsAsFactors = FALSE)
novel$group <- as.character(novel$group)

## Use the discovery table to recover chromosome and the best target trait.
get_target_trait <- function(x){
  qcols <- c('q_beta','q_area.F','q_SB','q_SF','q_FB','q_FS','q_BS','q_BF')
  qcols <- qcols[qcols %in% names(x)]
  if(!length(qcols)) return('beta')
  vals <- as.numeric(x[1, qcols, drop = TRUE])
  vals[!is.finite(vals)] <- Inf
  if(all(!is.finite(vals))) return('beta')
  sub('^q_', '', qcols[which.min(vals)])
}

extra_groups <- list()
extra_focal  <- list()
for(i in seq_len(nrow(seed_candidates))){
  g <- seed_candidates$group[i]
  z <- novel[novel$group == g, , drop = FALSE]
  if(!nrow(z)){
    warning('Novel group ', g, ' not found in discovery table; skipping')
    next
  }
  region_name <- seed_candidates$region[i]
  target_trait <- get_target_trait(z)
  extra_groups[[length(extra_groups) + 1L]] <- data.frame(
    region = region_name,
    CHROM = as.character(z$CHROM[1]),
    target_trait = target_trait,
    group = g,
    include = TRUE,
    stringsAsFactors = FALSE
  )
  ## Use the cluster centroid as a focal descriptive site so script 28 can draw
  ## exact trajectories in addition to the representative R5 tag SNPs.
  focal_pos <- suppressWarnings(as.integer(round(as.numeric(z$R5_centroid[1]))))
  if(is.finite(focal_pos)){
    extra_focal[[length(extra_focal) + 1L]] <- data.frame(
      region = region_name,
      CHROM = as.character(z$CHROM[1]),
      POS = focal_pos,
      hypothesis = 'R5_centroid_marker',
      stringsAsFactors = FALSE
    )
  }
}

extra_groups <- if(length(extra_groups)) do.call(rbind, extra_groups) else old_groups[0, ]
extra_focal  <- if(length(extra_focal)) do.call(rbind, extra_focal) else old_focal[0, ]

all_groups <- unique(rbind(old_groups, extra_groups))
all_groups$group <- as.character(all_groups$group)
all_groups$include <- ifelse(isTRUE(all_groups$include) | as.character(all_groups$include) %in% c('TRUE','True','true'), TRUE, FALSE)

all_focal <- unique(rbind(old_focal, extra_focal))
all_focal$POS <- as.integer(all_focal$POS)

write.csv(all_groups, out_group_file, row.names = FALSE)
write.csv(all_focal, out_focal_file, row.names = FALSE)
manifest <- merge(extra_groups, novel[, intersect(c('group','start','end','effect_beta','q_beta','effect_area.F','q_area.F','effect_SB','q_SB','score','reason'), names(novel))], by = 'group', all.x = TRUE, sort = FALSE)
write.csv(manifest, out_manifest, row.names = FALSE)
cat('Wrote ', out_group_file, ', ', out_focal_file, ' and ', out_manifest, '\n', sep = '')
cat('If you want additional regions, edit seed_candidates in this script and re-run it.\n')
