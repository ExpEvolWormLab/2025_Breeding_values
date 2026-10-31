## 48_Run_unified_main_analysis.R
## Run from PROJECT ROOT after corrected 08, 09f, 11, 13, 33.
## Deterministic, checkpointed, expensive primary analyses.
## GA3 excluded from PRIMARY; included only in the supplement.
source('43_Unified_config_and_helpers.R')
source('43_Prepare_unified_R5_catalogue.R')
calibration_gate <- file.path(uf_config$prefix,'43b_Calibration_interpretation_gate.txt')
calibration_signature_file <- file.path(uf_config$prefix,'43b_calibration_input_signature.txt')
if(!file.exists(calibration_gate) || !file.exists(calibration_signature_file) ||
   !identical(readLines(calibration_signature_file,warn=FALSE)[1L],
              uf_calibration_signature(uf_config)))
  source('43b_Empirical_WF_null_validation.R')
source('44_Calibrated_WF_selection_all_R5_clusters.R')
source('45_Unified_GWAS_evolution_LOCO_prediction.R')
source('46_Unified_main_figures_and_focal_annex.R')
cat('\nUnified final analysis COMPLETE. Results: ',uf_config$prefix,'\n',sep='')
