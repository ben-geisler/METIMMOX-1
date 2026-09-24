# Aggregate-only distribution audit before downstream regeneration.
source("analysis/02_setup_and_global_variables.R")
source("analysis/03_biomarker_strategies.R")
stopifnot(nrow(data_complete) == 68, sum(data_complete$Progression) == 50)
old <- economic_prediction_population(derive_pfs_endpoint(data, Inf, pfs_last_assessment))
stopifnot(nrow(old) == 68, sum(old$Progression) == 63,
          setequal(old$ID, data_complete$ID))
source("analysis/04_parametric_survival_analysis.R")
print(models$best_fit[c("os_distribution", "pfs_distribution", "ordering_aic_penalty")])
print(head(models$ordered_selection$pairs, 10))
cat("Complete cases:", nrow(data_complete), "PFS events:", sum(data_complete$Progression), "\n")
