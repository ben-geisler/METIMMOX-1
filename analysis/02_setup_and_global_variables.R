# load faster binary format
rds_path <- here::here("data", "tidy", "METIMMOX.rds")
if (file.exists(rds_path)) {
  data <- readRDS(rds_path)
} else {
  stop("File not found: ", rds_path)
}
rm(rds_path)

data$PFSwk <- data$`Days until progression`/7
data$OSwk <- data$`Days until death/last follow up`/7

# First on-treatment CT (issue #155): weeks from inclusion (the OS/PFS clock
# origin) to the scan at which tumour lesion reduction (TLR) is read. The trial
# export holds the scan date in the positional column `Date...122`, the date
# column that follows the baseline target-lesion block (`TL LD...121`), just as
# 01_data_prep.R reads the first on-treatment target-lesion sum from
# `TL LD...130`. CT1wk is NA exactly when TLR is NA (no first on-treatment
# scan). It defines the per-patient landmark in tlr_landmark.R and the
# first-scan diagnostics of the week-9 landmark; the economic model ignores it.
data$CT1wk <- as.numeric(as.Date(data$`Date...122`) -
                           as.Date(data$`Date of inclusion`)) / 7

# Rename "Progression exit" to "Progression" for simplicity
names(data)[names(data) == "Progression exit"] <- "Progression"
names(data)[names(data) == "Sex 0female"] <- "sex"

# PFS endpoint (issue #149): the trial export flags exit-for-progression only and
# censors deaths without progression. Recode PFS as progression OR death, with
# time to death for patients who died progression-free, matching the trial
# definition (Ree et al. 2024) and the partitioned-survival-model state
# definition. Raw values are kept as ProgressionExit / TTPwk.
source(here::here("R/pfs_endpoint.R"))
data <- derive_pfs_endpoint(data)

# global variables
## time parameters
cl <- 1/52         # 1 week cycle (not accounting for leap years)
time_horizon = 520 # 10 years in weeks

## analysis parameters
WTP <- 51000       # CE threshold (in Euros)
DSA_mult <- 0.2    # +/- 20% variations as standard for DSA
n_samples <- 5000  # number of resampled models
n_sim <- n_samples # number of PSA simulations
analysis_seed <- 123L  # default passed explicitly to each stochastic block

## global discount rate
dr = 0.04

## population parameters for scaling EVPPI to population level
## Source: Norwegian Cancer Registry (Kreftregisteret)
## - Total CRC incidence: ~4,000 cases/year
## - MSS/pMMR proportion: ~70-75% of metastatic CRC (international literature)
## - Metastatic proportion: ~50% present with metastatic disease
## - Conservative estimate: 1,500 annual eligible patients accounts for
##   diagnostic testing rates and treatment eligibility criteria
annual_incidence_norway <- 1500  # Annual MSS/pMMR mCRC cases in Norway
research_horizon_years <- 10     # Time horizon for research value (years)
discount_rate_research <- 0.035  # Discount rate for future research benefits (3.5%)
                                 # Standard rate for public health research in Norway
                                 # (slightly lower than 4% used for costs/QALYs)

## Utility source switch
## Controls which health state utility values are used as base case
## Options:
##   0 = "ipd"      - IPD-derived from METIMMOX trial (u_np = 0.9077, u_p = 0.9005)
##   1 = "correct"  - EQ-5D from CORRECT trial (u_np = 0.73, u_p = 0.59)
##                     (Grothey et al. Lancet 2013, cited by Gourzoulidis et al. J Comp Effect Res 2018)
##
## IMPORTANT: When changed, PSA and EVPPI caches must be regenerated.
## Sampling cache (survival models) is NOT affected.
UTILITY_SOURCE <- 1  # Default: CORRECT trial utilities (Gourzoulidis et al. 2018)
utility_source_label <- c("ipd", "correct")[UTILITY_SOURCE + 1]

# ===============================================================================
# SOURCE MODEL CONFIGURATIONS (Single Source of Truth)
# ===============================================================================
# This sources the central economic model configuration file which defines all
# formulas and provides helper functions: get_strategies(), get_biomarkers(), etc.
source(here::here("R/model_configs.R"))
source(here::here("R/prediction_population.R"))

# Source cache-path helpers (utility label + cached-object file locations).
source(here::here("R/cache_paths.R"))

# Source the shared PSA/DSA/EVPPI parameter specification.
source(here::here("R/parameter_distributions.R"))

# Set global strategy and biomarker vectors for backward compatibility.
# Economic analyses include control, CRP, and TMB/BRAF only.
strategies <- get_strategies()
biomarkers <- get_biomarkers()
