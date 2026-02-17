# clear all objects from the work space
rm(list = ls())
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

# Rename "Progression exit" to "Progression" for simplicity
names(data)[names(data) == "Progression exit"] <- "Progression"
names(data)[names(data) == "Sex 0female"] <- "sex"

# global variables
## time parameters
cl <- 1/52         # 1 week cycle (not accounting for leap years)
time_horizon = 520 # 10 years in weeks

## analysis parameters
WTP <- 51000       # CE threshold (in Euros)
DSA_mult <- 0.2    # +/- 20% variations as standard for DSA
n_samples <- 5000  # number of resampled models
n_sim <- n_samples # number of PSA simulations
set.seed(123)      #set seed for reproducibility

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

## switch to use either just the full model (0) or also the age- and sex-adjusted model
USE_BOTH_MODELS <- 0

## Model structure switch for survival analysis
## Controls which model formula structure is used for BOTH base case and PSA
## Options:
##   0 = "joint"    - All biomarkers + all treatment interactions in ONE model (Model A)
##   1 = "focused"  - Per-strategy formulas (Model B):
##                    - CRP/TMB_BRAF: ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx
##                    - TLR: ~ Age + sex + Rx + crp + tmb_braf + tlr:Rx
##   2 = "separate" - Only ONE biomarker + its interaction per model (Model C)
##
## All three biomarker strategies (CRP, TLR, TMB/BRAF) are available in ALL models.
## See scripts/R/functions/model_configs.R for complete formula definitions.
##
## IMPORTANT: When changed, BOTH sampling cache AND PSA cache must be regenerated.
## Cache filenames will include this setting for safety.
MODEL_STRUCTURE <- 0  # Default: joint model (current base case behavior)

# ===============================================================================
# SOURCE MODEL CONFIGURATIONS (Single Source of Truth)
# ===============================================================================
# This sources the central model configuration file which defines all formulas
# and provides helper functions: get_strategies(), get_biomarkers(), etc.
source(here::here("scripts/R/functions/model_configs.R"))

# Set global strategy and biomarker vectors for backward compatibility
# These are now dynamically determined based on MODEL_STRUCTURE
strategies <- get_strategies()
biomarkers <- get_biomarkers()