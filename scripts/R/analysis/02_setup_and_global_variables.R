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

## switch to use either just the full model (0) or also the age- and sex-adjusted model
USE_BOTH_MODELS <- 0

## Model structure switch for survival analysis (matches clinical_effectiveness.qmd)
## Controls which model formula structure is used for BOTH base case and PSA
## Options:
##   0 = "joint"    - All biomarkers + all treatment interactions in ONE model (Model A)
##                    Formula: ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
##   1 = "focused"  - CRP and TMB/BRAF with BOTH interaction terms (Model B)
##                    Formula: ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx
##                    NOTE: TLR is NOT included in Model B
##   2 = "separate" - Only ONE biomarker + its interaction per model (Model C)
##                    Formula: ~ Age + sex + Rx + [biomarker]:Rx
##
## IMPORTANT: When changed, BOTH sampling cache AND PSA cache must be regenerated.
## Cache filenames will include this setting for safety.
MODEL_STRUCTURE <- 0  # Default: joint model (current base case behavior)