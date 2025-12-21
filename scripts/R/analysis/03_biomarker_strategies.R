#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(dplyr)

# Create binary biomarker variables with simplified names
data$crp <- as.numeric(data$CRP1cat == 1)
data$tlr <- as.numeric(data$TLRcat == 1)
data$tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))

# Convert sex to factor BEFORE creating subsets
# This ensures all subset datasets (data_control, data_crp, etc.) have sex as factor
data$sex <- as.factor(data$sex)

# Create limited dataset with only the variables we need
data <- data %>%
  select(ID, PFSwk, Progression, OSwk, Death, Rx, crp, tlr, tmb_braf, Age, sex)

# Extract control and experimental groups
data_control <- subset(data, Rx == levels(data$Rx)[1])
data_experimental <- subset(data, Rx == levels(data$Rx)[2])

# Calculate prevalence rates
p_crp <- mean(data$crp, na.rm = TRUE)
p_tlr <- mean(data$tlr, na.rm = TRUE)
p_tmb_braf <- mean(data$tmb_braf, na.rm = TRUE)

# Create biomarker-guided datasets
# CRP combined dataset: CRP+ from experimental, CRP- from control
data_crp <- rbind(
  transform(subset(data_experimental, crp == 1), biomarker_status = "positive"),
  transform(subset(data_control, crp == 0), biomarker_status = "negative")
)

# TLR combined dataset: TLR+ from experimental, TLR- from control
data_tlr <- rbind(
  transform(subset(data_experimental, tlr == 1), biomarker_status = "positive"),
  transform(subset(data_control, tlr == 0), biomarker_status = "negative")
)

# TMB/BRAF combined dataset: TMB/BRAF+ from experimental, TMB/BRAF- from control
data_tmb_braf <- rbind(
  transform(subset(data_experimental, tmb_braf == 1), biomarker_status = "positive"),
  transform(subset(data_control, tmb_braf == 0), biomarker_status = "negative")
)

# Add strategy labels for model identification
data_control$strategy <- "control"
data_crp$strategy <- "crp"
data_tlr$strategy <- "tlr"
data_tmb_braf$strategy <- "tmb_braf"

# Define key vectors
strategies <- c("control", "crp", "tlr", "tmb_braf")
biomarkers <- c("crp", "tlr", "tmb_braf")

# Create a comprehensive strategy dataframe with all relevant information
strategies_df <- data.frame(
  id = strategies,
  name = c("Standard of care: FLOX chemotherapy only", 
           "Biomarker-guided: C-reactive protein", 
           "Biomarker-guided: tumor lesion reduction", 
           "Biomarker-guided: tumor mutation burden or BRAF mutation"),
  description = c(
    "Standard of care - All patients receive only FLOX chemotherapy",
    "C-reactive protein with cut-off of <5 for biomarker-positive status. If CRP-positive: alternating two cycles each of FLOX (chemotherapy) and nivolumab (anti-PD1 immunotherapy); if CRP-negative: chemotherapy only",
    "Tumor Lesion Reduction with cut-off of >=10% for biomarker-positive status. If TLR-positive: alternating two cycles each of FLOX (chemotherapy) and nivolumab (anti-PD1 immunotherapy); if TLR-negative: chemotherapy only", 
    "Combined biomarker: either Tumor Mutation Burden >= 9 or BRAF V600 mutation positive (both from next-generation sequencing). If TMB/BRAF-positive: alternating two cycles each of FLOX (chemotherapy) and nivolumab (anti-PD1 immunotherapy); if TMB/BRAF-negative: chemotherapy only"
  ),
  # Store prevalence rates in dataframe
  prevalence = c(1.0, p_crp, p_tlr, p_tmb_braf),
  n_patients = c(nrow(data_control), 
                 nrow(data_crp), 
                 nrow(data_tlr), 
                 nrow(data_tmb_braf)),
  stringsAsFactors = FALSE
)

# Print the final dataframe
print(strategies_df)

# Clean up intermediate variables
rm(list = setdiff(ls(), c(
  # Main datasets
  "data", "data_control", "data_crp", "data_tlr", "data_tmb_braf",
  # Model parameters and structure
  "strategies", "biomarkers", "strategies_df",
  # Biomarker prevalence
  "p_crp", "p_tlr", "p_tmb_braf",
  # Other essential variables from 02_setup_and_global_variables.R
  "time_horizon", "cl", "WTP", "DSA_mult", "n_samples", "n_sim", "dr",
  # Model configuration switches
  "USE_BOTH_MODELS", "MODEL_STRUCTURE"
)))