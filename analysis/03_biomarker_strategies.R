#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(dplyr)

# Create binary biomarker variables with simplified names
# crp uses CRP1cat: the week-4 (cycle 3 day 1, visit 3) CRP < 5 mg/L indicator,
# measured after two FLOX cycles common to both arms and before the first
# nivolumab dose. It is the value available when the immunotherapy decision is
# made (issue #150). CRP0cat (cycle 1 day 1) is the true baseline and is unused.
data$crp <- as.numeric(data$CRP1cat == 1)
data$tlr <- as.numeric(data$TLRcat == 1)
data$tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))

# Convert sex to factor BEFORE creating subsets
# This ensures all subset datasets (data_control, data_crp, etc.) have sex as factor
data$sex <- as.factor(data$sex)

# Create limited dataset with only the variables we need.
# ProgressionExit / TTPwk (raw progression-exit flag and time, issue #149) and
# CT1wk (first on-treatment CT week, issue #155) are carried for the clinical
# and DAG reports: the time-dependent PFS -> OS test needs the progression-only
# event, and the TLR landmark analyses need the scan date. No economic script
# reads them. LastEvalwk and PFS_rule also support clinical sensitivity and
# landmark diagnostics; every complete-case subset names its columns explicitly.
data <- data %>%
  select(ID, PFSwk, Progression, OSwk, Death, Rx, crp, tlr, tmb_braf, Age, sex,
         ProgressionExit, TTPwk, CT1wk, LastEvalwk, PFS_rule)

# Extract control and experimental groups
data_control <- subset(data, Rx == levels(data$Rx)[1])
data_experimental <- subset(data, Rx == levels(data$Rx)[2])

# One complete-case economic target population; clinical TLR uses its own cohort.
data_complete <- economic_prediction_population(data)
p_crp <- mean(data_complete$crp)
p_tlr <- mean(data$tlr, na.rm = TRUE)
p_tmb_braf <- mean(data_complete$tmb_braf)

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

# Define key vectors from the central economic-model configuration.
strategies <- get_strategies()
biomarkers <- get_biomarkers()

# Join keyed metadata and calculated values by strategy ID. This remains correct
# if the configured strategy order changes.
strategy_metadata <- get_strategy_metadata(strategies)
prevalence_by_strategy <- c(
  setNames(1, get_control_strategy()),
  setNames(vapply(biomarkers, function(x) get(paste0("p_", x)), numeric(1)),
           biomarkers)
)
n_by_strategy <- setNames(
  rep(nrow(data_complete), length(strategies)),
  strategies
)
strategies_df <- transform(
  strategy_metadata,
  prevalence = unname(prevalence_by_strategy[id]),
  n_patients = unname(n_by_strategy[id])
)

# Print the final dataframe
print(strategies_df)
