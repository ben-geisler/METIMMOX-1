# clear all objects from the work space
rm(list = ls())
# load faster binary format
data <- readRDS(file="data/tidy/METIMMOX.rds")

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
n_samples <- 5000  # number of bootstrap samples
n_sim <- n_samples # number of PSA simulations
set.seed(123)      #set seed for reproducibility

## global discount rate
dr = 0.04