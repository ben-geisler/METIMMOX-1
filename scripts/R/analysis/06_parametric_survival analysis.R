# Create time points sequence for prediction
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

# Initialize lists to store models and predictions
models <- list()
predictions <- list()

# Fit single models for each endpoint with all biomarkers and their interactions
# Note: Including Rx as main effect plus interactions
models$os <- flexsurvreg(Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx, 
                         data = data, dist = "weibull")
models$pfs <- flexsurvreg(Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx, 
                          data = data, dist = "weibull")

# Get reference values for treatments
exp_rx <- levels(data$Rx)[2]
ctrl_rx <- levels(data$Rx)[1]

# Calculate mean values for covariates (for prediction)
mean_age <- mean(data$Age, na.rm = TRUE)
# Assuming sex is binary (0/1), get the most common value or mean
mean_sex <- mean(as.numeric(data$sex), na.rm = TRUE)

# Predictions for control strategy 
# Setting biomarkers to 0 makes the interaction terms zero, effectively ignoring them
newdata_control <- data.frame(
  Age = mean_age,
  sex = mean_sex,
  Rx = factor(ctrl_rx, levels = levels(data$Rx)),
  crp = 0,
  tlr = 0,
  tmb_braf = 0
)

# Generate control predictions
control_os_pred <- predict(models$os, newdata = newdata_control, type = "survival", times = time_points)
control_pfs_pred <- predict(models$pfs, newdata = newdata_control, type = "survival", times = time_points)

# Unnest control predictions
control_os_unnested <- unnest(control_os_pred, .pred)
control_pfs_unnested <- unnest(control_pfs_pred, .pred)

# Store control predictions
predictions$control <- list(
  os = control_os_unnested$.pred_survival,
  pfs = control_pfs_unnested$.pred_survival
)

# Predictions for each biomarker strategy
for (biomarker in biomarkers) {
  # Create newdata for biomarker positive + experimental treatment
  newdata_pos_exp <- data.frame(
    Age = mean_age,
    sex = mean_sex,
    Rx = factor(exp_rx, levels = levels(data$Rx)),
    crp = 0,
    tlr = 0,
    tmb_braf = 0
  )
  newdata_pos_exp[[biomarker]] <- 1
  
  # Create newdata for biomarker negative + control treatment
  newdata_neg_ctrl <- data.frame(
    Age = mean_age,
    sex = mean_sex,
    Rx = factor(ctrl_rx, levels = levels(data$Rx)),
    crp = 0,
    tlr = 0,
    tmb_braf = 0
  )
  newdata_neg_ctrl[[biomarker]] <- 0
  
  # Generate predictions for all four scenarios
  # 1. OS - biomarker positive + experimental treatment
  os_pos_exp <- predict(models$os, 
                        newdata = newdata_pos_exp, 
                        type = "survival", 
                        times = time_points)
  
  # 2. PFS - biomarker positive + experimental treatment
  pfs_pos_exp <- predict(models$pfs, 
                         newdata = newdata_pos_exp, 
                         type = "survival", 
                         times = time_points)
  
  # 3. OS - biomarker negative + control treatment
  os_neg_ctrl <- predict(models$os, 
                         newdata = newdata_neg_ctrl, 
                         type = "survival", 
                         times = time_points)
  
  # 4. PFS - biomarker negative + control treatment
  pfs_neg_ctrl <- predict(models$pfs, 
                          newdata = newdata_neg_ctrl, 
                          type = "survival", 
                          times = time_points)
  
  # Unnest all predictions
  os_pos_exp_unnested <- unnest(os_pos_exp, .pred)
  pfs_pos_exp_unnested <- unnest(pfs_pos_exp, .pred)
  os_neg_ctrl_unnested <- unnest(os_neg_ctrl, .pred)
  pfs_neg_ctrl_unnested <- unnest(pfs_neg_ctrl, .pred)
  
  # Store all four predictions for this biomarker
  predictions[[biomarker]] <- list(
    os_pos_exp = os_pos_exp_unnested$.pred_survival,   # OS - biomarker+ & experimental
    pfs_pos_exp = pfs_pos_exp_unnested$.pred_survival, # PFS - biomarker+ & experimental
    os_neg_ctrl = os_neg_ctrl_unnested$.pred_survival, # OS - biomarker- & control
    pfs_neg_ctrl = pfs_neg_ctrl_unnested$.pred_survival # PFS - biomarker- & control
  )
  
  # Calculate weighted averages for the biomarker-guided strategy
  p_biomarker <- switch(biomarker,
                        "crp" = p_crp,
                        "tlr" = p_tlr,
                        "tmb_braf" = p_tmb_braf)
  
  predictions[[biomarker]]$os <- 
    p_biomarker * predictions[[biomarker]]$os_pos_exp + 
    (1 - p_biomarker) * predictions[[biomarker]]$os_neg_ctrl
  
  predictions[[biomarker]]$pfs <- 
    p_biomarker * predictions[[biomarker]]$pfs_pos_exp + 
    (1 - p_biomarker) * predictions[[biomarker]]$pfs_neg_ctrl
}

objects_to_remove <- c("biomarker", "time_points", "time_points_length", 
                       "exp_rx", "ctrl_rx", "mean_age", "mean_sex",
                       "newdata_control", "control_os_pred", "control_pfs_pred",
                       "control_os_unnested", "control_pfs_unnested",
                       "newdata_pos_exp", "newdata_neg_ctrl", 
                       "os_pos_exp", "pfs_pos_exp", "os_neg_ctrl", "pfs_neg_ctrl",
                       "os_pos_exp_unnested", "pfs_pos_exp_unnested", 
                       "os_neg_ctrl_unnested", "pfs_neg_ctrl_unnested",
                       "p_biomarker")
rm(list = objects_to_remove[objects_to_remove %in% ls()])
rm(objects_to_remove)
ls()
