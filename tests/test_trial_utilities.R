# Synthetic tests for the trial EQ-5D-5L utility derivation (issue #161); no trial data.
source("R/eq5d5l_utility.R")
source("R/trial_utilities.R")
source("R/scenario_analysis.R")

# Patient A progressed on 2020-03-01; patient B never progressed, and its
# "First progression" column holds the last assessment (2020-02-01).
df_trial <- data.frame(ID = c("A", "B"),
                       `First progression` = as.POSIXct(c("2020-03-01", "2020-02-01"), tz = "UTC"),
                       `Progression exit` = c(1, 0), check.names = FALSE)
df_eq5d <- data.frame(
  ID = c("A", "A", "A", "B", "B", "Z"),
  Date = as.Date(c("2020-01-01", "2020-03-01", "2020-04-01", "2020-01-15", "2020-05-01", "2020-01-01")),
  mo = c(1L, 2L, 3L, 1L, 1L, 5L), sc = 1L, ua = 1L, pd = c(1L, 1L, 2L, 1L, 2L, 5L), ad = 1L
)
v_dk <- eq5d5l_utility_dk(df_eq5d$mo, df_eq5d$sc, df_eq5d$ua, df_eq5d$pd, df_eq5d$ad)

l_u <- trial_eq5d_utilities(df_eq5d, df_trial, value_set = "dk")
# Progressed: only A's response after its progression date (row 3). A response
# on the progression date (row 2) and B's response after its last assessment
# (row 5) are progression-free; the unmatched patient Z is excluded.
stopifnot(
  abs(l_u$u_p - v_dk[3]) < 1e-12,
  abs(l_u$u_np - mean(v_dk[c(1, 2, 4, 5)])) < 1e-12,
  l_u$n_responses == 5L, l_u$n_patients == 2L,
  l_u$n_responses_p == 1L, l_u$n_patients_p == 1L, l_u$n_responses_np == 4L,
  l_u$n_unmatched_responses == 1L, l_u$value_set == "dk"
)
l_uk <- trial_eq5d_utilities(df_eq5d, df_trial, value_set = "uk")
stopifnot(abs(l_uk$u_p - eq5d5l_utility_uk(3L, 1L, 1L, 2L, 1L)) < 1e-12)
cat("PASS: progression requires a recorded progression and a later response date.\n")

# Danish value set: published worked example 13224 = 0.439 (Jensen et al. 2021)
# and full health = 1.
stopifnot(abs(eq5d5l_utility_dk(1L, 3L, 2L, 2L, 4L) - 0.439) < 1e-12,
          eq5d5l_utility_dk(1L, 1L, 1L, 1L, 1L) == 1)
cat("PASS: Danish value set reproduces the published worked example.\n")

# Invalid input stops.
#' TRUE when the expression raises an error.
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
df_bad_level <- df_eq5d; df_bad_level$pd[1] <- 6L
df_bad_date <- df_eq5d; df_bad_date$Date[1] <- NA
df_no_flag <- df_trial; df_no_flag[["Progression exit"]] <- NULL
df_no_progression <- df_trial; df_no_progression[["Progression exit"]] <- 0
stopifnot(
  fails(trial_eq5d_utilities(df_bad_level, df_trial)),
  fails(trial_eq5d_utilities(df_bad_date, df_trial)),
  fails(trial_eq5d_utilities(df_eq5d, df_no_flag)),
  fails(trial_eq5d_utilities(df_eq5d, df_no_progression)),  # no progressed response
  fails(trial_eq5d_utilities(df_eq5d, rbind(df_trial, df_trial)))
)
cat("PASS: invalid levels, dates, flags, empty states and duplicate IDs stop.\n")

# The Table S8 scenario carries the derived pair and a citation, never a placeholder.
df_scen <- define_utility_scenarios(l_u)
stopifnot(nrow(df_scen) == 1L, df_scen$u_np == l_u$u_np, df_scen$u_p == l_u$u_p,
          grepl("Jensen et al. 2021", df_scen$source, fixed = TRUE),
          !grepl("PLACEHOLDER", df_scen$source, fixed = TRUE))
cat("PASS: utility scenario built from the trial-derived pair.\n")

# Patient-clustered bootstrap (issue #188): reproducible, leaves the caller's
# RNG untouched, point estimates equal the full-data values, SEs positive.
set.seed(5)
n_pat <- 40L
df_trial_b <- data.frame(ID = sprintf("P%02d", seq_len(n_pat)),
                         `First progression` = as.POSIXct("2020-06-01", tz = "UTC"),
                         `Progression exit` = rep(c(1, 0), length.out = n_pat), check.names = FALSE)
df_eq5d_b <- do.call(rbind, lapply(seq_len(n_pat), function(i) data.frame(
  ID = sprintf("P%02d", i), Date = as.Date(c("2020-01-01", "2020-09-01")),
  mo = sample(1:3, 2, TRUE), sc = 1L, ua = sample(1:2, 2, TRUE), pd = sample(1:3, 2, TRUE), ad = 1L)))
invisible(runif(1)); rng_before <- .Random.seed
l_b1 <- trial_eq5d_utility_uncertainty(df_eq5d_b, df_trial_b, "dk", n_boot = 200L, seed = 7L)
stopifnot(identical(.Random.seed, rng_before))
l_b2 <- trial_eq5d_utility_uncertainty(df_eq5d_b, df_trial_b, "dk", n_boot = 200L, seed = 7L)
l_b3 <- trial_eq5d_utility_uncertainty(df_eq5d_b, df_trial_b, "dk", n_boot = 200L, seed = 8L)
l_full <- trial_eq5d_utilities(df_eq5d_b, df_trial_b, "dk")
stopifnot(identical(l_b1$se, l_b2$se), !identical(l_b1$se, l_b3$se),
          l_b1$u_np == l_full$u_np, l_b1$u_p == l_full$u_p,
          abs(l_b1$u_decrement - (l_full$u_np - l_full$u_p)) < 1e-15,
          all(l_b1$se > 0), identical(names(l_b1$se), c("u_np", "u_p", "u_decrement")),
          l_b1$n_boot_discarded == 0L)
# Clustering: duplicating every response within its patient adds no
# information, so the patient-level SE stays (about) the same.
df_dup <- rbind(df_eq5d_b, df_eq5d_b)
l_dup <- trial_eq5d_utility_uncertainty(df_dup, df_trial_b, "dk", n_boot = 200L, seed = 7L)
stopifnot(abs(l_dup$se["u_np"] / l_b1$se["u_np"] - 1) < 0.25)
cat("PASS: patient-clustered bootstrap is reproducible, RNG-neutral and cluster-aware.\n")

# The SEs replace the assumed CV in the PSA distributions (issue #188).
source("R/model_configs.R")
source("R/parameter_distributions.R")
v_spec <- parameter_distribution_spec()
v_spec <- v_spec[v_spec$parameter %in% c("u_np", "u_decrement"), ]
l_params <- list(u_np = 0.9, u_decrement = 0.02)
l_cv <- create_parameter_distributions(l_params, v_spec)
l_params$utility_se <- c(u_np = 0.01, u_decrement = 0.02)
l_se <- create_parameter_distributions(l_params, v_spec)
v_beta_var <- with(l_se$u_np, shape1 * shape2 / ((shape1 + shape2)^2 * (shape1 + shape2 + 1)))
stopifnot(abs(l_cv$u_decrement$shape - 1 / 0.15^2) < 1e-12,
          abs(sqrt(v_beta_var) - 0.01) < 1e-12,
          abs(l_se$u_decrement$shape - 1) < 1e-12,           # CV = 0.02 / 0.02 = 1
          abs(l_se$u_decrement$shape / l_se$u_decrement$rate - 0.02) < 1e-12)
l_params$utility_se <- c(u_np = -1)
stopifnot(fails(create_parameter_distributions(l_params, v_spec)))
cat("PASS: utility_se replaces the assumed CV of u_np and u_decrement.\n")

cat("Trial utility tests passed.\n")
