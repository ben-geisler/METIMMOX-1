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

cat("Trial utility tests passed.\n")
