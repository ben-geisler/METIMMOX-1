# Trial-based EQ-5D-5L health-state utilities (issue #161) -----------------------
#
# Derives progression-free and progressed utilities from the METIMMOX EQ-5D-5L
# responses for the Table S8 utility scenario. Requires R/eq5d5l_utility.R.
#
# A response counts as progressed when the patient has a recorded progression
# ("Progression exit" = 1) and the response is dated strictly after the
# "First progression" date. For patients without a progression that column
# holds the last assessment date, so the date alone does not identify
# progression (the archived notebook, archive/05_QALYs.Rmd, used the date alone
# and gave 0.910 / 0.910 instead of 0.912 / 0.897 with the Danish value set).
# The utilities are means over responses, so patients with more responses weigh
# more; they are scenario inputs, not a fitted longitudinal model.

TRIAL_EQ5D_PATH <- "data/sensitive/EQ5Dmanual.xlsx"
TRIAL_EQ5D_DIMENSIONS <- c("mo", "sc", "ua", "pd", "ad")

#' Read the METIMMOX EQ-5D-5L responses
#'
#' @param path Path of the EQ-5D export (one row per response with ID, Date and
#'   the five dimension levels mo, sc, ua, pd, ad).
#' @return Data frame with ID (character), Date (Date) and integer levels.
read_trial_eq5d <- function(path = here::here(TRIAL_EQ5D_PATH)) {
  if (!file.exists(path)) stop("EQ-5D-5L export not found: ", path)
  df_raw <- readxl::read_excel(path, sheet = 1, col_types = "text")
  v_missing <- setdiff(c("ID", "Date", TRIAL_EQ5D_DIMENSIONS), names(df_raw))
  if (length(v_missing) > 0) stop("EQ-5D-5L export lacks columns: ", paste(v_missing, collapse = ", "))
  df_eq5d <- data.frame(ID = as.character(df_raw$ID),
                        Date = as.Date(df_raw$Date, format = "%Y-%m-%d"),
                        stringsAsFactors = FALSE)
  for (v in TRIAL_EQ5D_DIMENSIONS) df_eq5d[[v]] <- suppressWarnings(as.integer(df_raw[[v]]))
  df_eq5d
}

#' Progression-free and progressed EQ-5D-5L utilities from trial responses
#'
#' @param eq5d Responses as returned by read_trial_eq5d().
#' @param trial_data Trial data with ID, "First progression" (date) and
#'   "Progression exit" (0/1), e.g. data/tidy/METIMMOX.rds.
#' @param value_set "dk" (Danish EQ-5D-5L value set, Jensen et al. 2021) or "uk".
#' @return List with u_np, u_p, value_set, the numbers of responses and patients
#'   overall and per state, and the number of responses not matched to the trial.
trial_eq5d_utilities <- function(eq5d, trial_data, value_set = c("dk", "uk")) {
  value_set <- match.arg(value_set)
  v_needed <- c("ID", "First progression", "Progression exit")
  v_missing <- setdiff(v_needed, names(trial_data))
  if (length(v_missing) > 0) stop("Trial data lack columns: ", paste(v_missing, collapse = ", "))
  if (any(is.na(eq5d$Date))) stop("EQ-5D-5L responses with a missing or unparseable date.")
  m_levels <- as.matrix(eq5d[TRIAL_EQ5D_DIMENSIONS])
  if (any(is.na(m_levels)) || any(!m_levels %in% 1:5)) {
    stop("EQ-5D-5L dimension levels must be integers 1 to 5.")
  }
  df_trial <- data.frame(ID = as.character(trial_data$ID),
                         progression_date = as.Date(trial_data[["First progression"]]),
                         progressed = trial_data[["Progression exit"]] == 1,
                         stringsAsFactors = FALSE)
  if (anyDuplicated(df_trial$ID)) stop("Trial data IDs are not unique.")
  v_matched <- eq5d$ID %in% df_trial$ID
  df_resp <- merge(eq5d[v_matched, , drop = FALSE], df_trial, by = "ID")
  if (any(is.na(df_resp$progressed)) ||
      any(df_resp$progressed & is.na(df_resp$progression_date))) {
    stop("Matched patients need a progression flag and, if progressed, a progression date.")
  }
  scorer <- switch(value_set, dk = eq5d5l_utility_dk, uk = eq5d5l_utility_uk)
  v_utility <- scorer(df_resp$mo, df_resp$sc, df_resp$ua, df_resp$pd, df_resp$ad)
  v_after <- df_resp$progressed & df_resp$Date > df_resp$progression_date
  if (!any(v_after) || all(v_after)) stop("Both health states need at least one response.")
  list(u_np = mean(v_utility[!v_after]), u_p = mean(v_utility[v_after]),
       value_set = value_set,
       n_responses = nrow(df_resp), n_patients = length(unique(df_resp$ID)),
       n_responses_np = sum(!v_after), n_responses_p = sum(v_after),
       n_patients_p = length(unique(df_resp$ID[v_after])),
       n_unmatched_responses = sum(!v_matched))
}

#' Trial utilities with patient-clustered bootstrap standard errors (issue #188)
#'
#' Point estimates are those of trial_eq5d_utilities() on all responses. Patients
#' answer several times, so uncertainty comes from a patient-level bootstrap:
#' each replicate resamples matched patients with replacement (with all their
#' responses) and recomputes u_np, u_p and their difference u_decrement on the
#' same resample, which keeps the dependence between the two states. Replicates
#' in which a state has no response are discarded and counted.
#'
#' @param eq5d Responses as returned by read_trial_eq5d().
#' @param trial_data Trial data (see trial_eq5d_utilities()).
#' @param value_set "dk" or "uk".
#' @param n_boot Number of bootstrap replicates.
#' @param seed RNG seed; the caller's RNG state is restored afterwards.
#' @return The list of trial_eq5d_utilities() plus u_decrement, se (named
#'   numeric: u_np, u_p, u_decrement), n_boot and n_boot_discarded.
trial_eq5d_utility_uncertainty <- function(eq5d, trial_data, value_set = c("dk", "uk"),
                                           n_boot = 2000L, seed = 2026L) {
  value_set <- match.arg(value_set)
  l_point <- trial_eq5d_utilities(eq5d, trial_data, value_set)
  v_ids <- unique(eq5d$ID[eq5d$ID %in% as.character(trial_data$ID)])
  l_rows <- split(which(eq5d$ID %in% v_ids), eq5d$ID[eq5d$ID %in% v_ids])[v_ids]

  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv())
          else assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(seed)
  m_boot <- matrix(NA_real_, n_boot, 2L, dimnames = list(NULL, c("u_np", "u_p")))
  for (b in seq_len(n_boot)) {
    v_pick <- sample(v_ids, length(v_ids), replace = TRUE)
    df_b <- eq5d[unlist(l_rows[v_pick], use.names = FALSE), , drop = FALSE]
    # Resampled patients must stay distinct units: relabel duplicates.
    df_b$ID <- paste0(df_b$ID, "#", rep(seq_along(v_pick), lengths(l_rows[v_pick])))
    df_trial_b <- trial_data[match(sub("#.*$", "", unique(df_b$ID)), as.character(trial_data$ID)), , drop = FALSE]
    df_trial_b$ID <- unique(df_b$ID)
    l_b <- tryCatch(trial_eq5d_utilities(df_b, df_trial_b, value_set), error = function(e) NULL)
    if (!is.null(l_b)) m_boot[b, ] <- c(l_b$u_np, l_b$u_p)
  }
  v_ok <- stats::complete.cases(m_boot)
  if (sum(v_ok) < 0.9 * n_boot) stop("More than 10% of utility bootstrap replicates failed.")
  m_boot <- m_boot[v_ok, , drop = FALSE]
  v_decrement <- m_boot[, "u_np"] - m_boot[, "u_p"]
  c(l_point, list(
    u_decrement = l_point$u_np - l_point$u_p,
    se = c(u_np = stats::sd(m_boot[, "u_np"]), u_p = stats::sd(m_boot[, "u_p"]),
           u_decrement = stats::sd(v_decrement)),
    n_boot = n_boot, n_boot_discarded = sum(!v_ok), seed = seed))
}

#' Source note for trial-derived utilities, for table footnotes
#'
#' @param utilities List returned by trial_eq5d_utilities().
#' @return A single character string.
trial_eq5d_source <- function(utilities) {
  value_set <- switch(utilities$value_set,
    dk = "Danish EQ-5D-5L value set (Jensen et al. 2021, doi:10.1007/s40258-021-00639-3)",
    uk = "UK EQ-5D-5L value set")
  sprintf(paste0("METIMMOX EQ-5D-5L responses (%d from %d patients), %s; mean over ",
                 "responses, progressed = dated after a recorded progression ",
                 "(%d responses from %d patients)"),
          utilities$n_responses, utilities$n_patients, value_set,
          utilities$n_responses_p, utilities$n_patients_p)
}
