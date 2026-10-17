# Third and fourth treatment sequences (issue #208)
# -----------------------------------------------------------------------------
# Trial team, 15 October 2026 (docs/open_questions_2026-10-15-answers.csv):
# SEQ1/SEQ2 - new sequences were also started on progression during a break,
# so their timing is not a fixed schedule and progression during a break is a
# PFS event (the endpoint is unchanged); SEQ3 - take the data as they are, as an
# indication of the share of patients receiving a third and a fourth sequence,
# nothing beyond; SEQ5 - later sequences keep the alternating 2 FLOX + 2
# nivolumab pattern. SEQ4 (actual exposure) was not answered, so every
# administration of a sequence is costed in full.
#
# The per-patient sequence data come from the export's administration block
# (derive_treatment_sequences() in R/data_contract.R, saved by script 01 as
# Seq3, Seq3wk, Seq4, Seq4wk). Script 05 turns them into the model inputs
# below for the economic prediction population; calculate_outcomes() costs
# them (later_sequence_resource_use()).

#' Administration weeks of one sequence by arm, from the first sequence
#'
#' The first sequence of TREATMENT_SCHEDULE_POSITIONS (positions 1-16, weeks
#' 0-15) is the pattern of every later sequence: FLOX every two weeks in the
#' control arm (8 administrations), FLOX at weeks 0, 2, 8, 10 and nivolumab at
#' weeks 4, 6, 12, 14 in the experimental arm (trial team, SEQ5).
#'
#' @param positions TREATMENT_SCHEDULE_POSITIONS.
#' @param sequence_weeks Length of one sequence in weeks.
#' @return list(control = list(FLOX, nivolumab), experimental = list(FLOX,
#'   nivolumab)) of integer week offsets from the sequence start.
later_sequence_regimen <- function(positions = TREATMENT_SCHEDULE_POSITIONS,
                                   sequence_weeks = 16L) {
  first <- function(v) as.integer(v[v <= sequence_weeks] - 1L)
  list(control = list(FLOX = first(positions$FLOX_control), nivolumab = integer(0)),
       experimental = list(FLOX = first(positions$FLOX_exp),
                           nivolumab = first(positions$nivo)))
}

#' Later-sequence model inputs for a population (issue #208)
#'
#' Counts by arm the patients of `population` who started a third and a
#' fourth sequence, and the observed start weeks (rounded to the weekly model
#' grid) of those who did. The shares are the base-case values of
#' LATER_SEQUENCE_SHARE_PARAMETERS: p_seq3_<arm> = n_seq3 / n_patients and
#' p_seq4_given_seq3_<arm> = n_seq4 / n_seq3 (0 without a third sequence).
#'
#' @param population Patients of the analysis (data_complete), with ID and Rx
#'   (first level control, second experimental).
#' @param sequences Per-patient sequence data with ID, Seq3, Seq3wk, Seq4,
#'   Seq4wk (data/tidy/METIMMOX.rds).
#' @param regimen later_sequence_regimen().
#' @return list(later_sequences = list(counts, start_weeks, regimen, source),
#'   shares = named list of the four shares).
later_sequence_inputs <- function(population, sequences, regimen = later_sequence_regimen()) {
  v_required <- c("ID", "Seq3", "Seq3wk", "Seq4", "Seq4wk")
  v_missing <- setdiff(v_required, names(sequences))
  if (length(v_missing))
    stop("Missing later-sequence columns (re-run analysis/01_data_prep.R): ",
         paste(v_missing, collapse = ", "), call. = FALSE)
  v_row <- match(population$ID, sequences$ID)
  if (anyNA(v_row)) stop("Later-sequence data lack patients of the population.", call. = FALSE)
  df_seq <- as.data.frame(sequences)[v_row, v_required, drop = FALSE]
  if (any(!df_seq$Seq3 %in% 0:1) || any(!df_seq$Seq4 %in% 0:1) ||
      any(df_seq$Seq4 == 1 & df_seq$Seq3 == 0))
    stop("Later-sequence indicators must be 0/1 with a fourth sequence only after a third.",
         call. = FALSE)
  for (k in 3:4) {
    v_started <- df_seq[[paste0("Seq", k)]] == 1
    v_week <- df_seq[[paste0("Seq", k, "wk")]]
    if (any(v_started & (!is.finite(v_week) | v_week < 0)))
      stop("Sequence ", k, " starts must have a nonnegative start week.", call. = FALSE)
  }
  v_arm <- CELL_ARMS[as.integer(population$Rx)]
  if (anyNA(v_arm)) stop("Population arms must be the two levels of Rx.", call. = FALSE)

  df_counts <- data.frame(
    arm = CELL_ARMS,
    n_patients = vapply(CELL_ARMS, function(a) sum(v_arm == a), numeric(1)),
    n_seq3 = vapply(CELL_ARMS, function(a) sum(df_seq$Seq3[v_arm == a]), numeric(1)),
    n_seq4 = vapply(CELL_ARMS, function(a) sum(df_seq$Seq4[v_arm == a]), numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL)

  l_rows <- list()
  for (a in CELL_ARMS) for (k in 3:4) {
    v_started <- v_arm == a & df_seq[[paste0("Seq", k)]] == 1
    if (!any(v_started)) next
    v_week <- as.integer(round(df_seq[[paste0("Seq", k, "wk")]][v_started]))
    v_tab <- table(v_week)
    l_rows[[length(l_rows) + 1L]] <- data.frame(
      arm = a, sequence = as.integer(k), week = as.integer(names(v_tab)),
      n = as.numeric(v_tab), stringsAsFactors = FALSE)
  }
  df_start <- if (length(l_rows)) do.call(rbind, l_rows) else
    data.frame(arm = character(), sequence = integer(), week = integer(), n = numeric())
  df_start <- df_start[order(df_start$arm, df_start$sequence, df_start$week, method = "radix"), ,
                       drop = FALSE]
  rownames(df_start) <- NULL

  v_shares <- list()
  for (i in seq_len(nrow(df_counts))) {
    a <- df_counts$arm[i]
    v_shares[[paste0("p_seq3_", a)]] <- df_counts$n_seq3[i] / df_counts$n_patients[i]
    v_shares[[paste0("p_seq4_given_seq3_", a)]] <-
      if (df_counts$n_seq3[i] > 0) df_counts$n_seq4[i] / df_counts$n_seq3[i] else 0
  }
  list(later_sequences = list(
         counts = df_counts, start_weeks = df_start, regimen = regimen,
         source = paste("METIMMOX administration block (export columns 2640-2879),",
                        "economic prediction population")),
       shares = v_shares)
}

#' Median and range of the start weeks of a later sequence by arm
#'
#' @param later params$later_sequences.
#' @return Data frame: arm, sequence, n (starters), median, min, max (weeks
#'   from inclusion, on the weekly model grid).
summarise_later_sequence_starts <- function(later) {
  df_start <- later$start_weeks
  if (!nrow(df_start)) return(data.frame(arm = character(), sequence = integer(), n = numeric(),
                                         median = numeric(), min = numeric(), max = numeric()))
  l_rows <- lapply(split(df_start, list(df_start$arm, df_start$sequence), drop = TRUE), function(d) {
    v_week <- rep(d$week, d$n)
    data.frame(arm = d$arm[1], sequence = d$sequence[1], n = length(v_week),
               median = stats::median(v_week), min = min(v_week), max = max(v_week))
  })
  df_out <- do.call(rbind, l_rows)
  df_out <- df_out[order(df_out$arm, df_out$sequence, method = "radix"), , drop = FALSE]
  rownames(df_out) <- NULL
  df_out
}
