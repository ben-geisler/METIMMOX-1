# Raw trial-export contract, single source of truth (issues #174, #181).
# readxl repairs duplicated headers by appending the column position
# ("CRP...540"), so the biomarker and scan columns are addressed by position.
# A revised export that inserts or removes a column would silently change what
# those names mean. check_export_contract() stops analysis/01_data_prep.R
# before anything is derived unless headers, types, ranges, dates, identifiers
# and row count match; check_export_counts() then checks the derived biomarkers
# against the known trial counts. Messages report column names and aggregate
# counts only, never patient values.

#' Expected layout and counts of "METIMMOX w TMB 241101.xlsx"
#' @return List of header positions (second header row, as repaired by readxl),
#'   top-row block labels (first header row, unrepaired), value windows and
#'   expected counts.
export_contract <- function() {
  list(
    # Header row read with skip = 1; position -> exact repaired name.
    headers = c(
      `1` = "ID", `2` = "Study1/control0", `9` = "Date of inclusion",
      # Endpoints (issue #181).
      `11` = "Days until last evaluation", `13` = "Days until progression",
      `15` = "Progression exit", `18` = "Days until death/last follow up",
      `19` = "Death",
      `26` = "Mutation", `29` = "TMB", `32` = "Sex 0female", `34` = "Age",
      # Baseline radiology block: scan date, target-lesion sum (TLR denominator).
      `114` = "Date...114", `115` = "NTL baseline (n)", `120` = "TL5 LD...120",
      `121` = "TL LD...121",
      # First on-treatment CT: scan date and target-lesion sum (TLR numerator).
      `122` = "Date...122", `123` = "Overall Recist...123",
      `129` = "TL5 LD...129", `130` = "TL LD...130", `131` = "TL...131",
      # Cycle 1 day 1 lab block (baseline CRP).
      `468` = "Date...468", `477` = "Thrombocytes...477", `478` = "CRP...478",
      `479` = "Creatinine...479",
      # Cycle 3 day 1 lab block (week-4 CRP, issue #150).
      `530` = "Date...530", `539` = "Thrombocytes...539", `540` = "CRP...540",
      `541` = "Creatinine...541",
      # Treatment-administration block (issue #208): the column before it and
      # the per-patient totals after it (export_sequence_block()).
      `2638` = "END...2638", `2880` = "Total treatment",
      `2881` = "Cycles Nivolumab", `2882` = "Cycles FLOX",
      sequence_block_headers()),
    # First header row (merged block labels). Its labels sit ahead of the blocks
    # they name; at these positions they fingerprint the block layout, so a whole
    # radiology or lab block inserted or removed (which would keep "CRP...540"
    # but give it another visit's value) is caught.
    top_labels = c(`114` = "Radiology 3", `127` = "Radiology 4",
                   `434` = "SEQ1 V1", `496` = "SEQ1 V3",
                   sequence_block_top_labels()),
    n_rows = 74L,
    arms = c(0, 1),
    mutations = c("BRAF", "KRAS", "NRAS", "wt"),
    na_tokens = list(`TL LD...130` = "NA", TMB = "Missing"),
    trial_period = as.Date(c("2018-01-01", "2023-12-31")),
    # Days from inclusion: baseline scan per patient; lab-block medians. The
    # median windows separate adjacent lab blocks (observed medians: screening
    # 0, cycle 1 day 1 5.5, cycle 2 21, cycle 3 day 1 37, cycle 4 52 days).
    baseline_scan_days = c(-90, 60),
    lab_median_days = list(`Date...468` = c(2, 14), `Date...530` = c(28, 45)),
    # Weeks from inclusion to the first on-treatment CT (scheduled at week 8).
    first_ct_weeks = c(4, 20),
    # Positives / classified patients by arm (0 = control, 1 = experimental).
    counts = list(
      crp_week4 = list(pos = c(7, 17), n = c(35, 36)),
      crp_baseline = list(pos = c(7, 8), n = c(35, 38)),
      tmb_braf = list(pos = 31, n = 69),
      tlr = list(pos = 44, n = 68),
      # Patients with at least one administration in treatment sequence k of
      # the administration block (issue #208), by arm.
      sequence_2 = list(pos = c(16, 20), n = c(36, 38)),
      sequence_3 = list(pos = c(5, 13), n = c(36, 38)),
      sequence_4 = list(pos = c(0, 4), n = c(36, 38))),
    # Patients whose administration-block counts equal the export's totals
    # "Cycles FLOX" and "Cycles Nivolumab". The totals count only the exact
    # spellings: one "Oxalliplatin" and one "Nivoliumab" visit are
    # administrations in the block but not in the totals, and one further
    # patient has one nivolumab administration more in the block.
    cycles_agreement = list(FLOX = 73, nivolumab = 72),
    # Allowed non-date entries of the block's date columns and the drug tokens
    # of its first drug column (spelling variants as exported, lower case).
    sequence_date_tokens = c("Pause", "Not given"),
    sequence_treatments = list(FLOX = "^oxal+iplatin$", nivolumab = "^nivoli?umab$",
                               none = "^not given$"),
    # The "Dates" sheet's sequence columns, whose entry counts the trial team
    # was asked about (15 October 2026, SEQ3); see check_dates_sheet().
    dates_sheet = list(
      headers = c(`1` = "ID", `2` = "Study1/control0", `6` = "Date of inclusion",
                  `58` = "Start first cycle", `59` = "End first cycle",
                  `60` = "Start second cycle", `61` = "End second cycle",
                  `62` = "Start third cycle", `63` = "End third cycle",
                  `64` = "Start fourth cycle", `65` = "End fourth cycle"),
      entries = c(74, 37, 22, 17)))
}

# ---------------------------------------------------------------------------
# Treatment-administration block (issue #208)
# ---------------------------------------------------------------------------
# The Datasheet records every administration of up to six treatment sequences
# of eight visits each (first-row labels "S1V1" ... "S6V8"). Each visit has
# five columns: date, first drug (oxaliplatin or nivolumab, or "Not given"),
# dose, second drug (5-FU), dose. The block starts at column 2640; as
# elsewhere in the first header row, its labels sit 34 columns ahead of the
# columns they name. The "Dates" sheet's "Start/End first..fourth cycle"
# columns are formulas that read this block (check_dates_sheet()).

#' Layout of the treatment-administration block
#' @return List: first_column, columns_per_visit, visits_per_sequence,
#'   n_sequences, top_label_offset.
export_sequence_block <- function() {
  list(first_column = 2640L, columns_per_visit = 5L, visits_per_sequence = 8L,
       n_sequences = 6L, top_label_offset = 34L)
}

#' Column position of a visit's date in the administration block
#' @param sequence Treatment sequence (1-6).
#' @param visit Visit within the sequence (1-8).
#' @param block export_sequence_block().
#' @return Integer column position (the first drug is the next column).
sequence_block_column <- function(sequence, visit, block = export_sequence_block()) {
  as.integer(block$first_column + block$columns_per_visit *
               ((sequence - 1L) * block$visits_per_sequence + (visit - 1L)))
}

#' Pinned headers of the administration block (date and both drug columns)
#' @param block export_sequence_block().
#' @return Named character vector, position -> repaired header.
sequence_block_headers <- function(block = export_sequence_block()) {
  v_date <- as.vector(outer(seq_len(block$visits_per_sequence), seq_len(block$n_sequences),
                            function(v, k) sequence_block_column(k, v, block)))
  v_pos <- sort(c(v_date, v_date + 1L, v_date + 3L))
  v_name <- c("Date", "Treatment", NA, "Treatment")[(v_pos - block$first_column) %%
                                                     block$columns_per_visit + 1L]
  stats::setNames(paste0(v_name, "...", v_pos), v_pos)
}

#' Pinned first-row labels of the administration block ("S1V1" ... "S6V8")
#' @param block export_sequence_block().
#' @return Named character vector, position -> label.
sequence_block_top_labels <- function(block = export_sequence_block()) {
  df_visits <- expand.grid(visit = seq_len(block$visits_per_sequence),
                           sequence = seq_len(block$n_sequences))
  v_pos <- sequence_block_column(df_visits$sequence, df_visits$visit, block) -
    block$top_label_offset
  stats::setNames(paste0("S", df_visits$sequence, "V", df_visits$visit), v_pos)
}

#' Parse a date column of the administration block or the Dates sheet
#'
#' readxl returns POSIXct for an all-date column and character (Excel serial
#' numbers as text) for a column that also holds a token such as "Pause".
#' @param x Column as read by readxl.
#' @param name Column name, for the error message.
#' @param tokens Allowed non-date entries (read as no date).
#' @return Date vector; stops on any other entry.
parse_export_date <- function(x, name, tokens = c("Pause", "Not given")) {
  if (inherits(x, c("POSIXct", "Date"))) return(as.Date(x))
  if (all(is.na(x))) return(as.Date(rep(NA_real_, length(x))))
  if (!is.character(x)) stop("Export contract: ", name, " is neither a date nor a text column.",
                             call. = FALSE)
  x <- trimws(x)
  v_serial <- suppressWarnings(as.numeric(x))
  v_bad <- !is.na(x) & is.na(v_serial) & !x %in% tokens
  if (any(v_bad)) {
    stop("Export contract: ", name, " has ", sum(v_bad), " entr",
         if (sum(v_bad) == 1) "y" else "ies", " that are neither dates nor ",
         paste(dQuote(tokens, FALSE), collapse = "/"), ".", call. = FALSE)
  }
  as.Date(v_serial, origin = "1899-12-30")
}

#' Classify the first drug column of a visit
#' @param x Treatment column.
#' @param name Column name, for the error message.
#' @param patterns Named regular expressions (export_contract()$sequence_treatments),
#'   matched against the lower-case entry.
#' @return Character vector "FLOX", "nivolumab", "none" or NA (empty); stops
#'   on an unknown entry.
classify_export_treatment <- function(x, name, patterns = export_contract()$sequence_treatments) {
  x <- tolower(trimws(as.character(x)))
  v_out <- rep(NA_character_, length(x))
  for (nm in names(patterns)) v_out[!is.na(x) & grepl(patterns[[nm]], x)] <- nm
  v_bad <- !is.na(x) & is.na(v_out)
  if (any(v_bad)) {
    stop("Export contract: ", name, " has ", sum(v_bad), " unknown treatment entr",
         if (sum(v_bad) == 1) "y" else "ies", ".", call. = FALSE)
  }
  v_out
}

#' Administrations of every treatment sequence, per patient
#'
#' An administration is a visit with a date and oxaliplatin (FLOX) or
#' nivolumab in its first drug column; a visit marked "Not given", without a
#' drug or without a date is not one.
#' @param data Export read with skip = 1 (layout checked).
#' @param contract export_contract().
#' @param block export_sequence_block().
#' @return List of patient x sequence matrices: `n_FLOX`, `n_nivolumab`,
#'   `first` and `last` (first and last administration date as days since
#'   1970-01-01, NA without administration) and `drug_without_date` (visits
#'   with a drug but no date).
export_administrations <- function(data, contract = export_contract(),
                                   block = export_sequence_block()) {
  n <- nrow(data)
  k_max <- block$n_sequences
  m_zero <- matrix(0, n, k_max, dimnames = list(NULL, paste0("S", seq_len(k_max))))
  out <- list(n_FLOX = m_zero, n_nivolumab = m_zero, first = m_zero * NA_real_,
              last = m_zero * NA_real_, drug_without_date = m_zero)
  for (k in seq_len(k_max)) for (v in seq_len(block$visits_per_sequence)) {
    j <- sequence_block_column(k, v, block)
    v_date <- parse_export_date(data[[j]], names(data)[j], contract$sequence_date_tokens)
    v_drug <- classify_export_treatment(data[[j + 1L]], names(data)[j + 1L],
                                        contract$sequence_treatments)
    v_is_drug <- v_drug %in% c("FLOX", "nivolumab")
    v_given <- !is.na(v_date) & v_is_drug
    out$drug_without_date[, k] <- out$drug_without_date[, k] + (is.na(v_date) & v_is_drug)
    out$n_FLOX[, k] <- out$n_FLOX[, k] + (v_given & v_drug == "FLOX")
    out$n_nivolumab[, k] <- out$n_nivolumab[, k] + (v_given & v_drug == "nivolumab")
    v_day <- ifelse(v_given, as.numeric(v_date), NA_real_)
    out$first[, k] <- pmin(out$first[, k], v_day, na.rm = TRUE)
    out$last[, k] <- pmax(out$last[, k], v_day, na.rm = TRUE)
  }
  out
}

#' Check the administration block: dates, drugs and sequence order
#' @param data Export read with skip = 1 (layout checked).
#' @param contract export_contract().
#' @return Character vector of violations (empty when the block is valid).
sequence_block_problems <- function(data, contract = export_contract()) {
  adm <- tryCatch(export_administrations(data, contract), error = identity)
  if (inherits(adm, "error")) return(conditionMessage(adm))
  v_problems <- character()
  v_period <- as.numeric(contract$trial_period)
  v_days <- c(adm$first, adm$last)
  if (any(!is.na(v_days) & (v_days < v_period[1] | v_days > v_period[2])))
    v_problems <- c(v_problems, paste0("administration block: ",
                                       sum(!is.na(v_days) & (v_days < v_period[1] | v_days > v_period[2])),
                                       " first/last date(s) outside the trial period"))
  if (sum(adm$drug_without_date) > 0)
    v_problems <- c(v_problems, paste0("administration block: ", sum(adm$drug_without_date),
                                       " visit(s) with a drug but no date"))
  if (any(adm$n_FLOX[, 1] + adm$n_nivolumab[, 1] == 0))
    v_problems <- c(v_problems, "administration block: every patient must start sequence 1")
  for (k in seq_len(ncol(adm$first))[-1]) {
    v_started <- !is.na(adm$first[, k])
    if (any(v_started & is.na(adm$first[, k - 1L])))
      v_problems <- c(v_problems, paste0("administration block: sequence ", k,
                                         " started without sequence ", k - 1L))
    if (any(v_started & !is.na(adm$last[, k - 1L]) & adm$first[, k] <= adm$last[, k - 1L]))
      v_problems <- c(v_problems, paste0("administration block: sequence ", k,
                                         " starts before sequence ", k - 1L, " ends"))
  }
  v_inclusion <- as.numeric(as.Date(data$`Date of inclusion`))
  if (any(adm$first[, 1] < v_inclusion, na.rm = TRUE))
    v_problems <- c(v_problems, "administration block: sequence 1 starts before inclusion")
  for (nm in c("Cycles FLOX", "Cycles Nivolumab")) {
    x <- data[[nm]]
    if (!is.numeric(x) || anyNA(x) || any(x < 0 | x != round(x)))
      v_problems <- c(v_problems, paste0(nm, " must be complete nonnegative whole numbers"))
  }
  v_problems
}

#' Per-patient treatment sequences 2-4 for the model (issue #208)
#'
#' For sequence k: `SeqK` (1 if the patient had at least one administration
#' in sequence k), `SeqKwk` (weeks from inclusion to its first
#' administration), `SeqK_FLOX` and `SeqK_nivo` (administrations delivered).
#' Sequences 5 and 6 are not modelled (trial team, SEQ3: "nothing beyond");
#' `AdmFLOX` and `AdmNivo` count the administrations of all six sequences for
#' the cross-check against the export's "Cycles FLOX" and "Cycles Nivolumab".
#' @param data Export read with skip = 1 (layout checked).
#' @param sequences Sequences returned.
#' @return Data frame aligned with the rows of `data`.
derive_treatment_sequences <- function(data, sequences = 2:4) {
  adm <- export_administrations(data)
  v_inclusion <- as.numeric(as.Date(data$`Date of inclusion`))
  df_out <- data.frame(row.names = seq_len(nrow(data)))
  for (k in sequences) {
    v_n <- adm$n_FLOX[, k] + adm$n_nivolumab[, k]
    df_out[[paste0("Seq", k)]] <- as.numeric(v_n > 0)
    df_out[[paste0("Seq", k, "wk")]] <- unname((adm$first[, k] - v_inclusion) / 7)
    df_out[[paste0("Seq", k, "_FLOX")]] <- unname(adm$n_FLOX[, k])
    df_out[[paste0("Seq", k, "_nivo")]] <- unname(adm$n_nivolumab[, k])
  }
  df_out$AdmFLOX <- unname(rowSums(adm$n_FLOX))
  df_out$AdmNivo <- unname(rowSums(adm$n_nivolumab))
  rownames(df_out) <- NULL
  df_out
}

#' Check the "Dates" sheet against the administration block (issue #208)
#'
#' The sheet's "Start/End first..fourth cycle" cells are formulas that read
#' the administration block, so each date they hold must equal the block's
#' date of the referenced visit; for the start columns that is the sequence's
#' first visit. A cell that refers to an empty block cell reads as Excel day
#' zero and is an entry without a date. The check pins the headers, the
#' patients, arms and inclusion dates and the entry counts the trial team was
#' asked about, and returns aggregate counts of how the sheet's entries relate
#' to the administrations (no patient values).
#' @param dates "Dates" sheet (header in the first row), analysis filter
#'   applied.
#' @param data Export read with skip = 1 (layout checked), same patients.
#' @param contract export_contract().
#' @return Data frame per sequence 1-4: `entries` (non-empty start cells),
#'   `real_dates` (start cells holding a date), `with_administration`
#'   (entries whose patient has an administration in that sequence) and
#'   `administered` (patients with an administration in that sequence); stops
#'   listing every violation.
check_dates_sheet <- function(dates, data, contract = export_contract()) {
  spec <- contract$dates_sheet
  v_problems <- character()
  flag <- function(...) v_problems <<- c(v_problems, paste0(...))
  v_pos <- as.integer(names(spec$headers))
  v_found <- names(dates)[v_pos]
  for (i in which(is.na(v_found) | v_found != spec$headers))
    flag("Dates sheet column ", v_pos[i], " is '", v_found[i], "', expected '",
         spec$headers[i], "'")
  if (length(v_problems)) stop("Export contract (Dates sheet layout) violated:\n  - ",
                               paste(v_problems, collapse = "\n  - "), call. = FALSE)
  v_row <- match(data$ID, dates$ID)
  if (nrow(dates) != nrow(data) || anyNA(v_row) || anyDuplicated(dates$ID))
    stop("Export contract: the Dates sheet does not hold exactly the analysed patients.",
         call. = FALSE)
  dates <- dates[v_row, , drop = FALSE]
  if (!identical(as.numeric(dates$`Study1/control0`), as.numeric(data$`Study1/control0`)))
    flag("Dates sheet arms differ from the Datasheet")
  if (!identical(as.Date(dates$`Date of inclusion`), as.Date(data$`Date of inclusion`)))
    flag("Dates sheet inclusion dates differ from the Datasheet")
  adm <- export_administrations(data, contract)
  v_start <- spec$headers[c("58", "60", "62", "64")]
  v_entries <- v_real <- v_with <- v_admin <- numeric(4)
  for (k in 1:4) {
    v_raw <- dates[[v_start[k]]]
    v_date <- parse_export_date(v_raw, v_start[k], contract$sequence_date_tokens)
    v_is_date <- !is.na(v_date) & v_date >= contract$trial_period[1]
    j <- sequence_block_column(k, 1L)
    v_block <- parse_export_date(data[[j]], names(data)[j], contract$sequence_date_tokens)
    v_differs <- v_is_date & (is.na(v_block) | v_date != v_block)
    if (any(v_differs))
      flag(v_start[k], ": ", sum(v_differs), " date(s) differ from the first visit of sequence ", k)
    v_has_admin <- (adm$n_FLOX[, k] + adm$n_nivolumab[, k]) > 0
    v_entries[k] <- sum(!is.na(v_raw))
    v_real[k] <- sum(v_is_date)
    v_with[k] <- sum(!is.na(v_raw) & v_has_admin)
    v_admin[k] <- sum(v_has_admin)
  }
  if (!identical(v_entries, as.numeric(spec$entries)))
    flag("Dates sheet start entries ", paste(v_entries, collapse = "/"), ", expected ",
         paste(spec$entries, collapse = "/"))
  if (length(v_problems)) stop("Export contract (Dates sheet) violated:\n  - ",
                               paste(v_problems, collapse = "\n  - "), call. = FALSE)
  data.frame(sequence = 1:4, entries = v_entries, real_dates = v_real,
             with_administration = v_with, administered = v_admin)
}

#' Convert an export column to numeric, allowing only listed missing tokens
#' @param x Column as read by readxl (numeric or character).
#' @param name Column name, for the error message.
#' @param na_tokens Strings that mean missing (e.g. "NA", "Missing").
#' @return Numeric vector; stops on any other non-numeric entry.
parse_export_numeric <- function(x, name, na_tokens = character()) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  v_missing <- is.na(x) | x %in% na_tokens
  v_out <- suppressWarnings(as.numeric(x))
  v_bad <- !v_missing & is.na(v_out)
  if (any(v_bad)) {
    stop("Export contract: ", name, " has ", sum(v_bad), " non-numeric entr",
         if (sum(v_bad) == 1) "y" else "ies", " other than ",
         paste(dQuote(na_tokens, FALSE), collapse = "/"), ".", call. = FALSE)
  }
  v_out[v_missing] <- NA_real_
  v_out
}

#' Check headers, types, ranges, dates, identifiers and row count
#' @param data Export read with skip = 1 and the analysis filter applied.
#' @param top_header Names of the first header row
#'   (read_excel(path, n_max = 0, .name_repair = "minimal")).
#' @param contract export_contract().
#' @return data, invisibly; stops listing every violation.
check_export_contract <- function(data, top_header, contract = export_contract()) {
  v_problems <- character()
  flag <- function(...) v_problems <<- c(v_problems, paste0(...))

  v_pos <- as.integer(names(contract$headers))
  v_found <- names(data)[v_pos]
  v_wrong <- is.na(v_found) | v_found != contract$headers
  for (i in which(v_wrong)) flag("column ", v_pos[i], " is '", v_found[i], "', expected '",
                                 contract$headers[i], "'")
  v_tpos <- as.integer(names(contract$top_labels))
  v_tfound <- top_header[v_tpos]
  v_twrong <- is.na(v_tfound) | v_tfound != contract$top_labels
  for (i in which(v_twrong)) flag("top header at column ", v_tpos[i], " is '", v_tfound[i],
                                  "', expected '", contract$top_labels[i], "'")
  # Values are meaningless under a shifted layout; report the layout first.
  if (length(v_problems)) {
    stop("Export contract (layout) violated:\n  - ",
         paste(v_problems, collapse = "\n  - "), call. = FALSE)
  }

  if (nrow(data) != contract$n_rows)
    flag(nrow(data), " rows after filtering, expected ", contract$n_rows)
  if (anyNA(data$ID) || anyDuplicated(data$ID))
    flag("ID must be non-missing and unique (", sum(is.na(data$ID)), " missing, ",
         sum(duplicated(data$ID)), " duplicated)")
  if (!all(data$`Study1/control0` %in% contract$arms)) flag("Study1/control0 must be 0/1")

  # Endpoints (issue #181).
  v_endpoints <- c("Days until last evaluation", "Days until progression",
                   "Progression exit", "Days until death/last follow up", "Death")
  for (nm in v_endpoints) {
    x <- data[[nm]]
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
      flag(nm, ": expected nonnegative finite numeric values")
    } else if (nm %in% c("Progression exit", "Death") && !all(x %in% c(0, 1))) {
      flag(nm, ": expected 0/1")
    }
  }
  if (all(vapply(data[v_endpoints], is.numeric, logical(1))) &&
      any(data$`Days until last evaluation` > data$`Days until death/last follow up`,
          na.rm = TRUE)) {
    flag("last evaluation exceeds death/last follow-up time")
  }

  # Dates.
  is_date <- function(x) inherits(x, c("POSIXct", "Date"))
  v_date_cols <- c("Date of inclusion", "Date...114", "Date...122",
                   names(contract$lab_median_days))
  for (nm in v_date_cols) if (!is_date(data[[nm]])) flag(nm, " is not a date column")
  if (length(v_problems)) {
    stop("Export contract violated:\n  - ", paste(v_problems, collapse = "\n  - "),
         call. = FALSE)
  }
  v_inclusion <- as.Date(data$`Date of inclusion`)
  if (anyNA(v_inclusion)) flag("Date of inclusion has missing values")
  for (nm in v_date_cols) {
    d <- as.Date(data[[nm]])
    v_out <- !is.na(d) & (d < contract$trial_period[1] | d > contract$trial_period[2])
    if (any(v_out)) flag(nm, ": ", sum(v_out), " date(s) outside ",
                         paste(format(contract$trial_period), collapse = " to "))
  }
  days <- function(nm) as.numeric(as.Date(data[[nm]]) - v_inclusion)
  v_base_days <- days("Date...114")
  if (anyNA(v_base_days) || any(v_base_days < contract$baseline_scan_days[1] |
                                v_base_days > contract$baseline_scan_days[2], na.rm = TRUE))
    flag("baseline scan (Date...114) missing or outside ",
         paste(contract$baseline_scan_days, collapse = " to "), " days of inclusion")
  v_ct1_weeks <- days("Date...122") / 7
  if (any(v_ct1_weeks <= contract$first_ct_weeks[1] |
          v_ct1_weeks > contract$first_ct_weeks[2], na.rm = TRUE))
    flag("first on-treatment CT (Date...122) outside (",
         paste(contract$first_ct_weeks, collapse = ", "), "] weeks after inclusion")
  for (nm in names(contract$lab_median_days)) {
    med <- stats::median(days(nm), na.rm = TRUE)
    v_win <- contract$lab_median_days[[nm]]
    if (is.na(med) || med < v_win[1] || med > v_win[2])
      flag(nm, ": median ", med, " days from inclusion, expected ",
           v_win[1], " to ", v_win[2], " (wrong visit block?)")
  }

  # Target-lesion sums (TLR) and CRP.
  v_tl_base <- data$`TL LD...121`
  if (!is.numeric(v_tl_base) || anyNA(v_tl_base) || any(v_tl_base <= 0))
    flag("TL LD...121 (baseline target-lesion sum) must be positive and complete")
  v_tl_ct1 <- tryCatch(parse_export_numeric(data$`TL LD...130`, "TL LD...130",
                                            contract$na_tokens$`TL LD...130`),
                       error = function(e) { flag(conditionMessage(e)); NULL })
  if (!is.null(v_tl_ct1)) {
    if (any(v_tl_ct1 < 0, na.rm = TRUE)) flag("TL LD...130 has negative values")
    if (!identical(is.na(v_tl_ct1), is.na(data$`Date...122`)))
      flag("TL LD...130 and Date...122 (first on-treatment CT) are not missing together")
  }
  for (nm in c("CRP...478", "CRP...540")) {
    x <- data[[nm]]
    if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < 0)))
      flag(nm, " must be nonnegative numeric or missing")
  }

  # TMB/BRAF.
  v_tmb <- tryCatch(parse_export_numeric(data$TMB, "TMB", contract$na_tokens$TMB),
                    error = function(e) { flag(conditionMessage(e)); NULL })
  if (!is.null(v_tmb) && any(v_tmb < 0, na.rm = TRUE)) flag("TMB has negative values")
  if (anyNA(data$Mutation) || !all(data$Mutation %in% contract$mutations))
    flag("Mutation must be one of ", paste(contract$mutations, collapse = ", "))

  # Treatment-administration block (issue #208).
  for (msg in sequence_block_problems(data, contract)) flag(msg)

  if (length(v_problems)) {
    stop("Export contract violated:\n  - ", paste(v_problems, collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(data)
}

#' Check derived biomarkers against the known trial counts
#' @param data Output of analysis/01_data_prep.R before saving (Rx, CRP0cat,
#'   CRP1cat, TMBcat, Mutation, TLRcat, Seq2-Seq4 and, for the cross-check of
#'   the administration block, AdmFLOX, AdmNivo, "Cycles FLOX" and
#'   "Cycles Nivolumab").
#' @param contract export_contract().
#' @return data, invisibly; stops listing every mismatch.
check_export_counts <- function(data, contract = export_contract()) {
  v_arm <- as.integer(data$Rx) - 1L  # 0 = control, 1 = experimental
  by_arm <- function(x) list(pos = vapply(0:1, function(a) sum(x[v_arm == a] == 1, na.rm = TRUE), numeric(1)),
                             n = vapply(0:1, function(a) sum(!is.na(x[v_arm == a])), numeric(1)))
  overall <- function(x) list(pos = sum(x == 1, na.rm = TRUE), n = sum(!is.na(x)))
  v_tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))
  observed <- list(
    crp_week4 = by_arm(as.numeric(as.character(data$CRP1cat))),
    crp_baseline = by_arm(as.numeric(as.character(data$CRP0cat))),
    tmb_braf = overall(v_tmb_braf),
    tlr = overall(as.numeric(as.character(data$TLRcat))),
    sequence_2 = by_arm(data$Seq2), sequence_3 = by_arm(data$Seq3),
    sequence_4 = by_arm(data$Seq4))
  fmt <- function(x) paste0(x$pos, "/", x$n, collapse = " vs ")
  v_problems <- character()
  for (nm in names(contract$counts)) {
    if (!identical(as.numeric(unlist(observed[[nm]])), as.numeric(unlist(contract$counts[[nm]]))))
      v_problems <- c(v_problems, paste0(nm, " positives ", fmt(observed[[nm]]),
                                         ", expected ", fmt(contract$counts[[nm]])))
  }
  # Cross-check: the administration block against the export's own
  # per-patient totals (patients whose counts agree).
  v_agree <- c(FLOX = sum(data$AdmFLOX == data$`Cycles FLOX`),
               nivolumab = sum(data$AdmNivo == data$`Cycles Nivolumab`))
  if (!identical(as.numeric(v_agree), as.numeric(unlist(contract$cycles_agreement))))
    v_problems <- c(v_problems, paste0(
      "administration block agrees with 'Cycles FLOX'/'Cycles Nivolumab' for ",
      paste(v_agree, collapse = "/"), " patients, expected ",
      paste(unlist(contract$cycles_agreement), collapse = "/")))
  if (length(v_problems)) {
    stop("Export contract (biomarker and sequence counts) violated:\n  - ",
         paste(v_problems, collapse = "\n  - "), call. = FALSE)
  }
  invisible(data)
}
