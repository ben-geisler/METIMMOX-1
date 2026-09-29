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
      `541` = "Creatinine...541"),
    # First header row (merged block labels). Its labels sit ahead of the blocks
    # they name; at these positions they fingerprint the block layout, so a whole
    # radiology or lab block inserted or removed (which would keep "CRP...540"
    # but give it another visit's value) is caught.
    top_labels = c(`114` = "Radiology 3", `127` = "Radiology 4",
                   `434` = "SEQ1 V1", `496` = "SEQ1 V3"),
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
      tlr = list(pos = 44, n = 68)))
}

#' Convert an export column to numeric, allowing only listed missing tokens
#' @param x Column as read by readxl (numeric or character).
#' @param name Column name, for the error message.
#' @param na_tokens Strings that mean missing (e.g. "NA", "Missing").
#' @return Numeric vector; stops on any other non-numeric entry.
parse_export_numeric <- function(x, name, na_tokens = character()) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- trimws(as.character(x))
  missing <- is.na(x) | x %in% na_tokens
  out <- suppressWarnings(as.numeric(x))
  bad <- !missing & is.na(out)
  if (any(bad)) {
    stop("Export contract: ", name, " has ", sum(bad), " non-numeric entr",
         if (sum(bad) == 1) "y" else "ies", " other than ",
         paste(dQuote(na_tokens, FALSE), collapse = "/"), ".", call. = FALSE)
  }
  out[missing] <- NA_real_
  out
}

#' Check headers, types, ranges, dates, identifiers and row count
#' @param data Export read with skip = 1 and the analysis filter applied.
#' @param top_header Names of the first header row
#'   (read_excel(path, n_max = 0, .name_repair = "minimal")).
#' @param contract export_contract().
#' @return data, invisibly; stops listing every violation.
check_export_contract <- function(data, top_header, contract = export_contract()) {
  problems <- character()
  flag <- function(...) problems <<- c(problems, paste0(...))

  pos <- as.integer(names(contract$headers))
  found <- names(data)[pos]
  wrong <- is.na(found) | found != contract$headers
  for (i in which(wrong)) flag("column ", pos[i], " is '", found[i], "', expected '",
                               contract$headers[i], "'")
  tpos <- as.integer(names(contract$top_labels))
  tfound <- top_header[tpos]
  twrong <- is.na(tfound) | tfound != contract$top_labels
  for (i in which(twrong)) flag("top header at column ", tpos[i], " is '", tfound[i],
                                "', expected '", contract$top_labels[i], "'")
  # Values are meaningless under a shifted layout; report the layout first.
  if (length(problems)) {
    stop("Export contract (layout) violated:\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }

  if (nrow(data) != contract$n_rows)
    flag(nrow(data), " rows after filtering, expected ", contract$n_rows)
  if (anyNA(data$ID) || anyDuplicated(data$ID))
    flag("ID must be non-missing and unique (", sum(is.na(data$ID)), " missing, ",
         sum(duplicated(data$ID)), " duplicated)")
  if (!all(data$`Study1/control0` %in% contract$arms)) flag("Study1/control0 must be 0/1")

  # Endpoints (issue #181).
  endpoints <- c("Days until last evaluation", "Days until progression",
                 "Progression exit", "Days until death/last follow up", "Death")
  for (nm in endpoints) {
    x <- data[[nm]]
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
      flag(nm, ": expected nonnegative finite numeric values")
    } else if (nm %in% c("Progression exit", "Death") && !all(x %in% c(0, 1))) {
      flag(nm, ": expected 0/1")
    }
  }
  if (all(vapply(data[endpoints], is.numeric, logical(1))) &&
      any(data$`Days until last evaluation` > data$`Days until death/last follow up`,
          na.rm = TRUE)) {
    flag("last evaluation exceeds death/last follow-up time")
  }

  # Dates.
  is_date <- function(x) inherits(x, c("POSIXct", "Date"))
  date_cols <- c("Date of inclusion", "Date...114", "Date...122",
                 names(contract$lab_median_days))
  for (nm in date_cols) if (!is_date(data[[nm]])) flag(nm, " is not a date column")
  if (length(problems)) {
    stop("Export contract violated:\n  - ", paste(problems, collapse = "\n  - "),
         call. = FALSE)
  }
  inclusion <- as.Date(data$`Date of inclusion`)
  if (anyNA(inclusion)) flag("Date of inclusion has missing values")
  for (nm in date_cols) {
    d <- as.Date(data[[nm]])
    out <- !is.na(d) & (d < contract$trial_period[1] | d > contract$trial_period[2])
    if (any(out)) flag(nm, ": ", sum(out), " date(s) outside ",
                       paste(format(contract$trial_period), collapse = " to "))
  }
  days <- function(nm) as.numeric(as.Date(data[[nm]]) - inclusion)
  base_days <- days("Date...114")
  if (anyNA(base_days) || any(base_days < contract$baseline_scan_days[1] |
                              base_days > contract$baseline_scan_days[2], na.rm = TRUE))
    flag("baseline scan (Date...114) missing or outside ",
         paste(contract$baseline_scan_days, collapse = " to "), " days of inclusion")
  ct1_weeks <- days("Date...122") / 7
  if (any(ct1_weeks <= contract$first_ct_weeks[1] |
          ct1_weeks > contract$first_ct_weeks[2], na.rm = TRUE))
    flag("first on-treatment CT (Date...122) outside (",
         paste(contract$first_ct_weeks, collapse = ", "), "] weeks after inclusion")
  for (nm in names(contract$lab_median_days)) {
    med <- stats::median(days(nm), na.rm = TRUE)
    win <- contract$lab_median_days[[nm]]
    if (is.na(med) || med < win[1] || med > win[2])
      flag(nm, ": median ", med, " days from inclusion, expected ",
           win[1], " to ", win[2], " (wrong visit block?)")
  }

  # Target-lesion sums (TLR) and CRP.
  tl_base <- data$`TL LD...121`
  if (!is.numeric(tl_base) || anyNA(tl_base) || any(tl_base <= 0))
    flag("TL LD...121 (baseline target-lesion sum) must be positive and complete")
  tl_ct1 <- tryCatch(parse_export_numeric(data$`TL LD...130`, "TL LD...130",
                                          contract$na_tokens$`TL LD...130`),
                     error = function(e) { flag(conditionMessage(e)); NULL })
  if (!is.null(tl_ct1)) {
    if (any(tl_ct1 < 0, na.rm = TRUE)) flag("TL LD...130 has negative values")
    if (!identical(is.na(tl_ct1), is.na(data$`Date...122`)))
      flag("TL LD...130 and Date...122 (first on-treatment CT) are not missing together")
  }
  for (nm in c("CRP...478", "CRP...540")) {
    x <- data[[nm]]
    if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < 0)))
      flag(nm, " must be nonnegative numeric or missing")
  }

  # TMB/BRAF.
  tmb <- tryCatch(parse_export_numeric(data$TMB, "TMB", contract$na_tokens$TMB),
                  error = function(e) { flag(conditionMessage(e)); NULL })
  if (!is.null(tmb) && any(tmb < 0, na.rm = TRUE)) flag("TMB has negative values")
  if (anyNA(data$Mutation) || !all(data$Mutation %in% contract$mutations))
    flag("Mutation must be one of ", paste(contract$mutations, collapse = ", "))

  if (length(problems)) {
    stop("Export contract violated:\n  - ", paste(problems, collapse = "\n  - "),
         call. = FALSE)
  }
  invisible(data)
}

#' Check derived biomarkers against the known trial counts
#' @param data Output of analysis/01_data_prep.R before saving (Rx, CRP0cat,
#'   CRP1cat, TMBcat, Mutation, TLRcat).
#' @param contract export_contract().
#' @return data, invisibly; stops listing every mismatch.
check_export_counts <- function(data, contract = export_contract()) {
  arm <- as.integer(data$Rx) - 1L  # 0 = control, 1 = experimental
  by_arm <- function(x) list(pos = vapply(0:1, function(a) sum(x[arm == a] == 1, na.rm = TRUE), numeric(1)),
                             n = vapply(0:1, function(a) sum(!is.na(x[arm == a])), numeric(1)))
  overall <- function(x) list(pos = sum(x == 1, na.rm = TRUE), n = sum(!is.na(x)))
  tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))
  observed <- list(
    crp_week4 = by_arm(as.numeric(as.character(data$CRP1cat))),
    crp_baseline = by_arm(as.numeric(as.character(data$CRP0cat))),
    tmb_braf = overall(tmb_braf),
    tlr = overall(as.numeric(as.character(data$TLRcat))))
  fmt <- function(x) paste0(x$pos, "/", x$n, collapse = " vs ")
  problems <- character()
  for (nm in names(contract$counts)) {
    if (!identical(as.numeric(unlist(observed[[nm]])), as.numeric(unlist(contract$counts[[nm]]))))
      problems <- c(problems, paste0(nm, " positives ", fmt(observed[[nm]]),
                                     ", expected ", fmt(contract$counts[[nm]])))
  }
  if (length(problems)) {
    stop("Export contract (biomarker counts) violated:\n  - ",
         paste(problems, collapse = "\n  - "), call. = FALSE)
  }
  invisible(data)
}
