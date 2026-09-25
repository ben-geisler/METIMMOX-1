# Independent export cross-check. Persist aggregates only, never patient rows.
source("R/pfs_endpoint.R")
raw <- suppressMessages(readxl::read_excel(
  "data/sensitive/METIMMOX w TMB 241101.xlsx", sheet = "Datasheet",
  range = cellranger::cell_limits(c(2, 1), c(NA, 21))))
raw <- raw[!is.na(raw$`Study1/control0`) & raw$`Study1/control0` != "", ]
# Match the two protocol exclusions already applied by analysis/01_data_prep.R.
raw <- raw[!raw$ID %in% c("1-030", "1-037"), ]
stopifnot(identical(names(raw)[c(18, 19, 20, 21)],
                   c("Days until death/last follow up", "Death", "SAM OS", "SAM PFS")))
parse_sam <- function(x, column) {
  pattern <- "^\\s*\\(\\s*([0-9]+(?:\\.[0-9]+)?)\\s*,\\s*([01])\\s*\\)\\s*$"
  ok <- !is.na(x) & grepl(pattern, x, perl = TRUE)
  if (!all(ok)) stop(column, ": ", sum(!ok), " missing or malformed tuples.")
  data.frame(days = as.numeric(sub(pattern, "\\1", x, perl = TRUE)),
             event = as.numeric(sub(pattern, "\\2", x, perl = TRUE)))
}
sam_os <- parse_sam(raw$`SAM OS`, "SAM OS")
sam_pfs <- parse_sam(raw$`SAM PFS`, "SAM PFS")
# Check the example supplied by the user without printing a patient record.
example <- which(raw$ID == "1-002")
stopifnot(length(example) == 1L, sam_os$days[example] == 320,
          sam_os$event[example] == 1, sam_pfs$days[example] == 112,
          sam_pfs$event[example] == 1)

source("analysis/02_setup_and_global_variables.R")
source("analysis/03_biomarker_strategies.R")
primary <- data_complete
example_analysis <- which(primary$ID == "1-002")
stopifnot(length(example_analysis) == 1L,
          abs(primary$OSwk[example_analysis] * 7 - 320) < 1e-7,
          primary$Death[example_analysis] == 1,
          abs(primary$PFSwk[example_analysis] * 7 -
                raw$`Days until progression`[example]) < 1e-7,
          primary$Progression[example_analysis] == 1)
cat("Supplied example: derived PFS date agrees with SAM:",
    abs(primary$PFSwk[example_analysis] * 7 - sam_pfs$days[example]) < 1e-7,
    "; SAM PFS date equals last evaluation:",
    abs(primary$LastEvalwk[example_analysis] * 7 - sam_pfs$days[example]) < 1e-7, "\n")
old <- derive_pfs_endpoint(primary, Inf, "TTPwk", verbose = FALSE)
alternative <- derive_pfs_endpoint(primary, 16, "LastEvalwk", verbose = FALSE)
index <- match(primary$ID, raw$ID)
stopifnot(!anyNA(index), !anyDuplicated(raw$ID))
compare <- function(label, reference, days, event) {
  same_time <- abs(reference$days - days) < 1e-7
  same_event <- reference$event == event
  stopifnot(!anyNA(same_time), !anyNA(same_event))
  data.frame(Comparison = label, N = length(days),
    Agree = sum(same_time & same_event),
    Time_only_differs = sum(!same_time & same_event),
    Event_only_differs = sum(same_time & !same_event),
    Both_differ = sum(!same_time & !same_event))
}
comparisons <- rbind(
  compare("SAM OS vs export OS (eligible export rows)", sam_os,
          raw$`Days until death/last follow up`, raw$Death),
  compare("SAM OS vs analysis OS", sam_os[index, ], primary$OSwk * 7, primary$Death),
  compare("SAM PFS vs primary (16 weeks, TTPwk)", sam_pfs[index, ],
          primary$PFSwk * 7, primary$Progression),
  compare("SAM PFS vs pre-181 (all deaths, TTPwk)", sam_pfs[index, ],
          old$PFSwk * 7, old$Progression),
  compare("SAM PFS vs alternative (16 weeks, LastEvalwk)", sam_pfs[index, ],
          alternative$PFSwk * 7, alternative$Progression))

sam <- sam_pfs[index, ]
no_progression_deaths <- primary$Death == 1 & primary$ProgressionExit == 0
progressors <- primary$ProgressionExit == 1
progression_date_differs <- progressors & abs(sam$days - primary$TTPwk * 7) > 1e-7
sam_censored_deaths <- raw$Death == 1 & sam_pfs$event == 0
details <- data.frame(
  Finding = c("Complete-case deaths without recorded progression",
    "Of those, SAM PFS censored", "Of those, primary rule censored",
    "Of those, both SAM and primary censored",
    "Of those, SAM event but primary censored",
    "Of those, SAM censored but primary event",
    "Recorded progressors with different SAM PFS dates",
    "Of those, SAM date equals last evaluation",
    "Eligible-export SAM PFS censored deaths",
    "Of those, SAM date equals last evaluation",
    "Eligible-export deaths censored by SAM OS",
    "Eligible-export non-deaths counted as SAM OS deaths",
    "SAM OS dates earlier than exported OS date",
    "SAM OS dates later than exported OS date"),
  Count = c(sum(no_progression_deaths),
    sum(no_progression_deaths & sam$event == 0),
    sum(no_progression_deaths & primary$Progression == 0),
    sum(no_progression_deaths & sam$event == 0 & primary$Progression == 0),
    sum(no_progression_deaths & sam$event == 1 & primary$Progression == 0),
    sum(no_progression_deaths & sam$event == 0 & primary$Progression == 1),
    sum(progression_date_differs),
    sum(progression_date_differs & abs(sam$days - primary$LastEvalwk * 7) < 1e-7),
    sum(sam_censored_deaths),
    sum(sam_censored_deaths &
      abs(sam_pfs$days - raw$`Days until last evaluation`) < 1e-7),
    sum(raw$Death == 1 & sam_os$event == 0),
    sum(raw$Death == 0 & sam_os$event == 1),
    sum(sam_os$days < raw$`Days until death/last follow up` - 1e-7),
    sum(sam_os$days > raw$`Days until death/last follow up` + 1e-7)))
print(comparisons, row.names = FALSE)
print(details, row.names = FALSE)
stopifnot(comparisons$Agree[1] == 71, comparisons$N[1] == 74,
          comparisons$Agree[4] == 42, comparisons$N[4] == 68)
path <- "validation/issue181_2026-09-24/"
write.csv(comparisons, paste0(path, "sam_endpoint_comparison.csv"), row.names = FALSE)
write.csv(details, paste0(path, "sam_reconciliation_counts.csv"), row.names = FALSE)
cat("PASS: SAM tuple parsing, supplied example, and aggregate endpoint cross-check.\n")
