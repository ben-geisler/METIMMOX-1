#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(readxl, dplyr)

#clear workspace
rm(list = ls())

#load EXCEL data
METIMMOXraw <- read_excel("data/sensitive/METIMMOX w TMB 241101.xlsx", skip = 1) %>%
  filter(!is.na(`Study1/control0`) & `Study1/control0` != "" & 
           !ID %in% c("1-030", "1-037"))

# rename some covariates and create categorical variables/factors
rm(list = setdiff(ls(), "METIMMOXraw"))

METIMMOX <- METIMMOXraw
rm(METIMMOXraw)

METIMMOX$Rx <- factor(METIMMOX$`Study1/control0`, levels = c(0, 1), labels = c("Control arm", "Experimental arm"))

METIMMOX$Female <- 1 - METIMMOX$"Sex 0female"
METIMMOX$Female <- factor(METIMMOX$Female, levels = c(0, 1), labels = c("Male", "Female"))

METIMMOX$ECOG0 <- 1- METIMMOX$S1V0
METIMMOX$ECOG0 <- factor(METIMMOX$ECOG0, levels = c(0, 1), labels = c("ECOG>1", "ECOG=0"))

METIMMOX$`RAS/RAF` <- factor(METIMMOX$`RAS/RAF`, levels = c(0, 1), labels = c("wildtype", "RAS/BRAF-mutant tumor"))

METIMMOX$RightLeft <- factor(METIMMOX$`Right0/Left`, levels = c(0, 1), labels = c("Right-sided primary tumor", "Left-sided primary tumor"))

METIMMOX$Liver <- as.factor(ifelse(grepl("^Liver", METIMMOX$`Metastatic organ`), 1, 0))

METIMMOX$InSitu <- factor (METIMMOX$`Primary tumor in situ`, levels = c(0, 1), labels = c("Primary tumor not in situ", "Primary tumour in situ"))

METIMMOX$nMets <- METIMMOX$`Number of metastatic sites`

METIMMOX$sMet <- 1- METIMMOX$`New0/Metachronous`
METIMMOX$sMet <- factor(METIMMOX$sMet, levels = c(0, 1), labels = c("Metachronous", "Synchronous"))

METIMMOX$CT1tl <- as.numeric(METIMMOX$`TL LD...130`)
METIMMOX$TLR <- METIMMOX$CT1tl / METIMMOX$`TL LD...121`
METIMMOX$TLRcat <- as.factor(ifelse(METIMMOX$TLR <= 0.9, 1, 0))

# CRP0: baseline CRP (lab block after "SEQ1 V1", cycle 1 day 1). Not used by the
# economic model or the clinical reports.
METIMMOX$CRP0 <- METIMMOX$`CRP...478`
METIMMOX$CRP0cat <- as.factor(ifelse(METIMMOX$CRP0 < 5, 1, 0))

# CRP1: week-4 CRP (lab block after "SEQ1 V3", cycle 3 day 1), measured after the
# two FLOX cycles given to both arms and before the first nivolumab dose. This is
# the CRP biomarker used throughout (data$crp in 03_biomarker_strategies.R);
# see issue #150.
METIMMOX$CRP1 <- METIMMOX$`CRP...540`
METIMMOX$CRP1cat <- as.factor(ifelse(METIMMOX$CRP1 < 5, 1, 0))

METIMMOX$TMB <- as.numeric(METIMMOX$TMB)
METIMMOX$TMBcat <- as.factor(ifelse(METIMMOX$TMB >= 9, 1, 0))

# Validate the export endpoint contract before conversion (issue #181).
endpoint_headers <- c("Days until last evaluation", "Days until progression",
                      "Progression exit", "Days until death/last follow up", "Death")
if (!identical(names(METIMMOX)[c(11, 13, 15, 18, 19)], endpoint_headers)) {
  stop("Endpoint export headers changed; check positions 11, 13, 15, 18, 19.")
}
for (nm in endpoint_headers) {
  x <- METIMMOX[[nm]]
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
    stop("Invalid endpoint column: ", nm, "; expected nonnegative finite numeric values.")
  }
  if (nm %in% c("Progression exit", "Death") && !all(x %in% c(0, 1))) {
    stop("Invalid endpoint flag: ", nm, "; expected 0/1.")
  }
}
if (any(METIMMOX$`Days until last evaluation` > METIMMOX$`Days until death/last follow up`)) {
  stop("Last evaluation exceeds death/last follow-up time.")
}
METIMMOX$LastEvalwk <- METIMMOX$`Days until last evaluation` / 7

METIMMOX$PFSmo <- METIMMOX$`Days until progression`*12/365

METIMMOX$OSmo <-METIMMOX$`Days until death/last follow up`*12/365

# NOTE (issue #149): "Days until progression" / "Progression exit" are the raw
# trial variables (time to progression, deaths censored). The analysis PFS
# endpoint (progression or death within the assessment window, issue #181) is derived in 02_setup_and_global_variables.R
# via derive_pfs_endpoint(); reports reading this RDS directly must call it too.
METIMMOX$PFSwk <- METIMMOX$`Days until progression`/7

METIMMOX$OSwk <-METIMMOX$`Days until death/last follow up`/7


saveRDS(METIMMOX, file="data/tidy/METIMMOX.rds")
