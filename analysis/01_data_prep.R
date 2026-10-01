#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(readxl, dplyr)

#clear workspace
rm(list = ls())

#load EXCEL data
export_path <- "data/sensitive/METIMMOX w TMB 241101.xlsx"
df_metimmox_raw <- read_excel(export_path, skip = 1) %>%
  filter(!is.na(`Study1/control0`) & `Study1/control0` != "" &
           !ID %in% c("1-030", "1-037"))

# rename some covariates and create categorical variables/factors
rm(list = setdiff(ls(), c("df_metimmox_raw", "export_path")))

df_metimmox <- df_metimmox_raw
rm(df_metimmox_raw)

# Validate the raw export contract before deriving anything (issues #174, #181):
# headers at the positional columns used below, top-row block labels, types,
# ranges, dates, unique ID and row count. Biomarker counts are checked before
# saving.
source("R/data_contract.R")
check_export_contract(df_metimmox,
                      top_header = names(read_excel(export_path, n_max = 0,
                                                    .name_repair = "minimal")))

df_metimmox$Rx <- factor(df_metimmox$`Study1/control0`, levels = c(0, 1), labels = c("Control arm", "Experimental arm"))

df_metimmox$Female <- 1 - df_metimmox$"Sex 0female"
df_metimmox$Female <- factor(df_metimmox$Female, levels = c(0, 1), labels = c("Male", "Female"))

df_metimmox$ECOG0 <- 1- df_metimmox$S1V0
df_metimmox$ECOG0 <- factor(df_metimmox$ECOG0, levels = c(0, 1), labels = c("ECOG>1", "ECOG=0"))

df_metimmox$`RAS/RAF` <- factor(df_metimmox$`RAS/RAF`, levels = c(0, 1), labels = c("wildtype", "RAS/BRAF-mutant tumor"))

df_metimmox$RightLeft <- factor(df_metimmox$`Right0/Left`, levels = c(0, 1), labels = c("Right-sided primary tumor", "Left-sided primary tumor"))

df_metimmox$Liver <- as.factor(ifelse(grepl("^Liver", df_metimmox$`Metastatic organ`), 1, 0))

df_metimmox$InSitu <- factor (df_metimmox$`Primary tumor in situ`, levels = c(0, 1), labels = c("Primary tumor not in situ", "Primary tumour in situ"))

df_metimmox$nMets <- df_metimmox$`Number of metastatic sites`

df_metimmox$sMet <- 1- df_metimmox$`New0/Metachronous`
df_metimmox$sMet <- factor(df_metimmox$sMet, levels = c(0, 1), labels = c("Metachronous", "Synchronous"))

# TLR: first on-treatment CT target-lesion sum over the baseline sum. CT1date is
# that scan's date, used for CT1wk in 02_setup_and_global_variables.R.
df_metimmox$CT1tl <- parse_export_numeric(df_metimmox$`TL LD...130`, "TL LD...130",
                                       export_contract()$na_tokens$`TL LD...130`)
df_metimmox$TLR <- df_metimmox$CT1tl / df_metimmox$`TL LD...121`
df_metimmox$TLRcat <- as.factor(ifelse(df_metimmox$TLR <= 0.9, 1, 0))
df_metimmox$CT1date <- df_metimmox$`Date...122`

# Both CRP columns are positional; their blocks are pinned by export_contract().
# CRP0: baseline CRP (lab block after "SEQ1 V1", cycle 1 day 1). Not used by the
# economic model or the clinical reports.
df_metimmox$CRP0 <- df_metimmox$`CRP...478`
df_metimmox$CRP0cat <- as.factor(ifelse(df_metimmox$CRP0 < 5, 1, 0))

# CRP1: week-4 CRP (lab block after "SEQ1 V3", cycle 3 day 1), measured after the
# two FLOX cycles given to both arms and before the first nivolumab dose. This is
# the CRP biomarker used throughout (data$crp in 03_biomarker_strategies.R);
# see issue #150.
df_metimmox$CRP1 <- df_metimmox$`CRP...540`
df_metimmox$CRP1cat <- as.factor(ifelse(df_metimmox$CRP1 < 5, 1, 0))

df_metimmox$TMB <- parse_export_numeric(df_metimmox$TMB, "TMB", export_contract()$na_tokens$TMB)
df_metimmox$TMBcat <- as.factor(ifelse(df_metimmox$TMB >= 9, 1, 0))

# Endpoint headers, types, ranges and 0/1 flags are checked by
# check_export_contract() above (issue #181).
df_metimmox$LastEvalwk <- df_metimmox$`Days until last evaluation` / 7

df_metimmox$PFSmo <- df_metimmox$`Days until progression`*12/365

df_metimmox$OSmo <-df_metimmox$`Days until death/last follow up`*12/365

# NOTE (issue #149): "Days until progression" / "Progression exit" are the raw
# trial variables (time to progression, deaths censored). The analysis PFS
# endpoint (progression or death within the assessment window, issue #181) is derived in 02_setup_and_global_variables.R
# via derive_pfs_endpoint(); reports reading this RDS directly must call it too.
df_metimmox$PFSwk <- df_metimmox$`Days until progression`/7

df_metimmox$OSwk <-df_metimmox$`Days until death/last follow up`/7

# Known biomarker counts (issue #174): week-4 CRP 7/35 vs 17/36, baseline CRP
# 7/35 vs 8/38 (control vs experimental), TMB/BRAF 31/69, TLR 44/68.
check_export_counts(df_metimmox)

saveRDS(df_metimmox, file="data/tidy/METIMMOX.rds")
