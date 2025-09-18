#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(readxl)

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

METIMMOX$CRP0 <- METIMMOX$`CRP...478`
METIMMOX$CRP0cat <- as.factor(ifelse(METIMMOX$CRP0 < 5, 1, 0))

METIMMOX$CRP1 <- METIMMOX$`CRP...540`
METIMMOX$CRP1cat <- as.factor(ifelse(METIMMOX$CRP1 < 5, 1, 0))

METIMMOX$TMB <- as.numeric(METIMMOX$TMB)
METIMMOX$TMBcat <- as.factor(ifelse(METIMMOX$TMB >= 9, 1, 0))

METIMMOX$PFSmo <- METIMMOX$`Days until progression`*12/365

METIMMOX$OSmo <-METIMMOX$`Days until death/last follow up`*12/365

METIMMOX$PFSwk <- METIMMOX$`Days until progression`/7

METIMMOX$OSwk <-METIMMOX$`Days until death/last follow up`/7


saveRDS(METIMMOX, file="data/tidy/METIMMOX.rds")

# load faster binary format
rm(list = ls())
METIMMOX <- readRDS(file="data/tidy/METIMMOX.rds")