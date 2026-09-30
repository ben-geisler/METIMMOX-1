# ===============================================================================
# TEST: cache fingerprints and manifest order do not depend on LC_COLLATE (#185)
# ===============================================================================
# Default sort()/order() collate by locale: "crp" sorts after "Death" in C but
# before it in English, so the same data hashed differently in RStudio and in a
# command-line Rscript session and every cache was reported stale. Each result
# below is computed under the C and an English collation locale and must be
# identical, and equal to the C (radix) order.
#
# Synthetic data and a temporary manifest only; no trial data or cache is read.
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_locale_independent_fingerprints.R
# ===============================================================================

source("R/cache_provenance.R")
source("R/cache_paths.R")
source("R/publication_artifacts.R")

original_collate <- Sys.getlocale("LC_COLLATE")
on.exit(Sys.setlocale("LC_COLLATE", original_collate), add = TRUE)

english <- c("English_United States.utf8", "English_United States.1252",
             "en_US.UTF-8", "en_US.utf8", "en_GB.UTF-8")
set_collate <- function(candidates) {
  for (locale in candidates) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) return(locale)
  }
  NA_character_
}
stopifnot(!is.na(set_collate("C")))
english_locale <- set_collate(english)
if (is.na(english_locale)) {
  cat("SKIP: no English collation locale available on this machine.\n")
  quit(save = "no", status = 0)
}
# The premise: the two locales order mixed-case names differently.
stopifnot(!identical(sort(c("Death", "crp")), { Sys.setlocale("LC_COLLATE", "C"); sort(c("Death", "crp")) }))

data <- data.frame(ID = as.character(1:6), Death = c(1, 0, 1, 1, 0, 1),
                   OSwk = c(10, 20, 30, 40, 50, 60), crp = factor(c(0, 1, 0, 1, 1, 0)),
                   sex = factor(c(0, 0, 1, 1, 0, 1)), stringsAsFactors = FALSE)
formulas <- list(os = Surv(OSwk, Death) ~ crp + sex)
mixed_list <- list(zeta = 1, Alpha = 2, beta = list(Death = 3, crp = 4))

manifest_root <- tempfile("locale-manifest-")
dir.create(file.path(manifest_root, "outputs", "tables"), recursive = TRUE)
dir.create(file.path(manifest_root, "outputs", "vignettes"), recursive = TRUE)
dir.create(file.path(manifest_root, "R"))
invisible(file.copy("R/publication_artifacts.R", file.path(manifest_root, "R")))
writeLines("---\n---", file.path(manifest_root, "outputs", "vignettes", "table_5.qmd"))
# Seed a manifest whose rows collate differently in C and English.
begin_artifact_render("table_5", manifest_root)
write_artifact_csv(data.frame(x = 1), file.path(manifest_root, "outputs/tables/table_5.csv"),
                   row.names = FALSE)
finish_artifact_render(new.env())
seed <- read.csv(artifact_manifest_path(manifest_root), stringsAsFactors = FALSE)
extra <- seed[rep(1, 3), ]
extra$file <- c("outputs/figs/figure_s4.png", "outputs/figs/figure1.png", "outputs/figs/Figure9.png")
write.csv(rbind(seed, extra), artifact_manifest_path(manifest_root), row.names = FALSE, na = "")

results <- lapply(c("C", english_locale), function(locale) {
  stopifnot(nzchar(Sys.setlocale("LC_COLLATE", locale)))
  sampling <- sampling_cache_fingerprint(formulas, data, list(os = "gamma"), 10, 1, "test")
  begin_artifact_render("table_5", manifest_root)
  write_artifact_csv(data.frame(x = 1), file.path(manifest_root, "outputs/tables/table_5.csv"),
                     row.names = FALSE)
  finish_artifact_render(new.env())
  list(sampling = sampling$fingerprint,
       data_hash = sampling$inputs$data_hash,
       data_columns = names(sampling$inputs$data_columns),
       list_hash = cache_fingerprint(mixed_list),
       frame_hash = cache_fingerprint(data),
       manifest_files = read.csv(artifact_manifest_path(manifest_root),
                                 stringsAsFactors = FALSE)$file)
})
names(results) <- c("C", english_locale)

for (key in names(results[[1]])) {
  if (!identical(results[[1]][[key]], results[[2]][[key]])) {
    stop("'", key, "' differs between LC_COLLATE = C and ", english_locale, ".")
  }
}
stopifnot(
  identical(results$C$data_columns, c("Death", "ID", "OSwk", "crp", "sex")),
  identical(results$C$manifest_files, sort(results$C$manifest_files, method = "radix"))
)

cat("PASS: sampling, data-frame and list fingerprints and manifest row order are",
    "identical under C and", english_locale, "collation (#185).\n")
