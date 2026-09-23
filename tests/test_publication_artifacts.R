# Synthetic zero/empty/failure and artifact-lifecycle regressions (#168).
source(here::here("R/scenario_analysis.R"))
source(here::here("R/evppi_functions.R"))
source(here::here("R/publication_artifacts.R"))

fixture <- function(id, effect, estimates = data.frame()) list(
  scenario_info = data.frame(scenario_id = id, scenario_name = id, wtp = 1),
  psa_obj = list(effect = effect, cost = matrix(0, nrow(effect), ncol(effect))),
  evppi_results = estimates)
groups <- list(utilities = c("u_np", "u_decrement"))
zero <- fixture("zero", cbind(c(2, 2), c(1, 1)))
small <- fixture("small", cbind(c(1, 1), c(1.01, .99)))
failed <- fixture("failed", cbind(c(1, 1), c(2, 0)), data.frame(
  parameter = "[GROUP] utilities", evppi = NA_real_, evppi_se = NA_real_,
  error = "Regression failed"))
result <- compile_evppi_results(list(zero, small, failed), groups)
stopifnot(setequal(result$scenario_id, c("zero", "small", "failed")),
  sum(result$parameter == "[EVPI]") == 3L)
zr <- subset(result, scenario_id == "zero" & parameter != "[EVPI]")
sr <- subset(result, scenario_id == "small" & parameter != "[EVPI]")
fr <- subset(result, scenario_id == "failed" & parameter != "[EVPI]")
stopifnot(zr$evppi == 0, zr$evppi_se == 0, zr$error == "",
  is.na(zr$evppi_percent_of_evpi), is.na(sr$evppi), nzchar(sr$error),
  is.na(fr$evppi), fr$error == "Regression failed")
stopifnot(nrow(compile_evppi_results(list(zero), list())) == 1L,
  nrow(compile_evppi_results(list(), groups)) == 0L)

root <- tempfile("artifact-contract-")
dir.create(file.path(root, "outputs/tables"), recursive = TRUE)
path <- file.path(root, "outputs/tables/table_5.csv")
writeLines("obsolete", path)
begin_artifact_render("table_5", root)
stopifnot(!file.exists(path))
value <- data.frame(label = c("comma, quote\"", "plain"), value = c(pi, NA_real_))
write_artifact_csv(value, path, row.names = FALSE)
env <- new.env(parent = emptyenv())
env$psa_obj <- list(fingerprint = "synthetic-psa")
finish_artifact_render(env)
m <- read.csv(artifact_manifest_path(root), stringsAsFactors = FALSE)
stopifnot(m$status == "csv_verified", m$comparison_with_predecessor == "changed",
  m$source_cache_fingerprints == "psa_obj=synthetic-psa",
  m$artifact_md5 == unname(tools::md5sum(path)))
# Empty successful render must remove both the predecessor and its old claim.
begin_artifact_render("table_5", root)
stopifnot(!file.exists(path), nrow(read.csv(artifact_manifest_path(root))) == 0L)
finish_artifact_render(env)
m <- read.csv(artifact_manifest_path(root))
stopifnot(m$status == "deleted_empty", !file.exists(path))
# A failed setup also cannot leave a previous output or provenance row behind.
begin_artifact_render("table_5", root)
stopifnot(nrow(read.csv(artifact_manifest_path(root))) == 0L,
  inherits(try(write_artifact_csv(value, file.path(root, "undeclared.csv")), silent = TRUE), "try-error"))
options(metimmox.artifact_context = NULL)

# Execute the real Figure 5 plotting chunk under all-zero, failed and entirely
# empty fixtures. This catches successful data compilation followed by a plot
# failure or a missing replacement image without needing confidential caches.
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
src <- readLines(here::here("outputs/vignettes/figure5.qmd"), warn = FALSE)
start <- grep("^```\\{r figure,", src) + 1L
end <- start + which(src[start:length(src)] == "```")[1] - 2L
plot_code <- parse(text = src[start:end])
plot_env <- new.env(parent = .GlobalEnv)
plot_env$evppi_grouped <- expand.grid(param_group = c("Costs", "Utilities"),
  Scenario = c("Base case", "Biosimilar", "Higher WTP"), stringsAsFactors = FALSE)
plot_env$evppi_grouped$evppi_pop <- 0
plot_env$evppi_grouped$se_pop <- 0
grDevices::pdf(NULL)
eval(plot_code, plot_env)
built <- ggplot_build(plot_env$p_evppi)
stopifnot(sum(built$data[[3]]$label == "0") == 6L)
plot_env$evppi_grouped$evppi_pop[1] <- NA_real_
plot_env$evppi_grouped$se_pop[1] <- NA_real_
eval(plot_code, plot_env)
stopifnot(sum(ggplot_build(plot_env$p_evppi)$data[[3]]$label == "NA") == 1L)
plot_env$evppi_grouped <- data.frame()
eval(plot_code, plot_env)
grDevices::dev.off()
begin_artifact_render("figure5", root)
image_path <- file.path(root, "outputs/figs/figure5.png")
save_artifact_plot(image_path, plot_env$p_evppi, width = 4, height = 3, dpi = 72)
finish_artifact_render(env)
stopifnot(file.exists(image_path), file.info(image_path)$size > 0)
begin_artifact_render("figure5", root)
finish_artifact_render(env)
stopifnot(!file.exists(image_path))
cat("PASS: zero, small-positive, failed, all-empty scenarios; real Figure 5 zero/NA/empty renders; CSV object agreement; manifest and stale-artifact lifecycle.\n")
