# Regression test for WB-005 / BB-PSM-WCHECK.
# Run from the repository root with Rscript. No trial data or caches are needed.

predict.mock_survival_model <- function(object, newdata, type, times, ...) {
  rx_survival <- ifelse(as.character(newdata$Rx) == "experimental", 0.8, 0.4)
  predictions <- lapply(rx_survival, function(probability) {
    data.frame(.pred_survival = c(1, rep(probability, length(times) - 1L)))
  })
  structure(
    list(.pred = predictions),
    class = "data.frame",
    row.names = seq_along(predictions)
  )
}

run_canonical_prevalence_test <- function() {
  source("R/model_configs.R", local = environment())
  source("R/prediction_functions.R", local = environment())

  data_complete <- data.frame(
    Rx = factor(
      c("control", "experimental", "control", "experimental"),
      levels = c("control", "experimental")
    ),
    crp = c(0, 0, 1, 1),
    tmb_braf = c(0, 1, 0, 1)
  )
  strategies_df <- data.frame(
    id = c("control", "crp", "tmb_braf"),
    prevalence = c(1, 0.25, 0.75)
  )
  canonical_prevalences <- setNames(
    strategies_df$prevalence,
    strategies_df$id
  )
  mock_model <- structure(list(), class = "mock_survival_model")

  predictions <- generate_population_averaged_predictions(
    models = list(os = mock_model, pfs = mock_model),
    strategies_df = strategies_df,
    data_complete = data_complete,
    time_points = c(0, 1),
    prevalences = canonical_prevalences,
    quiet = TRUE
  )

  # The complete-case prevalence is 0.5 for both biomarkers. The supplied
  # canonical values must instead determine the weighted curves.
  stopifnot(
    identical(predictions$crp$prevalence, 0.25),
    identical(predictions$tmb_braf$prevalence, 0.75),
    isTRUE(all.equal(predictions$crp$os[2], 0.5)),
    isTRUE(all.equal(predictions$tmb_braf$os[2], 0.7)),
    isTRUE(all.equal(predictions$crp$pfs[2], 0.5)),
    isTRUE(all.equal(predictions$tmb_braf$pfs[2], 0.7))
  )

  invalid_error <- tryCatch(
    {
      generate_population_averaged_predictions(
        models = list(os = mock_model, pfs = mock_model),
        strategies_df = strategies_df,
        data_complete = data_complete,
        time_points = c(0, 1),
        prevalences = c(crp = 0.25),
        quiet = TRUE
      )
      NULL
    },
    error = identity
  )
  stopifnot(
    inherits(invalid_error, "error"),
    grepl("prevalences must be a named vector", conditionMessage(invalid_error),
          fixed = TRUE)
  )

  cat("PASS: weighted survival curves use canonical full-cohort prevalences.\n")
}

run_canonical_prevalence_test()
