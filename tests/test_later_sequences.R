#' Third and fourth treatment sequences (issue #208), synthetic, no trial data:
#' 1. later_sequence_regimen() is the first sequence of the production
#'    schedules (8 FLOX in the control arm; 4 FLOX + 4 nivolumab otherwise).
#' 2. later_sequence_inputs(): counts, shares (p_seq4 given p_seq3) and start
#'    weeks by arm in the given population only; malformed data stop.
#' 3. later_sequence_resource_use() against a hand calculation: share x one
#'    sequence at each observed start week, two contacts per FLOX and one per
#'    nivolumab administration, one blood test per administration week,
#'    doses beyond the grid dropped; no inputs give no resource use; a
#'    positive share without start weeks stops.
#' 4. calculate_outcomes(): the later-sequence cost is the same whatever the
#'    occupancy, and QALYs and life-years do not change.
#' 5. The Second_sequence scenario removes the later sequences; the base list
#'    is unchanged when the scenario keeps the second sequence.
#'
#' Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_later_sequences.R

source("R/model_configs.R")
source("R/prediction_population.R")
source("R/parameter_distributions.R")
source("R/calculate_outcomes.R")
source("R/later_sequences.R")

#' Stop unless `expr` raises an error whose message contains `pattern`.
fails <- function(expr, pattern) {
  e <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e), fixed = TRUE))
}

# 1. Regimen of one sequence.
regimen <- later_sequence_regimen()
stopifnot(identical(regimen$control$FLOX, as.integer(seq(0, 14, by = 2))),
          length(regimen$control$nivolumab) == 0L,
          identical(regimen$experimental$FLOX, c(0L, 2L, 8L, 10L)),
          identical(regimen$experimental$nivolumab, c(4L, 6L, 12L, 14L)))

# 2. Inputs for a population of 10 (4 control, 6 experimental); patient 11 is
# outside the population and must not count.
df_population <- data.frame(ID = 1:10, Rx = factor(rep(c("Control arm", "Experimental arm"), c(4, 6)),
                                                   levels = c("Control arm", "Experimental arm")))
df_sequences <- data.frame(ID = c(11, 10:1),
  Seq3 = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0),
  Seq3wk = c(30, 50.4, 60, 49.6, NA, NA, NA, NA, NA, 70, NA),
  Seq4 = c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  Seq4wk = c(60, 80, NA, NA, NA, NA, NA, NA, NA, NA, NA))
inputs <- later_sequence_inputs(df_population, df_sequences)
later <- inputs$later_sequences
stopifnot(identical(later$counts$arm, c("control", "experimental")),
          identical(later$counts$n_patients, c(4, 6)),
          identical(later$counts$n_seq3, c(1, 3)),
          identical(later$counts$n_seq4, c(0, 1)),
          isTRUE(all.equal(unlist(inputs$shares),
                           c(p_seq3_control = 1 / 4, p_seq4_given_seq3_control = 0,
                             p_seq3_experimental = 3 / 6, p_seq4_given_seq3_experimental = 1 / 3))),
          setequal(names(inputs$shares), LATER_SEQUENCE_SHARE_PARAMETERS),
          identical(later$start_weeks$arm, c("control", "experimental", "experimental", "experimental")),
          identical(later$start_weeks$sequence, c(3L, 3L, 3L, 4L)),
          identical(later$start_weeks$week, c(70L, 50L, 60L, 80L)),
          identical(later$start_weeks$n, c(1, 2, 1, 1)))
df_summary <- summarise_later_sequence_starts(later)
stopifnot(all(df_summary$median == c(70, 50, 80)), all(df_summary$n == c(1, 3, 1)))
fails(later_sequence_inputs(df_population, transform(df_sequences, Seq3 = replace(Seq3, 2, 0))),
      "fourth sequence only after a third")
fails(later_sequence_inputs(df_population, transform(df_sequences, Seq3wk = replace(Seq3wk, 2, NA))),
      "Sequence 3 starts must have a nonnegative start week")
fails(later_sequence_inputs(df_population, df_sequences[-2, ]), "lack patients of the population")
fails(later_sequence_inputs(df_population, df_sequences[, -5]), "Missing later-sequence columns")

# 3. Resource use against a hand calculation on a 100-week grid.
n_points <- 101L
params <- c(inputs$shares, list(later_sequences = later))
use_exp <- later_sequence_resource_use(params, "experimental", n_points)
v_flox <- v_nivo <- rep(0, n_points)
for (start in c(50, 50, 60)) {  # third sequence, share 3/6 split over three starters
  v_flox[start + c(0, 2, 8, 10) + 1] <- v_flox[start + c(0, 2, 8, 10) + 1] + (3 / 6) / 3
  v_nivo[start + c(4, 6, 12, 14) + 1] <- v_nivo[start + c(4, 6, 12, 14) + 1] + (3 / 6) / 3
}
v_keep <- 80 + c(0, 2, 8, 10) <= 100  # fourth sequence at week 80, share 1/2 x 1/3
v_flox[80 + c(0, 2, 8, 10)[v_keep] + 1] <- v_flox[80 + c(0, 2, 8, 10)[v_keep] + 1] + 1 / 6
v_keep <- 80 + c(4, 6, 12, 14) <= 100
v_nivo[80 + c(4, 6, 12, 14)[v_keep] + 1] <- v_nivo[80 + c(4, 6, 12, 14)[v_keep] + 1] + 1 / 6
stopifnot(max(abs(use_exp$FLOX - v_flox)) < 1e-15, max(abs(use_exp$nivolumab - v_nivo)) < 1e-15,
          max(abs(use_exp$contacts - (2 * v_flox + v_nivo))) < 1e-15,
          max(abs(use_exp$blood - (v_flox + v_nivo))) < 1e-15,
          abs(sum(use_exp$FLOX) - (4 * 3 / 6 + 4 * 1 / 6)) < 1e-12)
use_ctl <- later_sequence_resource_use(params, "control", n_points)
stopifnot(abs(sum(use_ctl$FLOX) - 8 / 4) < 1e-12, all(use_ctl$nivolumab == 0),
          all(use_ctl$FLOX[71 + seq(0, 14, by = 2)] == 1 / 4))
use_none <- later_sequence_resource_use(list(), "control", n_points)
stopifnot(all(unlist(use_none) == 0))
fails(later_sequence_resource_use(modifyList(params, list(p_seq4_given_seq3_control = 0.5)),
                                  "control", n_points),
      "A positive share starts sequence 4 in the control arm but no start weeks are given")
fails(later_sequence_resource_use(modifyList(params, list(p_seq3_control = 1.2)), "control", n_points),
      "must be finite values in [0, 1]")

# 4. State independence in calculate_outcomes().
params <- c(params, list(u_np = 0.8, u_p = 0.6, c_drug_nivo = 1000, c_drug_FLOX = 400,
                         c_test_CT = 300, c_test_blood = 15, c_other_visit = 200,
                         c_other_baseline = 200, c_other_follow = 200, c_other_pp = 0,
                         c_other_last = 5000, dr_costs = 0.04, dr_effects = 0.04, cl = 1 / 52),
            build_treatment_schedules(n_points))
v_dw <- discount_weights(params, n_points)
run <- function(p, os, pfs, type) {
  s <- partitioned_survival_states(os, pfs)
  calculate_outcomes(p, s$p_pf, s$p_p, s$p_d, type, NULL, v_dw$cost, v_dw$effect, 1 / 52)
}
p_none <- params
p_none$later_sequences <- NULL
v_expected <- sum((use_exp$FLOX * 400 + use_exp$nivolumab * 1000 + use_exp$contacts * 200 +
                     use_exp$blood * 15) * v_dw$cost)
for (rate in c(0.005, 0.03)) {
  os <- exp(-rate * 0:100)
  pfs <- exp(-2 * rate * 0:100)
  with_later <- run(params, os, pfs, "experimental")
  without <- run(p_none, os, pfs, "experimental")
  stopifnot(abs(with_later$costs_total - without$costs_total - v_expected) < 1e-8,
            abs(sum(with_later$later_sequence_costs * v_dw$cost) - v_expected) < 1e-8,
            identical(with_later$qalys_total, without$qalys_total),
            identical(with_later$ly_disc_total, without$ly_disc_total),
            abs(with_later$nivo_administrations_total - without$nivo_administrations_total -
                  sum(use_exp$nivolumab * v_dw$cost)) < 1e-12)
}

# 5. The second-sequence scenario removes the later sequences.
dropped <- apply_structural_scenario("Second_sequence", 0, params, 100)$params
stopifnot(dropped$p_seq3_control == 0, dropped$p_seq3_experimental == 0,
          all(unlist(later_sequence_resource_use(dropped, "experimental", n_points)) == 0),
          identical(apply_structural_scenario("Second_sequence", 1, params, 100)$params, params),
          identical(structural_scenario_label("Second_sequence", 0),
                    "Second and later treatment sequences omitted"))
cat("PASS: later treatment sequences (issue #208): regimen, inputs by arm, resource use, state independence, Second_sequence scenario\n")
