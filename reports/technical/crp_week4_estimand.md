# Week-4 CRP Estimand
Ben Geisler
2026-09-25

- [Purpose](#purpose)
- [Estimand](#estimand)
- [Patients Without a Week-4 CRP](#patients-without-a-week-4-crp)
- [Missingness Sensitivity Analysis](#missingness-sensitivity-analysis)
  - [Methods](#methods)
  - [Scenario cohorts and fits](#scenario-cohorts-and-fits)
  - [Cost-effectiveness results](#cost-effectiveness-results)
- [Costs and QALYs Accrued Before the First Nivolumab
  Dose](#costs-and-qalys-accrued-before-the-first-nivolumab-dose)
  - [Methods](#methods-1)
  - [Results](#results)
- [Limitations](#limitations)

# Purpose

The CRP-guided strategy uses the week-4 C-reactive protein value (CRP,
cycle 3 day 1), measured after the two FLOX cycles that both arms
receive and before the first nivolumab dose. Issue \#175 asked which
estimand the economic model targets, and how much two consequences of
that choice matter:

- the patients without a week-4 CRP are excluded from the complete-case
  cohort, and
- the fitted treatment effects act from randomisation, although the
  treatment schedules of all strategies are identical until the first
  nivolumab dose.

The chosen approach is to keep the current model and quantify both
points with deterministic sensitivity analyses (option (c) of the
issue). This report states the estimand, describes the patients without
a week-4 CRP in aggregate, restores them under four biomarker
assignments, and measures the costs and QALYs the model accrues before
the first nivolumab dose. No model code and no cache is changed; all
results are deterministic.

# Estimand

The economic model estimates a **randomisation-time policy**. Each
strategy is assigned at randomisation (t = 0), and costs and QALYs are
counted from that point over the 520-week horizon. Under the CRP-guided
strategy, nivolumab is added from modeled week 4 for patients whose
week-4 CRP is below 5 mg/L; the other patients continue FLOX alone. The
week-4 value is available when the decision is made, because it is
measured before the first nivolumab dose. The target population is the
complete-case cohort of 68 patients with age, sex, arm, week-4 CRP,
TMB/BRAF and both endpoints recorded. Patients without a week-4 CRP are
excluded from this cohort and are not modeled as a separate decision
population. Under a randomisation-time estimand they belong to the
target population, so all 3 of them are missing data, including the
patient who died before week 4. The missingness sensitivity analysis
therefore restores all of them.

The treatment contrast is estimated from randomisation. The joint OS and
PFS models contain `Rx` and biomarker-by-`Rx` terms that act from t = 0,
so the strategies’ survival curves separate before week 4, although
every strategy receives the same FLOX, monitoring and visits until then.
Two alternatives were considered and not adopted:

- **(a)** force identical curves up to week 4 and start the treatment
  contrast at the decision point; or
- **(b)** a landmark estimand among patients alive at week 4. It would
  place the patient who died before week 4 outside the decision
  population and treat only the patients alive at week 4 as missing
  data.

The missingness sensitivity analysis and the pre-decision accrual below
quantify how much the results depend on the adopted estimand.

# Patients Without a Week-4 CRP

| Characteristic                                         | Value |
|:-------------------------------------------------------|------:|
| Patients without a week-4 CRP                          |     3 |
| Control arm / experimental arm                         | 1 / 2 |
| Died before week 4 (week-4 CRP could not be measured)  |     1 |
| Alive at week 4 without a week-4 CRP                   |     2 |
| Deaths during follow-up                                |     3 |
| Baseline CRP (cycle 1 day 1) available                 |     3 |
| TMB/BRAF known positive (BRAF-mutant)                  |     1 |
| TMB/BRAF unknown (KRAS-mutant, no TMB result)          |     2 |
| PFS events: recorded progression / death within window | 1 / 2 |

Patients without a week-4 CRP (aggregate)

3 patients have no week-4 CRP: 1 in the control arm and 2 in the
experimental arm. All 3 died, at weeks 2.4, 15.7, 20.9: the control-arm
patient at week 15.7, and the experimental-arm patients at weeks 2.4 and
20.9. Only 1 of them died before week 4; the other 2 were alive at week
4 but have no week-4 value. All 3 have a baseline CRP. The control-arm
patient is BRAF-mutant and therefore TMB/BRAF-positive. The 2
experimental-arm patients are KRAS-mutant with no TMB result, so their
TMB/BRAF status is unknown; they would be excluded from the
complete-case cohort even with a week-4 CRP. Their PFS endpoints follow
the issue \#181 rule already applied by script 02 (a 16-week death
window after the last assessment) and are used as they are.

The complete-case cohort also excludes 3 patients who have a week-4 CRP
but no TMB/BRAF result (all in the control arm). They are not part of
this analysis and stay excluded in every scenario.

# Missingness Sensitivity Analysis

## Methods

The 3 patients without a week-4 CRP are restored to the complete-case
cohort, giving 71 patients. Four scenarios cross two assignments:

- the CRP value assigned to all 3 patients (negative or positive), and
- the TMB/BRAF value assigned to the 2 patients with unknown TMB/BRAF
  status (negative or positive).

The control-arm patient keeps the known TMB/BRAF-positive status. For
each scenario:

1.  **Refit.** The joint OS and PFS models are refitted with the
    pipeline formulas (`get_model_formulas()`) and the distributions
    selected for the base case (gamma for OS, gamma for PFS). The
    distributions are not re-selected.
2.  **Rebuild the population.** The common prediction population, the
    biomarker prevalences and the joint biomarker cells are rebuilt from
    the 71-patient cohort with the pipeline helpers, among them
    `economic_prediction_population()` and
    `set_population_predictions()`.
3.  **Check ordering.** OS \>= PFS is checked for the control curve and
    the four biomarker-subgroup curves at every weekly time point. The
    deterministic model stops on a violation, so a scenario that
    violates the ordering is reported, not clamped.
4.  **Run the model.** The deterministic model is run with the unchanged
    schedules, prices and utilities, as script 05 builds them.

Pairwise ICERs and statuses versus standard of care come from
`calculate_pairwise_icers()`, and frontier statuses from dampack.
Applied to the unmodified 68-patient cohort, the same code reproduces
the deterministic base case (largest absolute difference in costs and
QALYs: 0, identical). Each scenario changes both the fitted treatment
effects and the target population (its size, prevalences and age/sex
mix), so a difference from the base case combines the two.

## Scenario cohorts and fits

| Cohort          |   n | Control/exp. |  CRP+ | TMB/BRAF+ | Converged | PFS \> OS points |
|:----------------|----:|-------------:|------:|----------:|:----------|-----------------:|
| Base case       |  68 |      32 / 36 | 33.8% |     44.1% | Yes       |                0 |
| CRP-, TMB/BRAF- |  71 |      33 / 38 | 32.4% |     43.7% | Yes       |                0 |
| CRP-, TMB/BRAF+ |  71 |      33 / 38 | 32.4% |     46.5% | Yes       |                0 |
| CRP+, TMB/BRAF- |  71 |      33 / 38 | 36.6% |     43.7% | Yes       |                0 |
| CRP+, TMB/BRAF+ |  71 |      33 / 38 | 36.6% |     46.5% | Yes       |                0 |

Cohorts, prevalences, convergence and OS \>= PFS ordering

*Note:* Converged: both refitted models converged with a finite
covariance matrix. PFS \> OS points: weekly time points at which PFS
exceeds OS, summed over the control curve and the four
biomarker-subgroup curves. The base-case row is the reproduction run on
the complete-case cohort.

All 4 scenario fits converged, and no refitted scenario puts PFS above
OS on any of the five curves at any of the 521 weekly time points, so
every scenario could be run through the deterministic model.

| Cohort          |    Rx | CRP x Rx | TMB/BRAF x Rx |
|:----------------|------:|---------:|--------------:|
| Base case       | 0.222 |   -0.171 |        -0.098 |
| CRP-, TMB/BRAF- | 0.325 |   -0.210 |        -0.225 |
| CRP-, TMB/BRAF+ | 0.268 |   -0.293 |         0.007 |
| CRP+, TMB/BRAF- | 0.306 |   -0.176 |        -0.344 |
| CRP+, TMB/BRAF+ | 0.248 |   -0.261 |        -0.096 |

OS model: treatment and biomarker-by-treatment coefficients

| Cohort          |    Rx | CRP x Rx | TMB/BRAF x Rx |
|:----------------|------:|---------:|--------------:|
| Base case       | 0.497 |   -0.957 |        -0.325 |
| CRP-, TMB/BRAF- | 0.569 |   -1.003 |        -0.450 |
| CRP-, TMB/BRAF+ | 0.535 |   -1.120 |        -0.219 |
| CRP+, TMB/BRAF- | 0.587 |   -0.981 |        -0.629 |
| CRP+, TMB/BRAF+ | 0.534 |   -1.089 |        -0.346 |

PFS model: treatment and biomarker-by-treatment coefficients

*Note:* Coefficients on the flexsurvreg rate scale of the gamma (OS) and
gamma (PFS) models; Rx is the experimental arm. A positive coefficient
raises the rate and shortens survival. Standard errors are not shown.

## Cost-effectiveness results

| Cohort          | Standard of Care | CRP-guided | TMB/BRAF-guided | Highest NMB      |
|:----------------|-----------------:|-----------:|----------------:|:-----------------|
| Base case       |           20,967 |     53,686 |          62,382 | Standard of Care |
| CRP-, TMB/BRAF- |           20,811 |     51,936 |          61,344 | Standard of Care |
| CRP-, TMB/BRAF+ |           20,818 |     52,308 |          61,226 | Standard of Care |
| CRP+, TMB/BRAF- |           20,668 |     54,435 |          62,734 | Standard of Care |
| CRP+, TMB/BRAF+ |           20,670 |     55,144 |          63,125 | Standard of Care |

Deterministic costs (EUR) by cohort

*Note:* Highest NMB: strategy with the highest net monetary benefit at
WTP = EUR 51,000 per QALY.

| Cohort          | Standard of Care | CRP-guided | TMB/BRAF-guided |
|:----------------|-----------------:|-----------:|----------------:|
| Base case       |           1.3400 |     1.3775 |          1.3389 |
| CRP-, TMB/BRAF- |           1.3176 |     1.3555 |          1.3409 |
| CRP-, TMB/BRAF+ |           1.3205 |     1.3695 |          1.2648 |
| CRP+, TMB/BRAF- |           1.2753 |     1.3349 |          1.3673 |
| CRP+, TMB/BRAF+ |           1.2788 |     1.3520 |          1.2942 |

Deterministic QALYs by cohort

| Cohort          | Inc. cost | Inc. QALYs |    ICER | Frontier             |
|:----------------|----------:|-----------:|--------:|:---------------------|
| Base case       |    32,720 |     0.0375 | 871,972 | On frontier          |
| CRP-, TMB/BRAF- |    31,124 |     0.0378 | 823,085 | On frontier          |
| CRP-, TMB/BRAF+ |    31,490 |     0.0490 | 642,793 | On frontier          |
| CRP+, TMB/BRAF- |    33,767 |     0.0596 | 566,110 | Extendedly dominated |
| CRP+, TMB/BRAF+ |    34,474 |     0.0732 | 471,086 | On frontier          |

CRP-guided: increments and ICER pairwise versus standard of care (EUR,
QALYs)

| Cohort          | Inc. cost | Inc. QALYs |      ICER | Frontier    |
|:----------------|----------:|-----------:|----------:|:------------|
| Base case       |    41,416 |    -0.0011 | Dominated | Dominated   |
| CRP-, TMB/BRAF- |    40,533 |     0.0232 | 1,744,687 | Dominated   |
| CRP-, TMB/BRAF+ |    40,408 |    -0.0557 | Dominated | Dominated   |
| CRP+, TMB/BRAF- |    42,065 |     0.0921 |   456,884 | On frontier |
| CRP+, TMB/BRAF+ |    42,455 |     0.0154 | 2,754,342 | Dominated   |

TMB/BRAF-guided: increments and ICER pairwise versus standard of care
(EUR, QALYs)

*Note:* Increments and ICERs are pairwise versus standard of care.
‘Dominated’ in the ICER column: pairwise status ‘Dominated by SoC’ (more
costly and no more effective than standard of care); a number: pairwise
status ‘Pairwise ICER vs SoC’; ‘Cost-saving’: ‘Cost-saving vs SoC’.
Frontier: dampack efficiency-frontier status across all three
strategies.

**Range across the four scenarios.**

- CRP-guided: pairwise ICERs from EUR 471,086 (CRP+, TMB/BRAF+) to EUR
  823,085 (CRP-, TMB/BRAF-) in 4 scenarios (base case: a pairwise ICER
  of EUR 871,972).
- TMB/BRAF-guided: pairwise ICERs from EUR 456,884 (CRP+, TMB/BRAF-) to
  EUR 2,754,342 (CRP+, TMB/BRAF+) in 3 scenarios; dominated by standard
  of care in 1 (CRP-, TMB/BRAF+) (base case: dominated by standard of
  care).
- CRP-guided was on the efficiency frontier in 3 of 4 scenarios;
  TMB/BRAF-guided was on the efficiency frontier in 1 of 4 scenarios.
- Standard of Care had the highest net monetary benefit at WTP = EUR
  51,000 in the base case and in every evaluated scenario, so the
  decision does not change.

Standard-of-care QALYs fall from 1.3400 in the base case to between
1.2753 and 1.3205, because every restored patient died within 21 weeks
and the population is standardised to the enlarged cohort. The
CRP-guided pairwise ICER is lower than in the base case in 4 of 4
evaluated scenarios, at 54% to 94% of its base-case value. Both
assignments enter the joint model, so each one moves both guided
strategies. Switching the CRP assignment from negative to positive
changes the incremental QALYs of the CRP-guided strategy by 0.0230 on
average and those of the TMB/BRAF-guided strategy by 0.0700. Switching
the TMB/BRAF assignment changes them by 0.0124 and 0.0778. These changes
are large relative to the base-case increments of 0.0375 (CRP-guided)
and -0.0011 (TMB/BRAF-guided). The interaction estimates rest on few
patients per arm and biomarker cell, so 3 patients with early deaths are
enough to move them.

# Costs and QALYs Accrued Before the First Nivolumab Dose

## Methods

All strategies receive the same treatment schedule until modeled week 4:
no nivolumab dose is scheduled before week 4, and the FLOX schedules of
the two arms are identical at grid points 0 to 3 (both checked when the
report is rendered). The first nivolumab dose is at week 4. Apart from
the one-time diagnostic test, any difference the model accrues between
strategies before week 4 therefore comes only from the fitted `Rx` and
biomarker-by-`Rx` effects acting from t = 0.

The deterministic base-case traces
(`model_fun(..., return_traces = TRUE)`) are passed through
`calculate_outcomes()` subgroup by subgroup, and each component is split
at week 4 using the model’s own discounting and trapezoidal interval
weights:

- **QALYs and ongoing progressed-state costs** (interval flows) are
  counted over the intervals from week 0 to week 4. This is the full
  trapezoidal weight of grid points 0 to 3 plus half the weight of grid
  point 4.
- **Scheduled and one-time charges** are counted at grid points 0 to 3:
  drugs, CT scans, blood tests, visits, the baseline visit and the
  diagnostic test (CRP or NGS). Every charge at grid point 4, including
  the first nivolumab dose, is excluded.
- **End-of-life costs** are counted for deaths in the intervals ending
  at weeks 1 to 4, charged at those grid points.

The recomputed full-horizon totals equal `model_fun()` (largest absolute
difference 0, identical). The diagnostic test is charged at t = 0 in
every biomarker-guided strategy, so it is shown separately from the
survival-dependent part of the cost increment.

## Results

| Measure                  | Standard of Care | CRP-guided | TMB/BRAF-guided |
|:-------------------------|-----------------:|-----------:|----------------:|
| QALYs, weeks 0-4         |          0.05600 |    0.05601 |         0.05600 |
| QALYs, full horizon      |           1.3400 |     1.3775 |          1.3389 |
| Cost, weeks 0-4          |        EUR 1,862 |  EUR 1,879 |       EUR 4,381 |
| of which diagnostic test |            EUR 0 |     EUR 16 |       EUR 2,518 |
| Cost, full horizon       |       EUR 20,967 | EUR 53,686 |      EUR 62,382 |
| OS at week 4             |           0.9991 |     0.9991 |          0.9990 |
| PFS at week 4            |           0.9845 |     0.9878 |          0.9847 |

Discounted QALYs and costs before week 4 and over the full horizon
(deterministic base case)

*Note:* weeks 0-4: accrued before the first nivolumab dose, as defined
in the Methods. OS and PFS: population-averaged survival probabilities.

| Measure                              | CRP-guided | TMB/BRAF-guided |
|:-------------------------------------|-----------:|----------------:|
| Inc. QALYs, weeks 0-4                |   0.000013 |     -0.00000065 |
| Inc. QALYs, full horizon             |     0.0375 |         -0.0011 |
| Share of full-horizon QALY increment |     0.035% |          0.061% |
| Inc. cost, weeks 0-4                 |  EUR 16.50 |    EUR 2,519.36 |
| of which diagnostic test             |  EUR 16.00 |    EUR 2,518.00 |
| of which survival-dependent          |   EUR 0.50 |        EUR 1.36 |
| Inc. cost, full horizon              | EUR 32,720 |      EUR 41,416 |
| Share of full-horizon cost increment |     0.050% |          6.083% |
| Share, survival-dependent part only  |     0.002% |          0.003% |

Increments versus standard of care: before week 4 and full horizon

*Note:* Survival-dependent: the pre-decision cost increment without the
diagnostic test, that is, the part driven by differences in state
occupancy.

Before week 4, the CRP-guided strategy accrues 0.000013 more QALYs than
standard of care, which is 0.035% of its full-horizon QALY increment of
0.0375. The survival-dependent part of its cost increment before week 4
is EUR 0.50 (0.002% of the full-horizon cost increment). The CRP test
adds EUR 16.00.

For the TMB/BRAF-guided strategy the corresponding values are
-0.00000065 QALYs (0.061% of -0.0011) and EUR 1.36. Its larger
pre-decision cost increment is the NGS test (EUR 2,518.00), which the
strategy genuinely incurs.

The costs and QALYs accrued before the first nivolumab dose therefore
contribute almost nothing directly to either strategy’s increments. The
curves have nevertheless separated by week 4: population-averaged OS
differs from standard of care by -0.00000512 (CRP-guided) and -0.0000975
(TMB/BRAF-guided), and PFS by 0.00333 and 0.000199. This separation
carries into the rest of the horizon. The split above does not attribute
its later consequences; that would need the curves re-anchored at the
decision point (option (a)), which is not estimated here.

# Limitations

- **Corner scenarios only.** The four scenarios give every restored
  patient the same CRP value, and every patient with unknown TMB/BRAF
  the same TMB/BRAF value. With 3 unknown CRP values and 2 unknown
  TMB/BRAF values there are 32 combinations. The model is not linear in
  the assignments, so intermediate combinations need not fall between
  the corner results; the reported range covers the corner scenarios
  only. No imputation model or missing-data mechanism is assumed.
- **Distributions not re-selected.** Each scenario keeps the base-case
  distributions and checks only the OS \>= PFS ordering. A new
  ordering-constrained selection on the enlarged cohort could choose a
  different pair.
- **Population and effects change together.** Restoring patients changes
  the target population (size, prevalences and age/sex mix) as well as
  the fitted effects, and the scenarios do not separate the two.
- **Deterministic only.** The scenarios quantify the sensitivity of the
  point estimates; parameter uncertainty is not propagated (no PSA).
- **Other exclusions.** The 3 patients with a week-4 CRP but no TMB/BRAF
  result remain excluded in every scenario.
- **Pre-decision split.** The split before week 4 measures the costs and
  QALYs accrued in that period and the survival difference reached by
  its end. It is not an estimate of option (a), which would refit or
  re-anchor the curves so that the strategies are identical up to the
  decision point.

------------------------------------------------------------------------

**Report completed on:** 2026-09-25  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.0
