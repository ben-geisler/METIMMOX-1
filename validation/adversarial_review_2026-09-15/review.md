# Adversarial codebase review

Date: 2026-09-15  
Reviewed commit: `89b2d040e45b613b407b07b38a6167889c881ea9`

## Assessment

The current code has several reproducible correctness problems despite substantial regression coverage. The most consequential are inconsistent population weights across comparators, a stale publication table, and incomplete cache provenance. These can produce plausible-looking results that represent different populations or different analysis versions.

This review identifies **17 actionable findings: five high-priority and twelve medium-priority**, followed by methodological and operational concerns. Severity reflects potential impact; conditional failure cases are distinguished from defects observed in current outputs. The evidence does **not** establish that fixing these issues would change the preferred strategy at EUR 51,000/QALY.

Production code and publication outputs were not edited. This directory contains the review and reproduction material only.

**Test side effect:** the existing `test_sim_idx_validation.R` sourced the full sampling pipeline, detected a changed `data_hash`, and regenerated `data/tidy/sampling_models_n5000_full.rds`. The PSA/EVPPI caches were not rewritten. A post-suite contract check now fails because the saved PSA references the previous sampling fingerprint. Strict economic report loaders will reject that mismatch until the corresponding results are regenerated or the pre-review sampling cache is restored from a backup. The original cache was not backed up by the test, and its prior bytes were not recoverable from the test run. No cache provenance was manually rewritten to conceal the mismatch. See F17.

## Scope and verification

- Examined model calculations, survival prediction/fitting/sampling, data preparation, PSA/DSA/EVPPI, enriched and scenario analyses, clinical/DAG helpers, report/vignette consumers, snapshots, publishing and regression tests.
- Executed all **18** existing `test_*.R` scripts: **17 passed; one failed in the initial sequential run**. That failure is the previously documented `test_psa_basecase_alignment.R` numerical criterion. The later sampling-cache side effect means a second test now fails on rerun; the initial count is not a clean final-suite result.
- The report/cache contract test initially passed **217 checks**, with two skips: zero-EVPI row checks and the absent current scenario cache. After the index-validation test regenerated sampling, rerunning the contract test produced one failure: the PSA/sampling fingerprint mismatch.
- Parsed **60 live R scripts** and extracted R code from **44 Quarto documents**: no syntax errors. Parsing is not a report-render test.
- Ran targeted adversarial probes using synthetic cases and aggregate results from the local trial data/caches. See [reproduce.R](reproduce.R), [probe_output.txt](probe_output.txt) and [test_results.csv](test_results.csv).
- Forced one assertion result to fail in a temporary copy of the index-validation test; it printed failure but exited successfully.
- Did not regenerate the 5,000-draw PSA/scenario analyses, render all reports, execute publishing, or audit remote deployment. The current `scenario_evppi_results_correct.rds` is absent locally, so current scenario calculations could not be revalidated from that cache. Archived code was treated as historical material, not part of the supported pipeline.

Reproduce the focused probes from the repository root:

```powershell
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' validation/adversarial_review_2026-09-15/reproduce.R
```

The probes require the local trial dataset and current sampling/PSA caches. They read those files without replacing them and print aggregate results only.

## Findings at a glance

P1 = address before relying on refreshed or published results. P2 = correctness issue requiring a fix, often triggered by a supported configuration change or an invalid draw.

| ID | Priority | Finding | Evidence/status |
|---|---|---|---|
| F01 | P1 | Guided and control strategies use different population weights | Same treatment still gives nonzero incremental QALYs |
| F02 | P1 | Tracked Table 5 contains obsolete PSA results | Direct disagreement with current cache and reports |
| F03 | P1 | PSA cache fingerprint ignores executable model logic | Doubling calculated QALYs leaves the fingerprint unchanged |
| F04 | P1 | EVPPI provenance omits WTP and other EVPPI inputs | Same accepted PSA fingerprint supports EVPI 0 or EUR 172.66 at different WTPs |
| F05 | P1 | PSA outcomes and parameter draws are not bound together | Reversed parameter rows pass the parameter-cache loader |
| F06 | P2 | Scenario and publication consumers bypass freshness checks | Static inspection of active loaders and scenario metadata |
| F07 | P2 | Two DAG independence tests test the wrong hypothesis | Synthetic test rejects an independence that holds by construction |
| F08 | P2 | Non-finite PSA results bypass the 2% failure limit | 3% invalid draws accepted with reported fallback rate 0% |
| F09 | P2 | QALY integration counts an extra grid point | 10-year no-mortality case yields 10.01923 QALYs; known issue |
| F10 | P2 | Base-case entry points ignore supplied time settings | Five-year input errors; changed cycle length silently ignored |
| F11 | P2 | Direct diagnostic-price changes can have no effect | NGS price +EUR 1,000 produces zero cost change without manual synchronization |
| F12 | P2 | PSA subgroup and base-case prediction cohorts can diverge | Endpoint-missing patient retained only in PSA subgroup prediction |
| F13 | P2 | Zero-EVPI scenarios disappear from publication outputs | Base scenario absent from Table S7; empty results skip artifact replacement |
| F14 | P2 | Invalid survival probabilities are accepted | OS = 1.1 produces negative dead-state occupancy and finite results |
| F15 | P2 | Short horizons fail inside the follow-up schedule | Fewer than 13 grid points raise a sequence error |
| F16 | P2 | An executable regression test exits successfully on failure | Forced failure prints 5/6 passed and exits 0 |
| F17 | P2 | Regression tests mutate live analysis caches | Observed sampling regeneration leaves downstream cache contracts failing |

## High-priority findings

### F01 — Standardize every strategy to the same population

Locations: [prediction_functions.R:184](../../R/prediction_functions.R#L184), [prediction_functions.R:292](../../R/prediction_functions.R#L292), [model_fun.R:498](../../R/model_fun.R#L498).

Control survival averages the complete-case cohort. Guided strategies average complete-case subgroups but combine them using the wider cohort's biomarker prevalence. Thus the comparator population changes with the strategy. PSA and DSA also vary these weights only for the guided strategies while leaving control's population fixed.

For TMB/BRAF, the canonical prevalence is **0.4492754**, versus **0.4411765** in the complete-case cohort. I assigned control treatment to both biomarker subgroups, equalized treatment schedules, and removed screening costs. The supposedly identical policies still differed:

| Strategy versus control | Incremental cost | Incremental QALYs |
|---|---:|---:|
| CRP | -EUR 0.14 | -0.0001106 |
| TMB/BRAF | EUR 5.14 | 0.0011760 |

Matching the subgroup weights to the control population makes both increments exactly zero. A supported +20% TMB/BRAF prevalence variation creates **0.0142229 artificial incremental QALYs** under identical treatment. This is material relative to the current base-case TMB/BRAF increment of -0.0105636 QALYs.

**Fix direction:** define one target population, including its joint biomarker distribution, and evaluate every policy on it. If prevalence changes, standardize control and all guided strategies to the changed population consistently. Using complete-case prevalence is a useful diagnostic, not an automatic substitute for deciding the intended target population. Add an identical-treatment invariant covering prevalence changes.

### F02 — Regenerate and reconcile the tracked PSA summary table

Locations: [table_5.csv:2](../../outputs/tables/table_5.csv#L2), [EVPPIs.md:31](../../reports/EVPPIs.md#L31), [table_5.qmd:25](../../outputs/vignettes/table_5.qmd#L25).

The committed publication table disagrees with the current PSA and both main economic reports:

| Control result | Tracked Table 5 | Current PSA / main reports |
|---|---:|---:|
| Mean cost | EUR 21,088 | EUR 21,467 |
| Mean QALYs | 1.29 | 1.416 |
| Probability cost-effective | 95.8% | 100.0% |

This is an observed publication defect, not a hypothetical stale-cache case. Table 5 could be used in a manuscript alongside current deterministic and EVPPI tables, producing mutually inconsistent claims.

**Fix direction:** regenerate Table 5 after validating its input cache; compare every tracked numerical publication table with its generating objects. Record source fingerprints in an artifact manifest so recently checked-out file modification times cannot be mistaken for analysis provenance.

### F03 — Include model implementation dependencies in PSA cache validity

Locations: [cache_paths.R:172](../../R/cache_paths.R#L172), [10_PSA.R:100](../../analysis/10_PSA.R#L100).

The PSA fingerprint hashes parameter values, distributions, sampling provenance, strategies, seed and time settings, but not the code that computes outcomes or predictions. Fixing `calculate_outcomes()`, `model_fun()` or a prediction helper can therefore leave the PSA cache accepted as current. This also affects reports using the otherwise stricter `expected_fingerprint` validation.

**Reproduction:** in the probe process only, I changed `calculate_outcomes()` to double QALYs. Control QALYs changed from **1.371420 to 2.742841**, while the expected PSA fingerprint remained identical.

**Fix direction:** hash the relevant implementation files/functions or maintain an explicit, enforced calculation version. Include dependency versions that can change model fitting/prediction behavior. Ensure uncommitted calculation changes invalidate results too; comparing a cache timestamp only with HEAD's commit time does not cover them.

### F04 — Fingerprint EVPPI inputs independently of PSA inputs

Locations: [11_EVPPIs.R:233](../../analysis/11_EVPPIs.R#L233), [EVPPIs.qmd:73](../../reports/EVPPIs.qmd#L73).

The saved EVPPI provenance is only `psa_obj$fingerprint`. The report checks that value, but WTP, parameter groups, estimator configuration and population-scaling inputs are not part of it. WTP correctly does not affect the underlying PSA cost/effect draws; it does affect EVPI and EVPPI.

**Reproduction:** on the same current PSA, EVPI is **EUR 0 at WTP 51,000**, **EUR 172.66 at WTP 100,000**, and **EUR 1,343.47 at WTP 150,000**. Changing WTP in setup without rerunning script 11 leaves the cached zero-EVPI result eligible for use. The report can then calculate a new decision table while repeating the old zero-EVPI narrative.

**Fix direction:** store and validate a separate EVPPI fingerprint containing the PSA result identity, WTP, complete parameter/group specification, estimator version/settings/seed and population parameters. Verify the cached EVPI against a cheap calculation from the loaded PSA and current WTP before rendering.

### F05 — Bind PSA outcome rows to their exact parameter rows

Locations: [10_PSA.R:38](../../analysis/10_PSA.R#L38), [10_PSA.R:209](../../analysis/10_PSA.R#L209), [cea_helpers.R:292](../../R/cea_helpers.R#L292).

`psa_obj` has an input fingerprint; `psa_params` has a seed but no corresponding fingerprint or content identity. Validation checks row count, column presence and seed, but does not establish that these parameter draws produced these outcomes. Presence of `model_idx` is checked without verifying the two files' row alignment. Separate, non-atomic saves allow an interrupted run or a partial cache restore to mix versions.

**Reproduction:** reversing all parameter rows while keeping their metadata passes `load_psa_params_cache(..., seed = analysis_seed, n_sim = psa$n_sim)`. The current cache's indices do match; the defect is the missing enforcement when files are replaced or restored.

This matters especially for EVPPI: an apparently valid regression can relate outcomes to the wrong sampled inputs and estimate a different value of information.

**Fix direction:** write outcomes, parameters and metadata as one atomic cache bundle, or give both files the same generation ID plus a hash of the exact ordered parameter table. Validate draw IDs and `model_idx` equality before any analysis.

## Other correctness findings

### F06 — Apply cache validation to scenarios and publication vignettes

Locations: [12_scenario_EVPPIs.R:175](../../analysis/12_scenario_EVPPIs.R#L175), [scenario_effect.qmd:49](../../reports/scenario_effect.qmd#L49), [figure3.qmd:27](../../outputs/vignettes/figure3.qmd#L27), [figure4.qmd:36](../../outputs/vignettes/figure4.qmd#L36), [figure5.qmd:33](../../outputs/vignettes/figure5.qmd#L33), [table_5.qmd:25](../../outputs/vignettes/table_5.qmd#L25).

The scenario report checks only the failed-draw-policy string. Several vignettes call `readRDS()` directly; Table 5 calls the shared loader without current-input expectations. The scenario cache stamps the sampling fingerprint but lacks a full identity for economic parameters and estimator inputs. Consequently, changing prices, utilities, survival inputs or estimator definitions can leave publication outputs based on stale scenarios or uncertainty draws while CEA.qmd rejects its own PSA.

**Fix direction:** one validated loading path for every consumer, with full scenario-specific fingerprints. Stop on missing or stale inputs. Current scenario-cache execution was unavailable locally; this finding follows directly from the saved fields and active loading conditions.

### F07 — Classify conditional independence involving a fully determined child correctly

Locations: [dag_association_tests.R:88](../../R/dag_association_tests.R#L88), [dag_association_tests.R:336](../../R/dag_association_tests.R#L336).

`TxCRP = T * CRP`. Given both `T` and `CRP`, `TxCRP` is fixed, so `TLR` is conditionally independent of `TxCRP` by construction. The same applies to `TxTMB` given `T` and `TMB_BRAF`. The current classification misses these cases and tests an interaction coefficient in a logistic model instead. A nonzero interaction coefficient is not a violation of the stated conditional independence.

**Reproduction:** a synthetic dataset has exactly one possible `TxCRP` value in every conditioning stratum, yet the repository's logistic test returns **p = 4.65e-11**. Current report p-values of 0.898 and 0.426 therefore do not validate those DAG statements either.

**Fix direction:** classify a statement as definitional whenever either tested variable is completely determined by the conditioning set. Retain interaction tests, if desired, under their actual hypothesis and outside the empirical-CI counts.

### F08 — Count non-finite initial PSA outcomes as failed draws

Locations: [psa_functions.R:232](../../R/psa_functions.R#L232), [psa_functions.R:283](../../R/psa_functions.R#L283), [psa_functions.R:403](../../R/psa_functions.R#L403).

The initial pass counts errors, flagged fallbacks and missing strategy rows. It does not count returned NA/NaN/Inf costs or effects. These are dropped only after the failure-rate check and replacement phase.

**Reproduction:** injecting NA costs into three of 100 initial draws produces `fallback_count = 0`, `fallback_rate = 0`, `dropped_count = 3`, and 97 retained draws despite a 2% limit. Those draws receive no replacement attempts. The same route can silently discard a much larger fraction.

**Fix direction:** validate one finite cost/effect result per strategy immediately after every initial model call, using the same completeness rule as replacements. Include all invalid initial draws in the threshold and replacement logic.

### F09 — Integrate QALYs over intervals, not all grid endpoints

Location: [calculate_outcomes.R:74](../../R/calculate_outcomes.R#L74).

There are 521 survival grid points but only 520 weekly intervals. Summing one week's utility at every point includes a full extra grid contribution. Under no mortality, no discounting and unit utilities, the stated 10-year model returns **10.01923077 QALYs**; a one-year grid returns **1.01923077**. This is already documented in AGENTS.md and remains present.

**Fix direction:** adopt an explicit interval-integration convention, such as trapezoidal integration for state utilities. Keep event and scheduled treatment costs on their appropriate time points; do not blindly apply the same correction to all costs. The report helper `restricted_mean_survival()` already uses trapezoidal integration, so survival totals and economic state rewards currently use different conventions.

### F10 — Honor time settings at base-case entry points

Locations: [model_fun.R:124](../../R/model_fun.R#L124), [cea_helpers.R:62](../../R/cea_helpers.R#L62), [08_basecase_analysis.R:11](../../analysis/08_basecase_analysis.R#L11).

`run_basecase(params)` and script 08 call `model_fun(params)` without passing time settings. `model_fun()` defaults to 520 cycles and `cl = 1/52`, regardless of the parameter list. This makes a globally changed horizon work in some explicitly parameterized callers and fail in the base case.

**Reproduction:** a valid five-year parameter list built by `build_horizon_params()` errors because its 261-point curves are checked against 521 points. Setting `params$cl = 1/12` leaves wrapper control QALYs at **1.371420**, whereas the explicit call returns **4.979967**. This is an API-consistency demonstration, not a clinically validated monthly model.

**Fix direction:** resolve horizon and cycle length once and pass them consistently. If arbitrary cycle lengths are unsupported because survival times and treatment schedules are weekly, reject them explicitly rather than silently ignoring or partially applying them.

### F11 — Remove the unsynchronized duplicate diagnostic-price input

Locations: [calculate_outcomes.R:94](../../R/calculate_outcomes.R#L94), [model_fun.R:45](../../R/model_fun.R#L45).

The parameter list stores diagnostic prices twice: scalar `c_test_NGS`/`c_test_CRP` and `c_test_biomarker`. The model validates the scalars but bills from the lookup. DSA and PSA explicitly synchronize the lookup; direct model/base-case callers do not.

**Reproduction:** raising `c_test_NGS` by EUR 1,000 changes no result. Calling `sync_biomarker_test_costs()` first produces the expected EUR 1,000 TMB/BRAF cost increase.

**Fix direction:** resolve prices from one canonical representation within the public calculation boundary, or assert that duplicates agree. Documented interactive parameter edits should not require knowing a hidden synchronization step.

### F12 — Use the same prediction cohort in deterministic and PSA paths

Locations: [model_fun.R:250](../../R/model_fun.R#L250), [04_parametric_survival_analysis.R:63](../../analysis/04_parametric_survival_analysis.R#L63).

PSA subgroup predictions receive `data`; control and deterministic predictions receive `data_complete`. Missing predictor values often conceal the mismatch because those predictions become NA. A patient missing an endpoint but having complete predictors can still be predicted in the PSA subgroup and is excluded from the other paths.

**Reproduction:** mark one CRP-positive patient's OS time missing. With identical fitted coefficients, including that patient changes the positive OS curve by up to **0.0038166** compared with the complete-case prediction set.

**Fix direction:** pass an explicitly defined common prediction population everywhere. If fitting and prediction populations intentionally differ, define both separately and apply the same prediction population to every strategy and uncertainty mode. This is conditional on endpoint missingness; it is not offered as the explanation for the current PSA/base-case mean gap.

### F13 — Preserve zero-EVPI scenarios in tables and clear superseded artifacts

Locations: [scenario_analysis.R:287](../../R/scenario_analysis.R#L287), [table_s7.qmd:77](../../outputs/vignettes/table_s7.qmd#L77), [table_s7.qmd:148](../../outputs/vignettes/table_s7.qmd#L148), [figure5.qmd:66](../../outputs/vignettes/figure5.qmd#L66).

Zero EVPI legitimately produces an empty EVPPI result. Compilation drops that scenario; Table S7 derives even its total-EVPI rows from the nonempty EVPPI table. The current tracked Table S7 therefore contains five scenarios and omits the base case. Figure 5 requests a base-case series but has no rows to draw for it.

If all estimates become empty, conditional output writes also leave any previous CSV/PNG in place. A successful render can consequently coexist with an obsolete publication artifact.

**Fix direction:** derive scenario presence and total EVPI from the scenario result objects, retaining explicit zero results. Distinguish legitimate zero from missing/failed estimates. Always replace the generated artifact with an explicit empty-state artifact or remove its obsolete predecessor according to a documented output contract.

### F14 — Validate probability ranges, finiteness and monotonicity before costing

Locations: [model_fun.R:316](../../R/model_fun.R#L316), [calculate_outcomes.R:10](../../R/calculate_outcomes.R#L10), [calculate_outcomes.R:129](../../R/calculate_outcomes.R#L129).

The model checks curve length and OS/PFS ordering but not whether survival probabilities are finite, within [0,1], or nonincreasing. Initial-state overrides and clipping negative death transitions do not repair invalid later probabilities.

**Reproduction:** setting `control_OS[2] = 1.1` is accepted. Dead-state occupancy becomes **-0.1**, while the model returns finite costs and QALYs. Clipping subsequent death transitions can then add spurious terminal-care costs.

**Fix direction:** validate the full survival contract before constructing states; route invalid PSA predictions through failed-draw handling. Also require finite scalar costs/rates and valid schedules, since checking only `< 0` permits some non-finite inputs.

### F15 — Handle horizons shorter than the first quarterly payment

Location: [calculate_outcomes.R:71](../../R/calculate_outcomes.R#L71); analogous unguarded sequences occur in [parameter_distributions.R:347](../../R/parameter_distributions.R#L347).

`seq(13, n_cycles, by = 13)` errors when `n_cycles < 13`; the later `length(quarterly_cycles)` guard is never reached. The calculation therefore cannot represent short valid horizons or early-horizon boundary tests.

**Reproduction:** 1, 2, 5 and 12 grid points all fail with `wrong sign in 'by' argument`.

**Fix direction:** return an empty schedule when the first due payment lies beyond the horizon, and use the same rule for rebuilt CT/blood schedules. Audit the quarterly offset too: positions 13, 26, ... correspond to weeks 12, 25, ..., not weeks 13, 26, ... on the declared zero-origin grid.

### F16 — Make regression assertion failures return a failing process status

Location: [test_sim_idx_validation.R:165](../../tests/test_sim_idx_validation.R#L165).

The test prints a failure summary but does not call `stop()` or `quit(status = 1)`. An automated runner relying on process status reports success even if the validation cases fail.

**Reproduction:** in a temporary copy, replace the valid-case result assignment with `FALSE`. The script prints **Passed: 5 / 6**, followed by **Some tests failed**, and exits **0**. The unmodified test did pass its current cases.

**Fix direction:** assert `all(results)` or explicitly exit nonzero after the failure summary. Keep expected skips distinct from passes.

### F17 — Isolate regression tests from production caches

Location: [test_sim_idx_validation.R:31](../../tests/test_sim_idx_validation.R#L31), calling the load/regenerate/save blocks in [06_sampling.R:277](../../analysis/06_sampling.R#L277).

The index-validation test sources the entire sampling pipeline against the real `data/tidy` directory. During this review it printed `Cache fingerprint does not match the current inputs`, identified `data_hash`, generated 5,000 draws, and overwrote the sampling cache. A test that previously passed the report/cache contract thus made that same contract fail later in the suite.

Observed file change: the sampling cache went from 700,411 bytes, modified 2026-09-08, to 700,408 bytes, modified 2026-09-15. The new fingerprint is `033e43f731c26d684837db6ddb68d5bd`. The cached PSA, EVPPI and violation diagnostics retain their earlier sampling provenance. This finding establishes the test side effect; it does not establish why the pre-existing data hash differed or that numerical draws changed.

A follow-up numerical check reproduced all five saved worst-case diagnostic example curves from the regenerated cache with **exactly zero OS and PFS differences**. That is evidence of numerical continuity for those examples, not proof of equality across all draws, and does not remove the provenance mismatch.

**Fix direction:** use a synthetic fixture or a temporary cache directory and extract/source only the functions required by the test. Snapshot and restore global configuration within tests. A validation run should not invalidate the analyst's live result files. Until isolated, run these tests on copied caches or an isolated workspace.

## Methodological concerns and accepted scope choices

These deserve attention but should not be confused with confirmed implementation errors above.

1. **The PSA mean and deterministic point estimate answer different questions.** The current failing test finds CRP incremental QALYs of **0.020418 deterministically versus 0.002040 in the PSA**, a 7.24-Monte-Carlo-SE gap. All three mean QALYs exceed their deterministic values. Nonlinear survival extrapolation means `E[f(theta)]` need not equal `f(E[theta])`; a discrepancy alone does not prove erroneous sampling. Preserve and investigate the known failure, validate equality when coefficient/parameter uncertainty is removed, and explain both result types. Do not recenter draws or loosen the test merely to obtain a pass.

2. **Independent OS/PFS draws and frequent survival crossings remain material assumptions.** The retained diagnostics report 2,422/5,000 draws with at least one crossing. This independence is an explicit user decision. The clamp changes the sampled survival distribution and should remain visible in uncertainty reporting. Existing diagnostics show its average incremental QALY impact is small relative to the full PSA/base-case gap; blaming that gap entirely on clipping would be unsupported.

3. **Prevalence uncertainty is more than two independent beta draws.** Both biomarkers are measured in the same small cohort and enter a shared treatment model. Independent marginal prevalence changes with fixed subgroup prognostic distributions do not describe a coherent change in their joint distribution. F01 is the directly demonstrated consequence for comparator standardization. Decide whether prevalence uncertainty represents sampling uncertainty or transport to a different target population and implement the corresponding common population weights.

4. **Delayed treatment and the week-4 CRP decision need an explicit estimand.** Treatment schedules are identical before nivolumab, but the fitted Rx effects operate from time zero. Equal schedules alone do not guarantee equal pre-week-4 expected costs or QALYs when state occupancies differ. Complete-case exclusion of patients without the week-4 measurement also deserves a missingness/selection sensitivity analysis. This review does not relabel CRP as baseline or assume randomization failed.

5. **The week-9 TLR landmark remains an exploratory analysis with residual timing issues.** Keep the user-selected primary landmark and its disclosed counts of late scans and retained first-scan progressors. Correct the claim that week 12 is after the latest scan: the same source reports a maximum of 12.1 weeks ([tlr_landmark.R:15](../../R/tlr_landmark.R#L15), [clinical_effectiveness.qmd:1264](../../reports/clinical_effectiveness.qmd#L1264)). Endpoint-specific OS and PFS landmark cohorts also support different event-free statements; avoiding progression before the landmark is a condition of the PFS cohort, not automatically of the OS cohort.

6. **The poster's single-versus-joint CEAC comparison changes covariance as well as model specification.** [poster_biomarker_correlation.qmd:385](../../outputs/vignettes/poster_biomarker_correlation.qmd#L385) pairs joint-model control draws with separately generated single-biomarker treatment draws in one panel, but shares joint draws across strategies in the other. Differences in cost-effectiveness probability therefore cannot be attributed solely to including both biomarkers in the mean survival model. Quantify this additional uncertainty-design difference before interpreting the panel contrast as the effect of biomarker adjustment.

7. **Firth is described too strongly.** The statements that it gives the best unbiased point estimate ([clinical_effectiveness.qmd:910](../../reports/clinical_effectiveness.qmd#L910), [clinical_effectiveness.qmd:1256](../../reports/clinical_effectiveness.qmd#L1256)) exceed what the method establishes. Describe it as penalized likelihood with bias reduction; the package's own documentation uses that characterization. [Official coxphf documentation](https://cran.r-project.org/web/packages/coxphf/coxphf.pdf). Check convergence warnings and returned iteration diagnostics rather than assuming penalization guarantees every fit and interval succeeded.

8. **Important uncertainty and cost omissions remain deliberate.** Fixed assumed nivolumab acquisition cost, omitted adverse-event effects, base-case omission of post-progression treatment costs and universal second-sequence costing are documented choices, not newly discovered bugs. They constrain what the reported zero EVPI means: zero under the represented strategies, distributions, prices and scope, not zero value of all possible future research. The literature-utility scenario still contains an explicit citation placeholder in [scenario_analysis.R:311](../../R/scenario_analysis.R#L311).

## Additional operational concerns

- **Publishing needs failure checks and a clear replacement policy.** [publish_reports.R:81](../../publish/publish_reports.R#L81) checks only for a local `gh-pages` branch, not an existing remote branch. Worktree creation, checkout, removal of old content, copying, staging and committing have unchecked return values; only the final push is checked. Copying new files over an existing worktree also retains PDFs removed from the new `_site`. These paths were reviewed statically; no publishing commands were executed.
- **`--no-render` can misstate PDF provenance.** The index is stamped with the current HEAD and current date even though existing PDFs may have been rendered earlier ([publish_reports.R:57](../../publish/publish_reports.R#L57)). Store actual render provenance with artifacts.
- **The environment is not locked.** Analysis scripts install/load packages dynamically with `pacman`; no checked-in lockfile pins the complete runtime. This limits clean-machine reproducibility and interacts with F03. The observed Windows R session also emitted locale/package-build-version warnings; those did not cause the reviewed test failure.
- **The raw export contract is fragile.** CRP, scan date and lesion values depend on repaired positional Excel names such as `CRP...540` and `Date...122`. Add explicit column/type/range and unique-patient checks before fitting; a revised spreadsheet can otherwise change the meaning of derived biomarkers. This review did not independently re-adjudicate the original clinical endpoints.
- **The PSA mode has an intentional but hazardous overload.** `determpsa = "psa", sim_idx = NULL` returns supplied-curve results, with `fallback_used = FALSE`; it is identical to the deterministic calculation when given base curves. The poster uses this deliberately. Make supplied-curve evaluation explicit and reject misspelled modes, so an accidentally omitted index cannot silently remove survival uncertainty.
- **Warnings can disappear from publication workflows.** Setup and document settings suppress warnings. Fit failures, missing estimate groups and skipped outputs should be persisted in a visible validation summary rather than depending on console warnings that readers never see.

## Suggested repair order

1. Resolve common population standardization (F01) and specify its uncertainty treatment.
2. Repair cache identities and validated loading (F03–F06) before generating replacement outputs.
3. Fix model invariants and time/price handling (F08–F12, F14–F15), with focused regression tests.
4. Correct the DAG hypotheses and isolate/enforce test behavior (F07, F16–F17).
5. Regenerate affected PSA/EVPPI/scenario results once, then regenerate and reconcile publication outputs, including explicit zero scenarios (F02, F13).

Retain the current documented methodological choices unless they are deliberately changed. Separate those decisions from mechanical corrections so impact comparisons remain interpretable.
