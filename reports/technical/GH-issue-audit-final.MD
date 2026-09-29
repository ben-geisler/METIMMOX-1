# GitHub issue audit: adjudication of Reports A and B

METIMMOX-1 Economic Evaluation, 29 September 2026. Adjudicated at HEAD `9cedd6b` (working tree unchanged except the three untracked audit inputs and the deleted `GH-issue-audit.md`).

## Method

Every one of the 23 issues in `reports/technical/issues.md` was checked against the repository before either report's reasoning was read for that issue. Evidence is file:line at HEAD, a commit hash, or a reproduction. The reproductions for #60 and #172 were run with a read-only synthetic script (no trial data, no cache, no repository write) against `R/model_configs.R`, `R/calculate_outcomes.R` and `R/model_fun.R`; its exact inputs and outputs are given under the verdicts. Issue bodies and comments were read with `gh issue view --json`. Nothing in the repository was modified, no cache was regenerated, no report was rendered.

The ranking rubric was not pasted into the request (the placeholder text was left in). The ranking therefore uses the three criteria named in the placeholder: user impact, blocking dependencies, effort.

## Verdict table

Verdicts: VALID = the issue still applies as written. PARTIALLY VALID = part of the reported behaviour is fixed, part remains. NOT VALID = resolved or correctly closed. UNVERIFIABLE = the repository cannot ground a verdict.

| Issue | A verdict | B verdict | Final verdict | Evidence | Confidence |
|---|---|---|---|---|---|
| #178 | Partly fixed, open | Partly resolved, open | **PARTIALLY VALID** | Render exit status is checked (`publish/publish_reports.R:41-42`). Still unchecked: worktree add/checkout/rm/copy/add/commit (`:83-91`); only the local branch is tested (`:81`); `file.copy(..., overwrite = TRUE)` onto an existing worktree never removes stale PDFs (`:89`); the index stamps `Sys.Date()` and current HEAD even with `--no-render` (`:57`, `:68`); `setup_report()` sets `warning = FALSE` and wraps sourcing in `suppressWarnings()` (`R/report_setup.R:66`, `:95`) with no warning summary. | High |
| #177 | Applies | Keep open | **VALID** | `outputs/vignettes/poster_SMDM.qmd:339-348`: single-biomarker models drawn with `sample_survival_coefficients()` and seeds 2027/2028. `:391-395`: control curve from `sm_joint` (sampling cache, `model_idx[r]`); the single panel's guided curves come from `sm_crp`/`sm_tmb` indexed by loop position `i`, unrelated to `r`. The joint panel uses `sm_joint` for all three strategies. Caption at `:423` says the panels "differ only in the guided-strategy curves", which the draw design contradicts. | High |
| #174 | Partly fixed, open | Partly resolved, open | **PARTIALLY VALID** | Commit `791443e` (#181) added header/type/range/0-1 checks for the five endpoint columns (`analysis/01_data_prep.R:59-77`). Unchecked positional columns remain: `TL LD...130`, `TL LD...121` (`:41-42`), `CRP...478` (`:47`), `CRP...540` (`:54`), `Date...122` (`analysis/02_setup_and_global_variables.R:24`). No `ID` uniqueness, row-count or biomarker-count assertion exists in script 01. | High |
| #172 | Partly fixed, open | Partly resolved, reproduced | **PARTIALLY VALID** | Fixed: commit `d3480ec` (#163) removed `quarterly_cycles <- seq(13, n_cycles, by = 13)` from `R/calculate_outcomes.R`; my 13-point probe returns 12/52 = 0.230769 QALYs, the value the issue said was wrong (0.25). Remaining, reproduced: `model_fun(params)` with a 261-point list fails with "length 261 but expected 521" (`R/model_fun.R:111` defaults, `:318`, `:326-331`); `run_basecase()` calls `model_fun(params)` bare (`R/cea_helpers.R:62`), as does `analysis/08_basecase_analysis.R:11`; `params$cl = 1/12` alone gives 10 QALYs, explicit `cl = 1/12` gives 43.33; `seq(13, n, by = 12)` and `seq(5, n, by = 4)` error for n < 13 / n < 5 (`analysis/05_basecase_input_parameters.R:29,33`; `R/parameter_distributions.R:361-362`). | High |
| #49 | Absorbed by #172 | Duplicate of #172, residual applies | **PARTIALLY VALID (duplicate of #172)** | `build_horizon_params()` truncates treatment positions to the horizon (`R/parameter_distributions.R:356-360`) and sets `p$time_horizon` (`:365`), so the "schedules ignore the horizon" part is done for scenarios. Residual: the same short-horizon `seq()` failures as #172, and script 05 assigns positions up to 39 directly (`:19,22,25`), which silently lengthens a vector shorter than 39 (probe: length 30 became 39). The issue's own 16 September comment defers closure to #172. | High |
| #60 | Applies in full | Reproduced | **VALID** | Only curve length (`R/model_fun.R:326-331`) and PFS <= OS (`:335-372`) are validated; `partitioned_survival_states()` overwrites the first point without a check (`R/calculate_outcomes.R:22-24`). Probe: `control_OS[2] = 1.1` accepted, dead occupancy -0.1, finite QALYs; curves starting at 0.8 accepted and overwritten to 1; an OS curve that rises from 0.5 to 0.9 accepted. Cost parameters now do have a finite-scalar check (`R/model_fun.R:77-79`); discount rates and utilities still only have sign/range tests (`:86-106`). | High |
| #36 | Applies; results change | Open input question; cannot verify EUR 33 is wrong | **VALID as an unresolved input; correctness UNVERIFIABLE** | `c_other_visit = 33` and `c_other_follow = 33` (`analysis/05_basecase_input_parameters.R:96,98`) with the comment "cost of standard outpatient visit"; no "DRG" string anywhere in `R/`, `analysis/`, `reports/` or `outputs/vignettes/`. The DRG weights are not in the repository, so whether 33 is wrong, and by how much, cannot be decided here. | Medium |
| #164 | Partly obsolete, open | Cadence half obsolete, open | **PARTIALLY VALID** | Point 1 (12- vs 13-week cadence) gone with `d3480ec`: `c_other_follow` is a quarterly rate integrated over progressed occupancy (`R/calculate_outcomes.R:119`). Point 2 stands: `l_CT` and `l_blood` run to the horizon (`05:29,33`), `l_visit` is 1 only at position 1 and at administrations (`05:35-37`), rebuilt the same way in `build_horizon_params()` (`:361-364`); progressed patients get no imaging (disclosed at `reports/input_parameters.qmd:270`). | High |
| #56 | Narrowed, open | Narrower, open | **PARTIALLY VALID (narrowed)** | CV 0.20 for the four resource-use costs (`R/parameter_distributions.R:45-46`), CV 0.15 for `u_np` and `u_decrement` (`:60`). `reports/input_parameters.qmd:178` states "CV 0.15" without a source; no citation for either CV in code or report. Unit prices are `psa = FALSE` (`:35`), so the original "costs" scope is narrower. | High |
| #161 | Applies | Keep open | **VALID** | `R/scenario_analysis.R:351` holds the literal placeholder; `outputs/tables/table_s8.csv:38-40` prints it. | High |
| #162 | Applies | Keep open; partial evidence exists | **VALID** | `docs/` has nothing; `reports/technical/psa_extrapolation_plausibility.md:70-71` and `:248-249` call the general-population comparator "too permissive to establish cancer-specific extrapolation validity". No file under `reports/`, `R/`, `analysis/` or `outputs/vignettes/` mentions NORDIC-VII, NORDCAN or registry survival (the only "Cancer Registry" hit is the incidence source comment at `analysis/02_setup_and_global_variables.R:55`). No fixed-year model-versus-KM table exists in `survival_model_specification.md` or `para_models.md`. | High |
| #160 | Applies | Keep open | **VALID** | `docs/manuscript/` and `docs/references/` contain only `.gitkeep`; commit `a3befd2` deleted the manuscript. The issue body's numbers are stale against `outputs/tables/table_5.csv` (SoC EUR 20,931.85 / 1.37753; CRP EUR 53,228.19 / 1.40075, ICER 1,390,875.92; TMB/BRAF dominated). B's extra finding is confirmed: `reports/CEA.qmd:88` says the randomisation-time clock "has no economic consequence" and that the three patients "died early" and the curves are "conditional on surviving to the week-4 measurement", whereas AGENTS.md and the estimand report record that only one of the three died before week 4 and that pre-week-4 increments are EUR 0.50 and 0.000013 QALYs (small, not zero). | High |
| #58 | Mostly fixed, open | Mostly reduced, open | **PARTIALLY VALID (mostly fixed)** | `strategy_label()` at `R/report_format.R:14`, used by 14 files. Literal label keys remain: `figure3.qmd:60-61`, `figure4.qmd:72-74`, `figure_s5.qmd:47-49`, `poster_SMDM.qmd:438-439`, `OWSA.qmd:216` (`short_strategy`), `figure_pfs_os_curves.qmd:113,181` ("Standard of care"). A change in `strategy_label()` would leave those colour maps silently unmatched. | High |
| #165 | Applies in full | Keep open; refresh inventory | **VALID** | No `CHANGELOG.md`, `tests/run_all.R`, `.Rprofile` or `renv/`; `renv.lock` and `DESCRIPTION` exist. 27 `tests/test_*.R` files (the body says 18). `tests/test_psa_basecase_alignment.R:1-5` and commit `2a48ae8` confirm the known-failure entry is obsolete. Package versions are already in the PSA fingerprint (`R/cache_provenance.R:28-33`; flexsurv, survival, dampack, mvtnorm only), so that sub-item of the 16 September comment is done. | High |
| #159 (closed) | Closed correctly | Keep closed | **NOT VALID (correctly closed)** | Commit `9b0049b` implemented joint OS/PFS sampling (`R/joint_survival_sampling.R`; `analysis/06_sampling.R:5`, `:79` `SAMPLING_METHOD <- "mvn_joint_v2"`). Commit `2a48ae8` (#180) made probabilistic results primary (`reports/CEA.qmd:140` "Primary Analysis: Probabilistic" precedes `:189` "Deterministic Point Evaluation") and retired the five-MCSE criterion (`tests/test_psa_basecase_alignment.R:1-5`). The 23 September comment set "closes when #180 lands" as the condition; #180 is closed. The title's premise (a centring defect) was not fixed but re-classified as a different estimand; that re-classification is the author's documented decision, recorded in AGENTS.md. | High |
| #124 (closed) | Superseded | Superseded / partly implemented | **NOT VALID (superseded)** | Models A/B/C removed by `70d672d` (#151). `figure1.qmd` has 500-draw 2.5/97.5% ribbons (`:124`, `:171-177`) and KM overlays (`:197-200`); `figure2.qmd` to `figure5.qmd` exist; `suppl_figure_pfs_plots.qmd` was deleted in `30d6629`. The ribbons are coefficient-draw percentiles, not the patient-bootstrap CIs the spec asked for; that follows from the later method change, not an omission. | High |
| #66 (closed) | Resolved | Resolved | **NOT VALID (resolved)** | 21 of 21 `.qmd` under `reports/` and `reports/technical/` carry `subtitle: "METIMMOX-1 Economic Evaluation"` and "Report version"; 0 contain `fancyhdr` or `include-in-header`. Source-level check only; no PDF was inspected visually. | High |
| #16 (closed) | Resolved | Resolved as a presentation request, with caveat | **NOT VALID (resolved as written)** | `R/report_setup.R:66` (`warning = FALSE, message = FALSE`) and `:95` (`suppressMessages(suppressWarnings(...))`). None of the 24 rendered `.md` files under `reports/` and none of the 26 vignette `.html` files contains a `## Warning` or `## Error` block (the only "Warning" hit is inside `GH-issue-audit-B.md`, which is not a report). Warnings are suppressed, not eliminated; that residual is #178. | High |
| #6 (closed) | Resolved | Resolved | **NOT VALID (resolved)** | `R/report_setup.R:70` sets `root.dir = here::here()`; no `setwd()` in `reports/`, `outputs/vignettes/`, `R/` or `analysis/`. | High |
| #29, #31, #33, #34 (closed) | Not planned | Not planned, correct | **NOT VALID (closed as not planned; author decision)** | `analysis/02_setup_and_global_variables.R:76` `UTILITY_SOURCE <- 1`; `analysis/05_basecase_input_parameters.R:72-73` selects 0.73/0.59. `archive/05_QALYs.Rmd` and `R/eq5d5l_utility.R` retained. A scope decision cannot be "verified" beyond confirming that the code is consistent with it; it is. | High (consistency), n/a (merit) |

## Verification of every NOT VALID / already-fixed claim

Each item below names the resolving change and states whether it addresses the reported behaviour.

- **#159**: `9b0049b` (joint covariance) plus `2a48ae8` (probabilistic-primary reporting, criterion retired). Addresses the closure condition the author set on 23 September. It does not make the PSA means equal the deterministic values; the author's position, recorded in AGENTS.md and the test header, is that they need not. Both reports state this correctly.
- **#124**: `70d672d` (single model), `30d6629` (deleted `suppl_figure_pfs_plots.qmd`), figure vignettes present with ribbons/KM/tornado/plane/CEAC/EVPPI. Addresses the specification in single-model form.
- **#66**: verified by scanning all 21 `.qmd` files. Addresses the request at source level.
- **#16**: `R/report_setup.R:66,95` and a clean scan of 24 `.md` and 26 `.html`. Addresses the literal request ("get rid of the messages"); does not address root causes, which the closure comment explicitly hands to #178.
- **#6**: `R/report_setup.R:70`. Addresses the request.
- **#29/#31/#33/#34**: not fixes; `UTILITY_SOURCE = 1` makes them out of scope. Consistent with the code.
- **#172, quarterly part**: `d3480ec` removed the `seq(13, n, by = 13)` schedule; my probe reproduces the corrected 13-point value (12/52). Addresses F15's follow-up-charge and QALY-overstatement claims; not the schedule-builder or entry-point claims.
- **#164, point 1**: same commit `d3480ec`; the 13-week cadence no longer exists. Addresses point 1 only.
- **#174, endpoint columns**: `791443e` added `01_data_prep.R:59-77`. Addresses the five endpoint columns only.
- **#178, render status**: `publish/publish_reports.R:41-42`. Addresses one of the five listed gaps.
- **#58, central labels**: `R/report_format.R:14` and 14 consumers. Addresses the report tables; not the colour maps.
- **#49, horizon truncation**: `R/parameter_distributions.R:356-365`. Addresses schedule truncation for the horizon scenarios only.
- **#165, package versions in fingerprint**: `R/cache_provenance.R:28-33`. Addresses the "record package versions" sub-item of the 16 September comment; the four main items remain.
- **#60, cost finiteness**: `R/model_fun.R:77-79`. Addresses the "a `< 0` check alone admits non-finite inputs" remark for costs only; not for rates or utilities, and not the curve contract.

## Disagreements resolved

A and B reach the same verdict on all 23 issues. The differences are in emphasis, and none required choosing one report over the other on evidence.

**#49.** A calls it "absorbed by #172"; B calls it "a dependent tracking issue or duplicate" and adds that script 05's direct position assignments can lengthen a short vector. Both are supported. B's extra point is real (probe: a 30-element vector became 39 elements), and it is an additional reason the fix belongs in #172's change set rather than a separate one.

**#16.** A closes it outright; B closes it "as a presentation request, not as proof that underlying warnings have disappeared". The repository supports B's framing: `setup_report()` suppresses warnings rather than removing their causes. A's closure note also points the residual to #178, so the two reports do not actually conflict.

**#162.** A says "applies"; B says the "no external comparison" wording is "too broad now" because #179 exists. Both hold. #179 is a general-population screen that its own text declares insufficient for the purpose #162 asks for, so the issue applies in full regardless of the wording.

**#165 and cache invalidation.** A argues that activating renv can change package versions and therefore invalidate the PSA caches, so the renv part should precede the single regeneration. The fingerprint does include package versions (`R/cache_provenance.R:28-33`), but only for flexsurv, survival, dampack and mvtnorm. A's ordering advice is sound only if `renv::restore()` changes one of those four; the lockfile documents 284 packages, and whether the installed versions of those four differ from the lockfile was not checked here. B does not address this ordering.

**#36.** A marks "Results: Yes"; B says the repository cannot show EUR 33 is wrong. Both are right at their scope: the result changes if the input changes, and the repository cannot say whether it should.

## Claims in A or B that the repository does not support

- **A, #6: "no `setwd()` anywhere."** `publish/publish_reports.R:29` calls `setwd(root)`. It is a command-line script, not a report, so #6 is still resolved, but the literal claim is false.
- **A, #165 (implicit) and issue #165's 16 September comment: "`renv.lock` (387 packages)".** `renv.lock` at HEAD has 284 `"Package":` entries. Neither report re-counted; the figure in the issue is stale or counted differently.
- **B, executable record: "593 checks passed" in `tests/test_report_contracts.R`.** Not rerun here; AGENTS.md records 576 on 25 September. Not contradicted (the test file may have grown), but unverified.
- **B, #16: "two Hessian/covariance warnings" during the zero-uncertainty test.** Not rerun here; unverified.
- **Issue #172's body, repeated in neither report as current: "13 points return 0.25 QALYs."** The probe now returns 0.230769; this specific number is fixed by `d3480ec`. Both reports say the overstatement is gone, which the probe confirms.
- **Issue #165's body: "the 18 tests."** There are 27 `test_*.R` files; B corrected this, A did not mention it.

Line references in both reports and in the 29 September issue comments were checked and hold at HEAD (`publish_reports.R:41-42,57,68,81,83-91`; `poster_SMDM.qmd:339-344,423,439`; `01_data_prep.R:59-77`; `02_setup_and_global_variables.R:24`; `cea_helpers.R:62`; `08_basecase_analysis.R:11`; `parameter_distributions.R:46,60,361-362`; `05_basecase_input_parameters.R:27-37`; `scenario_analysis.R:351`; `figure3.qmd:60`; `figure4.qmd:73`; `figure_s5.qmd:48`; `OWSA.qmd:216`).

## Top-7 ranking (user impact, blocking dependencies, effort)

Issues that share a root cause are ranked as one item.

1. **#160 Manuscript.** The deliverable; nothing else has user impact until it exists. Blocked by #161 and #162 and by any regeneration from items 2 to 4. Large effort. Also fix the `CEA.qmd:88` wording before copying prose from it.
2. **#172 with #49 (time-settings contract).** Reproducible silent inconsistency between entry points and hard failures below 13 grid points; blocks short-horizon boundary tests. Small effort (`model_fun()` reads `params$time_horizon`/`params$cl`; empty schedule when the first due point exceeds the horizon, in both `05_basecase_input_parameters.R` and `build_horizon_params()`). Touches `R/model_fun.R` and `R/parameter_distributions.R`, which are in the PSA fingerprint (`R/cache_provenance.R:22-30`), so it invalidates the PSA, EVPPI and scenario caches.
3. **#60 (survival-curve contract).** Reproducible silent acceptance of invalid curves (negative dead occupancy). Same file and same tests as item 2; fix together so the caches are invalidated once. Small effort.
4. **#36 + #164 + #56 (visit cost, post-treatment visits, cost/utility CVs).** The only group that changes headline numbers. Blocked on inputs outside the repository (DRG weights, protocol/practice, CV sources). #164's effect size is a multiple of #36's unit cost, and #56 sets the CV of the same parameter, so decide all three before one regeneration. Note that the CT/blood/visit schedules are built in both `05:27-37` and `build_horizon_params():361-364`; item 2 edits the same lines.
5. **#162 External validation.** Needed for the manuscript methods; independent of items 2 to 4 (depends only on the fitted models, which none of them changes), so it can run in parallel. Moderate effort; the external sources and tolerance must be fixed in advance.
6. **#161 Table S8 citation.** Trivial effort, blocks #160 and Table S8. Requires a literature decision, not code.
7. **#174 Raw export contract.** No current impact (the cohort is unchanged), but it guards every cache built after it, and the checks are cheap. Do it before the regeneration in item 4 so the final caches are built from asserted data.

Not in the top 7, in order: #178 (publishing robustness, before the next external publish), #58 (colour maps; re-render happens after item 4 anyway), #165 (hygiene, after the code settles), #177 (poster only; close as not planned if the SMDM poster is not shown again). #58 and #177 both edit `poster_SMDM.qmd` and can be done in one pass.

## Items not verified here and what a human should check

- **#36**: whether EUR 33 is wrong and what Eline's DRG weights imply. Not in the repository.
- **#56**: whether a citable source exists for CV 0.20 and 0.15. Literature question.
- **#161**: the reference for 0.80/0.65. Literature question.
- **#162**: which external mCRC source and tolerance to use. Not decidable from the repository.
- **B's test counts** ("593 checks", two Hessian warnings) and whether `renv::restore()` would change the four fingerprinted package versions. Rerunning `tests/test_report_contracts.R` and `renv::status()` answers both; neither was run here because the instruction was read-only.
- **#66 and #16 at PDF level**: the checks were on `.qmd` sources, rendered `.md` and vignette `.html`. No PDF was opened.
- **#178 deployment behaviour**: static inspection only; the `--push` path was not exercised in either report or here.
- **#159's re-classification**: that the remaining PSA-versus-deterministic gap is acceptable is an author decision recorded in AGENTS.md, not something the code can confirm. A reviewer who disputes the decision should reopen it as a methods question, not as a defect in `9b0049b`.

---

**Report completed on:** 2026-09-29\
**Repository:** ben-geisler/METIMMOX-1\
**Report version:** 1.0
