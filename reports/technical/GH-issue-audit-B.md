# GitHub issue audit

**METIMMOX-1 Economic Evaluation — 29 September 2026**

## Assessment

All **14 open issues** retain some applicable scope, but they do not represent 14 independent tasks. **#49 duplicates the remaining horizon work in #172**. Several other descriptions mix outstanding work with defects already resolved. I would keep 13 substantive issue scopes, with #49 treated as a dependent tracking issue or duplicate rather than a separate job.

I agree with closing all nine requested closed issues, with qualifications: **#16 is resolved as a report-presentation request, not as proof that underlying warnings have disappeared**; #124 is superseded in part by the single-model design; and #29, #31, #33 and #34 are valid **not-planned** closures following the decision to use CORRECT utilities. None needs reopening on the evidence reviewed.

The most useful sequence is to establish the data and model input contracts first, resolve the remaining cost/utility assumptions and external validation next, and then finalize the manuscript and publication outputs. Manuscript drafting and external-source work can start immediately. Neither the separate poster comparison nor every repository-hygiene task needs to block the core paper.

## Scope and verification

This audit covers every open issue returned by GitHub on 29 September 2026: **#36, #49, #56, #58, #60, #160, #161, #162, #164, #165, #172, #174, #177 and #178**, plus closed **#159, #124, #66, #16, #6, #29, #31, #33 and #34**. Issue bodies and comments were retrieved using `gh`; repository evidence was inspected at HEAD **`9cedd6ba8e8492d7f3e7ba372037d04bdf65665b`**. Issue links below identify the discussion assessed; local links identify the implementation or output checked.

The conclusions were independently checked against current code and rendered artifacts, supplemented by fresh executable probes and focused tests. Earlier audit documents were not used as evidence, although GitHub comments referring to an earlier audit were read as part of issue context. Existing uncommitted audit files and the deletion of `GH-issue-audit.md` were left untouched. No GitHub issues were edited, no model fixes were implemented, and no publishing commands or production-cache regeneration were run. Subagents were not used because these closely overlapping issues benefit from one shared inspection.

The utility decision is treated as settled: [script 02](../../analysis/02_setup_and_global_variables.R) sets `UTILITY_SOURCE = 1`; [script 05](../../analysis/05_basecase_input_parameters.R) selects **0.73 progression-free / 0.59 progressed**. This is an audit of issue applicability and dependencies, not a fresh appraisal of the CORRECT utility estimates, Norwegian tariffs, or external clinical literature. Where an external input was unavailable, that limitation is distinguished from a verified code defect.

## Open issues: what still applies

### Input integrity and model execution

| Issue | Verdict | Independent evidence and remaining scope |
|---|---|---|
| [#174 — Raw export contract](https://github.com/ben-geisler/METIMMOX-1/issues/174) | **Keep open; partly resolved.** | [Script 01](../../analysis/01_data_prep.R) now checks positions, names, numeric validity, nonnegative values and binary flags for five endpoint fields, plus last-evaluation ordering. Thus “no checks” is outdated. However, lesion measurements still use `TL LD...121` / `TL LD...130`, CRP uses `CRP...478` / `CRP...540`, and [script 02](../../analysis/02_setup_and_global_variables.R) uses `Date...122` without an equivalent export-header contract. Script 01 has no patient-ID uniqueness or cohort-count assertion. The remaining issue concerns the trustworthiness of future imports; it does not establish that the current cohort is wrong. |
| [#60 — Survival and input validation](https://github.com/ben-geisler/METIMMOX-1/issues/60) | **Keep open; directly reproduced.** | [model_fun.R](../../R/model_fun.R) checks curve lengths and OS/PFS ordering, while [partitioned_survival_states()](../../R/calculate_outcomes.R) overwrites the initial state. A fresh synthetic probe with `control_OS[2] = 1.1` returned finite outcomes and dead occupancy **−0.1**. Curves starting at 0.8 were also accepted and their initial state overwritten. The finite/range/monotonicity/origin contract remains incomplete. The expanded comment is partly outdated: required unit costs now receive finite numeric scalar checks. Discount-rate checks remain only nonnegativity checks, and schedules lack a full validity contract. This is input validation, distinct from the intentional handling of OS/PFS crossings in otherwise valid PSA curves. |
| [#172 — Time settings and short horizons](https://github.com/ben-geisler/METIMMOX-1/issues/172) | **Keep open; partly resolved and directly reproduced.** | [model_fun()](../../R/model_fun.R) still defaults to 520 and 1/52; [run_basecase()](../../R/cea_helpers.R) and [script 08](../../analysis/08_basecase_analysis.R) omit explicit time arguments. A valid synthetic 261-point parameter set failed through `run_basecase()` with “length 261 but expected 521”; an explicit 260-week call worked. Changing only `params$cl` was ignored. Short-horizon CT/blood `seq()` calls remain in [script 05](../../analysis/05_basecase_input_parameters.R) and [build_horizon_params()](../../R/parameter_distributions.R). However, the old quarterly-charge failure and 13-point QALY overstatement are gone: `calculate_outcomes()` now integrates occupancy over intervals. |
| [#49 — Horizon-dependent schedules](https://github.com/ben-geisler/METIMMOX-1/issues/49) | **Applicable residual, but duplicate of #172.** | `build_horizon_params()` already truncates treatment positions to the horizon. An unchanged treatment schedule through week 38 is not inherently wrong for a five-year rather than ten-year analysis; treatment timing need not scale with the analysis horizon. The residual short-horizon and entry-point problems belong to #172. Script 05's direct assignments can also extend short schedule vectors beyond their initialized length. This is part of the same horizon contract, not an additional modeling workstream. The issue comment explicitly makes closure depend on #172. |

These issues concern different layers: #174 validates what enters the analysis, #60 validates what enters state construction, and #172/#49 validate the time grid and schedules. They overlap operationally but are not interchangeable.

### Economic inputs and clinical validation

| Issue | Verdict | Independent evidence and remaining scope |
|---|---|---|
| [#36 — Visit costs](https://github.com/ben-geisler/METIMMOX-1/issues/36) | **Keep open as an unresolved input question.** | [Script 05](../../analysis/05_basecase_input_parameters.R) still specifies EUR 33 for `c_other_visit` and `c_other_follow`, and EUR 530 for baseline costs; [the rendered input report](../input_parameters.md) reports those values without a reconciliation to Eline's DRG weights. The issue contains neither the weights nor an expected replacement value. I can verify the missing documented reconciliation, but cannot independently establish that EUR 33 is numerically wrong or what it should be. Its scope should be read together with #164 and #56. |
| [#164 — Monitoring and visits](https://github.com/ben-geisler/METIMMOX-1/issues/164) | **Keep open; the cadence half is obsolete.** | [calculate_outcomes.R](../../R/calculate_outcomes.R) now treats progressed-state follow-up as a quarterly rate integrated over occupancy, not a payment at positions 13, 26, 39. There is therefore no discrete 12-versus-13-week payment drift to reconcile. The remaining concern is real: [script 05](../../analysis/05_basecase_input_parameters.R) schedules CT and blood monitoring through the horizon, but `l_visit` ends with treatment administration at position 39. PF surveillance therefore continues without subsequent modeled visits. Whether that is the intended resource-use assumption is unresolved. The progressed-state omission of imaging is already disclosed in [input_parameters.qmd](../input_parameters.qmd); disclosure does not establish that the pathway is appropriate. |
| [#56 — Uncertainty CV justification](https://github.com/ben-geisler/METIMMOX-1/issues/56) | **Keep open; narrower than the original issue.** | [parameter_distribution_spec()](../../R/parameter_distributions.R) fixes drug/test prices in PSA and applies CV 0.20 to four resource-use costs, CV 0.15 to `u_np`, and CV 0.15 to `u_decrement`; `u_p` is derived. The current input report explains dependence and why some costs remain uncertain, but not why these particular CV magnitudes are defensible. Selecting CORRECT means does not supply the uncertainty of those means or of the decrement. This remains relevant despite closing IPD utility work. |
| [#161 — Table S8 utility citation](https://github.com/ben-geisler/METIMMOX-1/issues/161) | **Keep open.** | [define_utility_scenarios()](../../R/scenario_analysis.R) still defines 0.80/0.65 and literally stores `[CITATION PLACEHOLDER: reference for u_np 0.80 / u_p 0.65 to be added]`; [Table S8's generator](../../outputs/vignettes/table_s8.qmd) consumes the source field. These are a separate sensitivity scenario, not the CORRECT 0.73/0.59 base case. The utility-source decision does not resolve this publication gap. |
| [#162 — External survival validation](https://github.com/ben-geisler/METIMMOX-1/issues/162) | **Keep open; acknowledge subsequent partial evidence.** | [The survival specification report](survival_model_specification.md), [parametric report](../para_models.md) and [Figure 1 source](../../outputs/vignettes/figure1.qmd) contain within-trial fit/KM comparisons. [The extrapolation plausibility report](psa_extrapolation_plausibility.md) additionally compares draws with a Norwegian general-population reference and examines tail contributions. Thus “no external comparison” is too broad now. However, that report expressly says the general-population comparator is too permissive to establish cancer-specific validity. I found neither the requested time-specific model-versus-KM table with uncertainty nor an independent comparable-mCRC validation against trial/registry survival. #179's diagnostics do not close #162. |

The CORRECT choice removes the rationale for doing further patient-level EQ-5D estimation for this paper. It does **not** remove the need to justify uncertainty (#56), support the separate scenario (#161), or validate survival (#162).

### Communication, reproducibility and publication

| Issue | Verdict | Independent evidence and remaining scope |
|---|---|---|
| [#160 — Manuscript](https://github.com/ben-geisler/METIMMOX-1/issues/160) | **Keep open; a final integration task.** | [docs/manuscript](../../docs/manuscript) and [docs/references](../../docs/references) each contain only `.gitkeep`. The issue body's numerical example is stale. [Current Table 5](../../outputs/tables/table_5.csv) reports probabilistic means of EUR 20,931.85 / 1.377530 QALYs for SoC, EUR 53,228.19 / 1.400750 for CRP, and EUR 61,751.88 / 1.369621 for TMB/BRAF; the CRP ratio-of-mean-increments ICER is EUR 1,390,876/QALY and TMB/BRAF is dominated by SoC. The old stale-Table-5 concern is resolved in current contracts. Drafting can begin now, but numerical and interpretive finalization depends on the remaining substantive decisions. |
| [#177 — Poster CEAC comparison](https://github.com/ben-geisler/METIMMOX-1/issues/177) | **Keep open; location and explanation need updating.** | The current source is [poster_SMDM.qmd](../../outputs/vignettes/poster_SMDM.qmd), CEAC section. Both panels use the joint-model control draw. However, the single-biomarker panel draws CRP and TMB/BRAF models separately, with seeds 2027 and 2028, while the joint panel uses the same joint draw for control and both guided strategies. Thus the comparison changes dependence between strategies as well as survival specification. Joint OS/PFS sampling *within* each model does not fix that *between-model* comparison. The caption's attribution remains insufficient. This is a separate poster inference, not evidence that the main PSA is defective. |
| [#58 — Hard-coded labels](https://github.com/ben-geisler/METIMMOX-1/issues/58) | **Keep open; mostly reduced to residual presentation coupling.** | [report_format.R](../../R/report_format.R) derives canonical labels from model metadata. But [Figure 3](../../outputs/vignettes/figure3.qmd), [Figure 4](../../outputs/vignettes/figure4.qmd) and [Figure S5](../../outputs/vignettes/figure_s5.qmd) still key colors by literal display labels; [OWSA](../OWSA.qmd) has a separate short-label mapping. A central wording change can therefore cease to match a plot's scale. This is maintenance risk, not evidence that current economic strategies are mislabeled. It does not justify reopening #124. |
| [#178 — Publishing and warnings](https://github.com/ben-geisler/METIMMOX-1/issues/178) | **Keep open; partly resolved.** | [publish_reports.R](../../publish/publish_reports.R) checks Quarto render exit status, but still leaves worktree creation, checkout, staging/committing and copying unchecked, checks only local `gh-pages`, and copies onto an existing branch without removing obsolete PDFs. Its index stamps current date/HEAD even with `--no-render`. [setup_report()](../../R/report_setup.R) suppresses warnings during sourcing and sets `warning = FALSE`, without collecting a general warning summary. Existing publication-artifact contracts do not establish correct PDF deployment or render provenance. These findings are from static inspection; deployment was not exercised. |
| [#165 — Repository hygiene](https://github.com/ben-geisler/METIMMOX-1/issues/165) | **Keep open; refresh the inventory and separate priorities.** | There is no `CHANGELOG.md`, `tests/run_all.R`, `.Rprofile` or `renv/` activation directory. There are now **27** `test_*.R` scripts, not 18. No comprehensive parameter-to-consumer table was found in the stated documentation; helpers such as `format_results()` in [script 09](../../analysis/09_DSA.R) still have plain comments rather than roxygen blocks. Naming/documentation completeness remains a broad maintenance judgment, not a demonstrated result defect. The old expected-failure instruction for `test_psa_basecase_alignment.R` is obsolete: it now passes a structural check. Cache/version provenance and test isolation already have dedicated implementations and regression coverage; those should not be treated as wholly missing. |

One integration detail for #160: [CEA.qmd](../CEA.qmd)'s CRP-timing paragraph still says the randomization-time clock has “no economic consequence” and describes the subgroup curves as conditional on surviving to week 4. The more precise [week-4 estimand report](crp_week4_estimand.md) distinguishes the one pre-week-4 death from the two patients alive without a measurement, and quantifies early effects. This is a concrete reason to reconcile narrative sources during #160, rather than copying report prose without checking it. It does not warrant reopening the deliberately bounded missing-data/estimand analysis.

## Closed issues: was closure justified?

### #159 — PSA means and joint OS/PFS sampling

**Yes, keep closed.** [The discussion](https://github.com/ben-geisler/METIMMOX-1/issues/159) ultimately made closure contingent on joint sampling plus the reporting work split into #180, not exact equality between uncertain PSA means and fitted-parameter outcomes.

[joint_survival_sampling.R](../../R/joint_survival_sampling.R) estimates dependence from paired patient bootstrap refits, calibrates to the fitted marginal covariance blocks, and draws the combined coefficient vector. [Script 06](../../analysis/06_sampling.R) uses that implementation. [CEA](../CEA.qmd) explicitly distinguishes nonlinear PSA expectations from deterministic point evaluations and leads with probabilistic results. The focused covariance, structural-control and zero-uncertainty tests all passed during this audit; current report/cache contracts passed as well.

The issue body's attribution of the full mean shift to independent endpoint draws is too strong. With fixed marginal coefficient distributions, changing their dependence cannot change the expectation of an additive raw-curve outcome; ordering corrections and other nonlinear joint operations can be affected. The current explanation recognizes that distinction. Remaining crossings are not evidence that the agreed joint-covariance implementation was left unfinished. Cancer-specific extrapolation credibility remains #162; reopening #159 would conflate a resolved sampling/reporting issue with a different validation question.

### #124 — Figure specifications

**Yes, as implemented in part and superseded in part.** [The issue](https://github.com/ben-geisler/METIMMOX-1/issues/124) requests Models A/B/C, which no longer belong to the economic design. Current [Figure 1](../../outputs/vignettes/figure1.qmd) has joint-model CRP/SoC OS/PFS panels, sampling-draw quantile ribbons and KM overlays; [Figure 2](../../outputs/vignettes/figure2.qmd) is an incremental-NMB tornado; [Figure 3](../../outputs/vignettes/figure3.qmd) is the PSA plane with WTP line; [Figure 4](../../outputs/vignettes/figure4.qmd) has base/biosimilar CEACs; [Figure 5](../../outputs/vignettes/figure5.qmd) uses group/scenario EVPPI. The obsolete supplementary PFS source is absent from the active vignette directory.

The ribbons are from the current coefficient-draw procedure, not literal resampled-patient bootstrap percentiles as the old specification requested. That is consistent with the subsequent methodological change. Final figure-to-text correspondence belongs to #160; residual label coupling belongs to #58.

### #66 — Standard headers and footers

**Yes.** [Closure](https://github.com/ben-geisler/METIMMOX-1/issues/66) is supported by a fresh scan of all **21 report `.qmd` files**: each has the standard subtitle and completion/repository/version signature fields; none contains `fancyhdr` or `include-in-header`. This verifies source conventions, not a fresh visual inspection of every PDF. There is no evidence supporting reopening.

### #16 — QMD warnings

**Yes for the stated presentation request, with an explicit distinction.** [The original request](https://github.com/ben-geisler/METIMMOX-1/issues/16) was to remove QMD warning/error messages. Current [report setup](../../R/report_setup.R) suppresses warnings/messages, and scans of rendered report Markdown and available vignette HTML found no matching R warning/error blocks. Reports are therefore clean in the relevant presentation sense.

This is not evidence that the underlying computations emit no warnings. The live zero-uncertainty test passed but emitted two Hessian/covariance warnings during candidate fitting, in addition to local R startup/package warnings. The general visibility of suppressed warnings remains legitimately open in **#178**. Keep that concern there rather than duplicating it by reopening #16. If #16 had been intended to require elimination of every underlying warning cause, closure would not be justified on suppression alone; its terse wording does not establish that stronger requirement.

### #6 — Working-directory warning

**Yes.** [The requested change](https://github.com/ben-geisler/METIMMOX-1/issues/6) is reflected in [setup_report()](../../R/report_setup.R), which uses `knitr::opts_knit$set(root.dir = here::here())`. The report-source scan found no `setwd()` calls, and the rendered-report scan found no corresponding working-directory warning. The closure comment's broader phrase “all file paths use here” need not be literally universal for this specific issue to be resolved.

### #29, #31, #33 and #34 — IPD utility work

| Closed issue | Opinion | Reason |
|---|---|---|
| [#29 — Parametric utility functions](https://github.com/ben-geisler/METIMMOX-1/issues/29) | **Correct to close as not planned.** | Further IPD utility-function exploration is not needed for the selected external health-state utilities. The sparse original issue provides no separate requirement that survives that scope decision. |
| [#31 — Norwegian EQ-5D tariffs](https://github.com/ben-geisler/METIMMOX-1/issues/31) | **Correct to close as not planned.** | Rescoring the trial's patient-level EQ-5D profiles under a Norwegian value set is a different project from using published CORRECT health-state values. Their applicability to this evaluation can be discussed without undertaking IPD rescoring. |
| [#33 — Utility table by time point](https://github.com/ben-geisler/METIMMOX-1/issues/33) | **Correct to close as not planned.** | A longitudinal IPD descriptive table is no longer an input-development prerequisite. It would become relevant only if a separate clinical HRQoL objective were introduced. |
| [#34 — Time-point-adjusted utilities](https://github.com/ben-geisler/METIMMOX-1/issues/34) | **Correct to close as not planned.** | The requested age/sex/time-point regression explicitly depended on #33 and belongs to the abandoned IPD estimation path. |

[The archived QALY notebook](../../archive/05_QALYs.Rmd) and [Danish/UK utility functions](../../R/eq5d5l_utility.R) remain available, but the economic base case selects the external constants. Retained historical code and an optional IPD switch do not mean these tasks remain required. These four closures are scope decisions, not claims that the requested analyses were completed. **#56 and #161 remain open for different reasons.**

## Overlap and order

### Dependencies that matter

| Connected issues | Nature of overlap | Ordering implication |
|---|---|---|
| #49 / #172 | Same remaining horizon and schedule contract. | One work item; assess both together when completed. |
| #174 / #60 / #172 | Raw data, probability inputs and time-grid validity are successive analysis boundaries. | Prioritize these before relying on a final regenerated result set. They can be investigated independently; not every code change has a strict dependency on the preceding one. |
| #36 / #164 / #56 | Cost per contact, which contacts are represented, and uncertainty in resource-use cost refer to the same quantities. | Settle the scope of contacts (#164) together with their valuation (#36), then finalize the associated uncertainty rationale (#56). Avoid treating the unit-cost number in isolation. |
| #56 / #161 | Both concern utilities, but one concerns uncertainty around the chosen inputs and the other a distinct alternative scenario. | Coordinate source review; neither requires reviving IPD estimation or necessarily waiting for the other. |
| #162 / #159 / #160 | Extrapolation validation informs the interpretation of the now-correctly reported PSA. | Keep #159 closed; carry #162 evidence and limitations into #160. General-population diagnostics can be reused, but are not a substitute for #162. |
| #177 / #159 / #124 | Shared sampling/figure terminology obscures distinct scopes. | #177 is an independent comparison between model specifications; it does not block accepting #159 or #124 as closed. |
| #58 / #124 / #160 | Central labels, existing figures and manuscript wording meet at publication. | Finish relevant label work before the final figure/text consistency pass. Do not reopen the superseded figure specification. |
| #178 / #16 / #6 / #165 | Warning visibility, clean rendering, directory handling and reproducibility overlap operationally. | #178 owns remaining warning/deployment concerns; #16/#6 stay closed. Useful test coordination from #165 can happen early, while cosmetic documentation can follow substantive work. |

### Recommended sequence

1. **Establish trusted inputs and execution: #174, #60, #172 with #49.** These have the broadest downstream reach. #60 and #172 have reproducible failures, while #174 guards the meaning of future data refreshes. Group this work before a final regeneration, without interpreting it as evidence that every current published number is wrong.

2. **Resolve substantive assumptions and validation in two concurrent tracks.** One track is **#164 + #36, followed by the related part of #56**, with **#161** alongside the utility-source discussion. The other is **#162**. External validation should start early because identifying a genuinely comparable population may take longer than repository work. Its final assessment should use the settled data definition and model behavior from step 1. These tracks need not wait for each other.

3. **Complete the relevant presentation work: #58 and, separately, #177.** #58 should precede final output regeneration and manuscript consistency checking. #177 can proceed independently, after any shared execution-contract changes it needs; it is a blocker only for using that poster's model-comparison claim. It need not hold up a manuscript that does not rely on that claim.

4. **Finalize #160 from settled evidence and current outputs.** Begin outline, methods and bibliography work immediately. Final numerical claims should follow steps 1–2 and any resulting regeneration; interpretation should incorporate #162 and accurately state the selected utilities, probabilistic estimand and deliberate scope choices. The manuscript should not wait for unrelated poster work or all low-priority hygiene. This ordering makes #160 a continuing integration task, not a reason to postpone all writing.

5. **Complete #178 before the next external publication.** Work can start earlier; warning visibility is useful during the final rendering pass. Correct provenance and deployment behavior are a release dependency, not a prerequisite for drafting text or validating model assumptions. A successful report/cache contract test does not remove this dependency.

6. **Handle #165 selectively throughout, then close the remaining documentation work.** Test coordination and environment/reproducibility clarity are useful early. The final changelog and input-to-output map should reflect the settled model. Roxygen/naming completeness should not become a serial blocker ahead of the substantive issues or manuscript.

This is a dependency-based order, not a restatement of GitHub's P1–P4 labels. In particular, #36 and #160 are both labeled P1, but resolving a price without its contact definition invites repeated work, and finalizing a manuscript before its inputs are settled invites repeated rewriting. Conversely, a P4 label does not make the useful test-coordination part of #165 inappropriate to do early.

## Executable verification record

The following checks were run afresh with the repository's R 4.3.2 executable. No results below are inferred solely from a historical “pass” record.

| Check | Result |
|---|---|
| `tests/test_joint_survival_sampling.R` | Exit 0. Paired refits, calibrated covariance, joint draws, RNG reproducibility and invalid-covariance guards passed. |
| `tests/test_psa_zero_uncertainty.R` | Exit 0. Trial-data fit selected gamma/gamma; every tested deterministic/zero-uncertainty PSA cost and QALY result agreed within the test's tolerances. Candidate-fitting Hessian warnings were emitted; this was not a warning-free run. |
| `tests/test_psa_basecase_alignment.R` | Exit 0. Current structural joint-model control prediction check passed. This no longer tests uncertain outcome means against deterministic outcomes. |
| `tests/test_report_contracts.R` | Exit 0. **593 checks passed; one empty-base-EVPPI row block skipped because EVPI is zero.** |
| Synthetic #60 probe | OS 1.1 accepted; dead occupancy −0.1; finite results. Initial OS/PFS 0.8 also accepted and overwritten in state construction. |
| Synthetic #172 entry-point probe | A 261-point parameter set failed in `run_basecase()` but produced 5 QALYs per strategy with explicit 260-week horizon, under unit utility/no mortality/no discount. |
| Synthetic #172 cycle-length probe | With 521 points, `params$cl = 1/12` alone still produced 10 QALYs; explicit `cl = 1/12` produced 43.33333. This demonstrates inconsistent argument handling, not a validated monthly model. |
| Short-grid probe | `seq(13, 12, by = 12)` and `seq(5, 4, by = 4)`, the current schedule expressions, failed. Direct model evaluation with supplied schedules at 1, 2, 5, 12 and 13 points succeeded, producing respectively 0, 1/52, 4/52, 11/52 and 12/52 QALYs. Thus the residual failure is in schedule construction, not the repaired interval integration. |
| Report-source scan | 21 `.qmd` reports had the expected subtitle/signature fields; no custom-header directives or `setwd()` calls were found. |

The synthetic probes used a temporary script and zero-cost, unit-utility curves; they do not estimate effects on the trial's cost-effectiveness results. No full PSA rerun, full report render, visual PDF review, external-data validation, or deployment exercise was performed. Those limits do not prevent the issue-scope judgments above, but they do preclude treating this audit as a complete model revalidation.

---

**Report completed on:** 2026-09-29  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.0
