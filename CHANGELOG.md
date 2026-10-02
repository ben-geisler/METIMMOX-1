# Changelog

One line per GitHub issue, from the validateHE run of July 2026 (#145) onwards, newest first. Each entry gives what changed, which caches had to be regenerated, and the headline result where it moved. Full details are in the issue, in the commit, in the bug-fix impact report (`reports/technical/bug_fix_impact.qmd`, from the snapshot pairs in `data/output/snapshots/`) and in AGENTS.md. Monetary values are EUR per patient; "PSA" means the probabilistic means of the 5,000-draw cache; WTP is EUR 51,000/QALY throughout.

Caches: **S** sampling (script 06), **P** PSA pair (10), **E** EVPPI (11), **Sc** scenario EVPPI (12), **O** PFS/OS ordering diagnostics, **R** reports and publication vignettes re-rendered.

## October 2026

- **#188** `UTILITY_SOURCE = 0` no longer uses the irreproducible constants 0.9077/0.9005: script 05 derives the IPD utilities with `trial_eq5d_utility_uncertainty()` (new `EQ5D_VALUE_SET` setting, Danish default: 0.912/0.897) and the PSA uses patient-clustered bootstrap SEs instead of CV 0.15. Both settings can be set with `options()`; IPD caches are labelled `ipd_<value set>`. Base case (`UTILITY_SOURCE = 1`) unchanged: P, E, Sc regenerated (code fingerprint only), snapshot pair `188` identical.
- **#161** Table S8 utility scenario: the uncited literature pair (0.80/0.65) replaced by utilities derived from the METIMMOX EQ-5D-5L responses with the Danish value set, `trial_eq5d_utilities()` (u_np 0.912, u_p 0.897; a response is progressed only after a recorded progression). The archived derivation and the hard-coded `UTILITY_SOURCE = 0` pair (0.9077/0.9005) could not be reproduced. Manifest records the EQ-5D file digest. Sc regenerated (code fingerprint only); Table S8 re-rendered. Standard of care remains optimal in the scenario.
- **#187** `analysis/07_traces.R` no longer overwrites the global `strategy_labels` of `R/report_format.R` (local renamed `v_trace_labels`); `test_model_config_coupling.R` now stops if any analysis script assigns a name that `R/` defines. No rendered output sourced script 07, so none changed; no cache change.
- **#186** `tests/test_sampling_failure_fallback.R` rewritten for the #159 joint sampler (real fits on simulated data: joint draws, covariance-unavailable stop, more than 10% failed bootstrap pairs stop) and removed from `KNOWN_FAILURES`. Test only; no cache change.
- **#165** Repository hygiene. Adds this changelog and `tests/run_all.R` (every test in its own process on a temporary cache copy; fails if a live cache changes). `renv.lock` now records the 9 packages it lacked (voi among them) and the README documents `renv::restore()`. Roxygen added for the helpers in `analysis/` and `tests/`. `input_parameters.qmd` (v4.7) gains a computed table that maps each `l_params_base` field to the step that consumes it. DARTH prefixes (`v_`, `m_`, `df_`) extended to function-local vectors, matrices and data frames in `R/`, `analysis/` and `tests/`. Caches S, P, E, Sc, O regenerated (code fingerprints changed); results unchanged (snapshot pair `165`).
- **#162** External validation of the extrapolated survival (`reports/technical/survival_external_validation.qmd`). Benchmarks and tolerance pre-specified (commit `90b093a`). All 17 published values lie within the parametric model's 95% interval; 11 of 17 agree on point estimates; tail point estimates are on the low side. No tail fix (user decision). No cache change.
- **#58, #178** Strategy labels centralised in `R/report_format.R`; `publish/publish_reports.R` hardened (step statuses, stale PDFs; #178 remains open). R.

## September 2026

- **#36, #164, #56** Visit costs from Norwegian DRG weights (EUR 213.46 per visit), a surveillance visit with every CT, a quarterly CT in the progressed state, PSA CVs documented as assumptions. P, E, Sc, R. PSA costs rise to 25,732 (control), 57,740 (CRP), 66,302 (TMB/BRAF); QALYs unchanged; CRP ICER EUR 1,378,448/QALY; TMB/BRAF dominated; EVPI 0.
- **#185** Fingerprints and manifest sorted in C order, independent of the collation locale. Fingerprints only; no result change.
- **#172, #49, #60** `model_fun()` takes horizon and cycle length from the parameter list; schedules built for any horizon (`build_treatment_schedules()`); survival curves validated (finite, in [0, 1], start at 1, nonincreasing). S, P, E, Sc regenerated; every result identical.
- **#174** Raw trial-export contract (`R/data_contract.R`) checked by script 01 before any derivation. Tidy data unchanged apart from `CT1date`.
- **#183** Ridge sensitivity penalises the biomarker main effects and interactions on an explicit grid; boundary solutions flagged. Clinical report v3.11.
- **#182, #176** Clinical report corrections: Rx row of the TLR tables, week-12 landmark statement, Firth wording, convergence table, exit reasons of the 3 patients without TLR.
- **#184, #170** TLR landmark base cohort shared by the clinical report, Figure 1 and the DAG tests (65 patients); 4 of 19 DAG independencies classified as definitional; marginal edge tests labelled as such.
- **#181** PFS counts a death without progression only within 16 weeks of the last assessment; PFS events 63 to 50 in the 68-patient cohort; gamma/gamma still selected. S, P, E, Sc, O, R. PSA: control 20,932 / 1.3775 QALYs, CRP 53,228 / 1.4008, TMB/BRAF 61,752 / 1.3696; CRP ICER EUR 1,390,876/QALY.
- **#175** Estimand of the week-4 CRP decision stated; missingness sensitivity analysis for the 3 patients without a week-4 CRP (`reports/technical/crp_week4_estimand.qmd`). No cache change.
- **#180** Probabilistic means, percentile intervals and ratio-of-mean-increments ICERs are the primary results (Table 5); the five-MCSE PSA-versus-base-case criterion retired as comparing different estimands. R.
- **#179** Technical report on the plausibility of the PSA survival extrapolations; general-population mortality floor recorded as out of scope.
- **#168** Publication artifacts written or removed under a manifest (`outputs/manifest.csv`); Table 5 rebuilt from the current PSA; zero-EVPI base case kept in Table S7 and Figure 5. P, E regenerated; results identical.
- **#173, #171, #163** Diagnostic-test prices read from one canonical scalar; non-finite PSA outcomes count as failures; trapezoidal integration of QALYs and progressed-state costs (the 521-point grid spans exactly 520 weeks). P, E, Sc.
- **#159** Joint OS/PFS multivariate-normal coefficient draws with a paired-bootstrap cross-covariance (996/1,000 pairs retained). S, P, E, Sc, O. EVPI EUR 0; crossing draws fall from 48% to 24%.
- **#166** One complete-case target population (68 patients) for control and all guided strategies; joint Dirichlet prevalence draws. S, P, E, Sc. Identical treatment now gives zero increments.
- **#167, #169** Cache fingerprints cover calculation code, package versions and the data; PSA outcome and parameter files bound as a pair; tests use temporary caches. S, P, E, Sc.
- **#157** Vignette and documentation housekeeping: subset-matched KM overlays, Table S1/S8 built from the pipeline, enriched-population screening cost per identified positive (enriched CRP ICER 1,587,067 to 1,578,523), restricted mean by the trapezoidal rule. R.
- **#156** Multivariate-normal survival draws replace the bootstrap; fingerprinted caches; `model_idx` recorded per PSA row; structural DSA scenarios disclosed; snapshots record PSA provenance. S, P, E, Sc. EVPI 22.13 to 0.
- **#155** DAG tests derived from the canonical graph; non-tautological `PFS -> OS` and `TLR -> PFS` tests; week-9 landmark diagnostics and sensitivity.
- **#154** Unit prices fixed in the PSA (DSA only); utilities drawn jointly (`u_p = u_np - u_decrement`); post-progression cost and second-sequence structural scenarios. P, E, Sc. EVPI 45.20 to 22.13.
- **#153** Reporting defects: arm-by-biomarker stratum labels, clinical cohort, selected-distribution reporting, pairwise versus frontier dominance labels. R.
- **#152** Regression EVPPI (`voi::evppi()`, GAM) with Monte Carlo SEs replaces the kNN estimator; group rows asserted at least their largest member. E, Sc.
- **#151** PSA control arm predicted from the joint model with Rx = control, as in the base case. P, E.
- **#150** CRP described as the week-4 value measured before the first nivolumab dose, not baseline.
- **#149** PFS recoded as progression or death. S, P, E, Sc.

## July-August 2026

- **#147** Three correctness fixes from the adversarial review (one changed EVPPI results); discounting on the model's cycle length (#82).
- **#146** A/B/C model-structure and TLR economic scaffolding retired; duplicated code consolidated. No result change.
- **#145** Critical findings of the validateHE run: OS/PFS crossings masked by a clamp, ordering-constrained distribution selection (gamma/gamma), documentation brought in line with the code.
