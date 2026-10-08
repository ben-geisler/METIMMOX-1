# METIMMOX-1

Secondary analyses of the randomised METIMMOX trial ([NCT03388190](https://clinicaltrials.gov/study/NCT03388190)) in metastatic microsatellite-stable (MSS) / mismatch repair-proficient (pMMR) colorectal cancer.

METIMMOX (*Colorectal Cancer METastasis: Shaping Anti-tumor IMMunity by OXaliplatin*) compared first-line alternating short-course oxaliplatin-based chemotherapy (Nordic FLOX) plus nivolumab with FLOX alone. This repository holds the code of three papers that share one data pipeline:

1. **Clinical effectiveness:** *Biomarker signals in MSS/pMMR mCRC are predominantly prognostic, not predictive* (working title). Do CRP, TMB/BRAF and target lesion reduction identify patients who benefit from nivolumab, or only patients with a better prognosis?
2. **Cost-effectiveness:** is giving nivolumab only to biomarker-selected patients cost-effective in Norway?
3. **Value of information:** what would further research to resolve the remaining uncertainty be worth?

Papers 2 and 3 use the same economic model. The code is written for clinical researchers, health economists and anyone building decision-analytic models in R.

## What this repository contains

The code that computes every number, table and figure of the three papers: model and helper functions (`R/`), the numbered pipeline (`analysis/`), the regression tests (`tests/`) and the Quarto sources of every report and publication figure or table (`reports/`, `reports/technical/`, `reports/vignettes/`).

It is a mirror of a private working repository, updated from it after every change. Three things are not here:

- **The trial data.** They are confidential (see [Data](#data)), so the pipeline cannot run from a fresh clone.
- **Rendered outputs** (figures, tables, rendered reports). The scripts regenerate them from the data into `outputs/`; the published versions are in the papers and their supplements.
- **Manuscript drafts.**

## Biomarkers

| Biomarker | Definition | Used in |
|---|---|---|
| CRP | C-reactive protein < 5 mg/L at week 4 (cycle 3 day 1), measured after the two FLOX cycles that both arms receive and before the first nivolumab dose | All papers |
| TMB/BRAF | Tumour mutational burden >= 9 mut/Mb or a BRAF mutation (baseline next-generation sequencing) | All papers |
| TLR | Target lesion reduction: at least 10% shrinkage in the sum of target-lesion diameters at the first on-treatment CT | Paper 1 only |

Week-4 CRP is not a baseline measurement, but it is known before the decision to add nivolumab; see the [CRP estimand report](reports/technical/crp_week4_estimand.qmd). TLR is a post-randomisation mediator measured on treatment, so it cannot select patients for treatment and is excluded from every economic analysis. Its clinical analyses use landmark cohorts (week 9 primary; scan date and week 12 as sensitivity analyses) to limit guarantee-time bias.

## Paper 1: Clinical effectiveness

**Methods.** Cox models of overall and progression-free survival with biomarker-by-treatment interactions (`Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx`) on the 68-patient complete-case cohort. Firth's penalised likelihood handles the small sample and few events; cross-validated ridge regression (`glmnet`) is a shrinkage sensitivity analysis. TLR is analysed on week-9 landmark cohorts. A directed acyclic graph (DAG) sets out the assumed causal structure, and every edge and implied conditional independence that can be tested against the trial data is tested. Sensitivity analyses cover the definition of the progression-free survival endpoint.

| Item | Generator | Output |
|---|---|---|
| Figure 1: Kaplan-Meier curves by biomarker and arm | [clin_effect_figure1.qmd](reports/vignettes/clin_effect_figure1.qmd) | `outputs/figs/clin_effect_figure1.png`, `.eps` |
| Tables 1-4: baseline characteristics; Firth Cox models; ridge versus Firth; TLR landmark analysis | [clinical_effectiveness.qmd](reports/clinical_effectiveness.qmd) | rendered report |
| Figure S2: directed acyclic graph | [clin_effect_figure_s2.qmd](reports/vignettes/clin_effect_figure_s2.qmd) | `outputs/figs/clin_effect_figure_s2.png` |
| Table S1: DAG node definitions | [dag.qmd](reports/dag.qmd) | rendered report |
| Table S2: bivariate association tests | [clin_effect_table_s2.qmd](reports/vignettes/clin_effect_table_s2.qmd) | `outputs/tables/clin_effect_table_s2.csv`, `_raw.csv` |
| Table S3: conditional-independence assessments | [clin_effect_table_s3.qmd](reports/vignettes/clin_effect_table_s3.qmd) | `outputs/tables/clin_effect_table_s3.csv`, `_raw.csv` |
| Tables S4-S11: Schoenfeld tests; landmark cohorts; standard versus Firth Cox; PFS endpoint sensitivity | [clinical_effectiveness.qmd](reports/clinical_effectiveness.qmd) | rendered report |

The DAG tests are also discussed in [dag_associations.qmd](reports/dag_associations.qmd); [biomarker_distributions.qmd](reports/biomarker_distributions.qmd) describes prevalence and overlap by arm. `tests/check_output_consistency.R` cross-checks the rendered clinical report against Tables S2 and S3. Pipeline: scripts 01-03.

## Paper 2: Cost-effectiveness

Three strategies are compared: standard of care (FLOX alone for everyone), CRP-guided and TMB/BRAF-guided treatment (FLOX + nivolumab for biomarker-positive patients, FLOX alone for the rest).

- **Structure.** Partitioned survival model with three states (progression-free, progressed, dead), weekly cycles, a 10-year horizon, Norwegian healthcare perspective, 4% discounting of costs and QALYs, willingness to pay EUR 51,000 per QALY.
- **Survival.** One joint parametric model per endpoint on the 68-patient complete-case cohort:
  ```r
  Surv(OSwk,  Death)       ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
  Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
  ```
  The distributions are chosen jointly from nine candidate families as the minimum-combined-AIC pair that keeps OS >= PFS in every modelled subgroup and week. Standard of care is the same model predicted with `Rx = control`. Curves are averaged over every patient's own covariates (population averaging).
- **Uncertainty.** The PSA (5,000 draws) samples survival coefficients from a joint multivariate normal whose OS-PFS dependence comes from a paired patient bootstrap, the joint biomarker distribution from a Dirichlet, and utilities and resource-use costs. Drug and test unit prices are fixed in the PSA and varied in the one-way analysis. Results are probabilistic means with percentile intervals.
- **Sensitivity and scenarios.** One-way DSA, structural scenarios (discount rate, horizon, post-progression cost, no second treatment sequence, alternative survival families), a biosimilar nivolumab price and biomarker-enriched populations.
- **Scope.** Deliberately excluded: separate adverse-event costs or disutilities, post-progression therapy cost in the base case, and a general-population mortality floor on extrapolated hazards. The nivolumab price of EUR 13,923 per administration is an assumption, not a tendered price.

| Item (numbering of the current draft) | Generator | Output |
|---|---|---|
| Table 1: model inputs | [table_1.qmd](reports/vignettes/table_1.qmd) | `outputs/tables/table_1.csv` |
| Table 2: base-case results (probabilistic; deterministic) | [table_5.qmd](reports/vignettes/table_5.qmd); [table_4.qmd](reports/vignettes/table_4.qmd) | `outputs/tables/table_5.csv`; `table_4.csv`, `table_s5.csv` |
| Table 3: biomarker-positive populations | [table_6.qmd](reports/vignettes/table_6.qmd) | `outputs/tables/table_6.csv` |
| Figure 1: survival of the CRP-guided strategy and standard care | [figure1.qmd](reports/vignettes/figure1.qmd) | `outputs/figs/figure1.png` |
| Figure 2 and S3: cost-effectiveness plane | [figure2.qmd](reports/vignettes/figure2.qmd) | `outputs/figs/figure2.png`, `figure_s3.png` |
| Figure 3 and S4: acceptability curves | [figure3.qmd](reports/vignettes/figure3.qmd) | `outputs/figs/figure3.png`, `figure_s4.png` |
| Figure 4: one-way sensitivity analysis | [figure4.qmd](reports/vignettes/figure4.qmd) | `outputs/figs/figure4.png` |
| Table S1: patient characteristics | [table_s1.qmd](reports/vignettes/table_s1.qmd) | `outputs/tables/table_s1.csv` |
| Table S2: biomarker correlations | [table_s2.qmd](reports/vignettes/table_s2.qmd) | `outputs/tables/table_s2.csv` |
| Table S3: survival model selection | [table_s3.qmd](reports/vignettes/table_s3.qmd) | `outputs/tables/table_s3.csv`, `table_s3_pairs.csv` |
| Table S4: survival model coefficients | [table_s4.qmd](reports/vignettes/table_s4.qmd) | `outputs/tables/table_s4.csv` |
| Table S5: external validation of the extrapolation | [survival_external_validation.qmd](reports/technical/survival_external_validation.qmd) | rendered report |
| Tables S6, S7: full deterministic and enriched-population results | [table_4.qmd](reports/vignettes/table_4.qmd), [table_s6.qmd](reports/vignettes/table_s6.qmd) | `outputs/tables/table_s5.csv`, `table_s6.csv` |
| Table S8: scenario analyses | [table_s8.qmd](reports/vignettes/table_s8.qmd) | `outputs/tables/table_s8.csv` |
| Figures S1, S2, S5: Kaplan-Meier versus fits; extrapolations; frontier | [figure_s1.qmd](reports/vignettes/figure_s1.qmd), [figure_s2.qmd](reports/vignettes/figure_s2.qmd), [figure_s5.qmd](reports/vignettes/figure_s5.qmd) | `outputs/figs/figure_s1.png`, `figure_s2.png`, `figure_s5.png` |
| Numbers in the text only (threshold price, PSA probabilities by price, landmark survival, ...) | [cea_numbers.qmd](reports/vignettes/cea_numbers.qmd) | `outputs/tables/cea_*.csv` |

Method reports: [para_models](reports/para_models.qmd), [input_parameters](reports/input_parameters.qmd), [CEA](reports/CEA.qmd), [OWSA](reports/OWSA.qmd), [biosimilar_scenario](reports/biosimilar_scenario.qmd), [enriched_population](reports/enriched_population.qmd), [biomarker_decomposition](reports/biomarker_decomposition.qmd), and the technical reports in [reports/technical/](reports/technical/) (survival model specification, joint survival sampling, PFS/OS ordering in the PSA, extrapolation plausibility, CRP estimand, external validation). Pipeline: scripts 01-10 with 08b; Figure 3 also needs script 12.

## Paper 3: Value of information

EVPI and EVPPI from the Paper 2 PSA. EVPPI uses nonparametric regression (Strong, Oakley and Brennan 2014, via `voi::evppi()`) on single parameters and parameter groups (biomarker-by-treatment interaction coefficients, joint biomarker prevalence, utilities, resource-use costs), each with a Monte Carlo standard error. Per-patient values are scaled to the Norwegian population (1,500 eligible patients per year, 10-year research horizon). Scenarios repeat the PSA and EVPPI at willingness-to-pay thresholds of EUR 51,000, 100,000 and 150,000 per QALY, at list and biosimilar nivolumab prices.

| Item | Generator | Output |
|---|---|---|
| Population EVPPI by parameter group and scenario | [figure5.qmd](reports/vignettes/figure5.qmd) | `outputs/figs/figure5.png` |
| EVPPI by group and scenario | [table_s7.qmd](reports/vignettes/table_s7.qmd) | `outputs/tables/table_s7.csv` |

Reports: [EVPPIs](reports/EVPPIs.qmd), [scenario_effect](reports/scenario_effect.qmd). Pipeline: scripts 01-12.

## Data

The trial data are confidential and not part of this repository. Pseudonymised individual patient data may be made available to qualified researchers on reasonable request, subject to approval and a data sharing agreement; open an issue or contact [ben-geisler](https://github.com/ben-geisler).

With access, the pipeline expects:

| File | Content |
|---|---|
| `data/sensitive/METIMMOX w TMB 241101.xlsx` | trial export (checked against the contract in `R/data_contract.R` before anything is derived) |
| `data/sensitive/EQ5Dmanual.xlsx` | EQ-5D-5L responses (trial-utility scenario and IPD utility mode) |
| `data/sensitive/excluded_ids.csv` | trial IDs excluded from every analysis, with the reason (columns `ID`, `reason`) |
| `data/external/survival_benchmarks.csv` | published survival estimates for the external validation, transcribed from Tveit et al. 2012 (doi:10.1200/JCO.2011.38.0915), Guren et al. 2017 (doi:10.1038/bjc.2017.93), Aasebo et al. 2019 (doi:10.1002/cam4.2205), Acta Oncologica 2025 (doi:10.2340/1651-226X.2025.42985) and Sorbye et al. 2013 (doi:10.1093/annonc/mdt197) |
| `data/external/ssb_07902_2025.csv` | Norwegian 2025 period life table, probabilities of death by age and sex, from [Statistics Norway table 07902](https://www.ssb.no/en/statbank/table/07902) (per 1,000 divided by 1,000) |

The two `data/external/` files contain only published values; they live with the data so that `data/` is one folder and are available on request with it.

## Running the analysis

Requires R 4.3, Quarto with a LaTeX installation for PDF output, and the data above.

```bash
Rscript make.R --dry-run      # list the steps
Rscript make.R                # analysis 01-12, publication vignettes, reports, tests
Rscript make.R --from 09      # rerun the analysis from script 09 (02-06 re-sourced as setup)
Rscript make.R --no-analysis --no-vignettes --no-tests --reports CEA,OWSA
```

From empty caches the analysis takes about 1.5 hours (script 06 several minutes, 10 about 20 minutes, 12 about 40 minutes); the vignettes and reports take about another hour. Scripts 06, 10, 11 and 12 write caches to `data/tidy/`, each with a fingerprint of its inputs, code and package versions; reports stop rather than render from a stale cache. Everything rendered goes to `outputs/`: figures and tables to `outputs/figs/` and `outputs/tables/` (with their provenance in `outputs/manifest.csv`), reports to `outputs/reports/` (PDF and Markdown).

**Package versions.** `renv.lock` lists the exact package versions the published results were produced with. The project does not switch R to a project-specific package library (renv is not "activated"), so the scripts use your normal library. To match the recorded versions, install them once:

```r
renv::restore(lockfile = "renv.lock", prompt = FALSE)
renv::status(lockfile = "renv.lock")   # should report no inconsistencies
```

The caches record the versions of the packages that compute the results, so a cache built under other versions is rebuilt rather than reused.

## Tests

`tests/run_all.R` runs every `tests/test_*.R` in its own process against a temporary copy of the caches, skips the tests that need the trial export when it is absent, and fails if any live cache changed. Most tests use synthetic data, so they run in a fresh clone:

```bash
Rscript tests/run_all.R --no-data     # without the trial data
Rscript tests/run_all.R               # all tests
Rscript tests/test_pfs_endpoint.R     # a single test
```

`tests/test_public_safety.R` scans the files of this repository for anything that must not be public (data files, patient identifiers, local paths, personal e-mail addresses).

## Related publications

1. Ree AH, Šaltytė Benth J, Hamre HM, et al. First-line oxaliplatin-based chemotherapy and nivolumab for metastatic microsatellite-stable colorectal cancer: the randomised METIMMOX trial. *Br J Cancer*. 2024;130(12):1921-1928. [doi:10.1038/s41416-024-02696-6](https://doi.org/10.1038/s41416-024-02696-6)
2. Meltzer S, Negård A, Bakke KM, et al. Early radiologic signal of responsiveness to immune checkpoint blockade in microsatellite-stable/mismatch repair-proficient metastatic colorectal cancer. *Br J Cancer*. 2022;127(12):2227-2233. [doi:10.1038/s41416-022-02004-0](https://doi.org/10.1038/s41416-022-02004-0)
3. Meltzer S, Berg JP, Hamre HM, et al. 632P Predictive value of C-reactive protein (CRP) in microsatellite-stable (MSS) metastatic colorectal cancer (mCRC) patients given first-line alternating short-course oxaliplatin-based chemotherapy (FLOX) and nivolumab. *Ann Oncol*. 2023;34:S449. [doi:10.1016/j.annonc.2023.09.1822](https://doi.org/10.1016/j.annonc.2023.09.1822)
4. Ree AH, Bousquet PA, Nilsen HL, et al. 543P Tumor mutational burden (TMB), BRAF status, and C-reactive protein (CRP) predict response to first-line alternating oxaliplatin-based chemotherapy and nivolumab in metastatic microsatellite-stable (MSS) colorectal cancer (CRC). *Ann Oncol*. 2024;35:S453. [doi:10.1016/j.annonc.2024.08.612](https://doi.org/10.1016/j.annonc.2024.08.612)

## Citation

See [CITATION.cff](CITATION.cff) (GitHub shows it under "Cite this repository").

## License and acknowledgements

MIT License; see [LICENSE](LICENSE). The license covers the code only; the trial data are confidential.

We thank the patients who took part in METIMMOX, the trial investigators and clinical staff, and Frederick Thielen for contributions to the survival modelling code.
